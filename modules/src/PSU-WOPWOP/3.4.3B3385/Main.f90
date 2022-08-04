program main 
! PSU-WOPWOP 
! $Id: Main.f90 3384 2018-07-30 20:13:54Z brentner $ 
! $LastChangedDate: 2018-07-30 16:13:54 -0400 (Mon, 30 Jul 2018) $
! 
! This code has been developed through support from The Pennsylvania State  
! University and the following U.S. Government Contracts: 
!  
! 9/1/00 - 8/31/01  NASA Cooperative Agreement NCC-1-406   
! (NASA Langley Research Center; H. E Jones, Technical Officer)  
! 
! 2/1/01 - 1/31/06  NASA Ames Training Grant No. NGT-2-52275  
! (NASA Ames Research Center; Y. H. Yung, Technical Officer)  
! 
! 5/15/2003 - 5/14/2004  NASA LaRC Purchase Order L-71744D  
! (NASA Langley Research Center; C. L. Burley, Technical Officer)  
! 
! 3/1/2003 - 2/28/2006 NASA Langley Grant NAG-03025  
! (NASA Langley Research Center, D. P. Lockard, Technical Officer) 
!
! DARPA Helicopter Quieting Program
! Awarding Agency: Army Research Office
! Prime Award No: W911NF-04-1-0419 (to Georgia Tech)
! Subcontract Award No: E-16-X20-G1
! 
! 
! Written by Guillaume Bres, Guillaume Perez, Leonard Lopes, Hsuan-Nien Chen, 
! Christopher Hennes, Rui Cheng and Benjamin Goldman. 
! Faculty advisor Dr. Kenneth S. Brentner.
 
 
!* 
!      MAIN PROGRAM  
!REMARKS: The file must be compile with output.f90, observer.f90, 
!         aircraft.f90, rotor.f90, blade.f90, patch.f90, constants.f90 
!         COBdata.f90, matrixRotation.f90, freqAnalysis.f90, mathmod.f90 
!  - All the constants of the calculations are in constants.f90 
!  - The other input data is read in the namelist file 
!     - All output is done in output.f90 
!     - This main routine reads all the input via createAllInput() in this file 
!       it then allocates the noise arrays, calls computeNoise() in this file, 
!       and finaly writes all the ouput in output.f90.  And finally 
!       calls destroyEnvironment() in this file also. 
!- Renamed to main.f90 and modified by Leonard Lopes during the summer of 02 
!  created environment type, createEnvironment, writeEnvironment,  
!  destryoEnvironment.  Moved master do loop, a.k.a. integrate routine,  
!  to computeNoise routine.  Modified routines from Brentner's WOPWOP to  
!  create plot3DOutput routine which outputs .x and .fn files to be read  
!  in by fieldview or a program like it. 
!* 
 
  use observerObject
  use TimeHistoryObject
  use observerContainerObject
  use containerObject 
  use debugModule
  use wallObject
  use BroadbandModule
  use constantsModule
  use COBObject, only:appendCB, CBStructure
  use frequencyAnalysisModule, only:InitializeFreqAnalysisData, ConvertdBtoMSP
  use surfaceObject, only:GetSurfaceType
  use loadingObject, only:GetLoadType
  use mathModule 
  use SigmaModule
  use MPIModule
  use strings
  use IOModule, only:InitializeStreamNumberCounter, GetStreamNumber, MakeDirectory
  implicit none 
 
  type environment 
     type(observerContainer),dimension(:),pointer::obsContainer =>NULL()
     type(container), dimension(:), pointer::srcContArray =>NULL()
     integer::nbSourceContainers
     integer::nbObserverContainers
     logical::aperiodicData
     real(kind=8):: minObsTimeInData, maxObsTimeInData
  end type environment 
 
  type(environment)::env 
  real(kind=8)::timer, timer0, timer1, timer2, timer3
 
  integer :: caseUnit, unitNumber, i,stat, numCases
  character(len=4096)::caseNameFile
  character(len=4096)::pressureFileName, pressureFolderName, pressurePath 
  character(len=4096)::SPLFileName, SPLFolderName, SPLPath
  character(len=4096)::OASPLPath,OASPLFileName
  character(len=4096)::audioFolderName, audioFileName, audioPath 
  character(len=4096)::sigmaFolderName, sigmaFileName, sigmaPath
  character(len=4096)::phaseFileName,phasePath, &
                       complexPressureFileName,complexPressurePath
  character(len=128),external::GetBuildNumber 
   
  nameList / caseName / globalFolderName, caseNameFile 
   
  call InitializeProgram() ! This function only applies to MPI - change its name
  call InitializeStreamNumberCounter()
  call InitializeFreqAnalysisData()   

  timer0 = getTime() 
 !
  debugLevel=1  ! use this debugLevel until the enviroment namelist is read.
  if (IsMaster()) then
    write (*,'(A)') '********** PSU-WOPWOP Version 3.4.3 (PSU Internal Build '//&
                    trim(GetBuildNumber())//') ***********'
    if (.not.MasterDoesComputation()) then
      write (*,*) 'Using ', GetNbSlaves(), ' Slaves.' 
    end if
  end if
  caseUnit = GetStreamNumber ()
  open(caseUnit, file = 'cases.nam',status='old',iostat=stat)
  if (stat /= 0) then
    call Error ("Could not find 'cases.nam' in the current directory.")
  end if    
  numCases = 0
  do
    timer1 = getTime()
    !isAperiodic = .false.  !this array needs to be reset for each case.
    if (.not.newSegLHS .and. .not.newSegRHS) then
      dataLimitLHS=.false.
      dataLimitRHS=.false.
      read(caseUnit, nml=caseName, iostat=stat)
      if (stat == 0) then
        if (len_trim(globalFolderName) == 0 .or. len_trim(caseNameFile) == 0) then
          call Error ("Error in Main: globalFolderName and caseNameFile ",&
                       "must be specified in cases.nam.")
        end if
        if (CheckString(globalFolderName) /= 0 .or. CheckString(caseNameFile) /=0) then
          call Error ("Could not read globalFolderName or caseNameFile.", &
                      "Make sure they are in quotes and try again.")
        end if
        unitNumber = GetStreamNumber()
        ! If we've gotten to this point in the code then the 'stat' variable 
        ! has served its original purpose and can be recycled.  In this case it's being
        ! used to make sure that a case file opens properly and otherwise exits the code.
        call CreateEnvironment(env, trim(globalFolderName)//trim(caseNameFile), stat)
        if (ObserverParallelized().and.env%nbObserverContainers.gt.1) then
          call Error("Observer parallel must be run with only 1 observer container")
        end if
      else if(stat==-1) then
        exit
      else
        call Error('Could not read caseName', &
                   'Error came back as '//IntegerToString(stat))
      end if
      if (.not.MasterDoesComputation().and.env%nbSourceContainers.gt.1) then
        call Error('Parallel integrations must have environment which ', &
                   'contains only 1 source container')
      end if
      if (.not.MasterDoesComputation().and.env%nbObserverContainers.gt.1) then
        call Error('Parallel integrations must have environment which ', &
                   'contains only 1 observer container')
      end if
      if (.not.MasterDoesComputation().and.NotMaster())then
        call SetupSlaveProcessor(env) 
      end if
      if (debugLevel.ge.13) then
        call WriteEnvDebugInfo(env)
      end if
    
!   ## Not sure of what this check is supposed to accomplish.
!   ##
!    if (.not.MasterDoesComputation()) then
!      if (debugLevel.gt.3 .and. IsMaster()) then
!        write(*,*) 'Master is deallocating source'
!        do i=1, env%nbSourceContainers
!          call DestroyContainer(env%srcContArray(i))
!        end do
!      end if
!    end if
      call Barrier()
    else
      do i=1, env%nbObserverContainers
        env%obsContainer(i)%nbObsTimeExpansions = env%obsContainer(i)%nbObsTimeExpansions + 1
        call CreateNewObservers(env%obsContainer(i), .true.)
        call InitializeObservers(env%obsContainer(i),env%minObsTimeInData, env%maxObsTimeInData)
!        call CheckObserverTimeRange(env%obsContainer(i))
        if (newSegLHS .or. newSegRHS) then
          ! Check the segment parameters, we know the observer time range.
          ! If we dont know the time range we need to check this later:
          call CheckSegmentTimeSettings(env%obsContainer(i))
          ! Lastly, finish the observer container /observer setup
          call FinishObserverContainerSetup(env, env%obsContainer(i))
        else
          ! If we've determined that adding another observer segment
          ! would push us past the available data we must simply 
          ! destroy the observer we just created.
          call CreateNewObservers(env%obsContainer(i), .false.)
        end if
      end do
      ! We also need to reset the KeyCount as it currently sits at 4, 
      ! rather than 0 as it should at the beginning of a case.
      do i=1,size(env%srcContArray)
        call ResetContainerKeyCount(env%srcContArray(i))
      end do
    end if

    if (stat.eq.0) then
      if (.not.dataLimitRHS .or. newSegLHS .or. newSegRHS) then
        !ComputeNoise is the top level routine that runs through the source containers
        !  in Environment and integrates over all patches.  It stores
        !  acoustic pressure time histories in each observer's pPrime array.
        if (IsMaster()) then
          call UpdateTimer()
         ! write(*,*) 'Setup Time = ', timer 
         write(*,*) 'Case Setup Time = ',ExecutionTimer(timer1)
         timer2 = GetTime()
        end if
        close(unitNumber)
        ! Calculate the noise
        if(.not.readInPressureFlag)then
          call ComputeNoise(env) 
        else
          call ReadInPressure(env%obsContainer(1))
        end if

        call Barrier()

        if (isMaster()) then
          call UpdateTimer()
          if(debugLevel.ge.2)then
            !write (*,*) 'Computation Time = ', timer 
            write(*,*) 'Case Computation Time = ', ExecutionTimer(timer2)
            timer3 = GetTime()
          end if  
        end if
      end if
      !
      !This is where we have to calculate the SPL results. At this point in
      !  the code, the only information stored is pPrime, calculated in ComputeNoise.
      if (.not.ObserverParallelized()) then
        call CalculateResults(env)
      end if
      ! We either have enough information to capture PNLTM-10dB on both sides
      ! or we have enough to capture just the right side (having already checked the left)
      if (.not.(newSegLHS .or. newSegRHS) .or. dataLimitRHS) then
        if (isMaster()) then
          ! We need to hold off closing the data files until we are sure
          ! they don't need to be accessed again (if PNLTM-10dB fails).

          if( env%aperiodicData ) then
            do i=1,env%nbSourceContainers
              call CloseDataFiles(env%srcContArray(i))
            end do
            write(*,*) 'Closed Data Files.'
          end if
          call WriteResults(env)
          if(sigmaFlag)then    
            call WriteSigmaData (env%srcContArray, trim(sigmaPath), env%obsContainer(1))                               
          end if
          call UpdateTimer()
          if(debugLevel.ge.2.and.IsMaster())then
            !write (*,*) 'Total Time = ', timer 
            write(*,*) 'Case Postprocessing Time = ',ExecutionTimer(timer3)
            write(*,*) 'Total Case Execution Time = ',ExecutionTimer(timer1)
          end if  
          call Message ('Destroying Environment')
          numCases = numCases+1
        end if
        call DestroyEnvironment(env)  !ksb debug: turned this back on 3/2/2014
        !!ksb debug: 
        !newSegLHS = .false.
        !newSegRHS = .false.
      end if
    else
      exit
    end if
  end do

  close(caseUnit)

  if (IsMaster().and.debugLevel>=2) then
    if (numCases > 1) then
      !call Message ("Ran "//trim(IntegerToString(numCases))//" cases.")
      write (*,*) 'Ran '//trim(IntegerToString(numCases))//' cases.'
    end if
  end if
  if (IsMaster().and.debugLevel>0) then
    write (*,*) 'Total Execution Time = ', ExecutionTimer(timer0) 
  end if
  call FinalizeProgram()

  contains 
   
  subroutine UpdateTimer()
    timer=GetTime()-timer0
  end subroutine UpdateTimer
  
  function ExecutionTimer(start_time) result(execution_time)
    real(kind=8):: start_time, execution_time
    execution_time = GetTime() - start_time
  end function ExecutionTimer
    
  
  subroutine SetupSlaveProcessor(env) 
    type(environment), intent(inout)::env

    if (ObserverParallelized()) then
      call ResetSlaveObserver(env%obsContainer(1))
    else
      call Error("Unknown parallel integration method in SetupSlaveProcessor")
    end if
  end subroutine SetUpSlaveProcessor
      
  subroutine WriteEnvDebugInfo(env)
    type(environment), intent(in)::env
    integer::unitnum, stat
    unitnum = GetStreamNumber()
    open(unitnum, file='DebugInfoForProc_'//trim(IntegerToString(ProcID()))//'.dat', &
         form='formatted', status='unknown',iostat=stat)
    write(unitnum,*) '*** EnvironmentDebugInfo ***'
    write(unitnum,*) 'nbContainer = ', trim(integertostring(env%nbSourceContainers))
    write(unitnum,*) 'ContainerArray Associated? ', associated(env%srcContArray)
    if(associated(env%srcContArray)) then
      write(unitnum,*) 'Size of containerarray is ', & 
                       trim(integertostring(size(env%srcContArray)))
    end if
    write(unitnum,*) 'nbObserverContainers = ', env%nbObserverContainers
    write(unitnum,*) 'associated obsContainer ', associated(env%obsContainer)
    if(associated(env%obsContainer)) then
      write(unitnum,*) 'size of obsContainer ', &
                       trim(integertostring(size(env%obsContainer)))
    end if
    write(unitnum,*) '--- End WriteEnvDebugInfo ---'
    do i=1, size(env%srcContArray)
      call  WriteContainerDebugInfo(env%srcContArray(i),unitnum)
    end do 
    do i=1, size(env%obsContainer)
      call WriteObserverContainerDebugInfo(env%obsContainer(i),unitnum)
    end do
  end subroutine WriteEnvDebugInfo
  
  !** 
  !SUBROUTINE CreateEnvironment() 
  !  This subroutine reads all the input from the namelist file and initializes 
  !   the environment which contains all the information we need for integration,  
  !   including constants and physical data.  If observer is attached to an object  
  !   then this will create a COB list identical to the attached object and give it 
  !   to the observer. This also reads in all the file names for output including 
  !   tecfile names and 2D output file name.  Lastly this routine will read in all 
  !   Sigma flags for the fieldview files. 
  !ARGUMENTS: 
  !   - env: The environment to be created. 
  !   - nameFile: The name of the namelist. 
  !**
  subroutine CreateEnvironment(env, nameFile, stat)
    implicit none 
    type(environment),intent(inout)::env
    character(len=*),intent(in)::nameFile
    real::tauMin, tauMax, dTau
    real(kind=8):: temp
    character(len=32)::sigmaStructuredFormat, sigmaUnstructuredFormat
    character(len=1024):: FAACertFile
    logical:: computePegg, computeBPM
    integer::i, nbAircraft, nbWall, nbSourceContainers, nTau, stat, wallIteration,&
                  nbObserverContainers, nbContainer, nbObs, totalNbObs, FAAunitNumber, contNum, &
                  temp1, temp2
     
    namelist /EnvironmentIn/ nbAircraft, nbSourceContainers, nbWall, wallIteration,    &
                             nbObserverContainers,nbContainer,nbObs,                   & 
              !Output Flags 
              acousticPressureFlag, thicknessNoiseFlag,loadingNoiseFlag,totalNoiseFlag,&
              spectrumFlag,OASPLFlag,OASPLdBFlag,OASPLdBAFlag,SPLdBFlag,SPLFlag,       &
              SPLdBAFlag,broadbandFlag,narrowbandFlag,audioFlag,phaseFlag,sigmaFlag,   &
              debugLevel,PNLFlag,PNLTFlag,PartH36Flag,EPNLFlag,forceEPNL,              &
              complexPressureFlag, SELFLag, forceSEL,     &
              
              !Read in pressure flag
              readInPressureFlag, &

              !Folder and file names 
              pressureFolderName, pressureFileName, SPLFolderName, SPLFileName,        & 
              sigmaFolderName, sigmaFileName, audioFolderName, audioFileName,          & 
              OASPLFileName,phaseFileName,complexPressureFileName,                     &
              
              !Sigma Flags 
              sigmaStructuredFormat, sigmaUnstructuredFormat,                          &
              isomsThicknessNoiseFlag, iBlankFlag,                                     &
              loadingNoiseSigmaFlag, thicknessNoiseSigmaFlag,                          & 
              totalNoiseSigmaFlag, normalSigmaFlag, machSigmaFlag,                     & 
              observerSigmaFlag, observerContainerSigmaFlag,                           &
              velocitySigmaFlag, accelerationSigmaFlag,                                & 
              densitySigmaFlag, momentumSigmaFlag, pressureSigmaFlag,                  & 
              loadingSigmaFlag, areaSigmaFlag, MdotrSigmaFlag,iblankSigmaFlag,         &
              ASCIIOutputFlag,                                                         &
              
              !Pressure Gradient
              PressureGradientFlag,PressureGradient1AFlag,                             &
              ASCIIOutputFlag,                                                         &
              
              !FAA Noise Certification
              FAACertFlag, FAACertFile
 
    ! Set up the default values 
    nbAircraft              = 0          ; nbObserverContainers    = 0  
    nbWall                  = 0          ; nbObs                   = 0
    nbSourceContainers      = 0          ; nbContainer             = 0
    dTau                    = 0.0        ; wallIteration           = 0
    tauMin                  = 0.0        ; tauMax                  = 0.0
    nTau                    = 0          ; acousticPressureFlag    = .false. 
    OASPLdBFlag             = .false.    ; OASPLdBAFlag            = .false. 
    SPLdBFlag               = .false.    ; OASPLFlag               = .false.
    SPLFlag                 = .false.    ; SPLdBAFlag              = .false.
    broadbandFlag           = .false.    ; narrowBandFlag          = .true.
    spectrumFlag            = .false.    ; phaseFlag               = .false. 
    sigmaFlag               = .false.    ; audioFlag               = .false.
    pressureGradientFlag    = .false.    ; pressureGradient1AFlag  = .false.
    PNLFlag                 = .false.    ; EPNLFlag                = .false.
    complexPressureFlag     = .false.    ; PNLTFlag                = .false.
    readInPressureFlag      = .false.    ; SELFlag                 = .false.
    pressureFileName        = "pressure" ; pressureFolderName      = ""         
    SPLFileName             = "spl"      ; sigmaFolderName         = "" 
    sigmaFileName           = "sigma"    ; audioFolderName         = "" 
    OASPLFileName           = "OASPL"    ; SPLFolderName           = ""
    audioFileName           = "audio"    ; phaseFileName           = "phase"
    thicknessNoiseFlag      = .false.    ; complexPressureFileName = "complexPressure"
    loadingNoiseFlag        = .false.    ; totalNoiseFlag          = .false.
    iBlankFlag              = .true.     ; loadingNoiseSigmaFlag   = .false. 
    thicknessNoiseSigmaFlag = .false.    ; totalNoiseSigmaFlag     = .false. 
    normalSigmaFlag         = .false.    ; machSigmaFlag           = .false. 
    observerSigmaFlag       = .false.    ; velocitySigmaFlag       = .false. 
    observerContainerSigmaFlag=.false.
    accelerationSigmaFlag   = .false.    ; densitySigmaFlag        = .false. 
    momentumSigmaFlag       = .false.    ; pressureSigmaFlag       = .false. 
    loadingSigmaFlag        = .false.    ; debugLevel              = 1  
    isomsThicknessNoiseFlag = .false.    ; areaSigmaFlag           = .false.
    MdotrSigmaFlag          = .false.    ; iblankSigmaFlag         = .false.
    ASCIIOutputFlag         = .false.
    sigmaStructuredFormat   = 'PLOT3D'   ; FAACertFlag             = .false.
    sigmaUnstructuredFormat = 'FIELDVIEW'; FAACertFile             = ""
    newSegLHS               = .false.    ; newSegRHS               = .false.
    dataLimitLHS            = .false.    ; dataLimitRHS            = .false.
    PNLTdBDrop              = .true.     ; SELdBDrop               = .false.
    !isAperiodic             = .false.    ; 
    forceEPNL               = .false.    ; PartH36Flag             = .true.
    HitMaxExpansions        = .false.    ; globalPeggNoiseFlag     = .false.    
    globalBPMNoiseFlag      = .false.    ; EPNLPrime	       	   = .false.
    
    env%minObsTimeInData = -huge(env%minObsTimeInData)
    env%maxObsTimeInData = huge(env%maxObsTimeInData)
        
    ! This if statement checks the input specified in the environment 
    ! namelist, reads in the enviromental constants, creates the observer, and 
    ! reads in the container array only on the master node. 
    ! After this is done the master node braodcasts it all to the rest of the 
    ! nodes.
    
    open(unitNumber,file=namefile,status='old',iostat=stat)
    if (stat /= 0) then
      write(*,*) "Could not open '"//trim(namefile)//"'.  ", &
                 "Case will not be run."
      return
    end if
    write(*,*) ' '
    write(*,*) "*****************************************************************"
    write(*,*) "***  Reading case file: "//trim(namefile) 
    write(*,*) "*****************************************************************"

    read(unitNumber,nml=EnvironmentIn)
    
    if (FAACertFlag) then
      if (len_trim(FAACertFile).gt.0) then
        FAAunitNumber = GetStreamNumber()      
        open(FAAunitNumber,file=trim(trim(globalFolderName)//FAACertFile),status='old',iostat=stat)
        if (stat /= 0) then
          write(*,*) "Could not open '"//trim(FAACertFile)//"'.", &
                  "Case will not be run."
          return
        end if    
        read(FAAunitNumber, nml=EnvironmentIn)
      else
        FAAunitNumber        = 0
      
        nbObserverContainers = 1
        nbSourceContainers 	 = 1
        thicknessNoiseFlag 	 = .true.
        loadingNoiseFlag	 = .true.
        totalNoiseFlag	     = .true.
        PNLFlag		         = .true.
        PNLTFlag		     = .true.
        EPNLFlag		     = .true.
        spectrumFlag		 = .true.
        SELFlag		         = .true.        
      end if
    end if
    
    ! Rewind to the top of the file and reread EnvironmentIn
    ! in case the user wants to overwrite any of the FAACert defaults.
    ! ( This is quite necessary since we can't anticipate 
    rewind(unitNumber)
    read(unitNumber, nml=EnvironmentIn)
    
    if (debugLevel >=2.and.IsMaster())  then
      write(*,*) "Main namelist read."
    end if
    
    ! There is a bug were if SELFlag is .true. PNLTFlag must be true too.
    if( SELFlag ) then
      PNLTFlag = .true.
    end if
    
    if (acousticPressureFlag )then
      if (.not.(thicknessNoiseFlag .or. loadingNoiseFlag .or. totalNoiseFlag) ) then
        totalNoiseFlag=.true. 
      end if
    end if

    call SetDebugLevel(debuglevel)
    call CheckFileNames()   
    
    !Check the environment namelist for number of observer containers:
    if((nbObs.ne.0).and.(nbObserverContainers.ne.0))then
      call Error('Only specify nbObserverContainers in the environment namelist.',&
                 'nbObs is now called nbObserverContainers.', &
                 'Stopping.')
    end if    
    if(nbObs.gt.0)then
      call Notice('Use the variable nbObserverContainers in place of nbObs', &
                  'in the environment namelist in future runs.')
      nbObserverContainers = nbObs
      nbObs = 0
    end if
    if((nbObs.eq.0).and.(nbObserverContainers.eq.0)) then
      call Warning('No number of observer containers specified.' , &
                   'Defaulting nbObserverContainers in the Environment in ', &
                   'namelist to 1')
      nbObserverContainers = 1
    end if
    if(nbObserverContainers.lt.1) then
      call Error('A number of observer containers greater than 0',&
                 'must be specified in the environment namelist.', &
                 'Stopping.')
    end if
    env%nbObserverContainers = nbObserverContainers    
    
    ! Check the source bodies and make sure input makes sense 
    if((nbAircraft.ne.0).and.((nbSourceContainers.ne.0).or.(nbContainer.ne.0)) .or. & 
       (nbAircraft.lt.0) .or.(nbSourceContainers.lt.0)) then 
      call Error('Either a number of aircraft OR a number of source containers', &
                 'must be specified and must be greater than zero.', &
                 'Error in environment namelist. Stopping.')
    end if
    if(nbAircraft==0.and.nbSourceContainers==0.and.nbContainer==0) then
      call Error('A number of aircraft or a number of source', &
                 'containers must be specified.  Stopping.')
    end if
    if((nbContainer.ne.0).and.(nbSourceContainers.ne.0))then
      call Error('Only specify nbSourceContainers in the environment namelist.',&
                 'nbContainer is now called nbSourceContainers.', &
                 'Stopping.')
    end if
    if(nbContainer.ne.0)then
      call Notice('Use the variable nbSourceContainers in place of nbContainer', &
                  'in the environment namelist in future runs.') 
      env%nbSourceContainers=nbContainer
      else if(nbAircraft.ne.0)then
        env%nbSourceContainers = nbAircraft
      else
        env%nbSourceContainers = nbSourceContainers          
    end if
    
    ! Check to make sure the integration parameters make sense:
    !
    if (isomsThicknessNoiseFlag .and..not. loadingNoiseFlag) then
      call Error("Isom's thickness noise requires the loading noise flag ",&
                 "to be on. Edit the EnvironmentIn section of the namelist.")
    end if
    !This is maintaining backwards compatability:
    if(spectrumFlag.and.(.not.(SPLdBFlag.or.SPLdBAFlag.or.phaseFlag.or. &
                               complexPressureFlag)))then
      call Notice('spectrumFlag can now be used with SPLdBFlag, SPLdBAFlag, and phaseFlag.', &
                  'Turning SPLdBFlag on to output narrow band SPL.')
    end if
    if(phaseFlag.and.(.not.spectrumFlag))then
      call Notice('phaseFlag must be used with spectrumFlag. Turning spectrumFlag on.')
      spectrumFlag = .true.
    end if
    if(SPLFlag)then
      SPLdBFlag = .true.
    end if
    if(OASPLFlag)then
      OASPLdBFlag = .true.
    end if
   
    if(pressureGradientFlag)then
      call Notice('pressureGradient1AFlag is now used in place of pressureGradientFlag.', &
                  'turning on pressureGradient1AFlag.')
      pressureGradientFlag   = .false.
      pressureGradient1AFlag = .true.
    end if    

    ! This call was moved further down in this subroutine to allow for of a check of whether
    ! broadband data files exist.  (If not then their noise flags are turned off) BG 7/2/2012
!    call SetPPrimeIndices()
    
    ! Set up the contants and create the observer on all nodes 
    call SetEnvironmentalConstants (unitNumber)
    call InitializeParallelParameters()    
        
    if (debugLevel >= 2.and.IsMaster())then
      write (*,*) "Environment set."
    end if    
    nullify(env%obsContainer)
    allocate(env%obsContainer(env%nbObserverContainers)) 
    do i=1, env%nbObserverContainers
      call InitializeObserverContainer (env%obsContainer(i), unitNumber, FAAunitNumber, SPLPath, &
                                        OASPLPath,pressurePath,audioPath,sigmaPath,phasePath,&
                                        complexPressurePath, env%nbObserverContainers.gt.1)
    end do

    if (debugLevel >= 2.and.IsMaster())  then
      write (*,*) "Observer Container(s) and observer(s) created." 
    end if
    if (nbWall.lt.0) then 
      call Error('Invalid number of walls')
    end if
    mg_NbWall        = nbWall
    mg_WallIteration = wallIteration
    if (mg_NbWall.gt.0) then
      allocate(mg_WallArray(mg_NbWall)) 
      do i=1, mg_NbWall 
        call CreateWall(unitNumber, mg_WallArray(i)) 
      end do 
    end if
    if (nbWall.gt.0) write (*,*) "Wall(s) created"
       
    call Message ("Creating source objects:")
    ! Read in the source bodies, if a number of aircraft were specified the  
    ! backwards compatible routine is called. 
    ! This is only done on the master node and then broadcasted to all other 
    ! nodes
    allocate (env%srcContArray(env%nbSourceContainers))
    env%aperiodicData = .false.
    if(nbAircraft.gt.0) then       
      do i=1, env%nbSourceContainers 
        call CreateAircraftContainer(unitNumber, env%srcContArray(i), prefix="  ") 
      end do 
    else 
      do i=1, env%nbSourceContainers
        contNum = 1
        call CreateContainer(unitNumber, env%srcContArray(i), & 
                             dTau, tauMin, tauMax, nTau, &
                             temp1, temp2, prefix="  ") 
        if( AperiodicCont(env%srcContArray(i)) .or. env%aperiodicData ) env%aperiodicData = .true.
      end do
    end if
    if (globalPeggNoiseFlag .or. globalBPMNoiseFlag) then
      computePegg = .false.; computeBPM = .false.
      do i=1,env%nbSourceContainers
        call CheckBBFlags(env%srcContArray(i), computePegg, computeBPM)
        if (computePegg .and. computeBPM) exit
      end do
      if (.not.computePegg) globalPeggNoiseFlag = .false.
      if (.not.computeBPM)  globalBPMNoiseFlag  = .false.
    end if
    
    call SetPPrimeIndices()    
    
    call Message("Creating observer objects")
    totalNbObs = 0
    do i=1, env%nbObserverContainers
      call FinishObserverContainerSetup(env, env%obsContainer(i))
      call InitializeObservers(env%obsContainer(i),-huge(temp),huge(temp))
      totalNbObs = totalNbObs + GetNbObservers(env%obsContainer(i))
    end do
    
    do i=1, env%nbSourceContainers
      call CreateContainerObsTimeArrays(env%srcContArray(i), totalNbObs)
    end do
    
    call ResolveIntegrationMethod()
    
    if (sigmaFlag) then
      call SigmaCheck() 
    end if
    ! Parse the sigma file formats (these are variables in constantsModule):
    call StrToUpper (sigmaStructuredFormat, len_trim(sigmaStructuredFormat))
    if (index(sigmaStructuredFormat,"PLOT3D") /= 0) then
      sigmaStructuredFileType = FORMAT_PLOT3D
    else
      call Warning("Unrecognized structured file format: "//&
                   trim(sigmaStructuredFormat), "Defaulting to Plot3D.")
      sigmaStructuredFileType = FORMAT_PLOT3D
    end if
    
    call StrToUpper (sigmaUnstructuredFormat, len_trim(sigmaUnstructuredFormat))
    if (index(sigmaUnstructuredFormat,"FIELDVIEW") /= 0) then
      sigmaUnstructuredFileType = FORMAT_FIELDVIEW
    else
      call Warning("Unrecognized structured file format: "//&
                   trim(sigmaUnstructuredFormat), "Defaulting to Fieldview (FV-UNS).")
      sigmaUnstructuredFileType = FORMAT_FIELDVIEW
    end if    
    if (IsMaster()) then
      write(*,*) 'Environment Created'
    end if
  end subroutine CreateEnvironment
  
  
  subroutine ResolveIntegrationMethod()
    if (GetNbSlaves().eq.0) then
      call NotParallel()
    else 
      call ObserverParallel()
    end if
  end subroutine ResolveIntegrationMethod
  
  subroutine SigmaCheck()
    if(.not.NotParallelized()) then
      call Warning('Will not record sigma surface and run in parallel')
      sigmaFlag=.false.
      return
    end if
     if(env%nbObserverContainers.gt.1) then
       call Error('Cannot record sigma surface for more than one observer container.')
     end if
     if(env%obsContainer(1)%nbObs.gt.1)then
       call Warning('The sigma surface will be recorded but functions that are observer', & 
                    'point specific will only recorded for the first observer point.')
     end if
  end subroutine SigmaCheck
  
  function GetEnvObsContainer(env,i) result(obsCont)
    implicit none
    type(environment),intent(in)::env
    type(observerContainer)::obsCont
    integer::i
    obsCont=env%obsContainer(i)
  end function GetEnvObsContainer

  function GetNbObserverContainers(env) result(nbObsCont)
    implicit none
    type(environment),intent(in)::env
    integer::nbObsCont
    nbObsCont = env%nbObserverContainers
  end function GetNbObserverContainers
  
  function GetNbSourceContainers(env) result(nbSourceCont)
    implicit none
    type(environment),intent(in)::env
    integer::nbSourceCont
    nbSourceCont=env%nbSourceContainers
  end function GetNbSourceContainers
  
  function GetSourceContainerArray(env) result(contArray)
    implicit none
    type(environment),intent(in)::env
    type(container),dimension(:),pointer::contArray
    allocate(contArray(size(env%srcContArray)))
    contArray=env%srcContArray
  end function GetSourceContainerArray
  
  !**
  !SUBROUTINE ComputeNoise
  !  This subroutine computes the noise for each observer container
  !  by calling the routine ComputeObsContNoise for each observer container.
  !ARGUMENTS:
  ! - env: the environment we are calculating the noise in
  !**  
  subroutine ComputeNoise(env) 
    implicit none
    type(environment),intent(inout)::env 
    integer::i, obsCount
    if(debugLevel.gt.1.and.IsMaster())then
      write(*,*)'Starting ComputeNoise:'
    end if
    ! If we are in serial or there is only one processor allocated the master
    ! does the computations
    obsCount = 0
    if (NotParallelized()) then
      do i=1, env%nbObserverContainers
        call ComputeObsContNoise(env%srcContArray, env%obsContainer(i), obsCount, &
                                 env%minObsTimeInData, env%maxObsTimeInData)
        obsCount = obsCount + GetNbObservers(env%obsContainer(i))   
      end do
    else 
      if (ObserverParallelized()) then
        if (IsMaster()) then
          call Message('Running Observer Parallel Integration')
        end if
        do i=1, env%nbObserverContainers
          if (GetNbSlaves().gt.GetNbObservers(env%obsContainer(i)).and.IsMaster()) then
            call Error('Running observer parallel with more slave', &
                       'processors than observer points', &
                       'Number of slaves = '//IntegerToString(GetNbSlaves()), &
                       'Number of observers = '// &
                       IntegerToString(GetNbObservers(env%obsContainer(i))))
          end if
          !call SetupParallelObserverContainer(env%obsContainer(i))
          call ObserverParallelIntegration(env, env%obsContainer(i), obsCount)
          obsCount = obsCount + GetNbObservers(env%obsContainer(i))
        end do
      else if (SourceParallelized()) then
        if (IsMaster()) then
          call Message('Running Source Parallel Integration')
        end if
        do i=1,GetNbObserverContainers(env)
          call SourceParallelIntegration(env, env%obsContainer(i))
        end do
      else if (BothParallelized()) then
        if (IsMaster()) then
          call Message('Running Observer And Source Parallel Integration')
        end if
        do i=1,GetNbObserverContainers(env)
          call BothMethodsParallelIntegration(env, env%obsContainer(i))
        end do
      else
        call Error('Unknown parallel method type')
      end if
    end if
    if(debugLevel.gt.1.and.IsMaster())then
      write(*,*)'Finished Computing Pressure.'
    end if
  end subroutine ComputeNoise
  
  subroutine BothMethodsParallelIntegration(env, obsCont)
    type(environment)::env
    type(observerContainer)::obsCont
    env=env
    obsCont=obsCont
    call Error("Should never get to BothMethodsParallelIntegration() in Main.f90")
  end subroutine BothMethodsParallelIntegration  
  
  subroutine SourceParallelIntegration(env, obsCont)
    type(environment)::env
    type(observerContainer)::obsCont
    env=env
    obsCont=obsCont
    call Error("Source parallel not working")
  end subroutine SourceParallelIntegration
    
  subroutine ObserverParallelIntegration(env, obsCont, obsCount)
    type(environment)::env
    type(observerContainer)::obsCont
    integer:: obsCount
  
    integer::k, nbTotalPoints, rcvPointIndex, limit 
    integer::pointIndex, tag, sender, master 
    logical::jobDoneFlag 
    jobDoneFlag=.false.
    if (IsMaster()) then 
      nbTotalPoints = obsCont%nbObs
      pointIndex    = 0 
      call SetMasterValue(master)
      if (GetnbSlaves().le.nbTotalPoints) then
        limit = GetNbSlaves()
      else
        limit = nbTotalPoints
      end if
      do k=1, limit
        if (debugLevel.ge.11) then
          write(*,*) 'Master init observer sending'
        end if
        call ProceedToNextPoint(obsCont, k, master, pointIndex)  
      end do
      ! Over each observer point the results are
      ! returned, and the next point is sent out if the jobs isnt done 
      do k=1, nbTotalPoints
        if (debugLevel.ge.11) then
          write(*,*) 'Master awaiting integer from slave'
        end if
        call ReceiveInteger(rcvPointIndex, GetAnyTag(), GetAnySource()) 
        if (debugLevel.ge.11) then
          write(*,*) 'Master rcvPointIndex received ', rcvPointIndex
        end if
        call SetSenderValue(sender)
        if (debugLevel.ge.11) then
          write(*,*) 'Master Set Sender Value'
        end if
        call ReceiveResultsFromSlave(obsCont, master, sender, rcvPointIndex)
        if (debugLevel.ge.11) then
          write(*,*) 'Received results from slave'
        end if
        call CheckIfDone(pointIndex.lt.nbTotalPoints, sender, master, jobDoneFlag)
        if (debugLevel.ge.11) then
          write(*,*) 'Master job done flag', jobDoneFlag
        end if
        if(.not.jobDoneFlag) then
          call ProceedToNextPoint(obsCont, sender, master, pointIndex)  
          if (debugLevel.ge.11) then
            write(*,*) 'Master proceeded to next point'
          end if
        end if
      end do
      !ksb debug: 
      !call ResetPPrimeIndices()
    else 
      ! This is the slave part of the code, it receives an observer point and 
      ! then integates and returns the result. 
      ! When it returns the result it waits until the next point is seant or 
      ! ends.
      do while (.not.jobDoneFlag)
        call SetMasterValue(master)
        if (debugLevel.ge.11) then
          write(*,*) ProcID(), 'slave master value set'
        end if
        call ReceiveObserverPoint(obsCont, GetAnyTag(), master)
        if (debugLevel.ge.11) then
          write(*,*) ProcID(), 'observer point received'
        end if
        call SetTagValue(tag) 
        if (debugLevel.ge.11) then
          write(*,*) ProcID(), 'tag value set'
        end if
        call ReceiveInteger(rcvPointIndex, tag, master)
        call SetObserverIndex(obsCont%obs(1),rcvPointIndex)
        if (debugLevel.ge.11) then
          write(*,*) ProcID(), 'rcvPointIndex received ', rcvPointIndex
        end if
        obsCount=0 ! Only one observer goes to each slave processor,
        call ComputeObsContNoise(env%srcContArray, obsCont, obsCount, &
                                  env%minObsTimeInData, env%maxObsTimeInData)
        if (debugLevel.ge.11) then
          write(*,*) ProcID(), 'observer integrated'
        end if

        call CalculateResults(env) 
        if (debugLevel.ge.11) then
          write(*,*) ProcID(), 'observer results calculated'
        end if
        call SendInteger(rcvPointIndex, master, tag)
        if (debugLevel.ge.11) then
          write(*,*) ProcID(), 'sent point index'
        end if

        call SendResultsToMaster(obsCont, master, tag)
        if (debugLevel.ge.11) then
          write(*,*) ProcID(), 'results sent to master'
        end if
        call ReceiveLogical(jobDoneFlag, tag, master)
        if (debugLevel.ge.11) then
          write(*,*) ProcID(), 'job done flag', jobDoneFlag
        end if
        !ksb debug:
        call ResetObserver(obsCont%obs(1),obsCont%tMin,obsCont%tMax,&
          obsCont%nt,obsCont%dt)
      end do
    end if
    if (debugLevel.ge.11) then
      write(*,*) myid, 'leaving observer routine'
    end if
  end subroutine ObserverParallelIntegration
  
  ! If it is not done it just continues by sending the jobDoneFlag out.
  subroutine CheckIfDone(flag, proc, master, jobDoneFlag) 
    integer, intent(in)::proc, master 
    logical::jobDoneFlag, flag 
    jobDoneFlag=.not.flag
    if(debugLevel.ge.11) then
      write(*,*) 'Master sending jobDoneFlag', jobDoneFlag
    end if
    call SendLogical(jobDoneFlag, proc, master) 
  end subroutine CheckIfDone 
     
  ! This routine is done on the master, it takes in a point index and then sends
  ! the point to an available node.
  subroutine ProceedToNextPoint(obsCont, proc, master, pointIndex) 
    type(observerContainer), intent(in)::obsCont 
    integer, intent(in)::proc, master 
    integer, intent(inout)::pointIndex 
   
    pointIndex = pointIndex + 1
    call PrintPointToScreen(obsCont, pointIndex)  
    call SendObserverPoint(obsCont, pointIndex, proc, master)
    if (debugLevel.gt.11) then
      write(*,*) 'Master sent observer to slave'
    end if
    call SendInteger(pointIndex, proc, master) 
    if (debugLevel.gt.11) then
      write(*,*) 'Master sent pointIndex'
    end if
  end subroutine ProceedToNextPoint 
  
  !**
  !SUBROUTINE ComputeObsContNoise
  !  This subroutine computes the noise for each observer by calling the routine 
  !  ComputeObsThickAndLoadNoise for each observer in the container.
  !ARGUMENTS:
  ! - sourceArray: the array of noise sources
  ! - obsCont: the observer container
  !** 
  subroutine ComputeObsContNoise(sourceArray, obsCont, obsCount, minObsTimeInData, maxObsTimeInData)
    implicit none
    type(container), dimension(:), intent(inout)::sourceArray
    type(observerContainer), intent(inout)::obsCont
    integer::i, obsCount, oldSegments, j
    real(kind=8)::minObsTimeInData, maxObsTimeInData
    
    if (IsMaster()) then
      write(*,*)'Computing noise for '//trim(obsCont%title)
    end if
    do i=1,GetNbObservers(obsCont)
      if(IsMaster()) then
        call PrintPointToScreen(obsCont, i)
      end if        
      if (obsCont%newObs(i)) then
        obsCount = obsCount + 1
        if (newSegLHS .or. newSegRHS) then
          if (globalPeggNoiseFlag .or. globalBPMNoiseFlag) then
            oldSegments = size(obsCont%obs(1)%mspBB)
          elseif (AtmAtten) then
            oldSegments = size(obsCont%obs(1)%radDistance)
          end if
        end if
        call ComputeObsThickAndLoadNoise(sourceArray, obsCont%obs(i), obsCont%CBList, &
                                         obsCont%nbSegments, obsCont%segmentSize, oldSegments, &
					                     obsCount, minObsTimeInData, maxObsTimeInData)
      end if
    end do
  end subroutine ComputeObsContNoise
      
  !**
  !SUBROUTINE ComputeObsThickAndLoadNoise
  !  This subroutine computes the thickness and loading noise for an observer
  !  by integrating over the source containers.
  !ARGUMENTS:
  ! - sourceArray: the array of noise sources
  ! - obs: the observer
  ! - COBList: the motion of the observer container
  !**     
  subroutine ComputeObsThickAndLoadNoise(sourceArray, obs, obsFrameMotion, nbSegments, segmentSize, &
                                         oldSegments, obsCount, minObsTimeInData, maxObsTimeInData)
    implicit none
    type(container),dimension(:),intent(inout)::sourceArray
    type(CBStructure),pointer::obsFrameMotion   
    type(observer),intent(inout)::obs
    integer::nt,i,j, obsCount, iCont, nbCont, nbSrc, nSeg, BBarraySize
    integer, intent(inout):: nbSegments, oldSegments
    real:: tMin,tMax,dt, segmentSize
    real(kind=8)::dtKind8, minObsTimeInData, maxObsTimeInData
    real        , dimension(:), allocatable::temp, freqArray
    real(kind=8), dimension(:), pointer::thickness, loading
    type(vectorKind8), dimension(:), pointer::PTGradient , PL1Gradient
    real, dimension(:,:), pointer:: IntensityKind4
    real(kind=8), dimension(:,:), pointer:: PeggIntensity=>null()
    real(kind=8), dimension(:,:,:), pointer:: BPMintensity=>null()
    real(kind=8), dimension(:), pointer:: srcTMinArray, srcTMaxArray
    real(kind=8):: A
    real, dimension(:), pointer:: radDistance
    character(len=1024):: Title
    type(Pegg), pointer:: PeggData=>null()
    type(BPM), pointer:: BPMData=>null()
    
    A=0.
    nSeg = 0 !ksb debug
    nt   = GetObserverNt(obs)
    tMin = GetObserverTMin(obs)
    tMax = GetObserverTMax(obs)
    dtKind8 = dble(tMax - tMin)/dble(nt-1)
    dt = dtKind8

    !Here we are allocating the arrays to be filled in the integratePatch routine:
    if(thicknessNoiseFlag.or.totalNoiseFlag)then
      if(pressureGradientFlag.or.pressureGradient1AFlag)then
        allocate(PTGradient(nt))
        PTGradient = vectorSetKind8Coordinates(A,A,A)
      end if
      allocate(thickness(nt))
      thickness = A
    end if
    if(loadingNoiseFlag.or.totalNoiseFlag)then
      if(pressureGradientFlag.or.pressureGradient1AFlag)then
        allocate(PL1Gradient(nt))
        PL1Gradient = vectorSetKind8Coordinates(A,A,A)
      end if
      allocate(loading(nt))
      loading = A
    end if

    if (globalPeggNoiseFlag .or. globalBPMNoiseFlag .or. AtmAtten) then
      if (nbSegments.eq.0) then
        nSeg = 1
      else
        nSeg = nbSegments
        if (newSegLHS .or. newSegRHS) nSeg = nSeg-oldSegments
      end if

      if (AtmAtten) then
        allocate(radDistance(nSeg), obs%radDistance(nseg))
        radDistance = A
        do i=1,nseg
          nullify(obs%radDistance(i)%t, obs%radDistance(i)%f)
          allocate(obs%radDistance(i)%t(1), obs%radDistance(i)%f(1))
          obs%radDistance(i)%t = A
          obs%radDistance(i)%f = A
          obs%radDistance(i)%nt = 1
          obs%radDistance(i)%title = 'Time '//integertostring(i)
        end do
      end if

      if (globalPeggNoiseFlag .or. globalBPMNoiseFlag) then
        if (globalPeggNoiseFlag) then
          allocate(PeggIntensity(BB_THIRD_OCTAVE_BANDS, nSeg))
          PeggIntensity = A
        elseif (globalBPMNoiseFlag) then
          allocate(BPMintensity(BPM_TERMS, BB_THIRD_OCTAVE_BANDS, nSeg))
          BPMintensity = A
        end if
      end if
    end if
    
    if(.not.associated(obs%iBlankArray)) then
      allocate(obs%iBlankArray(obs%nt))
    end if
    obs%iBlankArray = 0
    !ksb debug: 10/10/2017
    firstObsSetup=.true.
    do i=1,size(sourceArray)
      iCont = -1
      nbCont = 0
      call IntegrateContainer(sourceArray(i), GetObserverCoordinates(obs),          &
                              obsFrameMotion, thickness, loading,                   &
                              PTGradient, PL1Gradient, PeggIntensity, BPMIntensity, &
                              obs%iBlankArray, dtKind8, obsCount,                   &
                              CreateEvenlySpacedArray(tMin,tMax,nt), iCont, nbCont, &
                              nSeg, segmentSize, radDistance, PeggData, BPMData)

      minObsTimeInData = max(minObsTimeInData, getContainerMinObsTimeInData(sourceArray(i)))
      maxObsTimeInData = min(maxObsTimeInData, getContainerMaxObsTimeInData(sourceArray(i)))
    end do
    !  The temp array used in this routine is REQUIRED for compiling with 
    !  gFortran. It is a work around to passing the pressure gradient
    !  arrays to subroutines, and SHOULD NOT BE REMOVED. -JPE 3/14/2008
    allocate(temp(nt))
    
    !!!!!!!!!!!!!!!!!!!!!!!
    if (.not.associated(obs%pPrime)) then
      allocate(obs%pPrime(GetSizeOfPPrime()))
      do i=1, GetSizeOfPPrime()
        nullify(obs%pPrime(i)%t)
        nullify(obs%pPrime(i)%f)
      end do
    end if
    do i=1,GetSizeOfPPrime()
      call SetupTimeHistory(obs%pPrime(i), tMin, tMax, obs%nt, PPRIMETITLEARRAY(i))
    end do
    if(thicknessNoiseFlag.or.totalNoiseFlag)then
      temp = thickness
      call SettHFArray(obs%pPrime(THICK_APTH), temp)
      deallocate(thickness)
      if(pressureGradientFlag.or.pressureGradient1AFlag)then
        !Do not remove the temp array! (see note above) -JPE 3/14/2008
        temp = PTGradient(:)%A(1)
        call SettHFArray(obs%pPrime(THICK_PGX), temp)
        temp = PTGradient(:)%A(2)
        call SettHFArray(obs%pPrime(THICK_PGY), temp)
        temp = PTGradient(:)%A(3)
        call SettHFArray(obs%pPrime(THICK_PGZ), temp)
        deallocate(PTGradient)
      end if
    end if
    if(loadingNoiseFlag.or.totalNoiseFlag)then
      temp=loading
      call SettHFArray(obs%pPrime(LOAD_APTH),temp)
      deallocate(loading)
      if(pressureGradientFlag.or.pressureGradient1AFlag)then
        temp = PL1Gradient(:)%A(1)
        call SettHFArray(obs%pPrime(LOAD_PGX), temp)
        temp = PL1Gradient(:)%A(2)
        call SettHFArray(obs%pPrime(LOAD_PGY), temp)
        temp = PL1Gradient(:)%A(3)
        call SettHFArray(obs%pPrime(LOAD_PGZ), temp)
        deallocate(PL1Gradient)
      end if
    end if

    deallocate(temp)
    allocate(temp(nSeg))
    temp = CreateEvenlySpacedArray(tMin,tMax,nSeg+1)
    if (globalPeggNoiseFlag .or. globalBPMNoiseFlag) then
      ! Currently this array is 2D to allow the different BPM
      ! components to be output, but this structure also allows 
      ! the use of multiple broadband methods in the future.
      if (globalPeggNoiseFlag) then
        BBarraySize = 1
      elseif (globalBPMNoiseFlag) then
        BBarraySize = BPM_TERMS
      end if
      allocate(obs%mspBB(BBarraySize, nSeg), IntensityKind4(BB_THIRD_OCTAVE_BANDS, nSeg), freqArray(BB_THIRD_OCTAVE_BANDS))
      do i=1,size(obs%mspBB,1)
        if (globalPeggNoiseFlag) then
          IntensityKind4 = PeggIntensity
        elseif (globalBPMNoiseFlag) then
          IntensityKind4 = BPMIntensity(i,:,:)
        end if
        do j=1,size(obs%mspBB,2)
          allocate(obs%mspBB(i,j)%freq(BB_THIRD_OCTAVE_BANDS), obs%mspBB(i,j)%f(BB_THIRD_OCTAVE_BANDS))
          if (globalPeggNoiseFlag) then
            title = 'Pegg-Total (Pa^2)'
          elseif (globalBPMNoiseFlag) then
            if (i.eq.1) then
              title = 'BPM-TBLTEp (Pa^2)'
            elseif (i.eq.2) then
	      title = 'BPM-TBLTEs (Pa^2)'
            elseif (i.eq.3) then
	      title = 'BPM-TBLTEalpha (Pa^2)'
            elseif (i.eq.4) then
	      title = 'BPM-LBLVS (Pa^2)'
            elseif (i.eq.5) then
              title = 'BPM-Blunt (Pa^2)'
            elseif (i.eq.6) then
              title = 'BPM-Tip (Pa^2)'
            elseif (i.eq.7) then
              title = 'BPM-Total (Pa^2)'
            end if
          end if
          call SetFreqDomainTitle(obs%mspBB(i,j), trim(Title))
          call SetFreqDomainFArray(obs%mspBB(i,j), IntensityKind4(:,j))
          freqArray = GetThirdOctFreqArray()
          call SetFreqArray(obs%mspBB(i,j), freqArray)
          call SetNf(obs%mspBB(i,j), BB_THIRD_OCTAVE_BANDS)
          call SetMinFreq(obs%mspBB(i,j), freqArray(1))
          call SetMaxFreq(obs%mspBB(i,j), freqArray(BB_THIRD_OCTAVE_BANDS))
        end do
      end do
      deallocate(freqArray, IntensityKind4)
    end if
    if (AtmAtten) then
      do i=1,nSeg
        call SetupTimeHistory(obs%radDistance(i), temp(i), temp(i), 1, 'Radiation Distance (m)')
        obs%radDistance(i)%t(1) = temp(i)
        obs%radDistance(i)%f(1) = radDistance(i)
      end do
      deallocate(radDistance)
    end if

    deallocate(temp)
    if (globalPeggNoiseFlag) then
      deallocate(PeggIntensity)
      nullify(PeggIntensity)
    end if
    if (globalBPMNoiseFlag) then
      deallocate(BPMintensity)
      nullify(BPMIntensity)
    end if
  end subroutine ComputeObsThickAndLoadNoise

  !**
  !SUBROUTINE CalculateResults
  !  This is the top level routine that calculates all the 
  !  SPL results for each observer container in the environment.
  !ARGUMENTS:
  ! - env: the environment we are calculating results for
  subroutine CalculateResults(env)
    implicit none
    type(environment),intent(inout)::env
    integer::i
    if(debugLevel.gt.1.and.IsMaster())then
      write(*,*)'Starting to calculate results.'
    end if
    !This loops through each observer container stored in environment:
    !  It passes the same file path names but a different observer container 
    !  each time.
    do i=1,env%nbObserverContainers
      call CalcObsContResults(env%obsContainer(i))
    end do

    do i=1, env%nbObserverContainers
      if (newSegLHS.or.newSegRHS) then
        if (env%obsContainer(i)%nbObsTimeExpansions .ge. env%obsContainer(i)%maxObsTimeExpansions) then
          call Warning('The maximum number of observer time-range expansions',&
                       'has been exceeded.')
          !ksb debug:
          !print*,'env%obsContainer(i)%nbObsTimeExpansions= ',env%obsContainer(i)%nbObsTimeExpansions
          !print*,'env%obsContainer(i)%maxObsTimeExpansions= ', env%obsContainer(i)%maxObsTimeExpansions
          !ksb debug: end debug statements
          if (EPNLFlag .and. forceEPNL .and. .not.all(PNLTdBDrop)) &
              call Notice('The EPNL value output in this case is not a valid measure',&
                 'of EPNL for FAA Noise Certification.')
          if (SELFlag .and. forceSEL .and. .not.all(SELdBDrop))  &
              call Notice('The SEL value output in this case is not a valid measure',&
                 'of SEL for FAA Noise Certification.')
          newSegLHS=.false.
          newSegRHS=.false.
          dataLimitLHS=.true.
          dataLimitRHS=.true.
          HitMaxExpansions=.true.
          exit
        end if
      end if
    end do

    ! If we are going to allow a single case to be run multiple times to capture 
    ! PNLT-10dB (or whatever) then we need to wait to filter the pPrime arrays 
    ! until the very end.
    if (.not.(newSegLHS .or. newSegRHS) .or. dataLimitRHS) then
      do i=1, env%nbObserverContainers
        call FilterObsContPressure(env%obsContainer(i))
      end do
    end if
    
    if(debugLevel.gt.1.and.IsMaster())then
      write(*,*)'Finished calculating results.'
    end if
  end subroutine CalculateResults 
  
  !**
  !SUBROUTINE WriteResults
  !  This is the top level routine that writes all the 
  !  SPL results for each observer container in the environment.
  !ARGUMENTS:
  ! - env: the environment we are calculating results for
  subroutine WriteResults(env)
    implicit none
    type(environment),intent(inout)::env
    integer::i
    if(debugLevel.gt.1)then
      write(*,*)'Starting to write results.'
    end if
    !This loops through each observer container stored in environment:
    !  It passes the same file path names but a different observer container 
    !  each time.
    do i=1,env%nbObserverContainers
      call WriteObsContResults(env%obsContainer(i))
    end do    
    if(debugLevel.gt.1)then
      write(*,*)'Finished writing results.'
    end if
  end subroutine WriteResults   
 
  function GetNBenvironmentSurfaces (env) result (nbSurfaces)
    type(environment)::env
    integer :: i, nbSurfaces, nbStructuredSurfaces, nbUnstructuredSurfaces
    nbStructuredSurfaces = 0
    nbUnstructuredSurfaces = 0
    do i=1,env%nbSourceContainers
      call GetNbContainerSurfaces (env%srcContArray(i), &
                                   nbStructuredSurfaces, &
                                   nbUnstructuredSurfaces)
    end do
    nbSurfaces = nbStructuredSurfaces + nbUnstructuredSurfaces
  end function GetNBenvironmentSurfaces
   
  ! SImple Print rountine
  subroutine PrintSourcePositionToScreen(env, sourceIndex)
   type(environment)::env
   integer::sourceIndex
   if (debugLevel.gt.1) then
     write(*,'(A,I5,A,I5)') ' Integrating source ', sourceIndex, &
                            ' of ', GetNBenvironmentSurfaces(env)
   end if
  end subroutine PrintSourcePositionToScreen
  
 ! SImple Print rountine
  subroutine PrintPointToScreen(obsCont, i)
    type(observerContainer), intent(in)::obsCont
    integer, intent(in)::i
   
    if (debugLevel.ge.1) then
      write(*, '(A,I5,A,2(I5),A,I5)') ' Grid point:', &
                                      i, ' of ', GetNbObservers(obsCont) 
      write(*, '(A,3(G15.6))') ' At Location:',  GetObserverCoordinates(obsCont%obs(i)) 
    end if
  end subroutine PrintPointToScreen
 
  !** 
  !SUBROUTINE destroyEnvironment() 
  !    This subroutine destroys all allocated arrays in the system by calling  
  !  the environments children. 
  !ARGUMENTS: 
  !  - env: The environment to be destroyed. 
  !** 
  subroutine DestroyEnvironment(env) 
    type(environment),intent(inout)::env 
    integer::i

    do i=1, env%nbSourceContainers 
      call DestroyContainer(env%srcContArray(i)) 
    end do 

    do i=1, env%nbObserverContainers
      call DestroyObserverContainer(env%obsContainer(i))
    end do
    
    if(associated(env%srcContArray))then
      deallocate(env%srcContArray) 
      nullify(env%srcContArray)
    end if
    if(associated(env%obsContainer))then
      deallocate(env%obsContainer) 
      nullify(env%obsContainer)
    end if
    if(allocated(mg_wallArray))then
      deallocate(mg_wallArray) 
    end if
    if(allocated(pPrimeTitleArray))then
      deallocate(pPrimeTitleArray)
    end if
    if(IsMaster()) then
      write(*,*)'Environment Destroyed.'
    end if
  end subroutine DestroyEnvironment 
 
  subroutine CheckFileNames()
    ! Ensure that the user has specified filenames that are at least 
    ! somewhat valid (this is not a robust check yet) 
    if (acousticPressureFlag.and.len_trim(pressureFileName) == 0) then
      call Warning('Pressure filename not properly specified. Defaulting.')
      pressureFileName = "pressure"
    end if     
     
    if ((OASPLdBFlag.or.OASPLdBAFlag.or.SPLdBFlag).and.len_trim(SPLFileName) == 0) then 
      call Warning('SPL filename not properly specified. Defaulting.')
      SPLFileName = "spl"
    end if 
    
    if (audioFlag .and. len_trim(audioFileName) == 0) then 
      call Warning('Audio filename not properly specified. Defaulting.')
      audioFileName = "audio" 
    end if 
       
    if (sigmaFlag .and. len_trim(sigmaFileName) == 0) then 
      call Warning('Sigma surface filename not properly specified. Defaulting.')
      sigmaFileName = "sigma_surface" 
    end if 

    if (phaseFlag .and. len_trim(phaseFileName) == 0) then 
      call Warning('Phase filename not properly specified. Defaulting.')
      phaseFileName = "phase" 
    end if     
    
    if (complexPressureFlag .and. len_trim(complexPressureFileName) == 0) then 
      call Warning('Complex pressure filename not properly specified. Defaulting.')
      complexPressureFileName = "complexPressure" 
    end if

    ! Check filename path lengths 
    if (acousticPressureFlag) then
      if (len_trim(trim(globalFolderName)//trim(pressureFolderName)// &
                   trim(pressureFileName))+4 > 4096) then
        call Error('Path to pressure file to long.')
      end if
      call MakeDirectory (trim(globalFolderName)//trim(pressureFolderName))
    end if

    if (OASPLdBFlag.or.OASPLdBAFlag.or.SPLdBFlag.or.PNLFlag.or.PNLTFlag.or. &
        SPLdBAFlag.or.complexPressureFlag.or.phaseFlag) then
      if (len_trim(trim(globalFolderName)//trim(SPLFolderName)//      &
                   trim(SPLFilename))+4>4096) then                        
        call Error('Path to SPL file to long.')
      end if
      call MakeDirectory (trim(globalFolderName)//trim(SPLFolderName))
      if(spectrumFlag .eqv. .false.)then
        call Warning('No spectral data will be output because spectrum flag is false.')
      end if
    end if
    
    if (audioFlag) then
      if (len_trim(trim(globalFolderName)//trim(audioFolderName)//    & 
                   trim(audioFileName))+4 > 4096) then 
        call Error('Path to audio file to long.')
      end if
      call MakeDirectory (trim(globalFolderName)//trim(audioFolderName))
    end if
    if (sigmaFlag) then
      if (len_trim(trim(globalFolderName)//trim(sigmaFolderName)//    & 
                   trim(sigmaFileName))+4 > 4096) then 
        call Error('Path to sigma surface file to long.')
      end if
      call MakeDirectory (trim(globalFolderName)//trim(sigmaFolderName))
    end if
 
    ! If the names are fine, save them off 
    ! NOTE: The write statement is ensuring that the most that is put into the 
    ! string is 4096. This is basically just to keep the compiler quiet - it 
    ! knows that each of the component strings could be up to 4096 individually, 
    ! so it thinks that the length of the strings we're storing into should be 
    ! 'three times that. The checks above ensure that that is not the case. 
    write (SPLPath,      *)       trim((globalFolderName))//  &
                                  trim((SPLFolderName))//     &
                                  trim((SPLFilename)) 
    write (OASPLPath,    *)       trim((globalFolderName))//  &
                                  trim((SPLFolderName))//     &
                                  trim((OASPLFileName))                                   
    write (pressurePath, *)       trim((globalFolderName))//  &
                                  trim((pressureFolderName))//&
                                  trim((pressureFileName)) 
    write (audioPath,    *)       trim((globalFolderName))//  &
                                  trim((audioFolderName))//   &
                                  trim((audioFileName)) 
    write (sigmaPath,    *)       trim((globalFolderName))//  &
                                  trim((sigmaFolderName))//   &
                                  trim((sigmaFileName))
    write (phasePath,    *)       trim((globalFolderName))//  &
                                  trim((SPLFolderName))//     &
                                  trim((phaseFileName))
    write (complexPressurePath,*) trim((globalFolderName))//  & 
                                  trim((SPLFolderName))//     &
                                  trim((complexPressureFileName))
                                  
    ! Move the paths to the left of the string buffer: 
    SPLPath            = adjustl(SPLPath)
    OASPLPath          = adjustl(OASPLPath)
    pressurePath       = adjustl(pressurePath)
    audioPath          = adjustl(audioPath) 
    sigmaPath          = adjustl(sigmaPath)
    phasePath          = adjustl(phasePath)
    complexPressurePath= adjustl(complexPressurePath)
  end subroutine CheckFileNames 
       
  ! A function to ensure that all characters in a string are valid ASCII
  ! Printable characters.
  function CheckString (string) result (pos)
    character(len=*),intent(in) :: string
    integer :: pos

    character(len=26),parameter:: uppercase="ABCDEFGHIJKLMNOPQRSTUVWXYZ", &
                                  lowercase="abcdefghijklmnopqrstuvwxyz"
    character(len=32),parameter:: punctuation="`~!@#$%^&*()-_=+[{]}\|;:',<.>/? "
    character(len=10),parameter:: numbers="0123456789"
    character(len=1),parameter:: quote='"'

    pos = VERIFY(string,uppercase//lowercase//punctuation//numbers//quote)
    
  end function CheckString
  
  function GetNbContainers(env) result(nbCont)
    type(environment),intent(in)::env
    integer::nbCont
    nbCont=env%nbSourceContainers
  end function GetNbContainers
  
  subroutine FinishObserverContainerSetup(env, obsCont)
    type(environment),intent(inout)::env 
    type(observerContainer)::obsCont
    integer::i
    logical::success

    ! This finalized the observer 
    if(trim(obsCont%attachedTo) /= '') then 
      success=.false.
      do i=1, env%nbSourceContainers      
        call FindContainerForAttachedObs(env%srcContArray(i),obsCont%attachedTo,&
                                       obsCont%CBList, success)
      end do
      if(.not.success) then
        call Warning('Observer could not find object to attach to,', &
                     'or there are no change of base to specified object.')
      end if
      if(associated(obsCont%CBListLocalFrame)) then
        call appendCB(obsCont%CBList, obsCont%CBListLocalFrame)
      end if
    end if
!!! 
    ! Finally, set up the time range (if neccessary) 
    if (obsCont%tMin == obsCont%tMax) then
      if (.not.SourceParallelized()) then
        call Message('Calculating observer container time ranges; this may take a few minutes.')
        call SetObserverContainerTimeRange (env, obsCont)
        call CheckSegmentTimeSettings(obsCont)  
      else
        call Error('Cannot calculate observer time range when ', &
                   'using source parallelization',               &
                   'Specify observer time range (tMin, tMax, nt)')
      end if
    end if
    if( debuglevel> 1 ) call WriteSegmentDataToScreen(obsCont) 
    !Now we can check all the frequency analysis settings because we have
    !  all the segment time lengths, number of time steps, etc
    call CheckFreqAnalysisParameters(obsCont)
          
  end subroutine FinishObserverContainerSetup
  
  subroutine SetObserverContainerTimeRange (env, obsCont) 
    type(environment),intent(inout)::env 
    type(observerContainer)::obsCont
    integer :: i,j,obsNum
    type(vector),dimension(:,:),allocatable::obscoordinates
    real, dimension (:,:,:), allocatable :: obsTimes   
    
    allocate (obsTimes(obsCont%NbObsDim1,obsCont%NbObsDim2,2))
    allocate (obsCoordinates(obsCont%NbObsDim1,obsCont%NbObsDim2))
    obsTimes(:,:,1) = -HUGE(obsTimes)/100. 
    obsTimes(:,:,2) =  HUGE(obsTimes)/100. 
    obsNum=0
    
    do i=1,obsCont%NbObsDim1
      do j=1,obsCont%NbObsDim2
        obsNum=obsNum+1
        obscoordinates(i,j)=obsCont%obs(obsNum)%coordinates
      end do
    end do
        
    do i=1,env%nbSourceContainers
      if (obsCont%dt .gt. 0) then
        if (obsCont%tMax .gt. obsCont%tMin) then
        ! This should occur if tMax (and likely tMin) are given.  If tMin is not given
        ! it is still assumed to be 0.0
          obsCont%nt = floor(obsCont%tMin-obsCont%tMax)/obsCont%dt+1
          obsCont%tMax = (obsCont%nt-1)*obsCont%dt+obsCont%tMin
        else 
        !This will occur if only tMin, or neither tMin nor tMax are given.
          obsCont%tMax = (obsCont%nt-1)*obsCont%dt+obsCont%tMin
        end if
      end if    
      
      call GetContainerImplicitTimeRange (env%srcContArray(i), & 
                                          obsCont%tMin, obsCont%tMax, & 
                                          obscoordinates, obsCont%CBList, &
                                          obsTimes)
    end do
    
    ! Set the starting time to the time at which the signal starts for 
    ! ALL observers
    obsCont%tMin = MAXVAL (obsTimes(:,:,1))
    ! And set the max to the time at which the signal isn't reaching 
    ! ALL observers
    obsCont%tMax = MINVAL (obsTimes(:,:,2))
!   
!ksb - reduce time range.  This is a workaround because I can't find why the end points are 
!coming out wrong.
!
!    obsCont%tMin = obsCont%tMin + .04*(obsCont%tMax - obsCont%tMin)
!    obsCont%tMax = obsCont%tMax - .04*(obsCont%tMax - obsCont%tMin)

    if(obsCont%tMax .lt. obsCont%tMin)then
      write(*,*)'tMin = ',obsCont%tMin
      write(*,*)'tMax = ',obsCont%tMax
      call Error('Max observer time is less than min observer time.')
    end if
    if(parallelized()) then
      call BroadcastReal(obsCont%tMax)
      call BroadcastReal(obsCont%tMin)
    end if
    deallocate (obsTimes, obsCoordinates)
  end subroutine SetObserverContainerTimeRange
      
  subroutine ReadInPressure(obsCont)
    type(observerContainer),intent(inout)::obsCont
    integer::unitNum, obsNum, i, ndat, datnum, nFreq, j
    real(kind=4)::time, freq, f, radDist
    real(kind=4),dimension(:),allocatable::pressure
    character:: temp*33, BBtitle*14, obsTimeChar*16, radDistChar*21
    unitNum=GetStreamNumber()
    ! This counts the number of items that will be in the output pressure file
    ! - not the number of items that are in the pprime array. (I.e., if total
    ! pressure is requested, thickness and loading are computed and summed to 
    ! provide the total, but only total is output.)
    ndat = 0
    if( thicknessNoiseFlag ) ndat=ndat+1
    if( loadingNoiseFlag )   ndat=ndat+1
    if( totalNoiseFlag )     ndat=ndat+1
    allocate(pressure(ndat))      
    do obsNum=1,obsCont%nbObs
      if (.not.associated(obsCont%obs(obsNum)%pPrime)) then
        allocate(obsCont%obs(obsNum)%pPrime(ndat))
          do i=1, ndat
            nullify(obsCont%obs(obsNum)%pPrime(i)%t)
            nullify(obsCont%obs(obsNum)%pPrime(i)%f)
          end do
      end if
      do i=1,ndat
        call SetupTimeHistory(obsCont%obs(obsNum)%pPrime(i), &
                              obsCont%obs(obsNum)%tMin,      &
                              obsCont%obs(obsNum)%tMax,      &
                              obsCont%obs(obsNum)%nt,        &
                              PPRIMETITLEARRAY(i))
        allocate(obsCont%obs(obsNum)%pPrime(i)%f(obsCont%obs(obsNum)%nt))
      end do
    end do
    if(obsCont%nbObs==1)then !This is tecplot format
      write(*,'(A$)')'   Reading in pressure data from pressure.in file...'
      open(unit=unitNum,file=trim(globalfoldername)//'pressure.in',form='formatted', &
           access='sequential')    
      !Read four lines of header then get to the data:
      read(unitNum,*)
      read(unitNum,*)
      read(unitNum,*)
      read(unitNum,*)
      if(ndat==3)then
        do i=1,obsCont%nt
          read(unitNum,*)time,pressure(1),pressure(2),pressure(3)
          do datnum=1,ndat
            obsCont%obs(1)%pPrime(datnum)%t(i) = time
            obsCont%obs(1)%pPrime(datnum)%f(i) = pressure(datnum)
          end do
        end do        
      else if (ndat==2)then
        do i=1,obsCont%nt
          read(unitNum,*)time,pressure(1),pressure(2)
          do datnum=1,ndat
            obsCont%obs(1)%pPrime(datnum)%t(i) = time
            obsCont%obs(1)%pPrime(datnum)%f(i) = pressure(datnum)
          end do
        end do 
      else if (ndat==1)then
        do i=1,obsCont%nt
          read(unitNum,*)time,pressure(1)
          do datnum=1,ndat
            obsCont%obs(1)%pPrime(datnum)%t(i) = time
            obsCont%obs(1)%pPrime(datnum)%f(i) = pressure(datnum)
          end do
        end do
      end if
      close(unitNum)  
      write(*,*)'file read successfully.'

      if (globalPeggNoiseFlag .or. globalBPMNoiseFlag) then

    	allocate(obsCont%obs(1)%mspBB(1,obsCont%nbSegments))
        do i=1,obsCont%nbSegments
          write(*,'(A$)')'   Reading in broadband data from broadband' //trim(integertostring(i))// '.in file...'
          open(unit=unitNum,file=trim(globalfoldername)//'broadband' //trim(integertostring(i))// '.in',form='formatted', &
               access='sequential')
          !Read four lines of header then get to the data:
          read(unitNum,*)
          read(unitNum,*)
          read(unitNum,*)
          read(unitNum,*)
          allocate(obsCont%obs(1)%mspBB(1,i)%freq(BB_THIRD_OCTAVE_BANDS), obsCont%obs(1)%mspBB(1,i)%f(BB_THIRD_OCTAVE_BANDS))
!          call SetFreqDomainTitle(obsCont%obs(1)%mspBB(1,i), trim(BBtitle))
          call SetNf(obsCont%obs(1)%mspBB(1,i), BB_THIRD_OCTAVE_BANDS)

!          read(unitNum,*) obsTimeChar, time
          do j=1,obsCont%obs(1)%mspBB(1,1)%nf
            read(unitNum,*) freq,f
!          obsCont%obs(1)%mspBB(i)%t(i) = time
            obsCont%obs(1)%mspBB(1,i)%freq(j) = freq
            obsCont%obs(1)%mspBB(1,i)%f(j) = f   !ConvertdBtoMSP(f)
            if (j.eq.1) then
              call SetMinFreq(obsCont%obs(1)%mspBB(1,i), freq)
            else if (j.eq.obsCont%obs(1)%mspBB(1,1)%nf) then
              call SetMaxFreq(obsCont%obs(1)%mspBB(1,i), freq)
            end if
          end do
          close(unitNum)
          write(*,*)'file read successfully.'
        end do
      end if

      if (AtmAtten) then
        write(*,'(A$)')'   Reading in radiation distance data from RadDistance.in file...'
        open(unit=unitNum,file=trim(globalfoldername)//'radDistance.in',form='formatted', &
             access='sequential')    

        !Read four lines of header then get to the data:
	    read(unitNum,*)
	    read(unitNum,*)
	    read(unitNum,*)
	    read(unitNum,*)

	    allocate(obsCont%obs(1)%radDistance(obsCont%nbSegments))
        do i=1,obsCont%nbSegments
	      allocate(obsCont%obs(1)%radDistance(i)%t(1))
          allocate(obsCont%obs(1)%radDistance(i)%f(1))
          obsCont%obs(1)%radDistance(obsCont%nbSegments)%nt = 1
          obsCont%obs(1)%radDistance(obsCont%nbSegments)%title = 'Time ' //integertostring(i)// ''

          read(unitNum,*) time, radDist
          obsCont%obs(1)%radDistance(i)%t(1) = time
          obsCont%obs(1)%radDistance(i)%f(1) = radDist
        end do
        close(unitNum)  
        write(*,*)'file read successfully.'
      end if
    else
      call Error('ReadInPressureFlag is only supported for 1 observer. Stopping.')
    end if 

    if((obsCont%tMin.ne.obsCont%obs(1)%pPrime(1)%t(1)).or. &
        (obsCont%tMax.ne.obsCont%obs(1)%pPrime(1)%t(obsCont%nt)))then
      call Error('tMin and/or tMax in namelist do not match input pressure.in file. Stopping.')
    end if
  end subroutine ReadInPressure
end program main 
