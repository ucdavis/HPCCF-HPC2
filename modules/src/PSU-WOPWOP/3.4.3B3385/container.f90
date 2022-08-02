module containerObject
! PSU-WOPWOP
! $Id: container.f90 3372 2017-08-17 12:52:59Z brentner $
! $LastChangedDate: 2017-08-17 08:52:59 -0400 (Thu, 17 Aug 2017) $
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
!
! Written by Guillaume Bres, Guillaume Perez, Leonard Lopes, Hsuan-Nien Chen, 
! Christopher Hennes, Rui Cheng and Benjamin Goldman. 
! Faculty advisor Dr. Kenneth S. Brentner.


!********************************************************************************
!TYPE: Container module
!   This module encapsulates the type container. It contains all container data 
!  and functions including createContainer, destroyContainer, writeContainer,
!  integrateContainer, as well as some others.
!  - This file was created by Leonard Lopes 2/2704 based on aircraft created
!   during the summer of 2002
!********************************************************************************
  use COBObject
  use mathModule
  use patchObject
  use IOModule
  use debugModule
  use MPIModule
  use interpolateModule, only: sourceTimeInRange
  use broadbandModule
  use strings
  implicit none  

  !***************************************************************************
  !TYPE: patchFile
  !  This data type stores the file information for a single patch file. A
  !  single file may contain one or more patches, and a container can have
  !  several patch files associated with it.
  !***************************************************************************
  type PatchFile
    type(Patch), dimension(:), pointer::patchArray=>null()
    character(len=32) :: pressureUnit
    logical :: parallelAllocatedFlag     
    ! Geometry file information
    character(len=4096):: patchGeometryFile
    integer:: geometryStreamNumber, startGeometryData, geomStepInBytes
    integer:: geometryMajorVerNum, geometryMinorVerNum
    integer,dimension(GEO_SIZE)::geoInfo
    logical,dimension(:),pointer:: compThickness=>null() 
    logical:: geoAperiodic, geoBigEndian
    ! Loading/Flow data Information    
    character(len=4096):: patchLoadingFile
    integer:: surfaceDataStreamNumber, startSurfaceData, loadStepInBytes
    integer:: surfaceDataMajorVerNum,surfaceDataMinorVerNum
    integer,dimension(LOAD_SIZE)::loadInfo
    logical,dimension(:),pointer::hasLoading=>null()
    logical:: loadAperiodic, LoadBigEndian
    ! Mutliple time file (MTF) aperiodic information
    character(len=4096):: MTFGeometryFile, MTFLoadingFile
    character(len=4096), dimension(:), allocatable:: MTFGeoFileArray, MTFLoadFileArray
    integer:: nbMTFGeoFiles, nbMTFLoadFiles, MTFGeoNKey, MTFLoadNKey
    integer:: iMTFGeoFile, iMTFLoadFile, MTFPrevGeoFileNKey, MTFPrevLoadFileNKey
    integer, dimension(:), allocatable:: MTFGeoStartKeyIndex, MTFGeoEndKeyIndex 
    integer, dimension(:), allocatable:: MTFLoadStartKeyIndex, MTFLoadEndKeyIndex
    real, dimension(:), allocatable:: MTFGeoRefTMin, MTFGeoRefTMax
    real, dimension(:), allocatable:: MTFLoadRefTMin, MTFLoadRefTMax
    real:: GeoRefTMin, GeoRefTMax, LoadRefTMin, LoadRefTMax
    logical:: MTFgeoAperiodic, MTFLoadAperiodic, MTFgeoQperiodic, MTFLoadQperiodic
    ! This may be included in MTF file information, but it really is a way to 
    ! subset the data for Dave Lockard (manual way to run source parallel)
    ! These will contain the first and last node to use in each patch.
    integer, dimension(:), allocatable:: PatchFirstNode, PatchLastNode
  end type PatchFile

  !**************************************************************************
  !TYPE: Container
  !  The type container contains the patchArray and containerArray for this container,
  !  the container's title, the CBList which is a change of base list from its
  !  parent to itself, the number of base changes, and the number of subcontainers
  !  the container contains.
  !***************************************************************************
  type Container
    private
    character(len=4096)::title
    type(container), dimension(:), pointer::containerArray=>null()
    type(patchFile), dimension(:), pointer:: PFArray=>null()
    type(broadband), pointer:: BBdata=>null()
    logical:: PeggNoiseFlag, BPMNoiseFlag
    integer:: PeggUnitNumber, BPMUnitNumber
    character(len=4096)::PeggNoiseFile,BPMNoiseFile
    type(CBStructure), pointer::CBlist
    integer::nbContainer, nbBase, nbFiles, nTau
    real::tauMin, tauMax, dTau
    real(kind=8)::minObsTimeInData, maxObsTimeInData
    logical:: periodicBlades ! For backwards compatibility
    logical:: parentRotation
    logical:: aperiodicData
    real :: periodicKeyOffset
    character(len=32) :: prefix
    integer :: prefixLength
  end type Container 

  ! The maximum number of patch files that can be read in by a single container
  integer, parameter :: MAX_NB_PATCH_FILES=100
  
  !****
  ! We can't be sure that the source time array will be available before we get into 
  ! the integrateContainer routine.  At that point we need to be able to keep track of a 
  ! 'global' tauMin to correctly track a flight path.
  !****
  real:: mg_tauMin = -huge(mg_tauMin)
  
  private
  public:: Container, CreateContainer, DestroyContainer
  public:: GetNbContainerSurfaces, IntegrateContainer
  public:: FindContainerForAttachedObs
  public:: GetContainerStructuredDims, GetContainerImplicitTimeRange
  public:: CreateAircraftContainer, AperiodicCont
  public:: WriteContainerStructuredSigma, WriteContainerFVUNSHeader
  public:: WriteContainerFVUNSData, WriteContainerDebugInfo, CloseDataFiles
  public:: ResetContainerKeyCount, CreateContainerObsTimeArrays
  public:: getContainerMinObsTimeInData, getContainerMaxObsTimeInData
  public:: CheckBBFlags

contains


  !**************************************************************************
  ! SUBROUTINE: CreateContainer
  ! USAGE: This is the main creation routine for container. Note that it DOES
  !   NOT handle BladeIn, RotorIn or AircraftIn containers. It reads in the
  !   container namelist from the main namelist file and based on that 
  !   information allocates the next level of the container tree.
  ! PARAMETERS:
  !   unitNumber - the open file unit of the main namelist file.
  !   C - the container being initialized by the routine
  !   dTau - source time step from parent
  !   tauMin - min source time from parent
  !   tauMax - max source time from parent
  !   nTau - number of source timesteps from parent
  !   parentCBList - parent's list of base changes
  !   prefix - a string corresponding to the current output indentation level
  recursive subroutine CreateContainer(unitNumber, C, dTau, tauMin, tauMax, ntau, &
                                       PeggUnitNumber, BPMUnitNumber, parentCBList, prefix)
    type(Container), intent(inout) :: C
    type(CBStructure), optional, pointer::parentCBlist
    integer, intent(in)::unitNumber
    character(len=*), intent(in) :: prefix
    integer:: contNum, PeggUnitNumber, BPMUnitNumber
    integer, intent(inout) :: nTau
    real, intent(inout) :: dtau, tauMin, tauMax

    integer :: i,j,nbFiles, nbContainer, nbBase, stat
    type(CBStructure) :: tempCB, previousCB
    type(PeggNamelist), pointer:: PeggNamelistData
    type(BPMNamelist), pointer:: BPMNamelistData
    type(Pegg), pointer:: PeggData=>null()
    type(BPM), pointer:: BPMData=>null()
    character(len=4096) :: title
    !character(len=1024):: PeggNoiseFile, BPMNoiseFile
    character(len=1024), dimension(MAX_NB_PATCH_FILES) :: patchGeometryFile, patchLoadingFile
    logical:: PeggNoiseFlag, BPMNoiseFlag
    logical:: isPermeable, parentRotation
    real :: periodicKeyOffset
    ! Patch file information:
    namelist / ContainerIn / nbContainer, nbBase, title, dtau, patchGeometryFile, &
                             patchLoadingFile, isPermeable, tauMin, tauMax, &
                             nTau, periodicKeyOffset, PeggNoiseFlag, BPMNoiseFlag
    ! Set up the defaults
    patchGeometryFile = ''
    patchLoadingFile  = ''
    PeggNoiseFlag     = .false.
    BPMNoiseFlag      = .false.
    nbFiles           = 0
    C%nbFiles         = 0
    nbBase            = 0
    C%nbBase          = 0
    nbContainer       = 0
    C%nbContainer     = 0
    title             = "Generic Container"
    isPermeable       = .false. ! Now ignored--left for backwards compatibility
    periodicKeyOffset = 0.0
    C%aperiodicData   = .false.
    nullify(C%PFArray)
    nullify(C%containerArray)
    
    ! Read the namelist
    read(unitNumber, nml=ContainerIn)
    C%title  = title
    C%dTau   = dTau
    C%tauMin = tauMin
    C%tauMax = tauMax
    C%nTau   = nTau
    C%prefix = prefix
    C%prefixLength = len(prefix)
    C%periodicKeyOffset = periodicKeyOffset
    C%PeggNoiseFlag = PeggNoiseFlag
    C%BPMNoiseFlag  = BPMNoiseFlag

    call Message ('Container: '//C%prefix(1:C%prefixLength)//trim(C%title))
    if (C%dtau.eq.0.0.and.C%ntau.gt.0) then
      C%dtau = (C%tauMax-C%tauMin)/C%ntau
    end if
    if (nbContainer < 0) then
      call Error('Number of source containers negative in '//trim(C%title))
    end if
    if(nbBase < 0) then
      call Error('Negative number specified for the number of base ', &
                 'changes in '//trim(C%title))
      stop
    end if
    C%nbContainer = nbContainer
    C%nbBase      = nbBase
    ! Begin by assuming that the rotation started with a CB 
    ! from this container.  If the rotation is appended then 
    ! parentRotation=.false. as below
    parentRotation = .true.
    nullify(C%CBList)
    if(present(parentCBList)) then
      if(associated(parentCBList)) then
        call appendCB(C%CBList,parentCBlist)
        if (TrueRotation(parentCBList)) parentRotation = .false.
      end if
    end if
    
!ksb16 debug: it doesn't seem like this should be here after the if block immediately above this.    
! maye it should be C%parentRotation = parentRotation
! need to check this in the future. It seems to only be related to broadbnd noise.  3/2/2014
    C%parentRotation = .false.
    if (C%nbContainer.gt.0) then 
      allocate(C%containerArray(C%nbContainer))
    else
      nullify(C%containerArray)
    end if
    do i=1, C%nbBase
       call createCBData(tempCB,unitNumber)  
       call appendCB(C%CBlist,tempCB)
       if(tempCB%windframe) call buildWindFrame(tempCB,previousCB)
       previousCB = tempCB
    end do
    ! If the rotation was not appended (hence parentRotation is still .true.)
    ! but there is NOW a rotation, then it means that the rotation DOES start
    ! with a CB in this container.
    if (associated(C%CBList)) then
      if (parentRotation .and.  TrueRotation(C%CBList)) C%parentRotation=.true.
     end if
    
    C%nbFiles = 0
    ! Count the number of patch files:
    do i=1,MAX_NB_PATCH_FILES
      if( len_trim(patchGeometryFile(i)) > 0) then
        C%nbFiles = C%nbFiles + 1
      else
        exit
      end if
    end do
    if (C%nbFiles.gt.0) then
      allocate(C%PFArray(C%nbFiles))
    else 
      nullify(C%PFArray)
    end if
    ! This loop is necessary because the user can specifiy multiple patch or loading files
    ! in the namelist (for one container).  This do loop reads each file in turn.
    do i=1,C%nbFiles
      !ksb debug: 8/4/2014 - I don't know why these names have the "./" in the middle - it does work though.
      !C%PFArray(i)%patchGeometryFile = trim(globalFolderName)//'./'//patchGeometryFile(i)
      C%PFArray(i)%patchGeometryFile = trim(globalFolderName)//patchGeometryFile(i)
      if (len_trim(patchLoadingFile(i)) > 0) then
        C%PFArray(i)%patchLoadingFile = trim(globalFolderName)//patchLoadingFile(i) 
        !C%PFArray(i)%patchLoadingFile = trim(globalFolderName)//'./'//patchLoadingFile(i)
      else
        C%PFArray(i)%patchLoadingFile = ''
      end if
      ! Initialize aperiodic logical values for new Patch File
      C%PFarray(i)%geoAperiodic = .false. 
      C%PFarray(i)%MTFgeoAperiodic=.false.
      C%PFarray(i)%loadAperiodic = .false.
      !
      call CreatePatchesFromFile (C, C%PFArray(i), C%periodicKeyOffset)
      
      if( C%PFarray(i)%geoAperiodic .or. C%PFarray(i)%loadAperiodic) C%aperiodicData=.true.
    end do
    
    if (C%PeggNoiseFlag .or. C%BPMNoiseFlag) then
      if (TrueRotation(C%CBList)) then
        allocate(C%BBdata)
        ! Currently we're limiting the prediction to one BB method per run.
        if (C%PeggNoiseFlag) then
          if (.not.globalBPMNoiseFlag) then
            globalPeggNoiseFlag = .true.
            allocate(PeggNamelistData)
            call ReadPeggNamelist(unitNumber, PeggNamelistdata, C%PeggNoiseFile)
            call readBBHeader(C%PeggNoiseFlag, C%PeggNoiseFile, C%PeggUnitNumber)
            call CreatePeggArrays(C%BBData, C%PeggUnitNumber, PeggNamelistData)        
            PeggData => GetPeggData(C%BBData)          
          else
            C%PeggNoiseFlag=.false.
            Call Notice(' Currently, only one broadband prediction method is supported',&
                        ' to run at a time.  Broadband prediction will only be performed',&
                        ' with BPM`s method.')
          end if
        end if
        if (C%BPMNoiseFlag) then
          if (.not.globalPeggNoiseFlag) then
            globalBPMNoiseFlag  = .true.
            allocate(BPMNamelistData)
            call ReadBPMNamelist(unitNumber, BPMNamelistData, C%BPMNoiseFile)
            call readBBHeader(C%BPMNoiseFlag, C%BPMNoiseFile, C%BPMUnitNumber)
            call CreateBPMArrays(C%BBdata, C%BPMUnitNumber, BPMNamelistData)
            BPMData => GetBPMData(C%BBdata)
          else
            C%BPMNoiseFlag=.false.
            Call Notice(' Currently, only one broadband prediction method is supported',&
                        ' to run at a time.  Broadband prediction will only be performed',&
                        ' with Pegg`s method.')
          end if
        end if
        if (C%PeggNoiseFlag) then
          if (getPeggTimeType(PeggData).ne.LOAD_TIMETYPE_APERIODIC) call CloseBinaryFile (C%PeggUnitNumber)
          deallocate(PeggNamelistData)
          nullify(PeggData)
        end if
        if (C%BPMNoiseFlag) then
          if (getBPMTimeType(BPMData).ne.LOAD_TIMETYPE_APERIODIC) call CloseBinaryFile (C%BPMUnitNumber)
          call DestroyBPMNamelistData(BPMNamelistData)
          deallocate(BPMNamelistData)
          nullify(BPMData)
        end if
      else
        call Notice(' Broadband noise flags must be placed in a container in which',&
                    ' rotation=.true. such that the geometry data in the',&
                    ' subsequent containers is constant, requiring a CB for rotation.',&
                    ' Broadband noise will not be computed for this container.')
        C%PeggNoiseFlag = .false.
        C%BPMNoiseFlag  = .false.
      end if
    else
      C%PeggUnitNumber = PeggUnitNumber
      C%BPMUnitNumber  = BPMUnitNumber
    end if
    
    !ksb debug call Message (C%prefix(1:C%prefixLength)//trim(C%title))
    if (C%dtau.eq.0.0.and.C%ntau.gt.0) then
      C%dtau = (C%tauMax-C%tauMin)/C%ntau
    end if
    ! If this container contains other containers, this do loop accounts
    ! for it and calls CreateContainer again as many times as there are
    ! subContainers
    do i=1, C%nbContainer
       contNum = i
       call CreateContainer(unitNumber, C%containerArray(i), &
                            dTau, tauMin, tauMax, nTau, C%PeggUnitNumber, C%BPMUnitNumber, &
                            C%CBList, prefix//"  ")
       if( C%containerArray(i)%aperiodicData .or. C%aperiodicData ) C%aperiodicData=.true.
    end do
  
  end subroutine CreateContainer


  !**************************************************************************
  ! SUBROUTINE: CreateAircraftContainer
  ! USAGE: This routine exists for backwards compatibility with the original
  !   PSU-WOPWOP namelist structure AircraftIn. An aircraft is a container that
  !   is specialized to contain rotors (and nothing else). Because it is the top
  !   level object it does not inherit any parameters from a parent container.
  ! PARAMETERS:
  !   unitNumber - the open file unit of the main namelist file.
  !   C - the container being initialized by the routine
  !   prefix - a string corresponding to the current output indentation level
  subroutine CreateAircraftContainer(unitNumber,C, prefix)
    type(container), intent(inout)::C
    integer, intent(in)::unitNumber
    character(len=*), intent(in) :: prefix

    integer::i,nbRotor,nbBase, nbContainer, nTau
    type(CBStructure)::tempCB, previousCB
    real::dTau, tauMin, tauMax
    character(len=4096)::title
    
    namelist / AircraftIn / nbRotor, nbBase, title, dTau, tauMin, tauMax, nTau
    
    ! Set up the defaults
    nbBase      = 0
    nbContainer = 0
    nbRotor     = 0
    nTau        = 0
    tauMin      = 0.0
    tauMax      = 0.0
    dTau        = 0.0
    title       = "Aircraft"

    ! Read the namelist
    read(unitNumber, nml=AircraftIn)


    ! Set the object variablesenv%nbContainer
    C%nbFiles     = 0
    C%nbContainer = nbRotor
    C%nbBase      = nbBase
    C%title       = title
    C%prefix = prefix
    C%prefixLength = len(prefix)
    C%dTau        = dTau
    C%tauMin      = tauMin
    C%tauMax      = tauMax
    C%nTau        = nTau
    C%PeggNoiseFlag = .false.
    C%BPMNoiseFlag  = .false.
    call Message (C%prefix(1:C%prefixLength)//trim(C%title))
    if (C%dtau.eq.0.0.and.C%ntau.gt.0) then
      C%dtau = (C%tauMax-C%tauMin)/C%ntau
    end if

    ! Nullify the container's CBList
    nullify(C%CBList)
    nullify(C%PFArray)
    ! Allocate the container array for the aircraft rotors
    allocate(C%containerArray(C%nbContainer))
    ! Read in the CB list
    do i=1, C%nbBase
       call createCBData (tempCB, unitNumber)  
       call appendCB     (C%CBlist, tempCB)
       if(tempCB%windframe) call buildWindFrame(tempCB,previousCB)
       previousCB = tempCB
    end do
    do i=1, C%nbContainer
      call CreateRotorContainer(unitNumber, C%containerArray(i), C%CBList,& 
                                 dTau, tauMin, tauMax, nTau, prefix//"  ")
    end do

  end subroutine CreateAircraftContainer


  !**************************************************************************
  ! SUBROUTINE: CreateRotorContainer
  ! USAGE: This routine exists for backwards compatibility with the original
  !   PSU-WOPWOP namelist structure RotorIn. A rotor is a subcontainer of an
  !   aircraft (defined by the AircraftIn namelist) and inherits its parameters.
  !   It cannot contain any patches, and must contain BladeIn subcontainers.
  ! PARAMETERS:
  !   unitNumber - the open file unit of the main namelist file.
  !   C - the container being initialized by the routine
  !   parentCBList - parent's list of base changes
  !   dTau - source time step from parent
  !   tauMin - min source time from parent
  !   tauMax - max source time from parent
  !   nTau - number of source timesteps from parent
  !   prefix - a string corresponding to the current output indentation level
  subroutine CreateRotorContainer(unitNumber, C, parentCBList, &
                                  dTau, tauMin, tauMax, nTau, prefix)
    integer, intent(in)::unitNumber
    type(container), intent(inout)::C
    type(CBStructure),pointer::parentCBList
    character(len=*), intent(in) :: prefix

    integer::i, nbBlades, nbBase, iRotor, nTau, PeggUnitNumber, BPMUnitNumber
    logical::periodicBlades, PeggNoiseFlag, BPMNoiseFlag, parentRotation
    type(CBStructure)::tempCB, previousCB
    type(PeggNamelist), pointer:: PeggNamelistData
    type(BPMNamelist), pointer:: BPMNamelistData
    type(PEGG), pointer:: PeggData=>null()
    type(BPM), pointer:: BPMData=>null()
    real, dimension(100)::bladeAzimuths ! The 100 is completely arbitrary
    character(len=4096)::title
    real::dTau, tauMin, tauMax

    namelist / RotorIn / iRotor, title, nbBlades, nbBase, periodicBlades, bladeAzimuths, &
                         dTau, tauMin, tauMax, nTau, PeggNoiseFlag, BPMNoiseFlag

    ! Initialize the data to its defaults
    title          = 'Rotor'
    iRotor         = 1
    nbBlades       = 4
    nbBase         = 0
    periodicBlades = .false.
    bladeAzimuths  = 0.
    PeggNoiseFlag  = .false.
    BPMNoiseFlag   = .false.

    ! Read the namelist
    read(unitNumber, nml=RotorIn)
    ! nullify the CB list
    nullify(C%CBList)
    nullify(C%PFArray)
    ! If the parent list exists start this CB list with that one
    if(associated(parentCBList)) call appendCB(C%CBList, parentCBlist)
    ! Set the object variables
    C%nbFiles        = 0
    C%title          = title
    C%nbContainer    = nbBlades
    C%nbBase         = nbBase
    C%prefix         = prefix
    C%prefixLength   = len(prefix)
    C%periodicBlades = periodicBlades
    C%dTau           = dTau
    C%tauMin         = tauMin
    C%tauMax         = tauMax
    C%nTau           = nTau
    C%PeggNoiseFlag  = PeggNoiseFlag
    C%BPMNoiseFlag   = BPMNoiseFlag
    C%prefix = prefix
    C%prefixLength = len(prefix)
    C%parentRotation=.false.
    call Message (C%prefix(1:C%prefixLength)//trim(C%title))
    if (C%dtau.eq.0.0.and.C%ntau.gt.0) then
      C%dtau = (C%tauMax-C%tauMin)/C%ntau
    end if
    ! Allocate the container array for the blades
    allocate(C%containerArray(C%nbContainer))
    ! Read in the container CB list
    do i=1, C%nbBase
      call CreateCBData (tempCB, unitNumber)       
      call AppendCB     (C%CBlist, tempCB)
      if(tempCB%windframe)then
        call BuildWindFrame(tempCB,previousCB)
      end if
      previousCB = tempCB
    end do

    if (C%PeggNoiseFlag .or. C%BPMNoiseFlag) then
      if (TrueRotation(C%CBList)) then
        allocate(C%BBdata)
        ! Currently we're limiting the prediction to one BB method per run.
        if (C%PeggNoiseFlag) then
          if (.not.globalBPMNoiseFlag) then
            globalPeggNoiseFlag = .true.
            allocate(PeggNamelistData)
            call ReadPeggNamelist(unitNumber, PeggNamelistdata, C%PeggNoiseFile)
            call readBBHeader(C%PeggNoiseFlag, C%PeggNoiseFile, C%PeggUnitNumber)
            call CreatePeggArrays(C%BBData, C%PeggUnitNumber, PeggNamelistData)        
            PeggData => GetPeggData(C%BBData)          
          else
            C%PeggNoiseFlag=.false.
            Call Notice(' Currently, only one broadband prediction method is supported',&
                        ' to run at a time.  Broadband prediction will only be performed',&
                        ' with BPM`s method.')
          end if
        end if
        if (C%BPMNoiseFlag) then
          if (.not.globalPeggNoiseFlag) then
            globalBPMNoiseFlag  = .true.
            allocate(BPMNamelistData)
            call ReadBPMNamelist(unitNumber, BPMNamelistData, C%BPMNoiseFile)
            call readBBHeader(C%BPMNoiseFlag, C%BPMNoiseFile, C%BPMUnitNumber)
            call CreateBPMArrays(C%BBdata, C%BPMUnitNumber, BPMNamelistData)
            BPMData => GetBPMData(C%BBdata)
          else
            C%BPMNoiseFlag=.false.
            Call Notice(' Currently, only one broadband prediction method is supported',&
                        ' to run at a time.  Broadband prediction will only be performed',&
                        ' with Pegg`s method.')
          end if
        end if
        if (C%PeggNoiseFlag) then
          if (getPeggTimeType(PeggData).ne.LOAD_TIMETYPE_APERIODIC) call CloseBinaryFile (C%PeggUnitNumber)
          deallocate(PeggNamelistData)
          nullify(PeggData)
        end if
        if (C%BPMNoiseFlag) then
          if (getBPMTimeType(BPMData).ne.LOAD_TIMETYPE_APERIODIC) call CloseBinaryFile (C%BPMUnitNumber)
          call DestroyBPMNamelistData(BPMNamelistData)
          deallocate(BPMNamelistData)
          nullify(BPMData)
        end if
      else
        call Notice(' Broadband noise flags must be placed in a container in which',&
                    ' rotation=.true. such that the geometry data in the',&
                    ' subsequent containers is constant, requiring a CB for rotation.',&
                    ' Broadband noise will not be computed for this container.')
        C%PeggNoiseFlag = .false.
        C%BPMNoiseFlag  = .false.
      end if
    !ksb debug - I think this is just copied from CreateContainer (which is recursive)
    !ksb debug:  the problem with this is the PeggUnitNumber and BPMUnitNumber are not defined here.
    !ksb debug:  which is only noticed in the debug version (release goes right on).
    !else
    !  C%PeggUnitNumber = PeggUnitNumber
    !  C%BPMUnitNumber  = BPMUnitNumber
    end if
    
    ! Check to see if the blades are periodic:
    if(periodicBlades) then
       ! Create the first blade normally
       call createBladeContainer (C%containerArray(1), unitNumber, i, C%CBList, &
                                  bladeAzimuths(1), dTau, tauMin, tauMax, nTau,&
                                  C%PeggUnitNumber, C%BPMUnitNumber, prefix//"  ",periodicBlades)
       do i=2, C%nbContainer
          ! The rest of the blades are copies of the first one
          call copyBladeContainer(C%containerArray(i), C%containerArray(1),  unitNumber, & 
                                  i, C%CBList, bladeAzimuths(i), &
                                  dTau, tauMin, tauMax, nTau, &
                                  C%PeggUnitNumber, C%BPMUnitNumber, prefix//"  ")
       end do
       if (debugLevel >=2.and.IsMaster()) write (*,*)
       close(60)
    else
       do i=1, C%nbContainer
         call createBladeContainer(C%containerArray(i), unitNumber, i, C%CBList, &
                                   bladeAzimuths(i), dTau, tauMin, tauMax, nTau, &
                                   C%PeggUnitNumber, C%BPMUnitNumber, prefix//"  ",periodicBlades)
       end do
    end if
   end subroutine createRotorContainer


  !**************************************************************************
  ! SUBROUTINE: CreateBladeContainer
  ! USAGE: This routine exists for backwards compatibility with the original
  !   PSU-WOPWOP namelist structure BladeIn. A blade is a subcontainer of an
  !   rotor (defined by the RotorIn namelist) and inherits its parameters.
  !   It cannot contain any subcontainers, and is the only object that can
  !   contain integration surfaces (i.e. patches).
  ! PARAMETERS:
  !   C - the container being initialized by the routine
  !   unitNumber - the open file unit of the main namelist file
  !   bladeNum - the index of this blade
  !   parentCBList - parent's list of base changes
  !   initialAzimuth - the offset azimuth of this blade
  !   dTau - source time step from parent
  !   tauMin - min source time from parent
  !   tauMax - max source time from parent
  !   nTau - number of source timesteps from parent
  !   prefix - a string corresponding to the current output indentation level
  !   periodicBlades - a logical indicating whether we are making periodic
  !                    blades or not
  subroutine CreateBladeContainer(C, unitNumber, bladeNum, parentCBList,&
                                  initialAzimuth,dTau, tauMin, tauMax, nTau,&
                                  PeggUnitNumber, BPMUnitNumber, prefix, periodicBlades)
    type(container), intent(inout) :: C
    integer, intent(in) :: unitNumber, bladeNum, PeggUnitNumber, BPMUnitNumber
    real, intent(in) :: initialAzimuth
    type(CBStructure),pointer::parentCBList
    character(len=*), intent(in) :: prefix

    type(CBStructure)::tempCB
    integer::nbBase, i, iBlade, nbContainer, nTau
    real :: dtau, tauMin, tauMax
    character(len=4096) :: title 
    logical :: isPermeable,periodicBlades
    character(len=4096), dimension(MAX_NB_PATCH_FILES)::patchGeometryFile, patchLoadingFile
    
    ! Standard namelist file information
    namelist /BladeIn/ title, nbBase, iBlade, dtau, tauMin, tauMax, nTau, &
                       patchGeometryFile, patchLoadingFile, isPermeable

    ! Set the defaults
    Title       = "Blade "//trim(IntegerToString(bladeNum))
    iBlade      = bladeNum
    nbContainer = 0
    nbBase      = 0
    isPermeable = .false.
    patchLoadingFile=''
    patchGeometryFile=''
    C%PeggNoiseFlag = .false.
    C%BPMNoiseFlag = .false.
    ! Read the blade namelist
    read(unitNumber, nml=BladeIn)
    ! Nullify the containers CB list
    nullify(C%CBlist)
    nullify(C%PFArray)
    nullify(C%containerArray)
    ! If the parent CB list is associated then start with that list
    if(associated(parentCBList)) call appendCB(C%CBList,parentCBlist)
    ! Apply the objects data
    C%title       = title
    C%nbBase      = nbBase
    C%nbContainer = nbContainer
    C%dTau        = dTau
    C%tauMin      = tauMin
    C%tauMax      = tauMax
    C%nTau        = nTau
    C%prefix = prefix
    C%prefixLength = len(prefix)
    C%periodicBlades = periodicBlades
    C%PeggUnitNumber = PeggUnitNumber
    C%BPMUnitNumber = BPMUnitNumber
    call Message (C%prefix(1:C%prefixLength)//trim(C%title))
    if (C%dtau.eq.0.0.and.C%ntau.gt.0) then
      C%dtau = (C%tauMax-C%tauMin)/C%ntau
    end if
    ! Set up the changes of base:
    ! TODO: Someday, automatically create the offset COB from the offset that
    ! is passed in, eliminating it from the namelist file.
    do i=1, C%nbBase
      call createCBData (tempCB,unitNumber)
      call appendCB     (C%CBlist,tempCB)
      if (tempCB%rotation) CBRotation=i      
    end do
    ! 'Open the patch files and read the headers (assume native binary format for
    ! now):    
    C%nbFiles = 0
    ! Count the number of patch files:
    do i=1,MAX_NB_PATCH_FILES
      if(len_trim(patchGeometryFile(i))>0) then
        C%nbFiles = C%nbFiles + 1
      else
        exit
      end if
    end do
    allocate(C%PFArray(C%nbFiles))
    do i=1,C%nbFiles
      C%PFArray(i)%patchGeometryFile = trim(globalFolderName)//'/'//patchGeometryFile(i)
      if (len_trim(patchLoadingFile(i)) > 0) then
        C%PFArray(i)%patchLoadingFile = trim(globalFolderName)//'/'//patchLoadingFile(i)
      else
        C%PFArray(i)%patchLoadingFile = ''
      end if
      ! Initialize aperiodic logical values for new Patch File
      C%PFarray(i)%geoAperiodic = .false.  
      C%PFarray(i)%loadAperiodic = .false.
      C%PFArray(i)%geometryStreamNumber=0
      C%PFArray(i)%surfaceDataStreamNumber=0
      call CreatePatchesFromFile (C, C%PFArray(i), initialAzimuth)
    end do
  end subroutine createBladeContainer
  
  
  !**************************************************************************
  ! SUBROUTINE: CopyBladeContainer
  ! USAGE: This routine functions only with the original PSU-WOPWOP hierarchical
  !   layout consisting of Aircraft->Rotor->Blade sequences. It makes a "copy"
  !   of a given blade. The copy does not actually copy any of the blade data,
  !   but simply makes pointers and adds offsets.
  ! PARAMETERS:
  !   newBlade - the blade being copied TO
  !   oldBlade - the blade being copied FROM
  !   unitNumber - the open file unit of the main namelist file
  !   bladeNum - the index of this blade
  !   parentCBList - parent's list of base changes
  !   initialAzimuth - the offset azimuth of this blade
  !   dTau - source time step from parent
  !   tauMin - min source time from parent
  !   tauMax - max source time from parent
  !   nTau - number of source timesteps from parent
  !   prefix - a string corresponding to the current output indentation level
  subroutine CopyBladeContainer (newBlade, oldBlade, unitNumber, bladeNum, &
                                 parentCBList, initialAzimuth, &
                                 dTau, tauMin, tauMax, nTau, &
                                 PeggUnitNumber, BPMUnitNumber, prefix)
    type(container), intent(inout)  :: newBlade
    type(container), intent(in)     :: oldBlade
    type(CBStructure),pointer       :: parentCBList
    real,intent(in)                 :: initialAzimuth
    integer, intent(in)             :: unitNumber, bladeNum, PeggUnitNumber, BPMUnitNumber
    character(len=*), intent(in) :: prefix

    type(CBStructure)::tempCB
    integer::i,nbBase, iBlade, nTau
    real::tauMin, tauMax, dTau
    character(len=4096)::title

    namelist /BladeIn/ title, iBlade, nbBase, tauMin, tauMax, dTau, nTau

    ! Read the blade data
    title       = "Blade "//trim(IntegerToString(bladeNum))
    iBlade      = bladeNum
    nbBase      = 0
    read (unitNumber, nml=BladeIn)

    nullify(newBlade%CBList)
    nullify(newBlade%containerArray)
    
    ! Put in the unique data for this blade
    newBlade%title       = title
    newBlade%nbBase      = nbBase
    newBlade%nbFiles     = oldBlade%nbFiles
    newBlade%nbContainer = 0
    newBlade%prefix = prefix
    newBlade%prefixLength = len(prefix)
    newBlade%PeggUnitNumber = PeggUnitNumber
    newBlade%BPMUnitNumber  = BPMUnitNumber
    call Message (newBlade%prefix(1:newBlade%prefixLength)//trim(newBlade%title))
    if (newBlade%dtau.eq.0.0.and.newBlade%ntau.gt.0) then
      newBlade%dtau = (newBlade%tauMax-newBlade%tauMin)/newBlade%ntau
    end if
    
    ! Copy the rest of the data over
    newBlade%nbFiles        = oldBlade%nbFiles
    newBlade%periodicBlades = oldBlade%periodicBlades
    newBlade%PeggNoiseFlag  = oldBlade%PeggNoiseFlag
    newBlade%BPMNoiseFlag   = oldBlade%BPMNoiseFlag

    ! Set up the changes of base:
    ! TODO: Someday, automatically create the offset COB from the offset that
    ! is passed in, eliminating it from the namelist file.
    if(associated(parentCBList)) call appendCB(newBlade%CBList,parentCBlist)
    do i=1, nbBase
      call createCBData (tempCB,unitNumber)
      call appendCB (newBlade%CBlist,tempCB)
      if (tempCB%rotation) then
        CBRotation=i
      end if
    end do
    ! Copy the patches
    allocate(newBlade%PFArray(newBlade%nbFiles))
    do i=1, newBlade%nbFiles
      call CopyPatchFile (newBlade%PFArray(i), oldBlade%PFArray(i), initialAzimuth)
    end do

  end subroutine copyBladeContainer


  !**************************************************************************
  ! SUBROUTINE: CopyPatchFile
  ! USAGE: This routine functions only with the original PSU-WOPWOP hierarchical
  !   layout consisting of Aircraft->Rotor->Blade sequences. It makes a "copy"
  !   of a given patch file. The copy does not actually copy any of the patch data,
  !   but simply makes pointers and adds offsets.
  ! PARAMETERS:
  !   newPatchFile - the blade being copied TO
  !   oldPatchFile - the blade being copied FROM
  !   unitNumber - the open file unit of the main namelist file
  !   bladeNum - the index of this blade
  !   parentCBList - parent's list of base changes
  !   initialAzimuth - the offset azimuth of this blade
  !   dTau - source time step from parent
  !   tauMin - min source time from parent
  !   tauMax - max source time from parent
  !   nTau - number of source timesteps from parent
  !   prefix - a string corresponding to the current output indentation level
  subroutine CopyPatchFile (newPatchFile, oldPatchFile, initialAzimuth)
    type(PatchFile), intent(inout) :: newPatchFile
    type(PatchFile), intent(in) :: oldPatchFile
    real,intent(in) :: initialAzimuth
    integer :: i

    ! Copy the patches and othe variables.  This is now laid out in the order
    ! of the data in the PatchFile type definition.
    allocate (newPatchFile%patchArray(size(oldPatchFile%patchArray)))
    do i=1,size(newPatchFile%patchArray)
      call CopyPatch (newPatchFile%patchArray(i), oldPatchFile%patchArray(i), &
                      initialAzimuth)
    end do
    newPatchFile%pressureUnit            = oldPatchFile%pressureUnit    
    newPatchFile%parallelAllocatedFlag   = oldPatchFile%parallelAllocatedFlag
    ! Geometry file information
    newPatchFile%patchGeometryFile       = oldPatchFile%patchGeometryFile
    newPatchFile%geometryStreamNumber    = oldPatchFile%geometryStreamNumber
    newPatchFile%startGeometryData       = oldPatchFile%startGeometryData
    newPatchFile%geomStepInBytes         = oldPatchFile%geomStepInBytes
    newPatchFile%geometryMajorVerNum     = oldPatchFile%geometryMajorVerNum
    newPatchFile%geometryMinorVerNum     = oldPatchFile%geometryMinorVerNum    
    newPatchFile%geoInfo                 = oldPatchFile%geoInfo
    allocate (newPatchFile%compThickness(size(oldPatchFile%compThickness)))
    newPatchFile%compThickness           = oldPatchFile%compThickness 
    newPatchFile%geoAperiodic            = oldPatchFile%geoAperiodic
    newPatchFile%geoBigEndian             = oldPatchFile%geoBigEndian
    ! Loading/Flow data Information
    newPatchFile%patchLoadingFile        = oldPatchFile%patchLoadingFile
    newPatchFile%surfaceDataStreamNumber = oldPatchFile%surfaceDataStreamNumber 
    newPatchFile%startSurfaceData        = oldPatchFile%startSurfaceData
    newPatchFile%loadStepInBytes         = oldPatchFile%loadStepInBytes
    newPatchFile%surfaceDataMajorVerNum  = oldPatchFile%surfaceDataMajorVerNum
    newPatchFile%surfaceDataMinorVerNum  = oldPatchFile%surfaceDataMinorVerNum    
    newPatchFile%loadInfo                = oldPatchFile%loadInfo    
    allocate (newPatchFile%hasLoading(size(oldPatchFile%hasLoading)))
    newPatchFile%hasLoading              = oldPatchFile%hasLoading
    newPatchFile%loadAperiodic           = oldPatchFile%loadAperiodic
    newPatchFile%loadBigEndian            = oldPatchFile%loadBigEndian    
    ! Mutliple time file (MTF) aperiodic information    
    newPatchFile%MTFGeometryFile         = oldPatchFile%MTFGeometryFile   
    newPatchFile%MTFLoadingFile          = oldPatchFile%MTFLoadingFile
    if( allocated(oldPatchFile%MTFGeoFileArray) ) then
      allocate (newPatchFile%MTFGeoFileArray(size(oldPatchFile%MTFGeoFileArray)))
      newPatchFile%MTFGeoFileArray         = oldPatchFile%MTFGeoFileArray
    end if
    if( allocated(oldPatchFile%MTFLoadFileArray) ) then
      allocate (newPatchFile%MTFLoadFileArray(size(oldPatchFile%MTFLoadFileArray)))
      newPatchFile%MTFLoadFileArray        = oldPatchFile%MTFLoadFileArray
    end if
    newPatchFile%nbMTFGeoFiles           = oldPatchFile%nbMTFGeoFiles
    newPatchFile%nbMTFLoadFiles          = oldPatchFile%nbMTFLoadFiles
    newPatchFile%MTFGeoNKey              = oldPatchFile%MTFGeoNKey 
    newPatchFile%MTFLoadNKey             = oldPatchFile%MTFLoadNKey
    if( allocated(oldPatchFile%MTFGeoStartKeyIndex) ) then
      allocate (newPatchFile%MTFGeoStartKeyIndex(size(oldPatchFile%MTFGeoStartKeyIndex)))
      newPatchFile%MTFLoadFileArray        = oldPatchFile%MTFLoadFileArray
    end if
    if( allocated(oldPatchFile%MTFGeoEndKeyIndex) ) then
      allocate (newPatchFile%MTFGeoEndKeyIndex(size(oldPatchFile%MTFGeoEndKeyIndex)))
      newPatchFile%MTFLoadFileArray        = oldPatchFile%MTFLoadFileArray
    end if
    if( allocated(oldPatchFile%MTFLoadStartKeyIndex) ) then
      allocate (newPatchFile%MTFLoadStartKeyIndex(size(oldPatchFile%MTFLoadStartKeyIndex)))
      newPatchFile%MTFLoadStartKeyIndex    = oldPatchFile%MTFLoadStartKeyIndex
    end if
    if( allocated(oldPatchFile%MTFLoadEndKeyIndex) ) then
      allocate (newPatchFile%MTFLoadEndKeyIndex(size(oldPatchFile%MTFLoadEndKeyIndex)))
      newPatchFile%MTFLoadEndKeyIndex      = oldPatchFile%MTFLoadEndKeyIndex
    end if
    if( allocated(oldPatchFile%MTFGeoRefTMin) ) then
      allocate (newPatchFile%MTFGeoRefTMin(size(oldPatchFile%MTFGeoRefTMin)))
      newPatchFile%MTFGeoRefTMin      = oldPatchFile%MTFGeoRefTMin
    end if
    if( allocated(oldPatchFile%MTFGeoRefTMax) ) then
      allocate (newPatchFile%MTFGeoRefTMax(size(oldPatchFile%MTFGeoRefTMax)))
      newPatchFile%MTFGeoRefTMax      = oldPatchFile%MTFGeoRefTMax
    end if
    if( allocated(oldPatchFile%MTFLoadRefTMin) ) then
      allocate (newPatchFile%MTFLoadRefTMin(size(oldPatchFile%MTFLoadRefTMin)))
      newPatchFile%MTFLoadRefTMin      = oldPatchFile%MTFLoadRefTMin
    end if
    if( allocated(oldPatchFile%MTFLoadRefTMax) ) then
      allocate (newPatchFile%MTFLoadRefTMax(size(oldPatchFile%MTFLoadRefTMax)))
      newPatchFile%MTFLoadRefTMax      = oldPatchFile%MTFLoadRefTMax
    end if
    newPatchFile%GeoRefTMin         = oldPatchFile%GeoRefTMin
    newPatchFile%GeoRefTMax         = oldPatchFile%GeoRefTMax
    newPatchFile%LoadRefTMin        = oldPatchFile%LoadRefTMin
    newPatchFile%LoadRefTMax        = oldPatchFile%LoadRefTMax    
    newPatchFile%MTFgeoAperiodic         = oldPatchFile%MTFgeoAperiodic
    newPatchFile%MTFLoadAperiodic        = oldPatchFile%MTFLoadAperiodic
    newPatchFile%MTFgeoQperiodic         = oldPatchFile%MTFgeoQperiodic
    newPatchFile%MTFLoadQperiodic        = oldPatchFile%MTFLoadQperiodic
    if( allocated(oldPatchFile%PatchFirstNode) ) then
      allocate (newPatchFile%PatchFirstNode(size(oldPatchFile%PatchFirstNode)))
      newPatchFile%PatchFirstNode       = oldPatchFile%PatchFirstNode
    end if
    if( allocated(oldPatchFile%PatchLastNode) ) then
      allocate (newPatchFile%PatchLastNode(size(oldPatchFile%PatchLastNode)))
      newPatchFile%PatchLastNode        = oldPatchFile%PatchLastNode
    end if
    !should be done copying the patchFile now    
    
  end subroutine CopyPatchFile

  !**************************************************************************
  ! SUBROUTINE: DestroyContainer
  ! USAGE: This routine deallocates the memory associated with this container.
  recursive subroutine destroyContainer(C)
    type(Container), intent(inout)::C
    integer::i

    ! First, recurse to subcontainers
    if (associated(C%containerArray)) then
      do i=1, C%nbContainer
        call DestroyContainer(C%containerArray(i))
      end do
      deallocate (C%containerArray)
      nullify(C%containerArray)
    end if

    ! Next, destroy the files associated with this one:
    if (associated(C%PFArray)) then
      do i=1, C%nbFiles
        call DestroyPatchFile(C%PFArray(i))
      end do
      deallocate (C%PFArray)
    end if
    nullify (C%PFArray) 

    if(associated(C%CBList)) then
      call destroyCBlist(C%CBList)
    end if
    nullify(C%CBList)

    if(associated(C%BBdata)) then
      call destroyBBdata(C%BBdata)
    end if
    nullify(C%BBdata)
  end subroutine destroyContainer
  
  
  !**************************************************************************
  ! SUBROUTINE: DestroyPatchFile
  ! USAGE: This routine deallocates the memory associated with this patch file.
  subroutine DestroyPatchFile (PF)
    type(PatchFile), intent (inout) :: PF
    integer :: i
    do i=1, size(PF%patchArray)
      call DestroyPatch(PF%patchArray(i))
    end do
    if (associated(PF%patchArray)) deallocate (PF%patchArray)
    if (associated(PF%hasLoading)) deallocate (PF%hasLoading)
    if (associated(PF%compThickness)) deallocate (PF%compThickness)
    nullify(PF%patchArray, PF%hasLoading, PF%compThickness)
    if( allocated(PF%MTFGeoStartKeyIndex) ) deallocate (PF%MTFGeoStartKeyIndex)
    if( allocated(PF%MTFGeoEndKeyIndex) ) deallocate (PF%MTFGeoEndKeyIndex)
    if( allocated(PF%MTFGeoFileArray) ) deallocate (PF%MTFGeoFileArray)
    if( allocated(PF%MTFLoadStartKeyIndex) ) deallocate (PF%MTFLoadStartKeyIndex)
    if( allocated(PF%MTFLoadEndKeyIndex) ) deallocate (PF%MTFLoadEndKeyIndex)
    if( allocated(PF%MTFLoadFileArray) ) deallocate (PF%MTFLoadFileArray)
    if( allocated(PF%PatchFirstNode) ) deallocate(PF%PatchFirstNode)
    if( allocated(PF%PatchLastNode) ) deallocate(PF%PatchLastNode)
  end subroutine DestroyPatchFile

  
  !**************************************************************************
  ! SUBROUTINE: CreatePatchesFromFile
  ! USAGE: This routine takes in a patchFile and creates the associated patches.
  subroutine CreatePatchesFromFile (C, PF, initialAzimuth)
    type(container), intent(inout) :: C
    type(patchFile), intent (inout) :: PF
    real, intent(in) :: initialAzimuth
    integer :: i, keyIndex, keyIndexMax, nvars, stat
    
    ! The series of .true. statements found in DeterminePatchFileFormat and ReadFileHeaderV1_0
    ! have been added to allow for quick modification of the existing subroutine
    ! by placing checks for 'truth' at lines that need only be accessed at the beginning
    ! of WOPWOP.  These two subroutines are accessed again in IntegrateContainer
    ! but the various print and allocate statements are unnecessary at that point.
    
    ! Determine the file format
    call DeterminePatchFileFormat(PF, .true., .true.)
    ! Call the correct header input routine based on the file format:
    select case (PF%geometryMajorVerNum)
    case (VERSION_0 ) ! File format version zero, i.e. PSU-WOPWOP 3.1.X format
      call ReadFileHeader_W31 (PF)
    case (VERSION_1 ) ! File format version 1
      call ReadFileHeaderV1_0 (PF)
    case default
      call Error ("Unrecognized file format version.")
    end select
    
    ! Now that the file header is read in the memory can be allocated and the
    ! patches created
    allocate (PF%patchArray(PF%geoInfo(GEO_NBZONES)))
    keyIndex = 0
    ! The patches are only read in at this point if we are not in
    ! source-parallel mode
    if (.not.SourceParallelized()) then
      select case (PF%geometryMajorVerNum)
        case (VERSION_0)
          if ( (PF%geoInfo(GEO_TIMETYPE) == GEO_TIMETYPE_APERIODIC).or. &
          &    (PF%loadInfo(LOAD_TIMETYPE) == LOAD_TIMETYPE_APERIODIC) )then
            
            call Error("This version of PSU-WOPWOP only supports aperiodic cases", &
                       "in the Version 1 format.")
          end if
                      
          do i=1, PF%geoInfo(GEO_NBZONES)
            if (debugLevel >=4) then
              call Message ("Creating patch "//trim(IntegerToString(i))//" of "//&
                            trim(IntegerToString(size(PF%patchArray))))
            end if
            call CreatePatch (PF%patchArray(i), &
              PF%geometryStreamNumber, VERSION_0, VERSION_0, PF%geoInfo, &
              PF%surfaceDataStreamNumber, VERSION_0, VERSION_0, PF%loadInfo, &
              PF%hasLoading(i), PF%compThickness(i), C%prefix(1:C%prefixLength)//"  ") 
            
    !ksb debug: 2/14/2015 - but it may have been wrong for a while.
    !print*,'Not sure this section of the code is still working - need Version 0 files to check it.'
    ! 3/15/2015 - seems like it might be still working - conclusive testing would be good.
          !  if (i .eq. 1) then
          !    allocate(endofdata(PF%geoInfo(GEO_NBZONES)))
          !    allocate(endofLoadData(PF%geoInfo(GEO_NBZONES)))
          !    endofdata = .false.
          !    endofLoaddata = .false.
          !    ! If the surface or loading data are aperiodic then we don't read everything in at this point.
          !    ! This same basic statement is used within IntegrateContainer for reading aperiodic files.
          !  end if
          !  !ksb debug: MTF testing. if (.not. endofdata(i) .or. .not. endofLoadData(i)) then
          !  if( PatchITauInRange(PF%patchArray(i)) ) then ! true if patch%itau <= patch%nTau   
          !ksb debug: If I remember correctly, in Version 0 file format, we read all the data in at once.
          !If that is true, this may still work.
              call ReadPatchAtTime(PF%patchArray(i),PF%geometryStreamNumber,&
                                   PF%surfaceDataStreamNumber,PF%geoInfo,PF%loadInfo,&
                                   initialAzimuth, PF%hasLoading(i), PF%geometryMajorVerNum)
                            
          !  end if 
          !  !ksb debug. MTF testing. 
          !  if (all(endofdata) .and. all(endofLoadData)) then  !ksb debug: MTF attention needed
          !  !if( PatchLastITau(PF%patchArray(i)) ) then  ! patch%iTau == patch%nTau
          !    exit
          !  end if
          end do                                 

                   
        case default !(VERSION_1)
          do i=1, PF%geoInfo(GEO_NBZONES)
            if (debugLevel >=4) then
              call Message ("Creating patch "//trim(IntegerToString(i))//" of "//&
                            trim(IntegerToString(size(PF%patchArray))))
            end if
            call CreatePatch (PF%patchArray(i), &
              PF%geometryStreamNumber, PF%geometryMajorVerNum, PF%geometryMinorVerNum, PF%geoInfo, &
              PF%surfaceDataStreamNumber, PF%SurfaceDataMajorVerNum, PF%SurfaceDataMinorVerNum, PF%loadInfo, &
              PF%hasLoading(i), PF%compThickness(i), C%prefix(1:C%prefixLength)//"  ") 
!ksb debug:
              PF%MTFGeoNKey=GetPatchGeoNKey(PF%patchArray(i))
              if( PF%hasLoading(i) ) then
                PF%MTFLoadNKey=GetPatchLoadNKey(PF%patchArray(i))
              else
                PF%MTFLoadNKey=0
              end if
              
            if( allocated(PF%PatchFirstNode) ) then
              call setPatchNodeRange(PF%patchArray(i),PF%PatchFirstNode(i),PF%PatchLastNode(i))
            end if
          end do
          ! If the data is aperiodic, save the location for the start of the data - 
          ! after the various header data is read.
          if(PF%geoAperiodic) then
            PF%startGeometryData = GetBinaryFilePosition(PF%geometryStreamNumber)
            PF%geomStepInBytes = 0
            do i=1,PF%geoInfo(GEO_NBZONES)
              PF%geomStepInBytes = PF%geomStepInBytes + 4*(1+6*GetPatchNbNodes(PF%patchArray(i)))
            end do
          end if
          if(PF%loadAperiodic) then
            PF%startSurfaceData = GetBinaryFilePosition(PF%surfaceDataStreamNumber)
            PF%loadStepInBytes = 0
            select case (PF%loadInfo(LOAD_DATATYPE))
              case (LOAD_DATATYPE_PRESSURE)
                nvars = 1
              case (LOAD_DATATYPE_LOADING) 
                nvars = 3
              case (LOAD_DATATYPE_FLOW) 
                nvars = 5
              end select
            PF%loadStepInBytes = 0
            do i=1,PF%loadInfo(GEO_NBZONES) 
              if( PF%hasloading(i) )then
                PF%loadStepInBytes = PF%loadStepInBytes + &
                  4*(1+nvars*GetPatchNbNodes(PF%patchArray(i)))
              end if
            end do
          end if

          do
            keyIndex = keyIndex+1
            do i=1,PF%geoInfo(GEO_NBZONES)
              if( .not.PF%geoAperiodic ) then
                if( (PF%geoinfo(GEO_TIMETYPE) == GEO_TIMETYPE_QPERIODIC) .and. keyIndex==1 ) then
                  call ReadBinaryReal(PF%GeometryStreamNumber, PF%GeoRefTMin, stat)
                  call ReadBinaryReal(PF%GeometryStreamNumber, PF%GeoRefTMax, stat)
                end if
                call ReadPatchGeomAtTime(PF%patchArray(i),PF%geometryStreamNumber, &
                  PF%geoInfo, InitialAzimuth, PF%geometryMajorVerNum)              
              end if
              if( .not.PF%loadAperiodic .and. PF%hasLoading(i) ) then
                if( (PF%loadinfo(LOAD_TIMETYPE) == LOAD_TIMETYPE_QPERIODIC) .and. keyIndex==1 ) then
                  call ReadBinaryReal(PF%SurfaceDataStreamNumber, PF%LoadRefTMin, stat)
                  call ReadBinaryReal(PF%SurfaceDataStreamNumber, PF%LoadRefTMax, stat)
                end if
                call ReadPatchLoadAtTime(PF%patchArray(i), PF%surfaceDataStreamNumber, &
                  PF%loadInfo, InitialAzimuth, PF%surfaceDataMajorVerNum)                
              end if
              keyIndexMax = max( keyIndex, GetPatchNTau(PF%patchArray(i)) ) 
            end do

            if (keyIndex.eq.2) then
              do i=1,size(PF%patchArray)
                call ResetPatchKeyCount(PF%patchArray(i))
              end do
            end if  
            if( PF%GeoAperiodic .or. PF%LoadAperiodic .or. (keyIndex == keyIndexMax) ) then
              exit
            end if
          end do
      end select      

      if (PF%geoInfo(GEO_TIMETYPE) /= GEO_TIMETYPE_APERIODIC) then
        if (any(PF%hasLoading)) then
          if (PF%loadInfo(LOAD_TIMETYPE) /= LOAD_TIMETYPE_APERIODIC) then
            !
            ! Calculate the source time step
            !
            do i=1,PF%geoInfo(GEO_NBZONES)
              call DeterminePatchDtau(PF%patchArray(i),C%tauMin, C%tauMax,&
                                   C%nTau, C%dTau)
            end do
          end if
        end if
      end if 
      ! Close the files:
!ksb debug:  need to change this so quasiperiodic files remain open too. 5/5/2016
      ! If either of the files is aperiodic, then the geometry file needs to be left open.
      !if (.not.(PF%geoAperiodic .or. PF%loadAperiodic)) then
      !  call CloseBinaryFile (PF%geometryStreamNumber)
      !end if
      select case( PF%geoInfo(GEO_TIMETYPE) )
      case (GEO_TIMETYPE_APERIODIC, GEO_TIMETYPE_QPERIODIC)
        ! Don't close the file for these types because we are not done reading from it
      case default
        if (.not.(PF%loadAperiodic)) then  !I wish this wasn't necessary keeping geometry openfor loading - 5/15/16 KSB
          !For every other time type go ahead and close the file.
          call CloseBinaryFile (PF%geometryStreamNumber)
        end if
      end select
      ! If the loading is aperiodic or qperiodic, then the loading file needs to be left open.
      select case( PF%loadInfo(LOAD_TIMETYPE) )
      case (LOAD_TIMETYPE_APERIODIC, LOAD_TIMETYPE_QPERIODIC)
      ! Don't close the file for these types because we are not done reading from it
      case default
      !For every other time type go ahead and close the file.
        call CloseBinaryFile (PF%surfaceDataStreamNumber)      
      end select
      ! if (.not.PF%loadAperiodic .and. PF%surfaceDataStreamNumber/=0 ) then
        !call CloseBinaryFile (PF%surfaceDataStreamNumber)              
      ! end if
    end if    

  end subroutine CreatePatchesFromFile

  subroutine DeterminePatchFileFormat(PF, readLoading, writeInfo) 
    type(patchFile), intent(inout)::PF
    integer::stat, startpos
    integer :: magicGeoNum, magicLoadNum
    logical:: readLoading, writeInfo, firstRead
  
    if( PF%geometryStreamNumber == 0 ) then
      PF%geometryStreamNumber     = GetStreamNumber()
        
      !Determine File Format Version of File
      call OpenBinaryFile(PF%geometryStreamNumber, PF%patchGeometryFile, .false., stat)
      if (stat/=0) then
        call Error('subroutine DeterminePatchFileFormat', &
                   'Could not open the file '//trim(PF%patchGeometryFile), &
                   'stat came back as '//IntegerToString(stat))
        stop
      end if
      call ReadBinaryInteger(PF%geometryStreamNumber,magicGeoNum,stat)
      !determine endian correctness
      if (writeInfo .and. debuglevel>10) then  
        call Message('Geometry file magic number is read in as ' &
                     //IntegerToString(magicGeoNum))
      end if
      PF%GeoBigEndian = .false. ! native or little endian
      if (magicGeoNum == 704643072) then
        !file needs endian converted
        call CloseBinaryFile(PF%geometryStreamNumber)
        ! open file as BigEndian
        call OpenBinaryFile (PF%geometryStreamNumber, PF%patchGeometryFile, .true., stat)
        call ReadBinaryInteger(PF%geometryStreamNumber,magicGeoNum,stat)
        PF%GeoBigEndian = .true.  ! big endian
      end if

      if (magicGeoNum == 42) then
        ! Now we must have V1.X file, opened with proper endiannes, so now we can 
        ! read the file version numbers.
        ! The next two numbers (after the magic number) are the version numbers
        call ReadBinaryInteger(PF%geometryStreamNumber, PF%geometryMajorVerNum, stat)
        call ReadBinaryInteger(PF%geometryStreamNumber, PF%geometryMinorVerNum, stat)      
      else 
        ! This is not a V1.0 (V1.x) file because it does not have magic number
        ! Assume it is an old W31 file.  If it isn't, ReadFileHeader_W31 will check 
        ! and figure that out, report an error, and stop.
        call CloseBinaryFile(PF%geometryStreamNumber)
        PF%geometryMajorVerNum = 0
        PF%geometryMinorVerNum = 0
      end if
    else
      ! The file is already open and has the correct endianness, Now only need to
      ! position the file to the same point as if we had just read the magic number
      ! and the major and minor version numbers
      startpos = 12
      call RestartBinaryFileAtByte(PF%geometryStreamNumber, startpos, stat)
    end if
    ! Now check out the functional data (aka SurfaceData, aka loading) file
    if (readLoading .and. len_trim(PF%patchLoadingFile) > 0) then
      if( PF%surfaceDataStreamNumber == 0 ) then
        PF%surfaceDataStreamNumber  = GetStreamNumber()
        call OpenBinaryFile (PF%surfaceDataStreamNumber, PF%patchLoadingFile, .false., stat)
        if (stat/=0) then
          call Error('subroutine DeterminePatchFileFormt', &
                     'Could not open the file '//trim(PF%patchLoadingFile), &
                     'stat came back as '//IntegerToString(stat))
          stop
        end if
        call ReadBinaryInteger(PF%surfaceDataStreamNumber,magicLoadNum,stat)
        PF%LoadBigEndian = .false. ! native or little endian
        if (magicLoadNum == 704643072) then
          !file needs endian converted
          call CloseBinaryFile(PF%surfaceDataStreamNumber)
          call OpenBinaryFile (PF%surfaceDataStreamNumber, PF%patchLoadingFile, .true.,     stat)
          call ReadBinaryInteger(PF%surfaceDataStreamNumber,magicLoadNum,stat)
          PF%LoadBigEndian = .true. ! native or little endian
        end if
        if (magicLoadNum == 42) then
          ! Now we must have V1.X file, opened with proper endiannes, so now we can 
          ! read the file version numbers.
          ! The next two numbers (after the magic number) are the version numbers
          call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%surfaceDataMajorVerNum, stat)
          call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%surfaceDataMinorVerNum, stat)
        else
          ! This is not a V1.0 (V1.x) file becuase it does not have magic number
          ! Assume it is an old W31 file.  If it isn't, ReadFileHeader_W31 will check 
          ! and figure that out, report an error, and stop.
          call CloseBinaryFile(PF%surfaceDataStreamNumber)
          PF%surfaceDataMajorVerNum = VERSION_0
          PF%surfaceDataMinorVerNum = VERSION_0
        end if 
      else
        ! The file is already open and has the correct endianness, Now only need to
        ! position the file to the same point as if we had just read the magic number
        ! and the major and minor version numbers        
        startpos = 12
        call RestartBinaryFileAtByte(PF%surfaceDataStreamNumber, startpos, stat)
      end if
    else
      ! No loading to read or not needed
      magicLoadNum = 42
      PF%surfaceDataStreamNumber = 0
    end if
    
  end subroutine DeterminePatchFileFormat

  subroutine ReadFileHeader_W31 (PF)
    type(patchFile), intent (inout) :: PF
    character(len=32) :: geometryFileFormat, surfaceDataFileFormat
    integer :: stat, patchNumberOfZones
    
    ! We already have stream numbers: re-open the file(s)
    call OpenBinaryFile (PF%geometryStreamNumber, PF%patchGeometryFile, .false., stat)
    if (stat/=0) then
      call Error('Could not open the file '//trim(PF%patchGeometryFile))
      stop
    end if
    call ReadBinaryString(PF%geometryStreamNumber, geometryFileFormat, stat)
    call strToUpper (geometryFileFormat, 32)
    if (index(geometryFileFormat, 'PATCH') == 0) then
      call Error('The file '//trim(PF%patchGeometryFile), ' is not a patch file.')
      stop
    end if        
    call ReadBinaryInteger(PF%geometryStreamNumber, patchNumberOfZones, stat)
    ! allocate  has loading and set to false to avoid seg faults.
    allocate(PF%hasloading(patchnumberOfZones),PF%compThickness(patchNumberOfZones))
    PF%hasloading = .false.
    PF%compThickness = .true.
    ! Set the appropriate geometry data for a version zero file:
    PF%geoInfo(GEO_NBZONES) = patchNumberOfZones
    PF%geoInfo(GEO_CHECK) = 1
    PF%geoInfo(GEO_GRIDTYPE) = GEO_GRIDTYPE_STRUCTURED
    PF%geoInfo(GEO_TIMETYPE) = GEO_TIMETYPE_CONSTANT
    PF%geoInfo(GEO_NORMTYPE) = GEO_NORMTYPE_NODES
    PF%geoInfo(GEO_PRECISION) = GEO_PRECISION_SINGLE

    if (len_trim(PF%patchLoadingFile) > 0) then
      call OpenBinaryFile (PF%surfaceDataStreamNumber, PF%patchLoadingFile, .false., stat)
      if (stat/=0) then
        call Error('Could not open the file '//trim(PF%patchLoadingFile))
        stop
      end if
      PF%hasloading = .true.
      call ReadBinaryString(PF%surfaceDataStreamNumber, surfaceDataFileFormat, stat)
      call strToUpper (surfaceDataFileFormat, 32)
      if (index(surfaceDataFileFormat, 'LOAD') /= 0) then
        PF%loadInfo(LOAD_DATATYPE) = LOAD_DATATYPE_LOADING
        PF%loadInfo(LOAD_REFFRAME) = LOAD_REFFRAME_BLADE  ! BAG mod
      else if (index(surfaceDataFileFormat, 'FLOW') /= 0) then
        PF%loadInfo(LOAD_DATATYPE) = LOAD_DATATYPE_FLOW
        PF%loadInfo(LOAD_REFFRAME) = LOAD_REFFRAME_MIXED  ! BAG mod
      else
        call Error(trim(surfaceDataFileFormat)//' is not a valid flow data type.', &
                   'File: '//trim(PF%patchLoadingFile))
      end if
      ! Read in the time format and save it off:
      call ReadBinaryString (PF%surfaceDataStreamNumber, surfaceDataFileFormat, stat)
      call strToUpper(surfaceDataFileFormat,32)
      if (index(surfaceDataFileFormat, 'CONSTANT') /= 0)  then
        PF%loadInfo(LOAD_TIMETYPE) = LOAD_TIMETYPE_CONSTANT
      else if (index(surfaceDataFileFormat, 'APERIODIC') /= 0)  then
        PF%loadInfo(LOAD_TIMETYPE) = LOAD_TIMETYPE_APERIODIC 
      else if (index(surfaceDataFileFormat, 'PERIODIC') /= 0)  then
        PF%loadInfo(LOAD_TIMETYPE) = LOAD_TIMETYPE_PERIODIC
      else
        call Error (trim(surfaceDataFileFormat)//' is not a valid time type.',&
                   'File: '//trim(PF%patchLoadingFile))
        stop
      end if

      call ReadBinaryInteger(PF%surfaceDataStreamNumber, patchNumberOfZones, stat)
      if (stat/=0) then
        call Error('An error occured while reading ', trim(PF%patchLoadingFile))
        stop
      end if
      PF%loadInfo(LOAD_CHECK) = 2
      PF%loadInfo(LOAD_NBZONES) = patchNumberOfZones
      PF%loadInfo(LOAD_GRIDTYPE) = LOAD_GRIDTYPE_STRUCTURED
      PF%loadInfo(LOAD_NORMTYPE) = LOAD_NORMTYPE_NODES
      PF%loadInfo(LOAD_PRECISION) = LOAD_PRECISION_SINGLE
!      if (index(surfaceDataFileFormat, 'LOAD') > 0) then
!        PF%loadInfo(LOAD_DATATYPE) = LOAD_DATATYPE_LOADING
!        PF%loadInfo(LOAD_REFFRAME) = LOAD_REFFRAME_BLADE
!      else
!        PF%loadInfo(LOAD_DATATYPE) = LOAD_DATATYPE_FLOW
!        PF%loadInfo(LOAD_REFFRAME) = LOAD_REFFRAME_MIXED
!      end if
      !PF%loadInfo(LOAD_TIMETYPE) can't be set yet...
    else
      PF%surfaceDataStreamNumber = 0
    end if
    if (patchNumberOfZones > 1000) then 
      call Warning('This seems like a lot of patches ('//IntegerToString(patchNumberOfZones)//')',&
                   'Consider checking your patch file format.')
    end if
  end subroutine ReadFileHeader_W31


  subroutine ReadFileHeaderV1_0 (PF)
 !   use constantsModule, only:globalFolderName
    type(patchFile), intent (inout) :: PF
    integer :: i, nbLoadedZones, patchNum, nfiles, stat, startpos, ipatch, nbPatches
    character (len=1024) :: comments
    character (len=4096) :: filename
    integer, parameter:: mxFiles=1000
    
    if (PF%geometryStreamNumber /=0) then
      call ReadBinaryString(PF%geometryStreamNumber,  PF%pressureUnit, stat)
      call ReadBinaryString(PF%geometryStreamNumber,  comments, stat)
      call ReadBinaryInteger(PF%geometryStreamNumber, PF%geoInfo(GEO_CHECK), stat)
      call ReadBinaryInteger(PF%geometryStreamNumber, PF%geoInfo(GEO_NBZONES), stat)
      call ReadBinaryInteger(PF%geometryStreamNumber, PF%geoInfo(GEO_GRIDTYPE), stat) 
      call ReadBinaryInteger(PF%geometryStreamNumber, PF%geoInfo(GEO_TIMETYPE), stat)
      call ReadBinaryInteger(PF%geometryStreamNumber, PF%geoInfo(GEO_NORMTYPE), stat)
      call ReadBinaryInteger(PF%geometryStreamNumber, PF%geoInfo(GEO_PRECISION), stat)
      call ReadBinaryInteger(PF%geometryStreamNumber, PF%geoInfo(GEO_IBLANK), stat) 
      call ReadBinaryInteger(PF%geometryStreamNumber, PF%geoInfo(GEO_FUTURE2), stat)

      if (abs(PF%geoInfo(GEO_CHECK)) /= 1 ) then
        call Error('The file '//trim(PF%patchGeometryFile), 'is not a patch file.',&
                   "Input was "//IntegerToString(PF%geoInfo(GEO_CHECK)))
        stop
      end if   
      ! Record the position in the file - it is where to start the first file
      PF%startGeometryData =GetBinaryFilePosition(PF%GeometryStreamNumber)
      nbPatches = PF%geoInfo(GEO_NBZONES)
      if (PF%geoInfo(GEO_CHECK) == -1 ) then
        allocate( PF%PatchFirstNode(nbPatches), PF%PatchLastNode(nbPatches) )
        do iPatch=1, nbPatches
          ! Get the first and last node to use for each patch (zone)
          call ReadBinaryInteger(PF%GeometryStreamNumber, PF%PatchFirstNode(iPatch), stat)
          call ReadBinaryInteger(PF%GeometryStreamNumber, PF%PatchLastNode(iPatch),  stat)
        end do
        if (PF%geoInfo(GEO_TIMETYPE) /= GEO_TIMETYPE_MTFAPERIODIC) then
          ! Get the real patch file name
          call ReadBinaryString (PF%GeometryStreamNumber, filename, stat)
          filename = adjustl(filename)
          if( filename(1:1)=='/' .or. filename(1:2)=='~/' ) then
            ! filename is an absolution path (or at test not relative to the globalFolderName)
            PF%patchGeometryFile = trim(filename)
          else
            ! fliename must be a relative patch - assume it is relative to globalFolderName
            PF%patchGeometryFile = trim(globalFolderName)//trim(filename)
          end if
          ! now close this geometry file and open the real one
          call CloseBinaryFile(PF%GeometryStreamNumber)
          call OpenBinaryFile (PF%GeometryStreamNumber, PF%patchGeometryFile, PF%GeoBigEndian, stat)
          startpos = PF%startGeometryData
          call RestartBinaryFileAtByte(PF%GeometryStreamNumber, startpos,stat)
          ! set the geometry check to 1 like it would be for this file
          PF%geoInfo(GEO_CHECK) = 1
        end if
      end if
    
      PF%geoAperiodic = .false.
      PF%MTFgeoAperiodic = .false.
      PF%MTFgeoQperiodic = .false.
      PF%geoRefTMin = -huge(PF%geoRefTMin)
      PF%geoRefTMax = huge(PF%geoRefTMax)
      PF%loadRefTMin = -huge(PF%loadRefTMin)
      PF%loadRefTMax = huge(PF%loadRefTMax)
      
      select case(PF%geoInfo(GEO_TIMETYPE))
      case(GEO_TIMETYPE_APERIODIC)
        PF%geoAperiodic = .true.
      case(GEO_TIMETYPE_MTFAPERIODIC)
        ! This is a Multile Time File Aperiodic case (MTFAperiodic), so now
        ! read In the list of Aperiodic Time Files (must be in sequential order)
        PF%geoAperiodic = .true.
        PF%MTFgeoAperiodic = .true. 
        ! Read how many files there are
        call ReadBinaryInteger(PF%geometryStreamNumber,nFiles,stat)
        PF%nbMTFGeoFiles = nFiles
        ! Open the arrays to hold the file names, startKeyIndex, and endKeyIndex for
        ! each of the nFiles files.
        allocate( PF%MTFGeoFileArray(nfiles) )
        allocate( PF%MTFGeostartKeyIndex(nfiles), PF%MTFGeoendKeyIndex(nfiles) )
        PF%MTFGeoNKey = 0   !MTFGeoNkey is the total number of keys for all the nFiles files.
        do i=1,nFiles
          ! read list of file names, starting keyIndex, ending keyIndex
          ! Normally the starting keyIndex will equal 1 and ending keyIndex
          ! will equal nKey for the individual file.  This enables the user
          ! to not read all the data in the files - in cases there is overlap
          ! the file writer is responsible to ensure there are no "gaps"
          ! in the data (whether it is time, or some other keyvalue).
          call ReadBinaryString (PF%GeometryStreamNumber, filename, stat)
          filename = adjustl(filename)
          if( filename(1:1)=='/' .or. filename(1:2)=='~/' ) then
            ! filename is an absolution path (or at test not relative to the globalFolderName)
            PF%MTFGeoFileArray(i) = trim(filename)
          else
            ! fliename must be a relative patch - assume it is relative to globalFolderName
            PF%MTFGeoFileArray(i) = trim(globalFolderName)//trim(filename)
          end if
          call ReadBinaryInteger(PF%GeometryStreamNumber, PF%MTFGeoStartKeyIndex(i), stat)
          call ReadBinaryInteger(PF%GeometryStreamNumber, PF%MTFGeoEndKeyIndex(i), stat)
          PF%MTFGeoNKey = PF%MTFGeoNKey + (PF%MTFGeoEndKeyIndex(i)-PF%MTFGeoStartKeyIndex(i)+1)
        end do
        ! Now close the MTF Aperiodic Geometry file and open the first file in the
        ! sequence.  Need to open with the correct endianness and skip down to the
        ! place we would be if this were just an aperiodic loading file.
        PF%MTFGeometryFile = PF%patchGeometryFile
        PF%patchGeometryFile=PF%MTFGeoFileArray(1)
        PF%iMTFGeoFile = 1
        PF%MTFPrevGeoFileNKey = 0
        call CloseBinaryFile(PF%GeometryStreamNumber)
        call OpenBinaryFile (PF%GeometryStreamNumber, PF%patchGeometryFile, PF%GeoBigEndian, stat)
        startpos = PF%startGeometryData + &
                       (PF%MTFGeoStartKeyIndex(PF%iMTFGeoFile)-1)*PF%geomStepInBytes
        call RestartBinaryFileAtByte(PF%GeometryStreamNumber, startpos,stat)
        ! Reset the timetype as if this we read the header of this file.
        PF%geoInfo(GEO_TIMETYPE)=GEO_TIMETYPE_APERIODIC
      case(GEO_TIMETYPE_QPERIODIC)
    
      case(GEO_TIMETYPE_MTFQPERIODIC)     
        ! This is a Multile File Qperiodic case (MTFQperiodic), so there will be a list
        ! of Periodic Time Files (must be in sequential order).  We won't know how many
        ! files there are until we reach the end of the MTFQperiodic file, so we have to 
        ! count them.  Assume a big number to start with.
        !ksb debug: PF%nbMTFGeoFiles = mxFiles ! maximum number of Periodic Files set to 1000.  If this is not 
                                ! enough, then we will increate the array size by another 1000.
        ! Open the arrays to hold the file names, startKeyIndex, and endKeyIndex for
        ! each of the nFiles files.
        allocate( PF%MTFGeoFileArray(mxFiles) )
        allocate( PF%MTFGeoRefTMin(mxFiles), PF%MTFGeoRefTMax(mxFiles) )

        do i=1,mxFiles
          ! read list of file names, starting RefTMin, ending RefTMax. The reference times 
          ! (RefTMin and RefTMax) indicate the range of times over which the periodic data is
          ! to be used.  
          call ReadBinaryString (PF%GeometryStreamNumber, filename, stat)
          if (stat < 0 ) exit  !end of file condition
          filename = adjustl(filename)
          if( filename(1:1)=='/' .or. filename(1:2)=='~/' ) then
            ! filename is an absolute path (or at least not relative to the globalFolderName)
            PF%MTFGeoFileArray(i) = trim(filename)
          else
            ! fliename must be a relative patch - assume it is relative to globalFolderName
            PF%MTFGeoFileArray(i) = trim(globalFolderName)//trim(filename)
          end if
          call ReadBinaryReal(PF%GeometryStreamNumber, PF%MTFGeoRefTMin(i), stat)
          if (stat < 0 ) exit  !end of file condition
          call ReadBinaryReal(PF%GeometryStreamNumber, PF%MTFGeoRefTMax(i), stat)
          if (stat < 0 ) exit  !end of file condition
          nfiles = i
          if( nfiles+1 >mxFiles ) then
            call Error('Number of patch files greater than the dimension.  Need to change mxFiles and recompile.')
            stop
          end if
        end do
        ! Now close the MTF Aperiodic Geometry file and open the first file in the
        ! sequence.  Need to open with the correct endianness and skip down to the
        ! place we would be if this were just an periodic loading file.
        PF%nbMTFGeoFiles = nfiles
        PF%MTFGeometryFile = PF%patchGeometryFile
        PF%patchGeometryFile=PF%MTFGeoFileArray(1)
        PF%iMTFGeoFile = 1
        PF%GeoRefTMin = PF%MTFGeoRefTMin(1)
        PF%GeoRefTmax = PF%MTFGeoRefTMax(1)
        call CloseBinaryFile(PF%GeometryStreamNumber)
        call OpenBinaryFile (PF%GeometryStreamNumber, PF%patchGeometryFile, PF%GeoBigEndian, stat)
        startpos = PF%startGeometryData 
        call RestartBinaryFileAtByte(PF%GeometryStreamNumber, startpos,stat)
        !Need to read the Structured or Unstructured Header
        ! Reset the timetype as if this we read the header of this file.
        PF%geoInfo(GEO_TIMETYPE)=GEO_TIMETYPE_PERIODIC
        PF%MTFgeoQperiodic = .true.
      case default  ! constant, periodic, or qperiodic geometry      

      end select
    end if          
    ! initialize arrays needed for function data (loading).  These need to be set 
    ! even if there is no function data (so that we know it).
    allocate(PF%hasLoading(PF%geoInfo(GEO_NBZONES)),PF%compThickness(PF%geoInfo(GEO_NBZONES)))
    ! default values (what is required if there is no function data
    PF%hasLoading=.false.
    PF%compThickness=.true.      
    if (PF%surfaceDataStreamNumber /=0 .and. len_trim(PF%patchLoadingFile) > 0) then
      call ReadBinaryString (PF%surfaceDataStreamNumber, comments,stat)
      call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%loadInfo(LOAD_CHECK), stat)
      call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%loadInfo(LOAD_NBZONES), stat)   
      call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%loadInfo(LOAD_GRIDTYPE), stat)
      call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%loadInfo(LOAD_TIMETYPE), stat)
      call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%loadInfo(LOAD_NORMTYPE), stat)
      call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%loadInfo(LOAD_DATATYPE), stat)
      call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%loadInfo(LOAD_REFFRAME), stat)
      call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%loadInfo(LOAD_PRECISION), stat)
      call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%loadInfo(LOAD_FUTURE1), stat)
      call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%loadInfo(LOAD_FUTURE2), stat)
      PF%startSurfaceData =GetBinaryFilePosition(PF%surfaceDataStreamNumber)
      
      if (PF%loadInfo(LOAD_CHECK) /= 2 ) then
         call Error('The file '//trim(PF%patchLoadingFile), ' is not a loading or flow data file.')
        stop
      end if
      if (PF%geoInfo(GEO_NBZONES) /= PF%loadInfo(LOAD_NBZONES)) then
        call Error('Number of zones in patch must equal number in loading or flow data.', &
                   'Geometry patches: '//&
                   trim(IntegerToString(PF%geoInfo(GEO_NBZONES)))//&
                   ', loading patches: '//&
                   trim(IntegerToString(PF%loadInfo(LOAD_NBZONES))))
        stop
      end if
      if (PF%geoInfo(GEO_NBZONES) > 1000) then 
        call Warning('This seems like a lot of patches ('//IntegerToString(PF%geoInfo(GEO_NBZONES))//')',&
                     'Consider checking your patch file format.')
      end if

      PF%loadAperiodic = .false.
      PF%MTFLoadAperiodic = .false.
      PF%MTFloadQperiodic = .false.
      PF%loadRefTMin = -huge(PF%loadRefTMin)
      PF%loadRefTMax = huge(PF%loadRefTMax)
      select case(PF%loadInfo(LOAD_TIMETYPE))
      case(LOAD_TIMETYPE_APERIODIC)
        PF%loadAperiodic = .true.
      case(LOAD_TIMETYPE_MTFAPERIODIC)
        ! If this is a Multile Time File Aperiodic case (MTFAperiodic) then
        ! read In the list of Aperiodic Time Files (must be in sequential order)          
        PF%loadAperiodic = .true.
        PF%MTFLoadAperiodic = .true. 
        ! Record how far we ar in the file - use this later to open up the first file a this location
        PF%startSurfaceData =GetBinaryFilePosition(PF%surfaceDataStreamNumber)
        ! Read how many files there are
        call ReadBinaryInteger(PF%surfaceDataStreamNumber,nFiles,stat)
        PF%nbMTFLoadFiles = nFiles
        ! Open the arrays to hold the file names, startKeyIndex, and endKeyIndex for
        ! each of the nFiles files.
        allocate( PF%MTFLoadFileArray(nfiles) )
        allocate( PF%MTFLoadStartKeyIndex(nfiles), PF%MTFLoadEndKeyIndex(nfiles) )
        PF%MTFLoadNKey = 0  !MTFLoadNkey is the total number of keys for all the nFiles files.
        do i=1,nFiles
          ! read list of file names, starting keyIndex, ending keyIndex
          ! Normally the starting keyIndex will equal 1 and ending keyIndex
          ! will equal nKey for the individual file.  This enables the user
          ! to not read all the data in the files - in cases there is overlap
          ! the file writer is responsible to ensure there are no "gaps"
          ! in the data (whether it is time, or some other keyvalue).
          call ReadBinaryString (PF%surfaceDataStreamNumber, filename, stat)
          filename = adjustl(filename)
          if( filename(1:1)=='/' .or. filename(1:2)=='~/' ) then
            ! filename is an absolution path (or at test not relative to the globalFolderName)
            PF%MTFLoadFileArray(i) = trim(filename)
          else
            ! fliename must be a relative patch - assume it is relative to globalFolderName
            PF%MTFLoadFileArray(i) = trim(globalFolderName)//trim(filename)
          end if         
          call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%MTFLoadStartKeyIndex(i), stat)
          call ReadBinaryInteger(PF%surfaceDataStreamNumber, PF%MTFLoadEndKeyIndex(i), stat)
          PF%MTFLoadNKey = PF%MTFLoadNKey + (PF%MTFLoadEndKeyIndex(i)-PF%MTFLoadStartKeyIndex(i)+1)
        end do
        ! Now close the MTF Aperiodic Loading file and open the first file in the
        ! sequence.  Need to open with the correct endianness and skip down to the
        ! place we would be if this were just an aperiodic loading file.
        PF%MTFLoadingFile = PF%patchLoadingFile
        PF%patchLoadingFile=PF%MTFLoadFileArray(1)
        PF%iMTFLoadFile = 1
        call CloseBinaryFile(PF%surfaceDataStreamNumber)
        call OpenBinaryFile (PF%surfaceDataStreamNumber, PF%patchLoadingFile, PF%LoadBigEndian, stat)
        startpos = PF%startSurfaceData + &
                   (PF%MTFLoadStartKeyIndex(PF%iMTFLoadFile)-1)*PF%loadStepInBytes
        call RestartBinaryFileAtByte(PF%SurfaceDataStreamNumber, startpos,stat)
        ! Reset the timetype as if this we read the header of this file.
        PF%LoadInfo(LOAD_TIMETYPE)=LOAD_TIMETYPE_APERIODIC
      case(LOAD_TIMETYPE_QPERIODIC)
    
      case(LOAD_TIMETYPE_MTFQPERIODIC)     
        ! This is a Multile File Qperiodic case (MTFQperiodic), so there will be a list
        ! of Periodic Time Files (must be in sequential order).  We won't know how many
        ! files there are until we reach the end of the MTFQperiodic file, so we have to 
        ! count them.  Assume a big number to start with.
        allocate( PF%MTFLoadFileArray(mxFiles) )
        allocate( PF%MTFLoadRefTMin(mxFiles), PF%MTFLoadRefTMax(mxFiles) )

        do i=1,mxFiles
          ! read list of file names, starting RefTMin, ending RefTMax. The reference times 
          ! (RefTMin and RefTMax) indicate the range of times over which the periodic data is
          ! to be used.  
          call ReadBinaryString (PF%SurfaceDataStreamNumber, filename, stat)
          if (stat < 0 ) exit  !end of file condition
          filename = adjustl(filename)
          if( filename(1:1)=='/' .or. filename(1:2)=='~/' ) then
            ! filename is an absolute path (or at least not relative to the globalFolderName)
            PF%MTFLoadFileArray(i) = trim(filename)
          else
            ! fliename must be a relative patch - assume it is relative to globalFolderName
            PF%MTFLoadFileArray(i) = trim(globalFolderName)//trim(filename)
          end if
          call ReadBinaryReal(PF%SurfaceDataStreamNumber, PF%MTFLoadRefTMin(i), stat)
          if (stat < 0 ) exit  !end of file condition
          call ReadBinaryReal(PF%SurfaceDataStreamNumber, PF%MTFLoadRefTMax(i), stat)
          if (stat < 0 ) exit  !end of file condition
          nfiles = i
          if( nfiles+1 >mxFiles ) then
            call Error('Number of load files greater than the dimension.  Need to change mxFiles and recompile.')
            stop
          end if
        end do
        ! Now close the MTF Aperiodic Loading file and open the first file in the
        ! sequence.  Need to open with the correct endianness and skip down to the
        ! place we would be if this were just an periodic loading file.
        PF%nbMTFLoadFiles = nfiles
        PF%MTFLoadingFile = PF%patchLoadingFile
        PF%patchLoadingFile=PF%MTFLoadFileArray(1)
        PF%iMTFLoadFile = 1
        PF%LoadRefTMin = PF%MTFLoadRefTMin(1)
        PF%LoadRefTmax = PF%MTFLoadRefTMax(1)
        call CloseBinaryFile(PF%SurfaceDataStreamNumber)
        call OpenBinaryFile (PF%SurfaceDataStreamNumber, PF%patchLoadingFile, PF%LoadBigEndian, stat)
        startpos = PF%startSurfaceData
        call RestartBinaryFileAtByte(PF%SurfaceDataStreamNumber, startpos,stat)
        ! Reset the timetype as if this we read the header of this file.
        PF%loadInfo(LOAD_TIMETYPE)=LOAD_TIMETYPE_PERIODIC
        PF%MTFloadQperiodic = .true.
      case default  ! constant, periodic, or qperiodic loading      

      end select

      ! Determine which patches have loading

      call ReadBinaryInteger(PF%surfaceDataStreamNumber,nbLoadedZones,stat)

      do i=1,nbLoadedZones
        call ReadBinaryInteger(PF%surfaceDataStreamNumber,patchNum,stat)
        if (abs(patchNum) > PF%geoInfo(GEO_NBZONES) ) then
          call Error ("Invalid zone specification - check loading file format.")
          stop
        end if
        PF%hasLoading(abs(patchNum)) = .true.
        ! if patchNum < 0, then don't calculate thickness noise for the patch (e.g., compact patches)
        if( patchNum < 0 ) PF%compThickness(abs(patchNum)) = .false.
      end do
    end if
    
    if(PF%surfaceDataStreamNumber/=0 ) then
      PF%startSurfaceData =GetBinaryFilePosition(PF%surfaceDataStreamNumber)  !ksb debug: 5/22/16
    end if
      
  end subroutine ReadFileHeaderV1_0
  
  
  !subroutine ReadTopOfFile(PF)
  !  type(PatchFile), intent(inout)::PF
  !  integer:: surface = 1, loading = 2
  !  integer:: geoStreamNumber, loadStreamNumber, stat, temp, i
  !  character(len=32):: geoTitle, loadTitle
  !  
  !  !!ksb debug:
  !  ! reposition to the start of the data
  !  !print*,'PF%startGeometryData=',PF%startGeometryData,' PF%startSurfaceData=',PF%startSurfaceData
  !  if (PF%geoAperiodic) then
  !    call RestartBinaryFileAtByte(PF%geometryStreamNumber,PF%startGeometryData,stat)
  !  end if
  !  if (PF%loadAperiodic) then
  !    call RestartBinaryFileAtByte(PF%surfaceDataStreamNumber,PF%startSurfaceData,stat)
  !  end if
  !  return
  !  
  !  call DeterminePatchFileFormat(PF, PF%loadAperiodic, .false.) 
  !  call ReadFileHeaderV1_0(PF, PF%geoAperiodic, PF%loadAperiodic, .false.)
  !
  !  do i =1,size(PF%patchArray)
  !    if (PF%geoAperiodic) then
  !      geoStreamNumber = PF%geometryStreamNumber
  !      geoTitle = GetGeoTitle(PF%patchArray(i))
  !      call ReadBinaryString(geoStreamNumber,geoTitle, stat)
  !      call ReadBinaryInteger(geoStreamNumber, temp, stat)
  !    end if
  !    if (PF%hasLoading(i) .and. PF%loadAperiodic) then
  !      loadStreamNumber = PF%surfaceDataStreamNumber
  !      loadTitle = GetLoadTitle(PF%patchArray(i)) 
  !      call ReadBinaryString(loadStreamNumber,loadTitle, stat)
  !    end if
  !    if( PF%geoAperiodic .or. PF%loadAperiodic ) then
  !      call GetDimensions(PF%patchArray(i), geoStreamNumber, loadStreamNumber, PF%loadInfo)
  !    end if
  !  end do
  !  !!ksb debug:    integer::startGeometryData,startSurfaceData
  !  if(PF%geoAperiodic) then
  !      PF%startGeometryData=GetBinaryFilePosition(PF%geometryStreamNumber)
  !  else
  !      PF%startGeometryData=-1
  !  end if
  !  if(PF%loadAperiodic) then
  !      PF%startSurfaceData =GetBinaryFilePosition(PF%surfaceDataStreamNumber)
  !  else
  !     PF%startSurfaceData = -1
  !  end if
  !  print*,'PF%startGeometryData1=',PF%startGeometryData,' PF%startSurfaceData1=',PF%startSurfaceData
  !end subroutine ReadTopOfFile  
        
    
  recursive subroutine FindContainerForAttachedObs(C, attachedTo, CBList, flag)
    type(container)::C
    character(len=*)::attachedTo
    type(CBStructure), pointer::CBList
    logical::Flag
    integer::i
    if (trim(C%title).eq.trim(attachedTo)) then
      if (.not. associated(C%CBList)) then
        call Warning('Cannot attach observer to container with no changes of Base', &
                     'Observer will have zero changes of base')
        return
      else
        call appendCB(CBList, C%CBList)
      end if
      flag=.true.
    else
      do i=1, C%nbContainer
        call FindContainerForAttachedObs(C%containerArray(i), &
                                         attachedTo, CBList, Flag)
      end do
    end if

  end subroutine FindContainerForAttachedObs

  recursive subroutine WriteContainerDebugInfo(C,unitnum)
    type(container), intent(in)::C
    integer::i,unitnum
    write(unitnum,*) '*** ContainerDebugInfo ***' 
    write(unitnum,*) 'Title= ', trim(C%title)
    write(unitnum,*) 'Prefix= ', trim(C%prefix)
    write(unitnum,*) 'PrefixLength= ', trim(integertostring(C%prefixLength))
    write(unitnum,*) 'Periodic blades? ', C%periodicBlades
    write(unitnum,*) 'NbContainer= ', trim(integertostring(C%nbContainer))
    write(unitnum,*) 'NbBase= ', trim(integertostring(C%nbBase))
    write(unitnum,*) 'NbFiles= ', trim(integertostring(C%nbFiles))
    write(unitnum,*) 'NTau= ', trim(integertostring(C%nTau))
    write(unitnum,*) 'CBList associated? ', associated(C%CBList)
    if (associated(C%CBList)) then
      call WriteCBListDebugInfo(C%CBList, unitnum)
    end if
    write(unitnum,*) 'Patch file array associated? ', associated(C%PFArray)
    if (associated(C%PFArray)) then
      write(unitnum,*) 'Size of patch file array is ', trim(integertostring(size(C%PFArray)))
      do i=1, size(C%PFArray)
        call WritePatchFileDebugInfo(C%PFArray(i), unitnum)
      end do
    end if
    write(unitnum,*) 'Container array associated? ',associated(C%containerArray)
    if (associated(C%containerArray)) then
      write(unitnum,*) 'Size of container array ', trim(integertostring(size(C%containerArray)))
      do i=1, size(C%containerArray)
        call WriteContainerDebugInfo(C%containerArray(i), unitnum)
      end do
    end if
    write(unitnum,*) '--- End container Debug Info ---'
    call Barrier()
  end subroutine WriteContainerDebugInfo

  subroutine WritePatchFileDebugInfo(Pfile,unitnum)
    type(PatchFile), intent(in)::Pfile
    integer::i,unitnum
 
    write(unitnum,*) '*** PatchFile Debug Info ***' 
    write(unitnum,*) 'PatchGeometryFile= ', trim(Pfile%patchGeometryFile)
    write(unitnum,*) 'PatchLoadingFile= ', trim(Pfile%patchLoadingFile)
    write(unitnum,*) 'GeometryStreamNumber= ', trim(integertostring(Pfile%geometryStreamNumber))
    write(unitnum,*) 'SurfaceDataSteamNumber= ', trim(integertostring(Pfile%surfaceDataStreamNumber))
    write(unitnum,*) 'GeometryMajorVerNum= ', trim(integertostring(Pfile%geometryMajorVerNum))
    write(unitnum,*) 'GeometryMinorVerNum= ', trim(integertostring(Pfile%geometryMinorVerNum))
    write(unitnum,*) 'SurfaceDataMajorVerNum= ', trim(integertostring(Pfile%surfaceDataMajorVerNum))
    write(unitnum,*) 'SurfaceDataMinorVerNum= ', trim(integertostring(Pfile%surfaceDataMinorVerNum))
    write(unitnum,*) 'ParallelAllocatedFlag= ',  Pfile%parallelAllocatedFlag
    write(unitnum,*) 'PressureUnit= ', trim(Pfile%pressureUnit)
    write(unitnum,*) 'Size of hasLoading ', trim(integertostring(size(Pfile%hasLoading)))
    if(associated(Pfile%hasLoading)) then
      do i=1, size(Pfile%hasLoading)
        write(unitnum,*) 'Patch Numnber=', i, 'has Function Data=', PFile%hasLoading(i)
      end do
    end if
    write(unitnum,*) 'Size of geoInfo ', trim(integertostring(size(Pfile%geoInfo)))
    do i=1, size(Pfile%geoInfo)
      write(unitnum,*) 'GeoInfo value ', trim(integertostring(i)), ' has value of ', &
                  trim(integertostring(PFile%geoInfo(i))) 
    end do
    write(unitnum,*) 'Size of loadInfo is ', trim(integertostring(size(Pfile%loadInfo)))
    do i=1, size(Pfile%loadInfo)
      write(unitnum,*) 'LoadInfo value ', trim(integertostring(i)), ' has value of ', &
                  trim(integertostring(PFile%loadInfo(i))) 
    end do
    write(unitnum,*)'Size of patchArray ', trim(integertostring(size(Pfile%patchArray)))
    do i=1, size(Pfile%patchArray) 
      call WritePatchDebugInfo(Pfile%patchArray(i), unitnum)
    end do
    write(unitnum,*) '--- End PatchFile Debug Info ---'
  end subroutine WritePatchFileDebugInfo
 

  !*
  !SUBROUTINE integrateContainer(C, obsVector, obsCBList, thicknessNoise, loadingNoise)
  !   This subroutine takes the container then integrates over each of the 
  !  patches it contains then integrates over each subobject it contains
  !ARGUMENTS:
  !  - C: The container we are integrating over.
  !  - obsVector: The vector we are integrating with respect to.
  !  - obsCBList: The CBList form the base frame to the observer.
  !  - thicknessNoise: the thickness noise that will be created and returned.
  !  - loadingNoise: the loading noise that will be created and returned.
  !*
  recursive subroutine IntegrateContainer(C, obsInitVector, obsCBList, thicknessNoise, loadingNoise, &
                                          PTGradient,PL1Gradient, PeggIntensity, BPMintensity, iBlankArray, &
                                          dt, obsCount, observerTimeArray, iCont, nbCont, &
                                          nbSegments, segmentSize, radDistance, PeggData, BPMData)
    use constantsModule, only: aperiodicArraySize
    type(Container), intent(inout)::C
    type(Vector), intent(in)::obsInitVector
    type(CBStructure),pointer::obsCBList
    real(kind=8), dimension(:), intent(inout)::thicknessNoise, loadingNoise
    type(vectorKind8), dimension(:), intent(inout)::PTGradient,PL1Gradient    
    real(kind=8), dimension(:,:), pointer:: PeggIntensity
    real(kind=8), dimension(:,:,:), pointer:: BPMintensity
    integer, dimension(:)::iBlankArray
    real(kind=8), dimension(:), intent(in)::observerTimeArray
    real(kind=8), intent(inout):: dt
    real:: segmentSize
    real, dimension(:), pointer:: radDistance
    integer:: obsCount, iCont, nbCont, stat, int1, int2
    integer, intent(inout):: nbSegments
    type(PEGG), pointer:: PeggData
    type(BPM), pointer:: BPMData
    type(PatchFile), pointer:: PF

    integer::i, j, k, iTau, nTau, nTauMax, obsNT, isCompact, nSect
    integer:: keyIndexSurf, keyIndexLoad, maxTauRange, startpos
    integer, dimension(:), allocatable::patchNTau
    integer, pointer:: BBiTau=>null()
    real:: tMax, tauMax, tauRange, initialR, dTau
    real, dimension(:), pointer:: tauArray=>null()
    logical:: lastGeoZone, lastBlade
    character (len=32):: stemp32
    integer:: itemp
    real:: rtemp
    
    ! Integrate over all container objects
    if(debugLevel.ge.1.and.IsMaster()) write(*,*) C%prefix(1:C%prefixLength)//&
                                       'Integrating ', trim(C%title)
    ! Reset the min/max observer time range in the container.                                       
    C%minObsTimeInData = -huge(C%minObsTimeInData); C%maxObsTimeInData = huge(C%maxObsTimeInData)

    ! With regard to semi-empirical broadband noise prediction methods:
    ! We are assuming constant geometry which means the user
    ! must use rotation=.true. to define the fundamental (rotor) rotation
!ksb debug: 2/23/16 - need to check to see if constant geometry is really needed.  It looks like
! we only need rotation=.true.; which means the geometry does not have to be constant if this is
! still true.
    if (C%parentRotation) then
      iCont  = 0
      nbCont = C%nbContainer
      nullify(PeggData)
      nullify(BPMData)
    else
      if (iCont.ge.0) iCont = iCont + 1
    end if
    if (iCont.ne.nbCont) then
      lastBlade = .false.
    else
      lastBlade = .true.
    end if

    if (C%PeggNoiseFlag) then
      PeggData => GetPeggData(C%BBData)
    elseif (C%BPMNoiseFlag) then
      BPMData  => GetBPMData(C%BBData)

      if (associated(BPMdata%obsTimeArray)) then
        deallocate(BPMdata%obsTimeArray)
        nullify(BPMdata%obsTimeArray)
      end if
      allocate(BPMdata%obsTimeArray(size(BPMIntensity,3)))
      call BuildBroadbandObsTimeArray(BPMData%obsTimeArray, observerTimeArray)
    end if

    if (associated(PeggData)) then
      PeggData%iTau = 1
      if (associated(PeggData%obsTimeArray)) then
        if(size(PeggData%obsTimeArray)>0) deallocate(PeggData%obsTimeArray) !ksb debug: don't know why size <0 sometimes
        nullify(PeggData%obsTimeArray)
      end if
      allocate(PeggData%obsTimeArray(size(PeggIntensity,2)))
      call BuildBroadbandObsTimeArray(PeggData%obsTimeArray, observerTimeArray)
      ! Grab the hub CB list for Pegg's method
      PeggData%CBList => C%CBList
    elseif( globalBPMNoiseFlag ) then
      BPMData%iTau = 1

      if (associated(BPMData%Intensity)) then
        deallocate(BPMData%Intensity)
        nullify(BPMData%intensity)
      end if
      allocate(BPMData%Intensity(BPM_TERMS, BB_THIRD_OCTAVE_BANDS, nbSegments))
      BPMData%Intensity = 0.
    end if

    do i=1, C%nbFiles
      PF => C%PFArray(i)  ! use a pointer to make the following code easier to follow
      !
      ! Get the files in proper state to start.
      call PatchFileSetup(PF)
                 
      allocate(patchNTau(PF%geoInfo(GEO_NBZONES)))
      nTau = 0
      do j=1, PF%geoInfo(GEO_NBZONES)
        if (PF%geoAperiodic) then
          nTau = GetPatchSurfacenKey(PF%patchArray(j))
        end if
        if (PF%loadAperiodic) then
          if (GetAssociatedLoad(PF%patchArray(j))) then
            nTau = GetPatchLoadnKey(PF%patchArray(j)) 
          end if      
        end if
      end do

      if( PF%MTFGeoAperiodic .or. PF%MTFloadAperiodic ) then
        ntau = max(PF%MTFGeoNKey,PF%MTFLoadNKey)
        call SetPatchFileNKey(PF,ntau,PF%MTFGeoNKey,PF%MTFLoadNKey)
      end if
      obsNT = size(observerTimeArray)

      ! aperiodicArraySize = 4 (defined as a parameter in constantsModule) 
      ! We need these (4) points [L C R RR] in time for our calculations
      
      ! With the new design for aperiodic cases we are forced to split up the read and calculation segments
      ! as we must first read in enough data to calculate the first couple of vel, accel, etc. terms.  
      do iTau=1, aperiodicArraySize
        do j=1,PF%geoInfo(GEO_NBZONES)
          if( (PF%geoAperiodic .or. PF%loadAperiodic) .and. &
               PatchITauInRange(PF%patchArray(j)) ) then ! true if patch%itau <= patch%nTau
            call ReadPatchAtTime(PF%patchArray(j), PF%geometryStreamNumber, &
                                 PF%surfaceDataStreamNumber, PF%geoInfo, &
                                 PF%loadInfo, C%periodicKeyOffset,PF%hasLoading(j), & 
                                 PF%geometryMajorVerNum)
          end if
           
          ! Get source time array and calculate dTau
          if (iTau == 2) then
            if (.not.(newSegLHS .or. newSegRHS)) call ResetGlobals()
            
                        
            call DeterminePatchDtau(PF%patchArray(j),C%tauMin, C%tauMax,&
                                    C%nTau, C%dTau)
            
            ! Generate the source time array.  If any part of this case is aperiodic 
            ! then this is taken as from the aperiodic data, if not then it is 
            ! determined using CreateSourceTimeArray which is within GetPatchTimeArrays
            call GetPatchTimeArrays(PF%patchArray(j), C%CBList, obsInitVector, ObsCBList, &
                                    observerTimeArray, nTau, C%tauMin, &
                                    PF%geoInfo, PF%loadInfo)
                              
            if (PF%geoAperiodic .or. PF%loadAperiodic) then
              ! If this is a new segment then these arrays already exist.
              if (.not.(newSegLHS .or. newSegRHS)) call CreatePatchTimeRangeArrays(PF%patchArray(j), obsCount)
            end if
            
            ! This routine allocates the storedDataArray and also points the current patch module global variables
            ! to similarly-named variables within storedDataArray.  This allows us to maintain the structure of 
            ! Patch.f90 while still working simultaneously with multiple patches.
            call InitializeStoredData(PF%patchArray(j))
          end if                   
        end do
      end do

      ! This allows us to create a so-called ntauMax that lets us loop through integratePatch as many times as
      ! necessary.  If the current patch has fewer source times than ntauMax it is simply skipped over.
      nTauMax = 0
      do j=1,PF%geoInfo(GEO_NBZONES)
        patchNTau(j) = GetPatchNTau(PF%patchArray(j))
        if( PF%MTFGeoAperiodic .or. PF%MTFLoadAperiodic ) then
          nTauMax = PF%MTFLoadNKey
        else
          nTauMax = max(nTauMax,patchNTau(j))
        end if
      end do

      tauArray => GetPatchTauArray(PF%patchArray(maxloc(patchNTau,1)))
      dTau = tauArray(2) - tauArray(1)
      ! This only needs to be done once.
      if (i.eq.1 .and. (globalPeggNoiseFlag .or. globalBPMNoiseFlag)) then
        tauRange = -huge(tauRange)
        ! Get the widest range of source time to capture everything.
        do j=1,PF%geoInfo(GEO_NBZONES)
          ! Compare the period and stepsize of the BB data files with any loading data          
          if (PF%hasLoading(j)) then
            if (associated(PeggData)) then
              if (PeggData%timeType.ne.LOAD_TIMETYPE_CONSTANT) call CompareDataTimeInfo(PeggData%nTau, &
                  PeggData%period, PF%patchArray(j))
            end if
            if( globalBPMNoiseFlag ) then
              if (BPMData%timeType.ne.LOAD_TIMETYPE_CONSTANT) call CompareDataTimeInfo(BPMData%nTau, &
                  BPMData%period, PF%patchArray(j))
            end if
          end if
        end do
        call BuildBroadbandData(C%PFArray, C%CBList, PeggData, BPMData, &
			        nTauMax, PF%hasLoading, C%periodicKeyOffset, lastBlade, nbSegments)
        initialR = vectorAbsolute(obsInitVector)
      end if

      do iTau=1, aperiodicArraySize-2
        if (globalPeggNoiseFlag) then
          if (getPeggTimetype(PeggData).eq.GEO_TIMETYPE_APERIODIC) call readPeggData(PeggData, C%PeggUnitNumber,1)
          if(ComputePeggTerm(PeggData, PEGG_TOTALTHRUST) .and. PF%loadAperiodic) then
            PeggData%shaftAxisThrust = 0.
            do j=1,PF%geoInfo(GEO_NBZONES)
              call ComputePatchThrust(PeggData, PF%patchArray(j), C%CBlist, &
                                      PF%loadInfo(LOAD_TIMETYPE), 1, iTau, iTau)
            end do
          end if
          if (lastBlade .and. i.eq.C%nbFiles .and. PeggData%iTau.le.nbSegments) then
            ! Get the current hub position
            call CalculateResultantChangeOfBase(PeggData%CBList, GetSize(PeggData%CBList), tauArray(iTau))
            call PredictPeggBB(PeggData, initialR, obsCBList, tauArray(iTau), dt, dTau, iTau+PeggData%iTauMin-1, segmentSize)
          end if
        elseif (globalBPMNoiseFlag) then
          if (getBPMTimetype(BPMdata).eq.GEO_TIMETYPE_APERIODIC) call readBPMData(BPMData, C%BPMUnitNumber, 1)
        end if

        do j=1,PF%geoInfo(GEO_NBZONES)
          if( globalPeggNoiseFlag .or. globalBPMNoiseFlag ) then
            if (j.ne.PF%geoInfo(GEO_NBZONES)) then
              lastGeoZone=.false.
            else
              lastGeoZone=.true.
            end if
          end if
          !
          !ksb debug: need to figure out correct entry for taumax, nbZones arguement (everything after iTau)
          call IntegratePatch(PF%patchArray(j), BPMData, obsInitVector, obsCBList, C%CBlist, &
          &   iBlankArray, thicknessNoise, loadingNoise, PTGradient, PL1Gradient, observerTimeArray, &                      
          &   dt, iTau, patchNTau(j),iTau, iTau, j, PF%geoInfo(GEO_NBZONES), obsNT, obsCount, &
          &   lastGeoZone, lastBlade, radDistance, segmentSize)
        end do
      end do

      ! This allows us to quickly move through CalculateDifferentialNoise using most of the pre-existing code
      ! without having to worry about the size of the arrays.  For aperiodic data we only save 4 points in time,
      ! but we can use the same basic code if we simply redo our definition of iTau as this new variable keyIndex,
      ! whose value is dependent on whether or not the case is aperiodic.      
      keyIndexSurf = aperiodicArraySize - 2
      keyIndexLoad = aperiodicArraySize - 2
      ! (aperiodicArraySize-1) is one ahead of the Center position (R).  The first set 
      ! of LL, L and R values require a bit of specific logic to deal with but after 
      ! this everything is simply shifted to the left.
      do iTau = aperiodicArraySize-1,nTauMax
        if( tauArray(iTau)>PF%geoRefTMax )then  
          !move to next Reference Period, and next set of periodic data.
          if( PF%MTFgeoQperiodic ) then  
            !this is Multi-File Quasi-Periodic case, so need to move to next file.
            PF%iMTFGeoFile=PF%iMTFGeoFile+1
            PF%geoRefTMin = PF%MTFGeoRefTMin(PF%iMTFGeoFile)
            PF%geoRefTMax = PF%MTFGeoRefTMax(PF%iMTFGeoFile)
            PF%patchGeometryFile=PF%MTFGeoFileArray(PF%iMTFGeoFile)
            call CloseBinaryFile(PF%GeometryStreamNumber)
            call OpenBinaryFile (PF%GeometryStreamNumber, PF%patchGeometryFile, PF%GeoBigEndian, stat)
            startpos = PF%startGeometryData 
            call RestartBinaryFileAtByte(PF%GeometryStreamNumber, startpos,stat)
            do j=1,PF%geoInfo(GEO_NBZONES)
              call ResetPatchKeyCount(PF%patchArray(j),forceGeo=.true.)
              select case (PF%geoInfo(GEO_GRIDTYPE))
              case (GEO_GRIDTYPE_STRUCTURED)
                call ReadBinaryString(PF%GeometryStreamNumber, stemp32, stat)
                call ReadBinaryReal(PF%GeometryStreamNumber, rtemp , stat)
                call ReadBinaryInteger(PF%GeometryStreamNumber, PF%MTFgeoNKey, stat)
                call ReadBinaryInteger(PF%GeometryStreamNumber, itemp, stat)
                call ReadBinaryInteger(PF%GeometryStreamNumber, itemp, stat)
              case (GEO_GRIDTYPE_UNSTRUCTURED)
                print*,'IntegrateContainer(container.f90) - qperiod unstructured geometry not implemented.'
              !  call ReadBinaryString(PF%GeometryStreamNumber, PF%geoTitle, stat)
              !  call ReadBinaryReal(PF%GeometryStreamNumber, PF%geoRefTMin, stat)
              !  call ReadBinaryInteger(PF%GeometryStreamNumber, PF%geoRefTMin, stat)
              !  call ReadBinaryInteger(PF%GeometryStreamNumber, PF%geoRefTMin, stat)
              !  call ReadBinaryInteger(PF%GeometryStreamNumber, PF%geoRefTMin, stat)
              !  ! need to read connectivity here
              end select
            end do
          else 
            !this is a Quasi-Periodic case, so read the next Reference Period.
            call ReadBinaryReal(PF%GeometryStreamNumber, PF%geoRefTMin, stat)
            call ReadBinaryReal(PF%GeometryStreamNumber, PF%geoRefTMax, stat)
          end if
          ! Now read the periodic geometry data
          rtemp=0.
          do j=1,PF%geoInfo(GEO_NBZONES)
            call ResetPatchKeyCount(PF%patchArray(j),forceGeo=.true.)
          end do
          do itemp=1, PF%MTFgeoNKey
            do j=1,PF%geoInfo(GEO_NBZONES)
              call ReadPatchGeomAtTime(PF%patchArray(j),PF%geometryStreamNumber,&
                PF%geoInfo, rtemp, PF%geometryMajorVerNum)              
              !if( PF%hasLoading(j) ) then
              !  call ReadPatchLoadAtTime(PF%patchArray(j), PF%surfaceDataStreamNumber, &
              !    PF%loadInfo, rtemp, PF%surfaceDataMajorVerNum)                
              !end if
            end do
          end do
        end if
        if( tauArray(iTau)>PF%loadRefTMax )then  
          !move to next Reference Period, and next set of periodic data.
          if( PF%MTFloadQperiodic ) then  
            !this is Multi-File Quasi-Periodic case, so need to move to next file.
            PF%iMTFLoadFile=PF%iMTFLoadFile+1
            PF%loadRefTMin = PF%MTFloadRefTMin(PF%iMTFLoadFile)
            PF%loadRefTMax = PF%MTFloadRefTMax(PF%iMTFLoadFile)
            PF%patchLoadingFile=PF%MTFLoadFileArray(PF%iMTFLoadFile)
            call CloseBinaryFile(PF%SurfaceDataStreamNumber)
            call OpenBinaryFile (PF%SurfaceDataStreamNumber, PF%patchLoadingFile, PF%LoadBigEndian, stat)
            startpos = PF%startSurfaceData 
            call RestartBinaryFileAtByte(PF%SurfaceDataStreamNumber, startpos,stat)
            do j=1,PF%loadInfo(LOAD_NBZONES)
              call ResetPatchKeyCount(PF%patchArray(j),forceLoad=.true.)
              select case (PF%loadInfo(LOAD_GRIDTYPE))
              case (LOAD_GRIDTYPE_STRUCTURED)
                call ReadBinaryString(PF%SurfaceDataStreamNumber, stemp32, stat)
                call ReadBinaryReal(PF%SurfaceDataStreamNumber, rtemp , stat)

                call ReadBinaryInteger(PF%SurfaceDataStreamNumber, PF%MTFloadNKey, stat)
                call ReadBinaryInteger(PF%SurfaceDataStreamNumber, itemp, stat)
                call ReadBinaryInteger(PF%SurfaceDataStreamNumber, itemp, stat)
              case (LOAD_GRIDTYPE_UNSTRUCTURED)
                print*,'IntegrateContainer(container.f90) - qperiod unstructured loading not implemented.'

              !  call ReadBinaryString(PF%GeometryStreamNumber, PF%geoTitle, stat)
              !  call ReadBinaryReal(PF%GeometryStreamNumber, PF%geoRefTMin, stat)
              !  call ReadBinaryInteger(PF%GeometryStreamNumber, PF%geoRefTMin, stat)
              !  call ReadBinaryInteger(PF%GeometryStreamNumber, PF%geoRefTMin, stat)
              !  call ReadBinaryInteger(PF%GeometryStreamNumber, PF%geoRefTMin, stat)
              !  ! need to read connectivity here
              end select
            end do
          else 
            !this is a Quasi-Periodic case, so read the next Reference Period.
            call ReadBinaryReal(PF%SurfaceDataStreamNumber, PF%loadRefTMin, stat)
            if( stat /= 0 ) then
              call Error("Error reading loading file, loadRefTMin")
            end if
            call ReadBinaryReal(PF%SurfaceDataStreamNumber, PF%loadRefTMax, stat)
            if( stat /= 0 ) then
              call Error("Error reading loading file, loadRefTMin")
            end if
          end if
          ! Now read the periodic geometry data
          rtemp=0.
          do j=1,PF%loadInfo(LOAD_NBZONES) !PF%geoInfo(GEO_NBZONES)
            call ResetPatchKeyCount(PF%patchArray(j),forceLoad=.true.) !forceGeo=.true.)
          end do
          do itemp=1, PF%MTFloadNKey
            do j=1,PF%loadInfo(LOAD_NBZONES)
              !call ReadPatchGeomAtTime(PF%patchArray(j),PF%geometryStreamNumber,&
              !  PF%geoInfo, rtemp, PF%geometryMajorVerNum)              
              if( PF%hasLoading(j) ) then
                call ReadPatchLoadAtTime(PF%patchArray(j), PF%surfaceDataStreamNumber, &
                  PF%loadInfo, rtemp, PF%surfaceDataMajorVerNum)                
              end if
            end do
          end do
        end if
        
        if (.not. PF%geoAperiodic) then
          keyIndexSurf = iTau
        end if
        if (.not. PF%loadAperiodic) then
          keyIndexLoad = iTau
        end if   

        !BBobsTimeArray = observerTimeArray(1+(j-1)*nTauPerSegment)
        if (globalPeggNoiseFlag) then
          if (getPeggTimetype(PeggData).eq.GEO_TIMETYPE_APERIODIC) call readPeggData(PeggData, C%PeggUnitNumber,1)
          if(ComputePeggTerm(PeggData, PEGG_TOTALTHRUST) .and. PF%loadAperiodic) then
            PeggData%shaftAxisThrust = 0.
            do j=1,PF%geoInfo(GEO_NBZONES)
              call ComputePatchThrust(PeggData, PF%patchArray(j), C%CBlist, &
                      PF%loadInfo(LOAD_TIMETYPE), 1, keyIndexSurf, keyIndexLoad)
            end do
          end if
          if (lastBlade .and. i.eq.C%nbFiles .and. PeggData%iTau.le.nbSegments) then
            call CalculateResultantChangeOfBase(PeggData%CBList, GetSize(PeggData%CBList), tauArray(iTau))
            call PredictPeggBB(PeggData, initialR, obsCBList, tauArray(iTau), dt, dTau,  &
                keyIndexLoad+PeggData%iTauMin-1, segmentSize)
          end if
        elseif (globalBPMNoiseFlag) then
          if (getBPMTimetype(BPMData).eq.GEO_TIMETYPE_APERIODIC) call readBPMData(BPMData, C%BPMUnitNumber, 1)
        end if

        do j=1,PF%geoInfo(GEO_NBZONES)
          if( globalPeggNoiseFlag .or. globalBPMNoiseFlag) then
            if (j.ne.PF%geoInfo(GEO_NBZONES)) then
              lastGeoZone=.false.
            else
              lastGeoZone=.true.
            end if
          end if
          ! For aperiodic data, we need to rotate the data pointers and then read the
          ! patch data at the enxt time step.
          if (PF%geoAperiodic .or. PF%loadAperiodic) then
            call RotatePatchTimeData(PF%patchArray(j))
            if( PatchITauInRange(PF%patchArray(j)) ) then ! true if patch%itau <= patch%nTau
              call ReadPatchAtTime(PF%patchArray(j), PF%geometryStreamNumber, &
                                   PF%surfaceDataStreamNumber, PF%geoInfo, &
                                   PF%loadInfo, C%periodicKeyOffSet,PF%hasLoading(j), & 
                                   PF%geometryMajorVerNum)
            end if
          end if
          
          ! This lets us easily include the read and integrate subroutines together,
          ! since we can simply skip over the integration step for a patch
          ! once it has reached its unique nTau.
          if (iTau <= patchNTau(j)) then
            !ksb debug: need to figure out correct entry for taumax, nbZones arguement (everything after iTau)
            call IntegratePatch(PF%patchArray(j), BPMData, obsInitVector, obsCBList, C%CBlist, &
            &   iBlankArray, thicknessNoise, loadingNoise, PTGradient, PL1Gradient, observerTimeArray, &                      
            &   dt, iTau, patchNTau(j), keyIndexSurf, keyIndexLoad, j, PF%geoInfo(GEO_NBZONES), &
            &   obsNT, obsCount, lastGeoZone, lastBlade, radDistance, segmentSize)
          end if
        end do
        !ksb debug: 10/14/2013 This doesn't work properly because the check at the end of ReadAperiodicLoadingAtTime
        !"if (newLoading%keycount == newLoading%nkey) then" doesn't work, and  eloadFlag doesn't
        ! get set to .true. and (more importantly?) newLoading%keycount doesn't get set to 0.
        !commenting this if statement out means it will read all the time steps even if it doesn't need to.
        !if (iTau .gt. maxval(patchNTau)) exit
        
        ! iTau is essentially the global index for the MTF case (multiple time file case).  We need to 
        ! keep track of which timestep we are in for the CURRENT file.
        if( PF%MTFGeoAperiodic ) then
          if( PF%iMTFGeoFile == 1 ) then
            int1 = (iTau-PF%MTFPrevGeoFileNKey)+2
          else
            int1 = (iTau-PF%MTFPrevGeoFileNKey)
          end if
          int2 = (PF%MTFGeoEndKeyIndex(PF%iMTFGeoFile)-PF%MTFGeoStartKeyIndex(PF%iMTFGeoFile)+1)
          if(  int1 == int2 ) then
            call CloseBinaryFile(PF%GeometryStreamNumber)
            if( iTau == PF%MTFGeoNKey ) exit  ! I don't think this works, we have actually read a couple of steps ahead
            PF%MTFPrevGeoFileNKey=iTau
            PF%iMTFGeoFile=PF%iMTFGeoFile+1
            if( PF%iMTFGeoFile>PF%nbMTFGeoFiles ) exit
            PF%patchGeometryFile=PF%MTFGeoFileArray(PF%iMTFGeoFile)
            call OpenBinaryFile (PF%GeometryStreamNumber, PF%patchGeometryFile, PF%GeoBigEndian, stat)
            startpos = PF%startGeometryData + &
                       (PF%MTFGeoStartKeyIndex(PF%iMTFGeoFile)-1)*PF%geomStepInBytes
            call RestartBinaryFileAtByte(PF%GeometryStreamNumber, startpos,stat)
          end if
        end if
        if( PF%MTFLoadAperiodic) then
          if( PF%iMTFLoadFile == 1 ) then
            int1 = (iTau-PF%MTFPrevLoadFileNKey)+2
          else
            int1 = (iTau-PF%MTFPrevLoadFileNKey)
          end if
          int2 = (PF%MTFLoadEndKeyIndex(PF%iMTFLoadFile)-PF%MTFLoadStartKeyIndex(PF%iMTFLoadFile)+1)  
          if( int1 == Int2 ) then
            call CloseBinaryFile(PF%SurfaceDataStreamNumber)
            If( iTau == PF%MTFLoadNKey ) exit  ! I don't think this works, we have actually read a couple of steps ahead
            PF%MTFPrevLoadFileNKey=iTau
            PF%iMTFLoadFile=PF%iMTFLoadFile+1
            if( PF%iMTFLoadFile > PF%nbMTFLoadFiles ) exit
            PF%patchLoadingFile=PF%MTFLoadFileArray(PF%iMTFLoadFile)
            call OpenBinaryFile (PF%SurfaceDataStreamNumber, PF%patchLoadingFile, PF%LoadBigEndian, stat)
            startpos = PF%startSurfaceData + &
                       (PF%MTFLoadStartKeyIndex(PF%iMTFLoadFile)-1)*PF%loadStepInBytes
            call RestartBinaryFileAtByte(PF%SurfaceDataStreamNumber, startpos,stat) 
          end if
        end if
      end do
      ! In case the observer time range needs to be expanded, the minimum and maximum observer times from the data
      ! must be stored at this point so they can be compared to the minimum and maximum of a new observer time range.
      if (PF%geoAperiodic .or. PF%loadAperiodic) then
        do j=1,PF%geoInfo(GEO_NBZONES)  
          C%minObsTimeInData = max( C%minObsTimeInData, &
            getPatchMinObsTimeInData(PF%patchArray(j),obsCount) )
          C%maxObsTimeInData = min( C%maxObsTimeInData, &
            getPatchMaxObsTimeInData(PF%patchArray(j),obsCount) )
        end do
        
        if (C%maxObsTimeInData.lt.observerTimeArray(obsNT)) then
          tMax = observerTimeArray(obsNT)
          tauMax = GetPatchMaxSourceTime(PF%patchArray(1))
          call Warning('There is insufficient data in the aperiodic data file(s) to',&
                    'obtain pressure information through the observer time range requested.')
          !ksb debug:  The previous warning is probably enough.  Don't want this write.
          !write(*,'(2(A,G15.6,/))') ' The maximum observer time tMax =', trim(RealToString(tMax)), &
          !       ' has a respective source time after tauMax =', trim(RealToString(tauMax))           
          dataLimitRHS = .true.                
        end if
      end if

      if (i.eq.C%nbFiles .and. lastBlade .and. globalPeggNoiseFlag) then
        do j=1,nbSegments
!ksb debug: 2/23/16
          if(PeggData%nbContributions(j)==0) cycle
          PeggIntensity(:,j) = PeggIntensity(:,j) + PeggData%Intensity(:,j)/PeggData%nbContributions(j)
        end do
      elseif (i.eq.C%nbFiles .and. globalBPMNoiseFlag ) then
        do j=1,nbSegments
          !   Compute total Intensity for the segment for all mechanisms
          !   --------------------------------------------------------
          do k=1,BB_THIRD_OCTAVE_BANDS
            BPMdata%Intensity(BPM_SPLTOTAL,k,j) = sum(BPMData%Intensity(1:(BPM_SPLTOTAL-1),k,j))
          end do
          BPMIntensity(:,:,j) = BPMIntensity(:,:,j) + BPMData%Intensity(:,:,j)/BPMData%nbContributions(j)
        end do
      end if
      deallocate(patchNTau) 
    end do  
    
    do i=1, C%nbContainer 
      call IntegrateContainer(C%containerArray(i), obsInitVector, obsCBList, &
                              thicknessNoise, loadingNoise, PTGradient, &
                              PL1Gradient, PeggIntensity, BPMintensity, iBlankArray, &
                              dt, obsCount, observerTimeArray, iCont, &
                              nbCont, nbSegments, segmentSize, radDistance, PeggData, BPMData)
    end do

  end subroutine IntegrateContainer 
  
  
  subroutine BuildBroadbandData(PFarray, CBList, PeggData, BPMData, &
                                nTauMax, hasLoading, periodicKeyOffset, lastBlade, nbSegments)
    implicit none
    type(CBStructure), pointer:: CBList
    type(patchFile), dimension(:), pointer:: PFArray
    type(PEGG), pointer:: PeggData
    type(BPM), pointer:: BPMData
    real, intent(inout):: periodicKeyOffset
    real:: omega
    integer:: nTauMax, nbFiles
    integer, intent(inout):: nbSegments
    logical:: lastBlade
    logical, dimension(:):: hasLoading
    
    integer:: i, j, nTau, isCompact

    nbFiles = size(PFArray)
    do i=1,nbFiles
      !ksb debug:  The original has a problem when the geoInfo is not constant geometry.  It 
      ! deallocates the data structure that is then later used.  It seems to work OK if the data
      ! structure is not deallocated and the total thrust is "UserValue'.  This needs to be checked more
      ! throughly.  1/30/2016
      if( (PFArray(i)%geoInfo(GEO_TIMETYPE).ne.GEO_TIMETYPE_CONSTANT) .and. &
          (PeggData%computeTerm(PEGG_TOTALTHRUST).ne.'USERVALUE') ) then
        call Warning ('Broadband noise cannot be predicted correctly with non-constant geometry',&
                      'and TotalThrustFlag /= "UserValue".   Broadband noise will not be', &
                      'predicted for this container.')
      end if
      !if (PFArray(i)%geoInfo(GEO_TIMETYPE).ne.GEO_TIMETYPE_CONSTANT) then
        !call Warning ('Broadband noise cannot be predicted with non-constant geometry.',&
        !             'Broadband noise will not be predicted for this container.')           
        !ksb debug: 
        !if (associated(PeggData)) then
        !  deallocate(PeggData)
        !  nullify(PeggData)
        !end if
        !if (associated(BPMData)) then
        !  deallocate(BPMData)
        !  nullify(BPMData)
        !end if
      !end if
    end do
    
    omega = GetOmega(CBList)
    ! Broadband calculations
    if (globalPeggNoiseFlag) then
      if (associated(PeggData%nbContributions)) then
        if (size(PeggData%nbContributions).ne.nbSegments) then
          deallocate(PeggData%nbContributions)
          nullify(PeggData%nbContributions)
        end if
      end if
      allocate(PeggData%nbContributions(nbSegments))
      PeggData%nbContributions = 0
      if (associated(PeggData%Intensity)) then
        if (size(PeggData%Intensity,2).ne.nbSegments) then
          deallocate(PeggData%Intensity)
          nullify(Peggdata%intensity)
        end if
      end if
      allocate(PeggData%Intensity(BB_THIRD_OCTAVE_BANDS, nbSegments))
      PeggData%Intensity = 0.
      do j=1,PFArray(1)%geoInfo(GEO_NBZONES)
        if (hasLoading(j) .and. PeggData%timeType.eq.LOAD_TIMETYPE_PERIODIC) then
          PeggData%timeoffset = (periodicKeyOffset - PeggData%tau(1))/(PeggData%tau(size(PeggData%tau))  &
              - PeggData%tau(1))*PeggData%period
          PeggData%iTauMin  = GetPatchITauMin(PFArray(1)%patchArray(j))
          exit
        end if
      end do    
    
      if(ComputePeggTerm(PeggData, PEGG_TOTALTHRUST)) then
        nTau = 0
        do i=1,nbFiles
          if (PFarray(i)%geoAperiodic .or. PFarray(i)%loadAperiodic) then
            do j=1,PFArray(i)%geoInfo(GEO_NBZONES)
              if (hasLoading(j)) then
                nTau = max(nTau,GetNKey(PFArray(i)%patchArray(j)))
                call ComputePatchThrust(PeggData, PFArray(i)%patchArray(j), CBlist, &
                &  PFArray(i)%loadInfo(LOAD_TIMETYPE), nTau, 1, 1)
              end if
            end do
          end if
        end do
        if (PFarray(i)%geoAperiodic .or. PFarray(i)%loadAperiodic) then
          if (lastBlade .and. PFArray(i)%loadInfo(LOAD_TIMETYPE).eq.LOAD_TIMETYPE_PERIODIC) then
            PeggData%shaftAxisThrust = PeggData%shaftAxisThrust/(1.0*nTau-1)
          end if
        end if
      end if

      if (ComputePeggTerm(PeggData, PEGG_TOTALBLADEAREA)) then
        do i=1,nbFiles
          do j=1,PFArray(i)%geoInfo(GEO_NBZONES)            
            call ComputePatchArea(PeggData, PFArray(i)%patchArray(j))
          end do
        end do
        if (PeggData%totalBladeArea.eq.0.) then
          call Error (' The total blade area must be greater than zero.',&
                      ' Either there is an error in the data or the blade',&
                      ' is defined only by a compact patch which has no area.')
        end if
      end if
      
      if (ComputePeggTerm(PeggData, PEGG_bladeRadius)) then
        call CalculateResultantChangeOfBase(CBList, GetSize(CBList), 0.0)
        do i=1,nbFiles
          do j=1,PFArray(i)%geoInfo(GEO_NBZONES)
            call ComputePatchSpan(Peggdata, PFArray(i)%patchArray(j), CBList, resultantBase%position)
          end do
        end do
      end if
      
      if (ComputePeggTerm(PeggData, PEGG_ROTSPEED)) then
        PeggData%RotSpeed = omega
      end if
      
      if (PeggData%period .eq. 0.) PeggData%period = 2*pi/omega
    end if
    
    if (associated(BPMData)) then
      if (associated(BPMData%nbContributions)) then
        deallocate(BPMData%nbContributions)
        nullify(BPMData%nbContributions)
      end if
      allocate(BPMData%nbContributions(nbSegments))
      BPMData%nbContributions = 0

      if (BPMData%period .eq. 0.) BPMData%period = 2*pi/omega
      
      do j=1,PFArray(1)%geoInfo(GEO_NBZONES)
        if (hasLoading(j) .and. BPMData%timeType.eq.LOAD_TIMETYPE_PERIODIC) then
          BPMData%timeoffset = (periodicKeyOffset - BPMData%tau(1))/(BPMData%tau(size(BPMData%tau))  &
             - BPMData%tau(1))*BPMData%period
          exit
        end if
      end do

      isCompact = 0
      do i=1,nbFiles
        do j=1,PFArray(i)%geoInfo(GEO_NBZONES)
          if (BPMData%computeTerm(BPM_SECTLENGTH).eq.'COMPUTE') then
            call CheckPatchNSect(PFArray(i)%patchArray(j), BPMData, isCompact)
            if (isCompact.eq.1) then
              call ComputePatchSectLength(BPMData, PFArray(i)%patchArray(j))
            elseif (isCompact.gt.1) then
              call Error(' There must be exactly 1 compact patch per container when using ',&
                         ' the BPM method to calculate broadband noise.')
              exit
            end if
          end if
        end do
        if (isCompact.gt.1) exit
      end do
    end if
  end subroutine BuildBroadbandData
  
  
  subroutine BuildBroadbandObsTimeArray(BBobsTimeArray, obsTimeArray)
    implicit none
    real, dimension(:), pointer:: BBobsTimeArray
    real(kind=8), dimension(:):: obsTimeArray
    integer:: i, nbObsTimes, nbBBObsTimes
    nbObsTimes   = size(obsTimeArray)
    nbBBObsTimes = size(BBobsTimeArray)
    do i=1,nbBBObsTimes
      BBobsTimeArray(i) = obsTimeArray(1) + (i-1)*(obsTimeArray(nbObsTimes)-obsTimeArray(1))/(1.0*nbBBObsTimes)
    end do
  end subroutine BuildBroadbandObsTimeArray
  
  recursive subroutine ResetContainerKeyCount(C)
    type(container):: C
    integer:: i, j
    
    do i=1, C%nbFiles
      do j=1,C%PFArray(i)%geoInfo(GEO_NBZONES)
        call ResetPatchKeyCount(C%PFArray(i)%patchArray(j))
      end do
    end do
   
    do i=1,C%nbContainer
      call ResetContainerKeyCount(C%containerArray(i))
    end do
  end subroutine ResetContainerKeyCount
  
  subroutine ResetPatchFileKeyCount(PF,geoKeyVal,loadKeyVal)
    type(patchFile):: PF
    integer, intent(in), optional:: geoKeyVal, loadKeyVal
    integer:: j
    
    do j=1,PF%geoInfo(GEO_NBZONES)
      if( present(geoKeyVal) .and. present(loadKeyVal) ) then
        call ResetPatchKeyCount(PF%patchArray(j), geoKeyVal, loadKeyVal)
      else
        call ResetPatchKeyCount(PF%patchArray(j))
      end if
    end do

  end subroutine ResetPatchFileKeyCount
  
  subroutine SetPatchFileNKey(PF,nKey,geoNKey,loadNKey)
    type(patchFile):: PF
    integer, intent(in):: nKey
    integer, intent(in), optional:: geoNKey, loadNKey
    integer:: j
    
    do j=1,PF%geoInfo(GEO_NBZONES)
      if( present(geoNKey) .and. present(loadNKey) ) then
        call SetPatchNKey(PF%patchArray(j), nKey, geoNKey, loadNKey)
      else
        call SetPatchNKey(PF%patchArray(j),nKey)
      end if
    end do

  end subroutine SetPatchFileNKey
  
  recursive subroutine CheckBBFlags(C, computePegg, computeBPM)
    implicit none
    type(container), intent(inout):: C
    logical:: computePegg, computeBPM
    integer:: i
    
    if (C%PeggNoiseFlag) computePegg = .true.
    if (C%BPMNoiseFlag)  computeBPM  = .true.
    if (computePegg .and. computeBPM) return
    
    do i=1,C%nbContainer
      call CheckBBFlags(C%containerArray(i), computePegg, computeBPM)
    end do
  end subroutine CheckBBFlags
  
  
  recursive subroutine CreateContainerObsTimeArrays(C, nbObs)
    type(Container), intent(inout)::C
    integer:: i,j, nbObs
 
    do i=1,C%nbFiles
      do j=1,C%PFArray(i)%geoInfo(GEO_NBZONES)
        call CreatePatchObsTimeArrays(C%PFArray(i)%patchArray(j), nbObs)
      end do
    end do
 
    do i=1,C%nbContainer
      call CreateContainerObsTimeArrays(C%containerArray(i), nbObs)
    end do
  end subroutine CreateContainerObsTimeArrays
 
  !
  !SUBROUTINE GetNbContainerSurfaces(nbPatches,containr)
  !    This routine return the number of patches that are contained in this
  !  container as well as all the patches of the children.
  !ARGUMENTS:
  !  - nbPatches: This is the returned value of the number of patches
  !  - containr: The container we are counting the patches on.
  !*
  recursive subroutine GetNbContainerSurfaces(C,nbStructuredPatches,&
                                              nbUnstructuredPatches)
    type(Container), intent(in)::C
    integer, intent(inout)::nbStructuredPatches, nbUnstructuredPatches

    integer::i,j

    ! Loop over the patches and figure out if they are structured or
    ! unstructured:
    do i=1,C%nbFiles
      do j = 1,size(C%PFArray(i)%patchArray)
        if (IsStructured(C%PFArray(i)%patchArray(j))) then
          nbStructuredPatches = nbStructuredPatches + 1
        else
          nbUnstructuredPatches = nbUnstructuredPatches + 1
        end if
      end do
    end do
   
    ! Now loop over the sub-containers
    do i=1,C%nbContainer
       call GetNbContainerSurfaces(C%containerArray(i), &
                                   nbStructuredPatches, &
                                   nbUnstructuredPatches)
    end do

  end subroutine GetNbContainerSurfaces


  recursive subroutine GetContainerStructuredDims (C, dimensionsList, &
                                                   counter)
    type(container)::C
    integer, dimension(:,:)::dimensionsList
    integer::counter, i,j

    do i=1,C%nbFiles
      do j = 1,size(C%PFArray(i)%patchArray)
        if (IsStructured(C%PFArray(i)%patchArray(j))) then
          counter = counter + 1
          call GetPatchDimensions(C%PFArray(i)%patchArray(j), &
                                  dimensionsList(counter,:))
        end if
      end do
    end do    
    do i=1, C%nbContainer
       call GetContainerStructuredDims (C%containerArray(i), &
                                        dimensionsList, counter)
    end do

  end subroutine GetContainerStructuredDims 


  recursive subroutine GetContainerImplicitTimeRange (C, obsTMin, obsTMax, &
                                                      obsCoordinates, &
                                                      obsCBList, &
                                                      obsTimes)
    type(Container), intent(inout)::C
    type(vector), dimension(:,:), intent(in)::obsCoordinates
    type(CBStructure), pointer::obsCBList
    real, intent(inout)::obsTMin, obsTMax
    real, dimension(:,:,:), intent(inout) :: obsTimes
    integer :: i,j
    if (debugLevel >= 2.and.IsMaster()) then
      write (*,'(A,A,A,$)') C%prefix, trim(C%title), ' patch '
    end if

    do i=1,C%nbFiles
      if (debugLevel >= 2.and.IsMaster()) then
        write (*,'("(",A,")",$)') trim(IntegerToString (i))
      end if
      do j = 1,size(C%PFArray(i)%patchArray)
        if (debugLevel >= 2.and.IsMaster()) then
          write (*,'("[",A,"]",$)') trim(IntegerToString (j))
        end if
        call GetPatchImplicitTimeRange (C%PFArray(i)%patchArray(j), C%CBList, C%tauMin, C%tauMax, &
                                        C%dTau, C%nTau, obsCoordinates, obsCBList, obsTimes)
      end do
    end do
    if(IsMaster().and.debugLevel >= 2) write(*,*)
    do i=1, C%nbContainer
      call GetContainerImplicitTimeRange (C%containerArray(i), obsTMin, obsTMax, &
                                          obsCoordinates, obsCBList, obsTimes)
    end do
    
  end subroutine GetContainerImplicitTimeRange


  recursive function getContainerMinObsTimeInData(C)  result(minObsTimeInData)
    type(container), intent(in)::C
    real(kind=8):: minObsTimeInData
    integer:: i

    minObsTimeInData = -huge(minObsTimeInData)
    minObsTimeInData = max(minObsTimeInData,C%minObsTimeInData)
    
    do i=1, C%nbContainer
      minObsTimeInData = max(minObsTimeInData, getContainerMinObsTimeInData (C%containerArray(i)))
    end do
  end function getContainerMinObsTimeInData
  
  recursive function getContainerMaxObsTimeInData(C)  result(maxObsTimeInData)
    type(container), intent(in)::C
    real(kind=8):: maxObsTimeInData
    integer:: i
    
    maxObsTimeInData = huge(maxObsTimeInData)
    maxObsTimeInData = min(maxObsTimeInData,C%maxObsTimeInData)
    
    do i=1, C%nbContainer
      maxObsTimeInData = min(maxObsTimeInData, getContainerMaxObsTimeInData (C%containerArray(i)))
    end do
  end function getContainerMaxObsTimeInData

   
  recursive subroutine WriteContainerStructuredSigma (C, gridUnit, &
                                                      dataUnit, prefix)
    type(container), intent(in)::C
    integer, intent(in) :: gridUnit, dataUnit
    character(len=*), intent(in) :: prefix

    integer :: i,j
    
    call Message (prefix//trim(C%title))

    do i=1,C%nbFiles
      do j = 1,size(C%PFArray(i)%patchArray)
        call WritePatchPlot3DStructuredSigma(C%PFArray(i)%patchArray(j), gridUnit, dataUnit, &
                                              prefix//"  + ")
      end do
    end do   
    do i=1, C%nbContainer      
      call WriteContainerStructuredSigma(C%containerArray(i),gridUnit,dataUnit, &
                                         prefix//"  ")
    end do
  end subroutine WriteContainerStructuredSigma


  recursive subroutine WriteContainerFVUNSHeader (C, dataUnit)
    type(container), intent(in)::C
    integer, intent(in) :: dataUnit

    integer :: i,j

    do i=1,C%nbFiles
      do j = 1,size(C%PFArray(i)%patchArray)
        call WritePatchFVUNSHeader (C%PFArray(i)%patchArray(j), dataUnit, C%title)
      end do
    end do
    do i=1, C%nbContainer
      call WriteContainerFVUNSHeader (C%containerArray(i), dataUnit)
    end do

  end subroutine WriteContainerFVUNSHeader

  
  recursive subroutine WriteContainerFVUNSData (C, dataUnit, patchNumber,prefix)
    type(container), intent(in)::C
    integer, intent(inout) :: dataUnit, patchNumber
    character(len=*), intent(in) :: prefix

    integer :: i,j

    call Message (prefix//trim(C%title))

    do i=1,C%nbFiles
      do j = 1,size(C%PFArray(i)%patchArray)
        call WritePatchFVUNSData (C%PFArray(i)%patchArray(j), dataUnit, patchNumber, &
                                  prefix//"  + ")
      end do
    end do   
    do i=1, C%nbContainer
      call WriteContainerFVUNSData (C%containerArray(i), dataUnit, &
                                    patchNumber, prefix//"  ")
    end do
  end subroutine WriteContainerFVUNSData

  recursive subroutine CloseDataFiles(C)
    type(container)::C
    integer:: i, loading = 2
    type(Pegg), pointer:: PeggData=>null()
    type(BPM), pointer:: BPMData=>null()
    
    ! Close the files:
    do i=1,C%nbFiles
      if (C%PFarray(i)%geoAperiodic .or. C%PFarray(i)%loadAperiodic) then          
        call CloseBinaryFile (C%PFarray(i)%geometryStreamNumber)
      end if
      if (C%PFarray(i)%loadAperiodic) then
        call CloseBinaryFile (C%PFarray(i)%surfaceDataStreamNumber)  
      end if
    end do

    if (associated(C%BBdata)) then
      PeggData => GetPeggData(C%BBdata)
      if (associated(PeggData)) then
        if (getPeggTimeType(PeggData).eq.GEO_TIMETYPE_APERIODIC) then
          call CloseBinaryFile (C%PeggUnitNumber)
        end if
      end if
      nullify(PeggData)
      BPMData => GetBPMData(C%BBdata)
      if (associated(BPMData)) then
        if (getBPMTimeType(BPMData).eq.GEO_TIMETYPE_APERIODIC) then
          call CloseBinaryFile (C%BPMUnitNumber)
        end if
      end if
      nullify(BPMData)
    end if
    do i=1, C%nbContainer
      call CloseDataFiles(C%containerArray(i))
    end do
  end subroutine CloseDataFiles
    
  function AperiodicCont(C) result(aperiodic)
    type(container), intent(in)::C
    logical:: aperiodic
    
    aperiodic = C%aperiodicData
    
  end function AperiodicCont
  
  subroutine PatchFileSetup(PF)
    type(PatchFile),intent(inout) :: PF
    integer:: j, stat, startpos
    
    ! If the data is aperiodic, repostion to the start of the data section
    ! of the files.  For other cases (constant and periodic), all the data
    ! was read in when the patches were created.
    
    if (PF%geoAperiodic) then
      if( PF%MTFgeoAperiodic ) then
        call CloseBinaryFile(PF%GeometryStreamNumber)
        PF%MTFPrevGeoFileNKey=0
        PF%iMTFGeoFile=1
        PF%patchGeometryFile=PF%MTFGeoFileArray(PF%iMTFGeoFile)
        call OpenBinaryFile (PF%GeometryStreamNumber, PF%patchGeometryFile, PF%GeoBigEndian, stat)
        startpos = PF%startGeometryData + &
                  (PF%MTFGeoStartKeyIndex(PF%iMTFGeoFile)-1)*PF%geomStepInBytes
      else
        startpos = PF%startGeometryData
      end if
      call RestartBinaryFileAtByte(PF%GeometryStreamNumber, startpos,stat)
      if(stat /= 0 ) then
        call Error( "Problem restarting "//trim(PF%patchGeometryFile) )
      end if
    end if
    if (PF%loadAperiodic) then
      if( PF%MTFLoadAperiodic ) then
        call CloseBinaryFile(PF%SurfaceDataStreamNumber)
        PF%MTFPrevLoadFileNKey= 0
        PF%iMTFLoadFile=1
        PF%patchLoadingFile=PF%MTFLoadFileArray(PF%iMTFLoadFile)
        call OpenBinaryFile (PF%SurfaceDataStreamNumber, PF%patchLoadingFile, PF%LoadBigEndian, stat)
        startpos = PF%startSurfaceData + &
                   (PF%MTFLoadStartKeyIndex(PF%iMTFLoadFile)-1)*PF%loadStepInBytes
      else
        startpos = PF%startSurfaceData
      end if
      call RestartBinaryFileAtByte(PF%surfaceDataStreamNumber,startpos,stat)
      if(stat /= 0 ) then
        call Error( "Problem restarting "//trim(PF%patchLoadingFile) )
      end if
    end if
    !!ksb debug:  I really should combine these with above, plus I need to reset everything for MTF Aperiodic files
    do j=1,PF%geoInfo(GEO_NBZONES)
      call ResetPatchKeyCount( PF%PatchArray(j) )  
    end do
  end subroutine PatchFileSetup
end module containerObject    
