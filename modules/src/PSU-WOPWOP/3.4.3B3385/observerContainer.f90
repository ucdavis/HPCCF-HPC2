module ObserverContainerObject

  use TimeHistoryObject
  use FrequencyDomainObject
  use ObserverUtility
  use ConstantsModule
  use ObserverObject
  use MathModule
  
  implicit none  
           
  type observerContainer
    logical, dimension(:), pointer::newObs=>NULL()
    character(len=4096)::title,attachedTo
    !This is the number of observers the container contains
    integer::nbObs
    !This is the array of observers within the container
    type(observer),dimension(:),pointer::obs
    !These parameters are for segmenting a time history.  They are applied
    !  to every timeHistory this observer container contains:
    real::tMin,tMax, dt
    integer::nt
    real::segmentSize,segmentStepSize
    ! nbNewSegments stores the number of segments to add to the observer time range
    ! with each observer time expansion. The total time expansion is thus
    ! nbNewSegments*segmentSize
    integer:: nbObsTimeExpansions, maxObsTimeExpansions, nbNewSegments
    ! nbSegments stores the number of segments in the current obs. time range
    integer::nbSegments,ntPerSegment
    real,dimension(:),pointer::segTimeArray=>null()
    !These parameters are for windowing a time history. They are applied       
    !  to every timeHistory this observer contains:
    character(len=4096)::windowFunction
    !These parameters set up which overall range of frequencies the observer
    !  hears.  highPassFrequency and lowPassFrequency are applied to each spectrum
    !  created for the observer.
    real::highPassFrequency,lowPassFrequency,octaveNumber
    integer::nbfreqRanges,segmentIncrement,nbHarmonics
    logical::octaveApproxFlag
    character(len=4096),dimension(:),pointer::fRangeTitle=>null()
    real,dimension(:),pointer::fRangeMin=>null(),fRangeMax=>null()
    !Observer motion
    type(CBStructure),pointer::CBList=>null(),CBListLocalFrame=>null()
    integer::NbObsDim1,NbObsDim2
    real::sigmaRadius
    character(len=4096)::FSCInputFileName,FSCOutputFileName
    character(len=4096)::SPLPath,pressurePath,phasePath,OASPLPath,audioPath,sigmaPath, &
                         complexPressurePath
  end type observerContainer 
  
  private::CreateNewObsArray
  public::observerContainer, InitializeObserverContainer, CalcObsContResults,    &
          DestroyObserverContainer, GetNbObservers, EnlargeObsContainer,         &
          WriteObsContSigmaSurface, ReceiveResultsFromSlave,                     &
          SendResultsToMaster, ResetSlaveObserver,  ReceiveObserverPoint, CopyTH,&
          SetupParallelObserverContainer, SendObserverPoint, InitializeObservers,&
          WriteObsContResults, WriteObserverContainerDebugInfo, CopyObserverData,&
          CheckSegmentTimeSettings, WriteSegmentDataToScreen, CreateNewObservers,&
          CheckFreqAnalysisParameters  
contains

subroutine ResetSlaveObserver(obsCont)
  type(observerContainer), intent(inout)::obsCont
  obsCont%nbObs = 1
  obsCont%nbObsDim1 = 1
  obsCont%nbObsDim2 = 1
  call ResetObserver(obsCont%obs(1), obsCont%tMin, obsCont%tMax, obsCont%nt, obsCont%dt)
end subroutine ResetSlaveObserver

subroutine ReceiveObserverPoint(obsCont, sender, master)
  type(observerContainer), intent(inout)::obsCont
  integer, intent(in)::master, sender
  call ReceiveVector(obsCont%obs(1)%Coordinates, sender, master)  
end subroutine ReceiveObserverPoint

subroutine SendObserverPoint(obsCont, pointIndex, proc, master) 
  type(observerContainer), intent(in)::obsCont
  integer, intent(in)::pointIndex, proc, master
  call SendVector(obsCont%obs(pointIndex)%coordinates, proc, master)
end subroutine SendObserverPoint

subroutine SetupParallelObserverContainer(obsCont)
  type(observerContainer), intent(inout)::obsCont
  if (.not.IsMaster()) then
    call ResetSlaveObserver(obsCont)
  end if
end subroutine SetupParallelObserverContainer

subroutine ReceiveResultsFromSlave(obsCont, master, sender, rcvPointIndex)
  type(observerContainer), intent(inout)::obsCont
  integer, intent(in)::master, sender, rcvPointIndex
  call ReceiveObserverResultsFromSlave(obsCont%obs(rcvPointIndex), master, sender)
end subroutine ReceiveResultsFromSlave

subroutine SendResultsToMaster(obsCont, master, tag)
  type(observerContainer), intent(inout)::obsCont
  integer, intent(in)::master, tag
  call SendObserverResultsToMaster(obsCont%obs(1), master, tag) 
end subroutine SendResultsToMaster

subroutine WriteObserverContainerDebugInfo(obsCont, unitnum)
  type(observerContainer), intent(in)::obsCont
  integer::unitnum
  write(unitnum,*) '*** ObserverContainerDebugInfo ***' 
  write(unitnum,*) 'Title= ', trim(obsCont%title)
  write(unitnum,*) 'attachedTo= ', trim(obsCont%attachedTo)
  write(unitnum,*) 'windowFunction= ', trim(obsCont%windowFunction)
  write(unitnum,*) 'tmin= ', obsCont%tMin
  write(unitnum,*) 'tmax= ', obsCont%tMax
  write(unitnum,*) 'nt= ', obsCont%nt
  write(unitnum,*) 'segmentSize= ', obsCont%segmentSize
  write(unitnum,*) 'segmentStepSize= ', obsCont%segmentStepSize
  write(unitnum,*) 'NbSegments= ', obsCont%nbSegments
  write(unitnum,*) 'highPassPrequency= ', obsCont%highPassFrequency
  write(unitnum,*) 'lowPassPrequency= ', obsCont%lowPassFrequency
  write(unitnum,*) 'octaveNumber ', obsCont%octaveNumber
  write(unitnum,*) 'segmentIncrement ', obsCont%segmentIncrement
  write(unitnum,*) 'nbHarmonics ', obsCont%nbHarmonics
  write(unitnum,*) 'octaveApproxFlag ', obsCont%octaveApproxFlag
  write(unitnum,*) 'sigmaRadius ', obsCont%sigmaRadius
  write(unitnum,*) 'FSCInputFileName ', trim(obsCont%FSCInputFileName)
  write(unitnum,*) 'FSCOuputFileName ', trim(obsCont%FSCOutputFileName)
  write(unitnum,*) 'SPLPath ',  trim(obsCont%SPLPath)
  write(unitnum,*) 'pressurePath ',  trim(obsCont%pressurePath)
  write(unitnum,*) 'phasePath ',  trim(obsCont%phasePath)
  write(unitnum,*) 'OASPLPath ',  trim(obsCont%OASPLPath)
  write(unitnum,*) 'audioPath ',  trim(obsCont%audioPath)
  write(unitnum,*) 'sigmaPath ',  trim(obsCont%sigmaPath)
  write(unitnum,*) 'complexPressurePath ',  trim(obsCont%complexPressurePath)
  if (associated(obsCont%CBList)) then
    call WriteCBListDebugInfo(obsCont%CBList, unitnum)
  end if
  if (associated(obsCont%CBListLocalFrame)) then
    call WriteCBListDebugInfo(obsCont%CBListLocalFrame, unitnum)
  end if
  write(unitnum,*) 'segTimeArray allocated? ', associated(obsCont%segTimeArray)
  if (associated(obsCont%segTimeArray)) then
    write(unitnum,*) 'Size of segTimeArray array is ', trim(integertostring(size(obsCont%segTimeArray)))
  end if
  write(unitnum,*) 'fRangeTitle allocated? ', associated(obsCont%fRangeTitle)
  if (associated(obsCont%fRangeTitle)) then
    write(unitnum,*) 'Size of fRangeTitle array is ', trim(integertostring(size(obsCont%fRangeTitle)))
  end if
  write(unitnum,*) 'fRangeMin allocated? ', associated(obsCont%fRangeMin)
  if (associated(obsCont%fRangeMin)) then
    write(unitnum,*) 'Size of fRangeMin array is ', trim(integertostring(size(obsCont%fRangeMin)))
  end if
  write(unitnum,*) 'fRangeMax allocated? ', associated(obsCont%fRangeMax)
  if (associated(obsCont%fRangeMax)) then
    write(unitnum,*) 'Size of fRangeMax array is ', trim(integertostring(size(obsCont%fRangeMax)))
  end if
  write(unitnum,*) '--- End observer container Debug Info ---'
  call Barrier()
end subroutine WriteObserverContainerDebugInfo


!**
!SUBROUTINE EnlargeObsContainer(obs)
!ARGUMENTS
! - obs:
!CALLS:
! - createSphericalSectionObserver  (observer.f90)
!**
subroutine EnlargeObsContainer(obsCont)
  type(observerContainer), intent(inout)::obsCont
  integer::i
  real, dimension(3)::oldCoordinates
  integer:: nt
  real:: tMin,tMax
  oldCoordinates=obsCont%obs(1)%coordinates%A(:)
  nt = obsCont%obs(1)%nt
  tMin = obsCont%obs(1)%tMin
  tMax = obsCont%obs(1)%tMax
  obsCont%NbObsDim1 = 21
  obsCont%NbObsDim2 = 21
  obsCont%nbObs     = 21**2    
  if(associated(obsCont%obs))then
    deallocate(obsCont%obs)
    nullify(obsCont%obs)
  end if
  allocate(obsCont%obs(obsCont%nbObs))
  call createSphericalSectionObsCont(obsCont,obsCont%NbObsDim1,obsCont%NbObsDim2, &
                                     obsCont%sigmaRadius,0.0, 2*PI, -PI, PI, 0, 0, .false.)
  do i=1, obsCont%nbObs
    call NullifyObserver(obsCont%obs(i))
    obsCont%obs(i)%nt = nt
    obsCont%obs(i)%tMin = tMin
    obsCont%obs(i)%tMax = tMax
    obsCont%obs(i)%coordinates=obsCont%obs(i)%coordinates+vectorSetValue(oldCoordinates)
  end do
end subroutine EnlargeObsContainer

!**
!SUBROUTINE CalcObsContResults(obsCont)
!  This routine first calculates total noise if requested in the namelist.
!  It then writes acoustic pressure files if acousticPressureFlag=.true.
!  If required, it windows the time histories and then calculates the requested
!  SPL data.  It performs these operations by calling subroutines which operate
!  on each observer point in the observer container.
!**
subroutine CalcObsContResults(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  if (IsMaster()) then
    write(*,*) 'Calculating results for '//trim(obsCont%title)
  end if
  !This routine calculates the total noise if it was specified
  !  in the namelist.  Also, if total noise was specified but thickness
  !  and loading were not, this destroys the thickness and loading noise:
  if(.not.ReadInPressureFlag)then
    call CalculateObsContTotalNoise(obsCont)
  end if 
  if (newSegLHS .or. newSegRHS) then
    ! Recombine our small observer segments with the 'original' observers
    ! before we redo our calculation of the results.
    call CreateNewObservers(obsCont, .false.)
    
    call CheckSegmentTimeSettings(obsCont)
    ! Reset the 'newSegment' logicals since the rest of the run
    ! will deal with the 'recombined' observer.
    newSegLHS = .false.
    newSegRHS = .false.    
  end if

  !This routine nullfies and allocates the pointers obsCont%obs(:)%pPrimeSegs%f
  !  and obsCont%obs(:)%pPrimeSegs%t, and obsCont%obs(:)%pPrimeSegs
  call InitializePPrimeSegments(obsCont)  
  if(spectrumFlag.or.OASPLdBFlag.or.OASPLdBAFlag.or.SPLdBFlag.or.SPLdBAFlag.or.&
     octaveFlag.or.PNLFlag.or.PNLTFlag.or.EPNLFlag.or.phaseFlag.or.complexPressureFlag.or. &
     FSCOutputFlag.or.SELFlag)then
    call SegmentAndWindowPPrime(obsCont)
    !Now calculate all the SPL data:
    call CalcObsContFreqData(obsCont) 
  end if

  if(audioFlag) call CreateObsWavFiles(obsCont%obs(1),obsCont%audioPath)

end subroutine CalcObsContResults

!**
!SUBROUTINE SegmentAndWindowPPrime
!  This routine calculates all the obsCont%obs(:)%pPrimeSegs data 
!  from obsCont%obs(:)%pPrime, and then windows them if necessary
!**
subroutine SegmentAndWindowPPrime(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  if(debugLevel.gt.5) then
    write(*,*) 'Filling pPrimeSegs Array.'
  end if
  call CreatePPrimeSegments(obsCont)
  !We must window all the time histories:
  if(trim(obsCont%windowFunction).ne.'')then
    if(debugLevel.gt.5) then
      write(*,*) 'Windowing pPrimeSegs Arrays.'
    end if
    call WindowPPrimeSegments(obsCont)
  else
    call Warning('No weighting window was applied to the time' , &
                 'histories before calculating spectrum data.')
  end if  
end subroutine SegmentAndWindowPPrime

!**
!SUBROUTINE CalculateObsContTotalNoise
!  This routine calculates all the total noise contained in obsCont%obs(:)%pPrime
!**
subroutine CalculateObsContTotalNoise(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  integer:: i
  if(totalNoiseFlag)then
    if(debugLevel.gt.5) then
      write(*,*) 'Calculating Total Noise'
    end if
    do i=1,obsCont%nbObs
      if (obsCont%newObs(i)) call CalcTotalNoise(obsCont%obs(i))
    end do
    !  If either the thickness or loading noise flags are off,
    !  now we can discard this pressure data.  However, 
    !  we still need to keep the pPrime array in observer
    !  indexed from 1,2,3,... with no gaps in between.
    !  The following routine resets the obs%pPrime array
    !  and resets the indices LOAD_APTH,THICK_APTH, and so on.
    if((.not.thicknessNoiseFlag).or.(.not.loadingNoiseFlag) )then
      if(debugLevel.gt.5) then
        write(*,*) 'Reseting p prime indices'
      end if
      call ResetPPrimeAndIndices(obsCont%obs)          
    end if  
  end if
end subroutine CalculateObsContTotalNoise

!**
!SUBROUTINE InitializePPrimeSegments
!  This routine nullfies and allocates the pointers obsCont%obs(:)%pPrimeSegs%f
!  and obsCont%obs(:)%pPrimeSegs%t, and obsCont%obs(:)%pPrimeSegs
!  If a segmentIncrement is specified in the namelist, or there is only one 
!  segment to be analyzed, it only allocates pPrimeSegs(nbSources,1)
!  If there are multiple segments and no segmentIncrement is specified,
!  it allocates pPrimeSegs(nbSources,nbSegments)
!**
subroutine InitializePPrimeSegments(obsCont)
  implicit none
  integer::i,j,k
  type(observerContainer),intent(inout)::obsCont
  if(debugLevel.gt.5) then
    write(*,*) 'Allocating pPrimeSegs array.'
  end if
  if(obsCont%nbSegments == 0)then
    obsCont%nbSegments = 1
  end if
  do i=1,obsCont%nbObs
    if((obsCont%nbSegments==1).or.(obsCont%segmentIncrement.ne.0))then
      if (.not.associated(obsCont%obs(i)%pPrimeSegs)) then
        allocate(obsCont%obs(i)%pPrimeSegs(GetNewSizeOfPPrime(),1))
        do j=1,GetNewSizeOfPPrime()
          nullify(obsCont%obs(i)%pPrimeSegs(j,1)%t, obsCont%obs(i)%pPrimeSegs(j,1)%f)
        end do
      end if
    else
      nullify(obsCont%obs(i)%pPrimeSegs)
!      if (.not.associated(obsCont%obs(i)%pPrimeSegs)) then
      allocate(obsCont%obs(i)%pPrimeSegs(GetNewSizeOfPPrime(),obsCont%nbSegments))
      do j=1,GetNewSizeOfPPrime()
        do k=1,obsCont%nbSegments
          nullify(obsCont%obs(i)%pPrimeSegs(j,k)%t , obsCont%obs(i)%pPrimeSegs(j,k)%f)
        end do
      end do
!      end if
    end if
  end do
end subroutine InitializePPrimeSegments  

!**
!SUBROUTINE FilterObsContPressure
!  This filters all the pressure in obsCont%obs(:)%pPrime if necessary
!**
subroutine FilterObsContPressure(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  integer::i,j
  real::minFreq,maxFreq
  if(acousticPressureFlag.or.pressureGradientFlag.or.pressureGradient1AFlag)then  
    !We have to filter the entire time histories and write them to file.
    !  This is done by the following:
    !  1.) FFTW each time history (get complexDFTResult)
    !  2.) Filter it according to high and low pass filter
    !  3.) iFFTW back to pressure time history
    minFreq = 1./(obsCont%tMax-obsCont%tMin)
    maxFreq = minFreq*(obsCont%nt/2+1)
    !Filter the pressure (if necessary) according to the 
    ! high and low pass filters specified in the namelist.
    !
    ! KSB change - if lowPassFrequency is < minFrequency, the filters is not applied.  This
    ! is a problem if you just want to eliminate the mean value of the pressure that the
    ! FW-H equation gives you (i.e., to mimic a microphone more accurately)
    !
    !if((minFreq.lt.obsCont%highPassFrequency).or.(maxFreq.gt.obsCont%lowPassFrequency))then
    if((obsCont%highPassFrequency.gt.0.).or.(maxFreq.gt.obsCont%lowPassFrequency))then

      !This makes sure that nbHarmonics, if specified, does not affect pressure output:
      if(obsCont%nbHarmonics.eq.0)then
        if(debugLevel.gt.5) then
          write(*,*) 'Filtering p prime with high and low pass filter frequencies.' 
        end if
        do i=1,obsCont%nbObs
          do j=1,GetNewSizeOfPPrime()
            call FilterPPrime(obsCont%obs(i)%pPrime(j), obsCont%highPassFrequency, &
                                                     obsCont%lowPassFrequency)
          end do
        end do
      end if  
    end if
  end if
end subroutine FilterObsContPressure  

!**
!SUBROUTINE WindowPPrimeSegments
!  This routine windows all the time histories contained 
!  in obsCont%obs(:)%pPrimeSegs
!**
subroutine WindowPPrimeSegments(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  integer::i,j,k
  do i=1,obsCont%nbObs
    do j=1,GetNewSizeOfPPrime()
      do k=1, obsCont%nbSegments
        call WindowTimeHistory(obsCont%obs(i)%pPrimeSegs(j,k), &
                               obsCont%windowFunction)
      end do
    end do
  end do
end subroutine WindowPPrimeSegments

!**
!SUBROUTINE CreatePPrimeSegments(obsCont,pPrime)
!ARGUMENTS:
! - obsCont: the observer container 
!**
subroutine CreatePPrimeSegments(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  integer::i,j,k
  if(obsCont%nbSegments==1)then
    do i=1,obsCont%nbObs
      do j=1,GetNewSizeOfPPrime()
        call CopyTimeHistory(obsCont%obs(i)%pPrime(j),      &
                             obsCont%obs(i)%pPrimeSegs(j,1))
      end do
    end do
  else if (obsCont%segmentIncrement.gt.0)then
    do i=1, obsCont%nbObs
      do j=1, GetNewSizeOfPPrime()
        call CreateTimeHistorySegment(obsCont%obs(i)%pPrime(j), &
                                      obsCont%segmentSize,      &
                                      obsCont%segmentStepSize,  &
                                      obsCont%segmentIncrement, &
                                      obsCont%obs(i)%pPrimeSegs(j,1))
      end do
    end do
  else 
    do i=1, obsCont%nbObs
      do j=1, GetNewSizeOfPPrime()
        do k=1, obsCont%nbSegments
          call CreateTimeHistorySegment(obsCont%obs(i)%pPrime(j),&
                                        obsCont%segmentSize,     &
                                        obsCont%segmentStepSize, &
                                        k, obsCont%obs(i)%pPrimeSegs(j,k))
        end do 
      end do
    end do  
  end if
end subroutine CreatePPrimeSegments

!**
!SUBROUTINE CalcObsContFreqData(obsCont)
!  This routine calculates and writes all the SPL data requested. It creates
!  temporary arrays of type freqDomain in order to write the data but not store
!  them in the observer container or individual observers.
!**
subroutine CalcObsContFreqData(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  integer::i,j,k
  type(observer), pointer:: obs !ksb debug
  if(debugLevel.gt.5) then
    write(*,*) 'Calculating Observer Container Frequency Dependent Data.'
  end if

  !  complexPressure is the first thing we need to calculate,
  !  This routine is the top level routine that calculates frequency domain
  !  data from the pressure time histories. Everything we need after this
  !  subroutine uses the complex pressure data.  
  if (debugLevel.gt.11) then
    write(*,*) 'Calculating complex pressure'
  end if
  call CalcObsContComplexPressure(obsCont)

  !  At this point we have the complexPressure for each observer, noise source, and segment.
  !  Now we can calculate both the mean-squared pressure (MSP) for SPL
  !  calculations, and the phase information (if necessary).
  !  This is where we calculate the phase data from the complexPressure data:
  if(phaseFlag)then
    if(debugLevel.gt.1.and.IsMaster())then
      write(*,*)'Calculating phase data for '//trim(obsCont%title)
    end if
    do i=1,obsCont%nbObs
      call CalcObsPhaseFromComplexPressure(obsCont%obs(i))
    end do
  end if

  !This is where we calculate the MSP data from the complexPressure data:
  if(OASPLdBFlag.or.OASPLdBAFlag.or. &
     SPLdBFlag.or.SPLdBAFlag.or.PNLFlag.or.EPNLFlag.or.&
     octaveFlag.or.PNLTFlag.or.SELFlag)then
    if(debugLevel.gt.11) then
      write(*,*) 'Calculating mean square pressure'
    end if
    do i=1,obsCont%nbObs
      call CalcObsMSPFromComplexPressure(obsCont%obs(i))
    end do
  end if

  ! Apply atmospheric attenuation to the discrete frequency MSP
  if (AtmAtten) then
    do i=1,obsCont%nbObs
      call AtmosAtten(obsCont%obs(i)%msp, obsCont%obs(i)%radDistance)
    end do
  end if

  !Here is where we convert the MSP to dB and/or dBA (narrow band spectrum):
  if(spectrumFlag.and.(SPLdBFlag.or.SPLdBAFlag))then
    if(SPLdBFlag)then
      if(debugLevel.gt.1.and.IsMaster())then
        write(*,*)'Calculating SPL (dB) spectral data for '//trim(obsCont%title)
      end if
      do i=1,obsCont%nbObs
        call CalcObsdBFromMSP(obsCont%obs(i)%dB,obsCont%obs(i)%msp)
      end do
    end if     
    if(SPLdBAFlag)then
      if(debugLevel.gt.1.and.IsMaster())then
        write(*,*)'Calculating SPL (dBA) spectral data for '//trim(obsCont%title)
      end if
      do i=1,obsCont%nbObs
        call CalcObsdBAFromMSP(obsCont%obs(i)%dBA,obsCont%obs(i)%msp)
      end do
    end if    
  end if
   
  if((octaveFlag.and.(SPLdBFlag.or.SPLdBAFlag)).or.PNLFlag.or.PNLTFlag.or.EPNLFlag.or.SELFlag)then
    do i=1,obsCont%nbObs
      if ((globalPeggNoiseFlag.or.globalBPMNoiseFlag) .and. AtmAtten) then    
        ! Apply AtmAtten to the broadband spectrum before the MSP results are put through
        ! the octave filter so that BB results can be properly added (with or without AtmAtten)
        ! to the total noise octaveFiltMSP array.
        allocate(obsCont%obs(i)%mspBBAtmAtten(size(obsCont%obs(i)%mspBB,1), size(obsCont%obs(i)%mspBB,2)))
        do j=1,size(obsCont%obs(i)%mspBB,1)
          do k=1,size(obsCont%obs(i)%mspBB,2)
            call CopyFreqDomain(obsCont%obs(i)%mspBB(j,k), obsCont%obs(i)%mspBBAtmAtten(j,k))
          end do
        end do
!         call AtmosAtten_old(obsCont%obs(i)%mspBBAtmAtten, obsCont%obs(i)%radDistance(:)%f(1))
        call AtmosAtten(obsCont%obs(i)%mspBBAtmAtten, obsCont%obs(i)%radDistance)
      end if

     !ksb debug: call CalcOctaveFiltMSP(obsCont%obs(i),obsCont%octaveNumber)
      call CalcOctaveFiltMSP(obsCont%obs(i),obsCont%octaveNumber,obsCont%highPassFrequency, &
                             obsCont%lowPassFrequency)
    end do

    if(SPLdBFlag)then
      if(debugLevel.gt.1.and.IsMaster())then
        write(*,*)'Calculating octave filtered SPL (dB) spectral data for '//trim(obsCont%title)
      end if
      do i=1,obsCont%nbObs
        call CalcObsdBFromMSP(obsCont%obs(i)%octaveFiltdB,obsCont%obs(i)%octaveFiltMSP)
      end do
    end if
    
    if(SPLdBAFlag)then
      if(debugLevel.gt.1.and.IsMaster())then
        write(*,*)'Calculating octave filtered SPL (dBA) spectral data for '//trim(obsCont%title)
      end if
      do i=1,obsCont%nbObs
        call CalcObsdBAFromMSP(obsCont%obs(i)%octaveFiltdBA,obsCont%obs(i)%octaveFiltMSP)
      end do
    end if
  end if  
  
  if(OASPLdBFlag)then
    if(debugLevel.gt.1.and.IsMaster())then
      write(*,*)'Calculating OASPL (dB) data for '//trim(obsCont%title)
    end if
    do i=1,obsCont%nbObs
      call CalcObsOASPLdBData(obsCont%obs(i),obsCont%segTimeArray)
    end do
  end if
  
  if(OASPLdBAFlag)then
    if(debugLevel.gt.1.and.IsMaster())then
      write(*,*)'Calculating OASPL (dBA) data for '//trim(obsCont%title)
    end if
    do i=1,obsCont%nbObs
      call CalcObsOASPLdBAData(obsCont%obs(i),obsCont%segTimeArray,obsCont%octaveNumber)
    end do
  end if  

  if((obsCont%nbFreqRanges.gt.0).and.(SPLdBFlag.or.SPLdBAFlag))then
    if(debugLevel.gt.1.and.IsMaster())then
      write(*,*)'Calculating frequency range SPL data for '//trim(obsCont%title)
    end if
    call CalcFreqRangeSPLData(obsCont)
  end if
  
  if(PNLFlag.or.PNLTFlag.or.EPNLFlag.or.SELFlag) then
    if (obsCont%octaveNumber.eq.3)then
      if (PNLFlag.or.PNLTFlag.or.EPNLFlag) then
        if(debugLevel.gt.1.and.IsMaster())then
          write(*,*)'Calculating PNL data for '//trim(obsCont%title)
        end if
      end if
    
      if (SELFlag) then
        if(debugLevel.gt.1.and.IsMaster())then
          write(*,*)'Calculating SEL data for '//trim(obsCont%title)
        end if
      end if
    
      do i=1,obsCont%nbObs
        obs => obsCont%obs(i)  !ksb debug
        if (PNLFlag.or.PNLTFlag.or.EPNLFlag) then
          call CalcObsPNLData(obs,obsCont%segTimeArray)  !ksb debug
          ! call CalcObsPNLData(obsCont%obs(i),obsCont%segTimeArray)
        end if
        
        if (SELFlag) then
          call CalcObsSELData(obsCont%obs(i))
        end if
        if (.not.(all(PNLTdBDrop).or.forceEPNL)) exit
      end do
    else
      call Notice (' PNL, PNLT, EPNL and SEL must be calculated using', &
                   ' 1/3rd octave bands.  This is done by setting octaveFlag=.true.',&
                   ' and octaveNumber=3 in the &ObserverIn namelist.')
    end if
  end if

end subroutine CalcObsContFreqData


subroutine CalcObsContComplexPressure(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obscont
  real,dimension(:),allocatable::freqArray
  real::period,df,hpf,lpf
  integer::i,j,nt,nf,startBin,nbObs,nbSources,nbSegs,obsNum
  nbObs     = obsCont%nbObs
  nbSources = GetNewSizeOfPPrime()
  nbSegs    = obsCont%nbSegments
  period    = GetPeriod(obsCont%obs(1)%pPrimeSegs(1,1))
  df = 1./period
  nt = GetNt(obsCont%obs(1)%pPrimeSegs(1,1))
  lpf = obsCont%lowPassFrequency
  hpf = obsCont%highPassFrequency
  !First we have to determine how big the filtered spectrum will be:
  nf = 0                            
  !This frequency array is the same for both realP and imagP:
  do i=0,nt/2
    if((i*df.ge.hpf-1.e-5).and.(i*df.le.(lpf+1.e-5)))then
      nf = nf+1
    end if
  end do  
  allocate(freqArray(nf))
  !Fill the frequency array:                                             
  do i=1,nt/2+1
    if((i-1)*df.ge.hpf)then
      startBin=i
      do j=1,nf
        freqArray(j)=(i-1)*df+(j-1)*df
      end do
      exit      
    end if
  end do
  !This calculates the complexDFTResult for each pPrime:
  do obsNum=1,obsCont%nbObs
    call CalcObsComplexPressure(obsCont%obs(obsNum),freqArray, &
                                startBin,df,hpf,lpf)
  end do
  deallocate(freqArray)    
end subroutine CalcObsContComplexPressure

!**
!SUBROUTINE CalcFreqRangeSPLData
!  This routine calculations SPL and SPLdBA versus
!  segment time for any frequency range specified in the 
!  namelist.
!**
subroutine CalcFreqRangeSPLData(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  !MSPRanges(obsNum,noiseSource,segNum,rangeNum)
  type(freqDomain),dimension(:,:,:),allocatable::MSPRanges
  integer::obsNum,i,j,k,nbSegs,nbSources
  nbSegs    = size(obsCont%obs(1)%msp,2)
  nbSources = size(obsCont%obs(1)%msp,1)
  allocate(MSPRanges(nbSources,nbSegs,obsCont%nbFreqRanges))
  !First calculate seperate freqDomains for each observer,
  !  noise source, segment, and frequency range:
  do obsNum=1,obsCont%nbObs
    do i=1,nbSources
      do j=1,nbSegs      
        do k=1,obsCont%nbFreqRanges
          nullify(MSPRanges(i,j,k)%f,MSPRanges(i,j,k)%freq)        
          call FilterFreqDomain(obsCont%obs(obsNum)%msp(i,j), obsCont%fRangeMin(k), &
                                obsCont%fRangeMax(k), MSPRanges(i,j,k))
        end do
      end do
    end do
    if(SPLdBFlag)then    
      call CalcObsFreqRangedBData(obsCont%obs(obsNum),MSPRanges, &
                                  obsCont%segTimeArray)
    end if
    if(SPLdBAFlag)then
      call CalcObsFreqRangedBAData(obsCont%obs(obsNum),MSPRanges, &
                                   obsCont%segTimeArray)
    end if
  end do  
  do i=1,nbSources
    do j=1,nbSegs
      do k=1,obsCont%nbFreqRanges
        call DestroyFreqDomain(MSPRanges(i,j,k))
      end do
    end do
  end do
  deallocate(MSPRanges)
end subroutine CalcFreqRangeSPLData

!**
!FUNCTION CalcSegmentTimeArray(obsCont)
!  This function calculates an array of time values
!  according to the segment parameters stored in obsCont.
!  It creates an evenly spaced array of values, with each 
!  value being the starting time of each consecutive segment.
!ARGUMENTS:
! - obsCont: the observer container
! - tMin: the starting time of the array.
!**
subroutine CalcSegmentTimeArray(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  integer::i
  !The first time value is tMin and each consecutive value
  !  is an increment of the segmentStepSize.
  do i=1, size(obsCont%segTimeArray)
    obsCont%segTimeArray(i)=obsCont%tMin+(i-1)*obsCont%segmentStepSize
  end do
end subroutine CalcSegmentTimeArray

subroutine CreateObsContPlane(obsCont,nbx,nby,nbz,xMin,xMax,yMin,yMax,zMin,zMax,&
                              nbBase, unitNumber)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  integer,intent(in)::nbx,nby,nbz
  real,intent(in)::xMin,xMax,yMin,yMax,zMin,zMax
  integer, intent(in)::nbBase, unitNumber !The number of bases and namelist stream number
  real::dx,dy,dz  !The differential values between each point for each axis
  real::xCoord,yCoord,zCoord, xMinToAve,yMinToAve,zMinToAve
  integer::i,j,k,obsNum
  type(CBStructure)::tempCB, previousCB !For the list change of bases if there are some
  
  if(IsMaster()) then
    xMinToAve = 0.
    if(nbx==1) then
       obsCont%NbObsDim1=nby
       obsCont%NbObsDim2=nbz
       dx=0.0
       if((xMin.ne.xMax).and.debugLevel.gt.1)then
         call Warning('In CreateObsContPlane, xMin and xMax are not equal and nbx=1.',&
                      'The x location of each observer will be set equal to the',     &
                      'average of xMin and xMax.')
         xMinToAve = (-xMin+xMax)/2
       end if       
    else
      dx=(xMax-xMin)/real(nbx-1)
    end if
    yMinToAve = 0.
    if(nby==1) then
       obsCont%NbObsDim1=nbx
       obsCont%NbObsDim2=nbz
       dy=0.0
       if(yMin.ne.yMax.and.debugLevel.gt.1)then
         call Warning('In CreateObsContPlane, yMin and yMax are not equal and nby=1.',&
                      'The y location of each observer will be set equal to the',     &
                      'average of yMin and yMax.')
         yMinToAve = (-yMin+yMax)/2                     
       end if    
    else
      dy=(yMax-yMin)/real(nby-1)
    end if
    zMinToAve = 0.   
    if(nbz==1) then
       obsCont%NbObsDim1=nbx
       obsCont%NbObsDim2=nby  
       dz=0.0
       if(zMin.ne.zMax.and.debugLevel.gt.1)then
         call Warning('In CreateObsContPlane, zMin and zMax are not equal and nbz=1.',&
                      'The z location of each observer will be set equal to the',     &
                      'average of zMin and zMax.')
         zMinToAve = (-zMin+zMax)/2
       end if      
    else
      dz=(zMax-zMin)/real(nbz-1)
    end if   
    obsNum = 1
    do k=1, nbz
      do j=1, nby
        do i=1, nbx
          xCoord=xMin+(i-1)*dx + xMinToAve
          yCoord=yMin+(j-1)*dy + yMinToAve
          zCoord=zMin+(k-1)*dz + zMinToAve
          obsCont%obs(obsNum)%coordinates=vectorSetCoordinates(xCoord,yCoord,zCoord)
          obsCont%obs(obsNum)%obsIndex=obsNum
          obsNum=obsNum+1
        end do
      end do
    end do
  else
    obsCont%nbObsDim1 = 1
    obsCont%nbObsDim2 = 1
    obsCont%nbObs     = 1
    allocate(obsCont%obs(1))
    obsCont%obs(1)%coordinates=vectorSetCoordinates(0.0, 0.0, 0.0)
  end if
  call CreateNewObsArray(obsCont)
  !This do loop reads in any change of bases for the grid motion
  do i=1,nbBase
    call createCBData(tempCB,unitNumber)
    call appendCB(obsCont%CBlist,tempCB)
    if(tempCB%windframe) call buildWindFrame(tempCB,previousCB)
    previousCB=tempCB
  end do
end subroutine CreateObsContPlane

subroutine CreateMovingPoint(obsCont, xLoc, yLoc, zLoc, nbBase, unitNumber)
  implicit none
  type(observerContainer), intent(inout)::obsCont  !The observer
  real, intent(in)::xLoc, yLoc, zLoc  !Initial position of the observer
  !The number of bases that define the motion and the stream number of the namelist
  integer, intent(in)::unitNumber, nbBase 
  type(CBStructure)::tempCB, previousCB  !Some values that are needed for reading in a CB
  integer::i
  call CreateNewObsArray(obsCont)
  obsCont%obs(1)%obsIndex = 1
  obsCont%nbObsDim1 = 1
  obsCont%nbObsDim2 = 1
  !First set the coordinates to zero
  obsCont%obs(1)%coordinates=vectorSetCoordinates(xLoc, yLoc, zLoc)
  !Read in the bases that define its motion
  do i=1, nbBase
    call createCBData(tempCB,unitNumber)  
    call appendCB(obsCont%CBlist,tempCB)
    if(tempCB%windframe) call buildWindFrame(tempCB,previousCB)
    previousCB=tempCB
  end do
end subroutine CreateMovingPoint

!**
!SUBROUTINE CreateExternalDefinitionSurface(fileName, obsCont, nbBase, unitNumber)
!  This file reads in a plot3D formatted structured plot file and creates the 
!  coordinates of the observers accordingly.
!ARGUMENTS
! fileName:  The name of the file that contains the data
! obs: The observer
! nbBase: The number of bases that defines the motion
! unitNumber: The number associated with the namelist
!**
subroutine CreateExternalDefinitionSurface(fileName, obsCont, nbBase, unitNumber)
  type(observerContainer), intent(inout)::obsCont       !The observer
  character(len=*), intent(in)::fileName   !The filename that contains the information
  integer, intent(in)::nbBase, unitNumber  !number of bases and stream number of the namelist
  integer::i, j, obsNum  !Some counters
  integer::nbx, nby, nbz !The size of array that is in the file
  integer::unit, stat
  real, dimension(:,:), allocatable:: temp2DArray !A temporary array
  type(CBStructure)::tempCB, previousCB  !Some CB needed or the CB stuff
  if(IsMaster()) then
    !Open the file
    unit=GetStreamNumber()
    open(unit=unit,file=trim(globalFolderName)//trim(fileName), access="sequential", &
       form="formatted", status="old",iostat=stat)
    if( stat /= 0 ) then
      print*,'************************************************************'
      print*,'An error occured while trying to open the observer grid file'
      print*, '"',trim(globalFolderName)//trim(fileName),'".   Please check the file name and location.'
      print*,'************************************************************'
      stop
    end if
    read(unit,*) nbx, nby, nbz
    !Allocating the observer coordinates and the temporary array
    if(nbx==1)then
      obsCont%NbObsDim1=nby
      obsCont%NbObsDim2=nbz
      obsCont%nbObs=nby*nbz             
    else if(nby==1)then
      obsCont%NbObsDim1=nbx
      obsCont%NbObsDim2=nbz
      obsCont%nbObs=nbx*nbz 
    else
      obsCont%NbObsDim1=nbx
      obsCont%NbObsDim2=nby
      obsCont%nbObs=nbx*nby 
    end if     
    allocate(obsCont%obs(obsCont%nbObs),  &
             temp2DArray(3, obsCont%nbObs))           
    read(unit,*) & !Reading it all in
         (temp2DArray(1,obsNum),obsNum=1,obsCont%nbObs), &
         (temp2DArray(2,obsNum),obsNum=1,obsCont%nbObs), &
         (temp2DArray(3,obsNum),obsNum=1,obsCont%nbObs)
    close(unit)
    !This do loop creates the vectors in the observer coordinate array
    obsNum=1
    do j=1, obsCont%nbObsDim2
      do i=1, obsCont%nbObsDim1
        obsCont%obs(obsNum)%obsIndex=obsNum
        obsCont%obs(obsNum)%coordinates=vectorSetValue(temp2DArray(:,obsNum))
        obsNum=obsNum+1
      end do
    end do
  else
    obsCont%nbObsDim1 = 1
    obsCont%nbObsDim2 = 1
    obsCont%nbObs     = 1
    allocate(obsCont%obs(1))
    obsCont%obs(1)%coordinates=vectorSetCoordinates(0.0, 0.0, 0.0)
  end if
  call CreateNewObsArray(obsCont)
  !Reading in the bases that define the observer surface motion
  do i=1,nbBase
    call createCBData(tempCB,unitNumber)
    call appendCB(obsCont%CBlist,tempCB)
    if(tempCB%windframe) call buildWindFrame(tempCB,previousCB)
    previousCB=tempCB
  end do
  !Deallocating the temporary array
  if (allocated (temp2DArray)) deallocate(temp2DArray)
end subroutine CreateExternalDefinitionSurface

!**
!SUBROUTINE CreateFSCObservers(obsCont)
!  This subroutine reads in the file obsCont%FSCInputFile
!  which is a FSC collocation point input file.  PSU-WOPWOP
!  treats each collocation point location as an individual
!  observer point inside an observer container.  
!**
subroutine CreateFSCObservers(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  integer::unit,nsurf,i,surfPoints,totalPoints,j,obsNum
  character(len=4096)::title
  real,dimension(3)::coord
  if(trim(obsCont%FSCInputFileName)=='')then
    call Error('CreateFSCObservers called with no collocation point', &
               'input file specified. Stopping.')
  end if
  unit=GetStreamNumber()
  totalPoints = 0
  surfPoints  = 0
  obsNum      = 0
  if(IsMaster()) then
    open(unit=unit,file=trim(globalFolderName)//trim(obsCont%FSCInputFileName), &
         access='sequential', form='formatted', status='old')
    read(unit,*)nsurf
    do i=1,nsurf
      read(unit,*)
      read(unit,*)surfPoints
      totalpoints=totalPoints+surfPoints
      do j=1,surfPoints
        read(unit,*)
      end do
    end do
    close(unit)
    obsCont%nbObs     = totalPoints
    obsCont%NbObsDim1 = totalPoints
    obsCont%NbObsDim2 = 1
    allocate(obsCont%newObs(obsCont%nbObs))
    obsCont%newObs = .true.
    open(unit=unit,file=trim(globalFolderName)//trim(obsCont%FSCInputFileName), &
         access='sequential', form='formatted', status='old')
    read(unit,*)nsurf
    do i=1,nsurf
      read(unit,*)title
      read(unit,*)surfPoints
      do j=1,surfPoints
        obsNum=obsNum+1
        read(unit,*)coord(1),coord(2),coord(3)
        call SetObserverTitle(obsCont%obs(obsNum),trim(title))
        obsCont%obs(obsNum)%obsIndex=obsNum
        obsCont%obs(obsNum)%coordinates=VectorSetValue(coord)
      end do
    end do
  else
    obsCont%nbObsDim1 = 1
    obsCont%nbObsDim2 = 1
    obsCont%nbObs     = 1
    allocate(obsCont%obs(1))
    obsCont%obs(1)%coordinates=vectorSetCoordinates(0.0, 0.0, 0.0)
  end if
  call CreateNewObsArray(obsCont)
  close(unit)
end subroutine CreateFSCObservers

!**
!SUBROUTINE CreateSphericalSectionObsCont(obsCont, nbTheta, nbPsi, radius, thetaMin, &
!                                         thetaMax, psiMin, psiMax, nbBase, unitNumber)
!  This routine creates a spherical section around zero and reads in an change of bases
! that describe its motion
!ARGUMENTS
! obs: The observer that is being created
! radius: The radius of the sphere section
! nbTheta, nbPsi: The number of points along each rotation (theta is XY plane)
! thetaMin, thetaMax, psiMin, psiMax: Min and max values of the section
! nbBase: The number of bases that define its motion
! unitNumber: The stream number associated with the namelist
!CALLS:
! - CreateCBData       (COB.f90)
! - AppendCB           (COB.f90)
! - BuildWindFrame     (COB.f90)
!**
subroutine CreateSphericalSectionObsCont(obsCont, nbTheta, nbPsi, radius, thetaMin, &
                                         thetaMax, psiMin, psiMax, nbBase, unitNumber, &
                                         indexSwap)
  type(observerContainer), intent(inout)::obsCont
  !The radius is the mean distance from zero.
  real,intent(in)::radius
  !This is the number of evenly spaced angle value.
  ! dTheta is along the XY plane and dPsi is along the YZ or XZ plane.
  integer, intent(in)::nbTheta, nbPsi
  logical, intent(in):: indexSwap
  !These are the angle min and max values
  real, intent(in)::thetaMin, thetaMax, psiMin, psiMax
  !The number of bases that define the motion of the observer and the namelist stream number
  integer, intent(in)::nbBase, unitNumber
  !dTheta and dPsi are the differential angle values
  real::dTheta, dPsi
  !These are some temp values that are needed
  real::tempX1, tempX2, tempX3
  !Some counters that are needed
  integer::iTheta, iPsi, i
  !Some CB that are needed to create in the CB list
  type(CBStructure)::tempCB, previousCB  
  !Fist compute the differential values of each angle
  if(nbTheta==1) then
    dTheta=(thetaMax-thetaMin)/2
  else
    dTheta=(thetaMax-thetaMin)/(nbTheta-1)
  end if
  if(nbPsi==1) then
    dPsi=(psiMax-psiMin)/2
  else
    dPsi=(psiMax-psiMin)/(nbPsi-1)
  end if
  !The following do loop compute the coordinates and sets them in the observer coordinate array
  i=1
  if (IsMaster()) then
    obsCont%nbObsDim1 = nbTheta
    obsCont%nbObsDim2 = nbPsi
    obsCont%nbObs=nbTheta*nbPsi
    allocate(obsCont%obs(obsCont%nbObs))
    if( indexSwap ) then     ! swaped index order
      do iTheta=1, nbTheta
        do iPsi=1, nbPsi
          tempX1=radius*COS((thetaMin+((iTheta-1)*dTheta)))* &
                        COS((psiMin+((iPsi-1)*dPsi)))
          tempX2=radius*SIN((thetaMin+((iTheta-1)*dTheta)))* &
                        COS((psiMin+((iPsi-1)*dPsi)))
          tempX3=radius*SIN((psiMin+((iPsi-1)*dPsi)))
          obsCont%obs(i)%coordinates=vectorSetCoordinates(tempX1, tempX2, tempX3)
          i=i+1
        end do
      end do
    else                                ! original index order
      do iPsi=1, nbPsi
        do iTheta=1, nbTheta
          tempX1=radius*COS((thetaMin+((iTheta-1)*dTheta)))* &
                        COS((psiMin+((iPsi-1)*dPsi)))
          tempX2=radius*SIN((thetaMin+((iTheta-1)*dTheta)))* &
                        COS((psiMin+((iPsi-1)*dPsi)))
          tempX3=radius*SIN((psiMin+((iPsi-1)*dPsi)))
          obsCont%obs(i)%coordinates=vectorSetCoordinates(tempX1, tempX2, tempX3)
          i=i+1
        end do
      end do
    end if
  else
    obsCont%nbObsDim1 = 1
    obsCont%nbObsDim2 = 1
    obsCont%nbObs     = 1
    allocate(obsCont%obs(1))
    obsCont%obs(1)%coordinates=vectorSetCoordinates(0.0, 0.0, 0.0)
  end if
  call CreateNewObsArray(obsCont)
  !Reads in and creates the observer change of base array
  do i=1,nbBase
    call createCBData(tempCB,unitNumber)
    call appendCB(obsCont%CBlist,tempCB)
    if(tempCB%windframe) call buildWindFrame(tempCB,previousCB)
    previousCB=tempCB
  end do
end subroutine CreateSphericalSectionObsCont

! This subroutine creates a logical array in each observer container
! so the code can check whether or not calculations have been 
! performed for each observer.  Currently this routine is only
! used for FAA Certification cases where PNLTM-10dB was not
! achieved and the time range must be extended.
subroutine CreateNewObsArray(obsCont)
  type(observerContainer), intent(inout)::obsCont
  if (associated(obsCont%newObs)) deallocate(obsCont%newObs)
  allocate(obsCont%newObs(obsCont%nbObs))
  obsCont%newObs = .true.

end subroutine CreateNewObsArray


!**
!SUBROUTINE AttachObsContToObject
!ARGUMENTS
! -obs: The observer whose coordinates we are creating 
! -nbBaseObsContFrame: The number of bases in the observer frame
! -nbBaseLocalFrame: The number of bases in the local frame
! -unitNumber: The stream number associated with the namelist
!CALLS:
! - createCBData   (COB.f90)
! - appendCB       (COB.f90)
! - BuildWindFrame (COB.f90)
!**
subroutine AttachObsContToObject(obsCont, nbBaseObsContFrame, nbBaseLocalFrame, unitNumber)
  type(observerContainer), intent(inout)::obsCont  !The observer container
  !This is the number of bases that will be added to the obsCont list before
  !  the object that the observer is attached to.
  integer, intent(in)::nbBaseObsContFrame
  !After the CBList of the object is added to the observers CBList, these will
  !  be added.
  integer, intent(in)::nbBaseLocalFrame
  integer, intent(in)::unitNumber !namelist stream number
  type(CBStructure)::tempCB, previousCB  !needed to readin a CB
  integer::i !counter
  !Reading in the bases in the observer frame
  nullify(obsCont%CBList,obsCont%CBListLocalFrame)
  do i=1, nbBaseObsContFrame
    call createCBData(tempCB,unitNumber) 
    call appendCB(obsCont%CBList,tempCB)
    if (tempCB%windframe) call buildWindFrame(tempCB, previousCB)
   previousCB=tempCB
  end do
  !Reading in the bases in the local frame
  do i=1, nbBaseLocalFrame
    call createCBData(tempCB,unitNumber)  
    call appendCB(obsCont%CBListLocalFrame,tempCB)
    if(tempCB%windframe) call buildWindFrame(tempCB,previousCB)
   previousCB=tempCB
  end do
end subroutine AttachObsContToObject

subroutine DestroyObserverContainer(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  integer::i   
  do i=1,obsCont%nbObs
    call DestroyObserver(obsCont%obs(i))
  end do  
  deallocate(obsCont%obs)
  nullify(obsCont%obs)
  if(associated(obsCont%segTimeArray))then
    deallocate(obsCont%segTimeArray)
    nullify(obsCont%segTimeArray)
  end if
  if(associated(obsCont%fRangeTitle))then
    deallocate(obsCont%fRangeTitle)
    nullify(obsCont%fRangeTitle)
  end if
  if(associated(obsCont%CBlist)) then
    call DestroyCBList(obsCont%CBlist)
    nullify(obsCont%CBlist)
  end if
  if(associated(obsCont%CBlistLocalFrame)) then
    call DestroyCBList(obsCont%CBlistLocalFrame)
    nullify(obsCont%CBlistLocalFrame)
  end if
  if(associated(obsCont%fRangeMin)) then
    deallocate(obsCont%fRangeMin)
    nullify(obsCont%fRangeMin)
  end if
  if(associated(obsCont%fRangeMax)) then
    deallocate(obsCont%fRangeMax)
    nullify(obsCont%fRangeMax)
  end if
end subroutine DestroyObserverContainer

function GetNbObservers(obsCont) result(nbObs)
  implicit none
  type(observerContainer),intent(in)::obsCont
  integer::nbObs
  nbObs=obsCont%nbObs
end function GetNbObservers

subroutine InitializeObserverContainer(obsCont, unitNumber, FAAunitNumber, SPLPath, &
                                       OASPLPath,pressurePath,audioPath,sigmaPath, &
                                       phasePath,complexPressurePath, &
                                       appendTitleFlag)
  implicit none
  !The observer that we are creating
  type(observerContainer), intent(inout)::obsCont
  !A unit number associated with the namelist
  integer, intent(in)::unitNumber, FAAunitNumber
  character(len=*),intent(in)::SPLPath,OASPLPath,pressurePath,audioPath,sigmaPath, &
                               phasePath,complexPressurePath
  logical, intent(in)::appendTitleFlag
  integer::nt, i !i is a counter
  real::tMin, tMax, dt
  real::xLoc, yLoc, zLoc,xMax,yMax,zMax,xMin,yMin,zMin
  integer::nbx,nby,nbz
  !nbBaseObsFrame is for backwards compatability only.
  !  if nbBaseObsFrame is set to some value in the namelist,
  !  nbBaseObsContFrame is automatically set equal to nbBaseObsFrame
  integer::nbBaseObsContFrame,nbBaseObsFrame
  integer::nbBaseLocalFrame
  integer::nbBase
  real::radius
  real::thetaMin, thetaMax, psiMin, psiMax
  integer::nbTheta, nbPsi
  logical:: indexSwap
  character(len=4096)::fileName,windowFunction,Title,attachedTo
  character(len=4096)::FSCInputFileName,FSCOutputFileName
  integer::nbfreqRanges,segmentIncrement,spectrumIncrement
  integer::nbHarmonics, nbNewSegments, nbObsTimeExpansions, maxObsTimeExpansions
  real::octaveNumber,highPassFrequency,lowPassFrequency,minFrequency,maxFrequency
  real::segmentSize,segmentStepSize, windowSize, increment
  real::sigmaRadius, octaveBand
      
  !This is reading all possible values of the observer.  If these values don't exist
  !  they are set to some default value and ignored.  This is declaring the
  !  variables the namelist contains but not writing anything yet.  
  namelist /ObserverIn/ Title,nbBase, attachedTo, xLoc,yLoc,zLoc,           &
                     &  nbBaseObsContFrame,nbBaseObsFrame,nbBaseLocalFrame, &
                     &  nbx,nby,nbz,xMin,yMin,zMin,xMax,yMax,zMax,          &
                     &  fileName,radius,nbTheta,nbPsi, thetaMin,            &
                     &  thetaMax, psiMin, psiMax, indexSwap,                &
                     &  nbfreqRanges,segmentSize,                           &
                     &  segmentStepSize,nbNewSegments, nbObsTimeExpansions, &
                     &  maxObsTimeExpansions, nbHarmonics,tMin,tMax,nt,     &
                     &  highPassFrequency,lowPassFrequency,                 &
                     &  octaveBand,windowFunction,segmentIncrement,         &
                     &  octaveFlag,octaveApproxFlag,octaveNumber,           &
                     &  sigmaRadius,increment,windowSize,spectrumIncrement, &
                     &  FSCInputFileName,FSCOutputFileName, FSCOutputFlag,  &
                     &  FSCThicknessFlag,FSCLoadingFlag, AtmAtten
                     
  namelist /RangeIn/ Title, minFrequency, maxFrequency


  ! Set the default values:
  nbx = 0 ; nby = 0 ; nbz = 0
  xLoc = 0.0 ; yLoc = 0.0 ; zLoc = 0.0
  xMin = 0.0 ; xMax = 0.0 ; yMin = 0.0 ; yMax = 0.0 ; zMin = 0.0 ; zMax = 0.0
  !If nbBaseObsFrame is ever greater than -1, it is because the user 
  !  specified it in the namelist.  Because nbBaseObsFrame is really
  !  nbBaseObsContFrame now, write out an error and stop the user
  !  if nbBaseObsFrame was defined in the namelist.
  attachedTo         = ''              ; nbBaseObsFrame     = -1
  nbBaseObsContFrame = 0               ; nbBaseLocalFrame   = 0
  fileName           = ''              ; radius             = 0
  nbTheta            = 0               ; nbPsi              = 0
  thetaMin           = 0.0             ; thetaMax           = 0.0
  psiMin             = 0.0             ; psiMax             = 0.0; indexSwap=.false.
  nbBase             = 0               ; Title              = 'ObserverContainer'
  nbfreqRanges       = 0               ; sigmaRadius        = 0.1
  highPassFrequency  = 0.0             ; increment          = 0
  lowPassFrequency   = huge(lowPassFrequency)
  segmentSize        = 0.0             ; windowSize         = 0
  segmentStepSize    = 0.0             ; nbHarmonics        = 0
  windowFunction     = ''              ; octaveFlag         = .false.
  octaveApproxFlag   = .true.          ; octaveNumber       = 3.0
  segmentIncrement   = 0               ; spectrumIncrement  = 0
  FSCInputFileName   = ''              ; FSCOutputFlag      = .false.
  FSCOutputFileName  = ''              ; octaveBand         = 0
  tMin               = 0.0             ; FSCThicknessFlag   = .false.   
  tMax               = 0.0             ; FSCLoadingFlag     = .false.
  nbNewSegments      = 1               ; maxObsTimeExpansions= 100
  AtmAtten	     = .false.
 
  if (FAACertFlag) then
    if (FAAunitNumber.gt.0) then
      read(FAAunitNumber, nml=ObserverIn)
      close (FAAunitNumber)
    else
      nt 			       = 35001
      tMin 			       = 0.0
      tMax 			       = 35.0
      nbx 			       = 1
      nby 			       = 3
      nbz 			       = 1
      xmin 			       = 0.0
      xmax 			       = 0.0
      ymin 			       = -150.0
      ymax 			       = 150.0
      zmin 			       = 1.22
      zmax 			       = 1.22
      segmentSize 		   = 0.50
      segmentStepSize 	   = 0.50
      nbNewSegments        = 1
      maxObsTimeExpansions = 100      
      highPassFrequency    = 1.0
      windowFunction 	   = 'Hanning Window'
      atmAtten             = .true.
    end if
  end if
 
  ! Read the values from the namelist file:
  read(unitNumber, nml=ObserverIn)
  
  obsCont%title                = Title
  obsCont%attachedTo           = attachedTo
  obsCont%windowFunction       = windowFunction
  obscont%octaveNumber         = octaveNumber
  obsCont%highPassFrequency    = highPassFrequency
  obsCont%lowPassFrequency     = lowPassFrequency
  obsCont%nbFreqRanges         = nbFreqRanges
  obsCont%nbHarmonics          = nbHarmonics
  obsCont%sigmaRadius          = sigmaRadius
  obsCont%nbSegments           = 0
  obsCont%segmentIncrement     = 0
  obsCont%nbNewSegments        = nbNewSegments 
  obsCont%nbObsTimeExpansions = 0
  obsCont%maxObsTimeExpansions = maxObsTimeExpansions
  obsCont%tMin                 = tmin
  obsCont%tMax                 = tmax  
  obsCont%nt                   = nt
  obsCont%octaveApproxFlag     = octaveApproxFlag
  
  nullify(obsCont%obs,obsCont%segTimeArray,obsCont%fRangeTitle,&
          obsCont%fRangeMax,obsCont%fRangeMin,obsCont%CBList,  &
          obsCont%CBListLocalFrame) 
          
  obsCont%SPLPath             = trim(SPLPath) 
  obsCont%pressurePath        = trim(pressurePath)
  obsCont%audioPath           = trim(audioPath)
  obsCont%sigmaPath           = trim(sigmaPath)
  obsCont%OASPLPath           = trim(OASPLPath)
  obsCont%phasePath           = trim(phasePath)
  obsCont%complexPressurePath = trim(complexPressurePath)
 
  if (appendTitleFlag) then
    obsCont%SPLPath = trim(obsCont%SPLPath)//trim(obsCont%title)
    obsCont%pressurePath = trim(obsCont%pressurePath)//trim(obsCont%title)
    obsCont%audioPath = trim(obsCont%audioPath)//trim(obsCont%title)
    obsCont%sigmaPath = trim(obsCont%sigmaPath)//trim(obsCont%title)
    obsCont%OASPLPath = trim(obsCont%OASPLPath)//trim(obsCont%title)
    obsCont%phasePath = trim(obsCont%phasePath)//trim(obsCont%title)
    obsCont%complexPressurePath=trim(obsCont%complexPressurePath)//trim(obsCont%title)
  end if
    
  obsCont%FSCInputFileName  = FSCInputFileName
  obsCont%FSCOutputFileName = FSCOutputFileName  
  if(FSCOutputFlag)then
    if(trim(FSCOutputFileName).eq.'')then
      call Warning('FSCOutputFlag requires FSCOutputFileName to be specified in the', &
                   'observer namelist.  Defaulting to FSC.out.')
      obsCont%FSCOutputFileName = 'FSC.out'
    end if
    if((.not.pressureGradientFlag).and.(.not.pressureGradient1AFlag))then
      call Error('Outputting data for FSC requires pressureGradient1AFlag to be on.', &
                 'Check observer namelist. Stopping.')
    end if
    if(FSCThicknessFlag.and.(.not.(thicknessNoiseFlag.or.totalNoiseFlag)))then
      call Error('FSCThicknessFlag is true without either thicknessNoiseFlag or', &
                 'totalNoiseFlag true.  Check namelist input. Stopping.')
    end if
    if(FSCLoadingFlag.and.(.not.(loadingNoiseFlag.or.totalNoiseFlag)))then
      call Error('FSCLoadingFlag is true without either loadingNoiseFlag or', &
                 'totalNoiseFlag true.  Check namelist input. Stopping.')
    end if
    if((.not.(FSCLoadingFlag.or.FSCThicknessFlag)).and.(.not.totalNoiseFlag))then
      call Error('FSCOutputFlag requires totalNoiseFlag be true unless either', &
                 'FSCLoadingFlag or FSCThicknessFlag are true. Stopping.')
    end if
  end if 
  if((trim(FSCOutputFileName).ne.'').and.(.not.FSCOutputFlag))then
    call Warning('FSCOutputFileName is specified in the namelist and FSCOutputFlag', &
                 'is set to false.  No FSC data will be output.')
  end if 
      !The inputs needed from FSC collocation point file are just the x,y, and z
      !  location of the observers.  For FSC acoustic input they do not need
      !  a grid file, and we just need to provide the real and imaginary components
      !  of pressure gradient for each point.
  if(trim(FSCInputFileName).ne.'')then
    call CreateFSCObservers(obsCont)
  else if(trim(fileName).ne.'')then
    call CreateExternalDefinitionSurface(fileName,obsCont,nbBase,unitNumber)
  else if((nbTheta.ne.0).and.(nbPsi.ne.0).and.(radius.gt.0.0))then
    call CreateSphericalSectionObsCont(obsCont, nbTheta, nbPsi, radius, thetaMin, &
                                       thetaMax, psiMin, psiMax, nbBase, unitNumber,&
                                       indexSwap)
  else if((nbx.ne.0).and.(nby.ne.0).and.(nbz.ne.0))then
    obsCont%nbObs=nbz*nby*nbx
    allocate(obsCont%obs(obsCont%nbObs))
    call CreateObsContPlane(obsCont, nbx, nby, nbz, xMin, xMax, &
                                                    yMin, yMax, &
                                                    zMin, zMax, &
                                                    nbBase, unitNumber)
  else
    obsCont%nbObs=1
    allocate(obsCont%obs(1))
    call CreateMovingPoint(obsCont, xLoc, yLoc, zLoc, nbBase, unitNumber)
  end if  
  !This is for backwards compatibility:
  if(nbBaseObsFrame.ge.0)then
    call Error('nbBaseObsFrame and nbBaseObsContFrame are the same variable.', &
               'Use nbBaseObsContFrame in place of nbBaseObsFrame in the observer', &
               'namelist for future runs (see user manual).  Stopping.')
  end if 
  if (trim(attachedTo) /= '') then
    call AttachObsContToObject(obsCont, nbBaseObsContFrame, nbBaseLocalFrame, unitNumber)
  end if  

  !The following is for backwards compatibility:
  if(windowSize.gt.0)then
    segmentSize = windowSize
    call Notice('Use the variable segmentSize in place of windowSize', &
                'in future runs.  These variables are the same thing.')
  end if
  if(increment.gt.0)then
    segmentStepSize = increment
    call Notice('Use the variable segmentStepSize in place of increment', &
                'in future runs.  These variables are the same thing.')
  end if
  if((spectrumIncrement.gt.0).and.(segmentIncrement.gt.0))then
    call Error('Only specify segmentIncrement, not spectrumIncrement.', &
               'These two variables mean the same thing.', &
               'Stopping.')
  else if(spectrumIncrement.gt.0)then
    call Notice('Use the variable segmentIncrement in place of spectrumIncrement', &
                'for future runs.  They mean the same thing.')
    segmentIncrement = spectrumIncrement
  end if
  
  !We are going to set the segment parameters now. However if 
  !  tMin = tMax = 0. (value not set in the namelist) then we must
  !  check them later in the code after we calculate the observer
  !  after the sources are setup.
  obsCont%segmentSize        = segmentSize
  obsCont%segmentStepSize    = segmentStepSize
  obsCont%segmentIncrement   = segmentIncrement
  
  if(nbHarmonics.gt.0)then
    call Warning('NbHarmonics was specified. This will filter the FFT so only the first', &
                 'nbHarmonics number of bins are shown. However, unless the input ', &
                 'time parameters leads to a frequency bin width equal to the blade passage frequency,',&
                 'nbHarmonics will not output the first nbHarmonics of blade passage frequency.')
  end if
  
  !Check the segment parameters, we know the observer time range.
  !  If we dont know the time range we need to check this later:
  if(tMin.ne.tMax)then
    call CheckSegmentTimeSettings(obsCont)  
  end if
  
  if (octaveBand.ne.0) then
    call Notice('Use the variable octaveBand is no longer used.', &
                'Use octaveNumber in future runs.', &
                'Setting octaveNumber equal to octaveBand.')
    obsCont%octaveNumber = octaveBand
  end if
  
  !This is a check to make sure we are not writing out a sound file or a spectrum over a grid, 
  !which would make no sense
  if(audioFlag) then
    if(obsCont%nbObs.gt.1) then
      call Error('Cannot create audio file for observer surface or line.', &
                 'Check namelist file.  Stopping.')
    end if
    if(NotMaster())then
      call Error('Cannot create audio files when running a multiple observer', &
                 'case in parallel.  Check namelist file. Stopping.')
    end if
    if(integrationType==FREQUENCY_DOMAIN) then
      call Warning('Cannot write out audio when integration by frequency.',  &
                   'Turning audioFlag to .false.')
      audioFlag = .false.
    end if
  end if
  
  !If the acoustic pressure is recorded for a grid the function file can get extremely large
  if((acousticPressureFlag.or.spectrumFlag).and.(obsCont%nbObs.gt.1)) then 
    call Warning('Recording pressure or spectrum over grid.','This may be a VERY large file.')
  end if
  if(acousticPressureFlag.and.integrationType==FREQUENCY_DOMAIN) then
    call Warning('Writing out acoustic pressure time history when',               &
                 'adding the signal by amplitude not time history.',                &
                 'A random phase will be used to compute time history.')
    spectrumFlag = .true.
  end if  
 
  if(highPassFrequency.ge.lowPassFrequency)then
    call Error('highPassFrequency is greater than lowPassFrequency. Stopping.')
  end if
  !
  ! Check to see if frequency range is too small for any harmonics
  !
  if( tmax.ne.tmin .and. (lowPassFrequency-highPassFrequency).le. 1./(tMax-tMin) ) then
    call Error('The diff between lowPassFrequency and highPassFrequency is less that&
    &  1/(tMax-tMin). Stopping.')
  end if
  
  if(nbFreqRanges.ne.0) then
      ! Read in the ranges and create the data:
    allocate (obsCont%fRangeTitle(nbFreqRanges), &
              obsCont%fRangeMin(nbFreqRanges),   &
              obsCont%fRangeMax(nbFreqRanges))
    do i=1,nbFreqRanges
      read (unitNumber, nml=RangeIn)
      obsCont%fRangeTitle(i) = Title
      obsCont%fRangeMin(i)   = minFrequency
      obsCont%fRangeMax(i)   = maxFrequency
      if(obsCont%nt/2/(obsCont%tMax-obsCont%tMin).le.minFrequency)then
        call Error('A frequency range was requested where the minimum frequency is', &
                   'above the Nyquist frequency. Check the namelist. Stopping.')
        else if(obsCont%nt/2/(obsCont%tMax-obsCont%tMin).lt.maxFrequency)then
          call Warning('A frequency range was requested where the maximum frequency',     &
                       'is above the Nyquist frequency.  Changing the maximum frequency', &
                       'to equal the Nyquist frequency.')
          maxFrequency = obsCont%nt/2/(obsCont%tMax-obsCont%tMin)
      end if
    end do
  end if
  do i=1, obsCont%nbObs
    obsCont%obs(i)%nbFreqRanges = obsCont%nbFreqRanges
  end do
end subroutine InitializeObserverContainer

subroutine CreateNewObservers(obsCont, newObs)
  implicit none
  type(observerContainer), intent(inout)::obsCont
  type(observer), dimension(:), allocatable:: tempObsArray
  logical:: newObs,  copyTF
  type(timeHistory), dimension(:,:), allocatable:: tH
  type(timeHistory), dimension(:,:), pointer:: tempRadDistance
  type(freqDomain), dimension(:,:), pointer:: freqArray
  character(len=4096):: title
  integer:: i, j, k, pPrimeSize, nbObs, copyNbObs, segmentSizeObs1, segmentSizeObs2
  integer:: oldTHICK_APTH, oldLOAD_APTH, oldBB_PEGG_APTH, oldBB_BPM_APTH, oldTOTAL_APTH, oldTHICK_PGX, &
            oldTHICK_PGY, oldTHICK_PGZ, oldLOAD_PGX, oldLOAD_PGY, oldLOAD_PGZ, oldTOTAL_PGX, oldTOTAL_PGY, oldTOTAL_PGZ
  !
  ! If newObs = .true. then this subroutine is going to set up a second observer
  ! which will be one time segment long, which later should be added to the original
  ! observer.  So when newObs = .true., so we double the number of observers.
  ! 
  ! If newObs = .false. then this subroutine is going to combine the new single 
  ! segment observers with the original observers.  Hence the number of observers
  ! will be reduced into half.
  !
  ! The logical array newObs keeps track which of the observers are "new" (.true.) and
  ! which observers are "old" (.false.)
  !
  ! Wipe the newObs logical array.
  !
  deallocate(obsCont%newObs)
  ! Copy the current observer information into a temporary array.
  pPrimeSize = GetSizeOfPPrime()
  copyNbObs = obsCont%nbObs
  if (dataLimitRHS) copyNbObs = obsCont%nbObs/2  
  if (newObs) then
    allocate(obsCont%newObs(2*obsCont%nbObs))
    obsCont%newObs(1:obsCont%nbObs) = .false.
    obsCont%newObs(obsCont%nbObs+1:2*obsCont%nbObs) = .true.
  else
    allocate(obsCont%newObs(obsCont%nbObs/2))
    obsCont%newObs = .false.
    if (associated(ObsCont%obs(1)%mspBB)) then !ksb debug: not sure what happens if neither mspBB AND radDistance associated
      segmentSizeObs2 = size(ObsCont%obs(1+ObsCont%NbObs/2)%mspBB,2)
    elseif (associated(ObsCont%obs(1)%radDistance)) then
      segmentSizeObs2 = size(ObsCont%obs(1+ObsCont%NbObs/2)%radDistance)
    end if
  end if

  allocate(tempObsArray(copyNbObs))
  allocate(tH(copyNbObs, pPrimeSize))
  if (associated(ObsCont%obs(1)%mspBB)) then  !ksb debug: not sure what happens if neither mspBB AND radDistance associated
    segmentSizeObs1 = size(ObsCont%obs(1)%mspBB,2)
  elseif (associated(ObsCont%obs(1)%radDistance)) then
    segmentSizeObs1 = size(ObsCont%obs(1)%radDistance)
  end if

  do i=1,copyNbObs
    if (.not. dataLimitRHS) then
      call CopyObserverData(obsCont%obs(i), tempObsArray(i), obsCont%tMin, &
                            obsCont%tMax, obsCont%nt, obsCont%dt)
    else
      ! If we have hit the data limit on the RHS then we are simply destroying 
      ! the current data set as it contains erroneous/incomplete results.
      call CopyObserverData(obsCont%obs(i), tempObsArray(i), obsCont%obs(i)%tMin, &
                            obsCont%obs(i)%tMax, obsCont%obs(i)%nt, obsCont%dt)    
    end if
  end do

  ! Copy the time history arrays.  If making new observers then
  ! simply copy the time historys to a temporary array.
  ! If combining observers then we need a little extra code
  ! to merge the sets of data, specifically the location of 
  ! each observer's pressure history arrays in the arrays of the
  ! new observer.
  ! GetNewSizeOfPPrime()
  ! If the time range of one of the observers we are combining
  ! is not within the observer time range when we simply 
  ! copy the 'original' observer as if it were the only one.
  
  ! We need to be careful about copying things back if any of the
  ! calculation flags are turned off since the ResetPPrimeAndIndices routine 
  ! rearranges things. We need to save over the 'old' indices and then
  ! re-reset them to get their correct placement since we are effectively
  ! rerunning the case.
  oldTHICK_APTH = THICK_APTH
  oldLOAD_APTH = LOAD_APTH
  oldBB_PEGG_APTH = BB_PEGG_APTH
  oldBB_BPM_APTH = BB_BPM_APTH
  oldTOTAL_APTH = TOTAL_APTH
  oldTHICK_PGX = THICK_PGX
  oldTHICK_PGY = THICK_PGY
  oldTHICK_PGZ = THICK_PGZ
  oldLOAD_PGX = LOAD_PGX
  oldLOAD_PGY = LOAD_PGY
  oldLOAD_PGZ = LOAD_PGZ
  oldTOTAL_PGX = TOTAL_PGX
  oldTOTAL_PGY = TOTAL_PGY
  oldTOTAL_PGZ = TOTAL_PGZ
  call SetPPrimeIndices()
  if (newObs .or. dataLimitRHS) then
    if (globalPeggNoiseFlag .or. globalBPMNoiseFlag) allocate(freqArray(copyNbObs, segmentSizeObs1))
    if (AtmAtten) allocate(tempRadDistance(copyNbObs, segmentSizeObs1))

    tH%nt = obsCont%obs(1)%nt
    do i=1,copyNbObs
      do j=1,pPrimeSize
        allocate(tH(i,j)%T(tH(1,1)%nt), tH(i,j)%F(tH(1,1)%nt))
        call WhatToCopy(j, k, copyTF, oldTHICK_APTH, oldLOAD_APTH, oldTOTAL_APTH, &
                oldBB_PEGG_APTH, oldBB_BPM_APTH, oldTHICK_PGX, oldTHICK_PGY, oldTHICK_PGZ, &
                oldLOAD_PGX, oldLOAD_PGY, oldLOAD_PGZ, oldTOTAL_PGX, oldTOTAL_PGY, oldTOTAL_PGZ)
        title = PPRIMETITLEARRAY(j)
        call CopyTH(tH(i,j), tH(i,j)%nt, title, copyTF, obsCont%obs(i)%pPrime(k))
      end do
      if (globalPeggNoiseFlag .or. globalBPMNoiseFlag) then
        do j=1,segmentSizeObs1
         call CopyFreqDomain(obsCont%obs(i)%mspBB(1,j), freqArray(i,j))
        end do
      end if
      if (AtmAtten) then
        do j=1,segmentSizeObs1
          call CopyTimeHistory(obsCont%obs(i)%radDistance(j), tempRadDistance(i,j))
        end do
      end if
    end do
  else
    ! Condense the observer array back to its original size
    if (globalPeggNoiseFlag .or. globalBPMNoiseFlag) allocate(freqArray(obsCont%nbObs/2,segmentSizeObs1+segmentSizeObs2))
    if (AtmAtten) allocate(tempRadDistance(obsCont%nbObs/2, segmentSizeObs1+segmentSizeObs2))

    tH%nt = obsCont%obs(1)%nt + obsCont%obs(1+obsCont%nbObs/2)%nt-1
    do i=1,obsCont%nbObs/2
      call MergeObserverData(tempObsArray(i), tempObsArray(i+obsCont%nbObs/2))
      do j=1,pPrimeSize
        allocate(tH(i,j)%T(tH(1,1)%nt), tH(i,j)%F(tH(1,1)%nt))
        call CopyTH(tH(i,j), tH(i,j)%nt, obsCont%obs(i)%pPrime(j)%title, .true., obsCont%obs(i)%pPrime(j), &
                          obsCont%obs(i+obsCont%nbObs/2)%pPrime(j))
      end do
      if (newSegLHS) then
        if (globalPeggNoiseFlag .or. globalBPMNoiseFlag) then
          do j=1,segmentSizeObs2
	        call copyFreqDomain(obsCont%obs(i+obsCont%nbObs/2)%mspBB(1,j), freqArray(i,j))
          end do
          do j=1,segmentSizeObs1
            call copyFreqDomain(obsCont%obs(i)%mspBB(1,j), freqArray(i,j+segmentSizeObs2))
          end do
        end if

        if (AtmAtten) then
          do j=1,segmentSizeObs2
            call CopyTimeHistory(obsCont%obs(i+obsCont%nbObs/2)%radDistance(j), tempRadDistance(i,j))
          end do
          do j=1,segmentSizeObs1
            call CopyTimeHistory(obsCont%obs(i)%radDistance(j), tempRadDistance(i,j+segmentSizeObs2))
          end do
        end if

      elseif (newSegRHS) then
        if (globalPeggNoiseFlag .or. globalBPMNoiseFlag) then
          do j=1,segmentSizeObs1
	        call copyFreqDomain(obsCont%obs(i)%mspBB(1,j), freqArray(i,j))
          end do
          do j=1,segmentSizeObs2
            call copyFreqDomain(obsCont%obs(i+obsCont%nbObs/2)%mspBB(1,j), freqArray(i,j+segmentSizeObs1))
          end do
        end if

        if (AtmAtten) then
          do j=1,segmentSizeObs1
            call CopyTimeHistory(obsCont%obs(i)%radDistance(j), tempRadDistance(i,j))
          end do
          do j=1, segmentSizeObs2
            call CopyTimeHistory(obsCont%obs(i+obsCont%nbObs/2)%radDistance(j), tempRadDistance(i,j+segmentSizeObs1))
          end do
        end if
      end if
    end do
  end if

  ! Now that everything is copied over, destroy the current observer arrays.
  do i=1,obsCont%nbObs
    call DestroyObserver(obsCont%obs(i))
  end do  
  deallocate(obsCont%obs)
  nullify(obsCont%obs)

  if (newObs) then
    copyNbObs = obsCont%nbObs
    obsCont%nbObs = 2*obsCont%nbObs
  else 
    obsCont%nbObs = obsCont%nbObs/2
    copyNbObs = obsCont%nbObs     
  end if

  ! Reallocate the observer array and copy the data back into it.
  allocate(obsCont%obs(obsCont%nbObs))
  do i=1,copyNbObs
    call CopyObserverData(tempObsArray(i), obsCont%obs(i), &
                          tempObsArray(i)%tMin, tempObsArray(i)%tMax, tempObsArray(i)%nt, &
                          tempObsArray(i)%dt)
    if (newObs) then
      obsCont%obs(i+obsCont%nbObs/2)%coordinates = obsCont%obs(i)%coordinates
      obsCont%obs(i+obsCont%nbObs/2)%obsIndex    = 2*i
    end if
    ! Now reallocate the pPrime arrays and move everything back
    allocate(obsCont%obs(i)%pPrime(pPrimeSize))
    do j=1,pPrimeSize
      allocate(obsCont%obs(i)%pPrime(j)%T(tH(i,j)%nt), obsCont%obs(i)%pPrime(j)%F(tH(i,j)%nt))
      call CopyTH(obsCont%obs(i)%pPrime(j), tH(i,j)%nt, tH(i,j)%title, .true., tH(i,j))
    end do
    if (globalPeggNoiseFlag .or. globalBPMNoiseFlag) then
      allocate(obsCont%obs(i)%mspBB(1,size(freqArray,2)))
      do j=1,size(freqArray,2)
	call CopyFreqDomain(freqArray(i,j), obsCont%obs(i)%mspBB(1,j))
      end do
    end if
    if (AtmAtten) then
      allocate(obsCont%obs(i)%radDistance(size(tempRadDistance,2)))
      do j=1,size(tempRadDistance,2)
        call CopyTimeHistory(tempRadDistance(i,j), obsCont%obs(i)%radDistance(j))
      end do
    end if
  end do

  ! Finally, destroy the temporary arrays.
  do i=1,size(tH,1)
    call DestroyObserver(tempObsArray(i))
    do j=1,size(tH,2)
      call DestroyTimeHistory(tH(i,j))
    end do
  end do
  deallocate(tH, tempObsArray)
  if (globalPeggNoiseFlag .or. globalBPMNoiseFlag) then
    call Destroy2DFDArray(freqArray)
    deallocate(freqArray)
    nullify(freqArray)
  end if
  if (AtmAtten) then
    call Destroy2DTHArray(tempRadDistance)
    deallocate(tempRadDistance)
    nullify(freqArray)
  end if

  ! Lastly, if we're at the point of calculating results
  ! we need to put the PPrime indices back the way they were.
  if (.not. newObs) call ResetPPrimeIndices()  !not sure if this is needed any longer
end subroutine CreateNewObservers

subroutine WhatToCopy(j, k, copyTF, oldTHICK_APTH, oldLOAD_APTH, oldTOTAL_APTH, &
                oldBB_PEGG_APTH, oldBB_BPM_APTH, oldTHICK_PGX, oldTHICK_PGY, oldTHICK_PGZ, &
                oldLOAD_PGX, oldLOAD_PGY, oldLOAD_PGZ, oldTOTAL_PGX, oldTOTAL_PGY, oldTOTAL_PGZ)
  implicit none
  integer:: j, k
  logical::  copyTF
  integer, intent(inout):: oldTHICK_APTH, oldLOAD_APTH, oldTOTAL_APTH, &
            oldBB_PEGG_APTH, oldBB_BPM_APTH, oldTHICK_PGX, oldTHICK_PGY, oldTHICK_PGZ, &
            oldLOAD_PGX, oldLOAD_PGY, oldLOAD_PGZ, oldTOTAL_PGX, oldTOTAL_PGY, oldTOTAL_PGZ
  
  copyTF = .false.
  if (j.eq.THICK_APTH) then
      k = oldTHICK_APTH
      if (ThicknessNoiseFlag) copyTF = .true.
  else if (j.eq.LOAD_APTH) then
      k = oldLOAD_APTH
      if (LoadingNoiseFlag) copyTF = .true.
!  else if (j.eq.BB_PEGG_APTH) then
!      k = oldBB_PEGG_APTH
!      if (globalPeggNoiseFlag) copyTF = .true.
!  else if (j.eq.BB_BPM_APTH) then
!      k = oldBB_BPM_APTH
!      if (globalBPMNoiseFlag) copyTF = .true.
  else if (j.eq.TOTAL_APTH) then
    k = oldTOTAL_APTH
    if (TotalNoiseFlag) copyTF = .true.
  else if (j.eq.THICK_PGX) then
    k = oldTHICK_PGX
    if (PressureGradient1AFlag.and.ThicknessNoiseFlag) copyTF = .true.
  else if (j.eq.THICK_PGY) then
    k = oldTHICK_PGY  
    if (PressureGradient1AFlag.and.ThicknessNoiseFlag) copyTF = .true.
  else if (j.eq.THICK_PGZ) then
    k = oldTHICK_PGZ
    if (PressureGradient1AFlag.and.ThicknessNoiseFlag) copyTF = .true.
  else if (j.eq.LOAD_PGX) then
    k = oldLOAD_PGX
    if (PressureGradient1AFlag.and.LoadingNoiseFlag) copyTF = .true.
  else if (j.eq.LOAD_PGY) then
    k = oldLOAD_PGY
    if (PressureGradient1AFlag.and.LoadingNoiseFlag) copyTF = .true.
  else if (j.eq.LOAD_PGZ) then
    k = oldLOAD_PGZ
    if (PressureGradient1AFlag.and.LoadingNoiseFlag) copyTF = .true.
  else if (j.eq.TOTAL_PGX) then
    k = oldTOTAL_PGX
    if (PressureGradient1AFlag.and.TotalNoiseFlag) copyTF = .true.
  else if (j.eq.TOTAL_PGY) then
    k = oldTOTAL_PGY
    if (PressureGradient1AFlag.and.TotalNoiseFlag) copyTF = .true.
  else if (j.eq.TOTAL_PGZ) then
    k = oldTOTAL_PGZ
    if (PressureGradient1AFlag.and.TotalNoiseFlag) copyTF = .true.
  end if

end subroutine WhatToCopy


subroutine CopyObserverData(oldObs, newObs, tMin, tMax, nt, dt)
  implicit none
  type(observer):: oldObs, newObs
  integer:: nt
  real:: tMin, tMax, dt
  
  newObs%tMin = tMin
  newObs%tMax = tMax
  newObs%nt   = nt
  newObs%dt   = dt  
  newObs%title = oldObs%title
  newObs%obsIndex = oldObs%obsIndex
  newObs%coordinates = oldObs%coordinates
  
end subroutine CopyObserverData


subroutine CopyTH(tHTemp, nt, title, copyTF, tH1, tH2)
  implicit none
  type(timeHistory)::tHTemp, tH1 
  type(timeHistory),optional::tH2 
  integer:: nt
  character(len=4096):: title
  logical:: copyTF

  tHTemp%title            = title
  tHTemp%nt               = nt
  if (present(tH2)) then
    if (associated(tH1%T)) then
      if (newSegLHS) then
        tHTemp%T(1:(tH2%nt-1))    = tH2%T(1:(tH2%nt-1))
        tHTemp%F(1:(tH2%nt-1))    = tH2%F(1:(tH2%nt-1))
        tHTemp%T(tH2%nt:thTemp%nt) = tH1%T
        tHTemp%F(tH2%nt:tHTemp%nt) = tH1%F 
      else if (newSegRHS) then
        tHTemp%T(1:tH1%nt)    = tH1%T
        tHTemp%F(1:tH1%nt)    = tH1%F
        tHTemp%T(tH1%nt+1:tHTemp%nt) = tH2%T(2:tH2%nt)
        tHTemp%F(tH1%nt+1:tHTemp%nt) = tH2%F(2:tH2%nt)
      end if
    end if
    tHTemp%tMin             = minval((/tH1%tMin,tH2%tMin/))
    tHTemp%tMax             = maxval((/tH1%tMax,tH2%tMax/))
  else
    if (copyTF) then
      tHTemp%T                = tH1%T
      tHTemp%F                = tH1%F
    end if
    tHTemp%tMin             = tH1%tMin
    tHTemp%tMax             = tH1%tMax
  end if
    
end subroutine CopyTH

subroutine MergeObserverData(obs1, obs2)
  type(observer), intent(inout):: obs1, obs2

  obs1%tMin = minval((/obs1%tMin, obs2%tMin/))
  obs1%tMax = maxval((/obs1%tMax, obs2%tMax/))
  
end subroutine MergeObserverData

subroutine CheckFreqAnalysisParameters(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  
  if(obsCont%nbHarmonics.gt.(obsCont%ntPerSegment/2-1))then
    call Warning('nbHarmonics is greater than the maximum allowable for',    &
                 'the time specified in the namelist.  The entire spectrum', &
                 'will be analyzed.')
    obsCont%nbHarmonics = obsCont%ntPerSegment/2-1
  end if  

  if(obsCont%nbHarmonics.gt.0)then
    call Warning('NbHarmonics was specified and will override the use of', &
                 'highPassFrequency and lowPassFrequency')
    obsCont%highPassFrequency = .5/obsCont%segmentSize  
! KSB: might want to chanage this; there is no reason we 
!shouldn't be able to use a high pass filter and a limited 
!number of harmonics. Need to check if the highPassFrequency 
!is above the max frequency using nbHarmonics or if 
!lowPassFrequency is below the min frequency
    obsCont%lowPassFrequency = (real(obsCont%nbHarmonics)+1.)/obsCont%segmentSize &
                               + 0.5/obsCont%segmentSize
  end if

  !Here we are going to check if the frequency data calculated
  !  from the pressure time histories spans the frequency range
  !  necessary for PNL calculations.  This range is 50 - 10000Hz.
  !  If PNLFlag is on and the data does not span the range, write  
  !  an warning and turn the flag off.
  
  if(PNLFlag.or.PNLTFlag.or.EPNLFlag)then  
    if(((obsCont%ntPerSegment/2+1)/(obsCont%segmentSize).lt.10000) .or. &
       (obsCont%lowPassFrequency.lt.10000.) .or.                &
       (obsCont%highPassFrequency.gt.50.)   .or.                &
       ((obsCont%nbHarmonics.gt.0.).and.                        &
       obsCont%nbHarmonics/obsCont%segmentSize.lt.10000.))then !ksb debug - need to change this for helicopter (PartH36Flag?)
      call Notice('FAA regulation PNL calculations require a spectral range of 50 Hz ',&
                 'through 10000Hz. The largest frequency resulting from a FFT of the',&
                 'segmented observer time does not fall within these limits.')
    end if
  end if  
end subroutine CheckFreqAnalysisParameters


subroutine WriteSegmentDataToScreen(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  real::period
  period = obsCont%tMax - obsCont%tMin
  if(IsMaster()) then
    write(*,'(A,G15.6,A)')' Total Time =',period,' seconds.'
    write(*,'(A,G15.6,A)')' tMin =',obsCont%tMin,' seconds.'
    write(*,'(A,G15.6,A)')' tMax =',obsCont%tMax,' seconds.'
    if(obsCont%nbSegments.gt.0.)then
      write(*,*)'Number of segments = ',obsCont%nbSegments
      write(*,*)'Nt per segment = ',obsCont%ntPerSegment
      write(*,'(A,G15.6,A)')' Segment size =',obsCont%segmentSize,' seconds.'
      if(obsCont%nbSegments.gt.1)then  
        write(*,'(A,G15.6,A)')' Segment step size =',obsCont%segmentStepSize,' seconds.'
      end if
      write(*,*)'Nyquist frequency = ',(obsCont%ntPerSegment/2.)/(obsCont%segmentSize),' Hz.'
    else
      write(*,*)'Analyzing the entire time history.'
    end if
  end if                                       
end subroutine WriteSegmentDataToScreen


subroutine CheckSegmentTimeSettings(obsCont)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  real::segmentSize,segmentStepSize
  integer::segmentIncrement,ntPerSegment,nbSegments
  !First copy the values over into dummy storage:
  segmentSize = obsCont%segmentSize
  segmentStepSize = obsCont%segmentStepSize
  segmentIncrement = obsCont%segmentIncrement  
  if((segmentSize.gt.0).and.(segmentStepSize.eq.0))then
    call Error('If segmentSize is specified, segmentStepSize must', &
               'also be specified. Stopping.')
  else if((segmentSize.eq.0).and.(segmentStepSize.gt.0))then
    call Error('If segmentStepSize is specified, segmentSize must', &
               'also be specified. Stopping.')
  else if((segmentSize.gt.0).and.(segmentStepSize.gt.0))then
    call CheckAndFinalizeSegParameters(obsCont%tMin,obsCont%tMax,   &
                                       obsCont%nt,                  &
                                       segmentSize,segmentStepSize, &
                                       nbSegments,segmentIncrement, &
                                       ntPerSegment)
    if(nbSegments.gt.0)then    
      if( dataLimitRHS ) nbSegments = nbSegments-1 !ksb debug:
      obsCont%segmentSize     = segmentSize
      obsCont%segmentStepSize = segmentStepSize
      obsCont%nbSegments      = nbSegments
      obsCont%ntPerSegment    = ntPerSegment    
      allocate(obsCont%segTimeArray(nbSegments))                                          
      call CalcSegmentTimeArray(obsCont)
      !This checks to make sure that segmentIncrement is valid:
      if(segmentIncrement.gt.nbSegments)then
        call Error('segmentIncrement specified in the namelist is greater than',&
                     'the number of created segments. Stopping.')
      else if(segmentIncrement.lt.0)then
        call Error('segmentIncrement specified in the namelist is less than zero.',&
                   'Stopping.')
      else
        obsCont%segmentIncrement=segmentIncrement
        if(segmentIncrement.ne.0)then
          write(*,*)'Performing frequency analysis on segment '// &
                    trim(IntegerToString(segmentIncrement))//' only.'
        end if
      end if
    else if(nbSegments.eq.0)then
      allocate(obsCont%segTimeArray(1))
      obsCont%segTimeArray(1) = obsCont%tMin
      obsCont%segmentSize     = obsCont%tMax-obscont%tMin
!      obsCont%segmentStepSize = 0.
      obsCont%ntPerSegment    = obsCont%nt
    end if      
  else
    allocate(obsCont%segTimeArray(1))
    obsCont%segTimeArray(1) = obsCont%tMin
    obsCont%segmentSize     = obsCont%tMax-obsCont%tMin
    obsCont%segmentStepSize = 0.
    obsCont%ntPerSegment    = obsCont%nt
  end if

end subroutine CheckSegmentTimeSettings 


!**
!SUBROUTINE InitializeObservers
!  This routine initializes all arrays in an observer and
!  defines parameters that were set in the namelist.
!**
subroutine InitializeObservers(obsCont, minObsTimeInData, maxObsTimeInData)
  implicit none
  type(observerContainer),intent(inout)::obsCont
  real(kind=8):: minObsTimeInData, maxObsTimeInData
  
  integer::i, nt, startIndex, EndIndex
  real:: tMin, tMax, dt, temp
  
  do i=1,obsCont%nbObs
    if (obsCont%newObs(i)) then
      ! This is the most common case, we're simply putting the observer
      ! container information with each observer
      tMin = obsCont%tMin
      tMax = obsCont%tMax      
      dt   = (obsCont%tMax-obsCont%tMin)/(obsCont%nt-1)      
      if (newSegLHS .and. .not.dataLimitLHS) then
        ! Create a new segment just prior to the min. observer time
        tMax = obsCont%tMin
        tMin = obsCont%tMin - obsCont%nbNewSegments*obsCont%SegmentSize
        if (tMin .lt. (minObsTimeInData)) then
          tMin = tMax - obsCont%segmentSize
          if (tMin .lt. (minObsTimeInData)) then
            dataLimitLHS = .true.
            newSegLHS = .false.
            call Warning('Failed to expand the beginning of the observer time array',&
                       ' with the time range available in the aperiodic data files.')
          end if
        end if
        if (newSegLHS) then
          call Notice('Adding a '//trim(realtostring(tMax-tMin))//' second segment of observer time to ',&
                      'the beginning of the current observer time range.')
        end if
      end if
      if (newSegRHS .and. .not.newSegLHS .and. .not.dataLimitRHS) then
        ! Create a new segment just after the max. observer time.
        tMin = obsCont%tMax
        tMax = obsCont%tMax + obsCont%nbNewSegments*obsCont%SegmentSize
        call Notice('Adding a '//trim(realtostring(tMax-tMin))//' second segment of observer time to ',&
                    'the end of the current observer time range.')
      end if
      ! This ensures an 'nt' with the same dt.
      nt = NINT((tMax-tMin)/(obsCont%tMax-obsCont%tMin)*(1.0*obsCont%nt-1.0))+1
      call SetObserverNt(obsCont%obs(i),nt)
      call SetObserverTMin(obsCont%obs(i),tMin)
      call SetObserverTMax(obsCont%obs(i),tMax)
      call SetObserverTitle(obsCont%obs(i),'Observer'//IntegerToString(i))
      call NullifyObserver(obsCont%obs(i))
    end if
  end do 
  if (newSegLHS .or. newSegRHS) then
    obsCont%nt = obsCont%obs(1)%nt+obsCont%obs(obsCont%nbObs)%nt - 1
    obsCont%tMin = minval(obsCont%obs%tmin)
    obsCont%tMax = maxval(obsCont%obs%tmax)
  end if
end subroutine InitializeObservers

!**
!SUBROUTINE WriteObsContResults
!  This is the top level routine that writes out all the
!  results that were calculated for the observer container
!**
subroutine WriteObsContResults(obsCont)
  implicit none
  type(observerContainer),intent(in)::obsCont
  integer:: i

  if(acousticPressureFlag)then
    if (debugLevel.gt.2) then
      write(*,*) 'Pressure path = ', trim(obsCont%pressurePath)
    end if
    call WriteObsContPressureData(obsCont)
  end if  
  !This writes out OASPLdB and OASPLdBA:
  if(OASPLdBFlag.or.OASPLdBAFlag)then  
    if (debugLevel.gt.2) then
      write(*,*) 'OASPL path = ', trim(obsCont%OASPLPath)
    end if
    call WriteObsContOASPLData(obsCont)
  end if  
  if((SPLdBFlag.or.SPLdBAFlag).and.((obsCont%nbFreqRanges.gt.0)))then
    if (debugLevel.gt.2) then
      write(*,*) 'Freq range SPL path = ', trim(obsCont%SPLPath)
    end if
    call WriteObsContFreqRangeData(obsCont)
  end if
  if(PNLFlag.or.PNLTFlag.or.EPNLFlag)then
    if (debugLevel.gt.2) then
      write(*,*) 'PNL path = ', trim(obsCont%SPLPath)
    end if
    call WriteObsContPNLData(obsCont)
  end if
  if (SELFlag .and. (all(SELdBDrop).or.forceSEL)) then
    if (debugLevel.gt.2) then
      write(*,*) 'SEL path = ', trim(obsCont%SPLPath)
    end if
    call WriteObsContSELData(obsCont)
  end if    
  if(spectrumFlag)then
    if (debugLevel.gt.2) then
      if(phaseFlag)then
        write(*,*) 'Phase path = ', trim(obsCont%phasePath)
      end if
      if(complexPressureFlag)then
        write(*,*) 'Complex pressure path = ', trim(obsCont%complexPressurePath)
      end if
      if(SPLdBFlag.or.SPLdBAFlag)then
        write(*,*) 'SPL path = ', trim(obsCont%SPLPath)
      end if
    end if
    if(narrowBandFlag) then !phaseFlag .or. complexPressureFlag .or. SPLdBFlag .or. SPLdBAFlag ) then  !ksb debug: 8/9/2017 - I only want the spectral data output if thesee are true - NOT when spectrumFlag is true.
      call WriteObsContSpectralData(obsCont)
    end if
  end if  

  if(FSCOutputFlag)then
    if (debugLevel.gt.2) then
      write(*,*) 'FSC output path = ', trim(obsCont%FSCOutputFileName)
    end if
    call WriteFSCOutput(obsCont)
  end if

  if(octaveFlag.and.(SPLdBFlag.or.SPLdBAFlag))then
    if (debugLevel.gt.2) then
      write(*,*) 'Octave path = ', trim(obsCont%SPLpath)
    end if
    call WriteObsContOctaveData(obsCont)
  end if

  if (AtmAtten) then
    if (debugLevel.gt.2) then
      write(*,*) 'RadDistance path = ', trim(obsCont%SPLPath)
    end if  
    call WriteObsContRadDist(obsCont)
  end if

  if( broadbandFlag ) then
    if (globalPeggNoiseFlag) then
      if (debugLevel.gt.2) then
        write(*,*) 'PeggSPL path = ', trim(obsCont%SPLPath)
      end if
      call WriteObsContBBmsp(obsCont)
    elseif (globalBPMNoiseFlag) then
      if (debugLevel.gt.2) then
        write(*,*) 'BPMSPL path = ', trim(obsCont%SPLPath)
      end if
      call WriteObsContBBmsp(obsCont)
    end if
  end if
end subroutine WriteObsContResults

!**
!SUBROUTINE WriteFSCOutput
!  Binary support: no
!**
subroutine WriteFSCOutput(obsCont)
  implicit none
  type(observerContainer),intent(in)::obsCont
  integer::unit,nsurf,i,j,obsNum,k,nbObs,nbSources,firstFreqIndex, nf
  character(len=4096),dimension(:),allocatable::titles
  integer,dimension(:),allocatable::ncolls
  real,dimension(:),allocatable::freqArray
  real,dimension(1)::dummyArray
  real::vxr,vyr,vzr,vxi,vyi,vzi,tempdenom
  nbSources = size(obsCont%obs(1)%realCompP,1)
  nbObs     = obsCont%nbObs
  allocate(freqArray(GetNf(obsCont%obs(1)%realCompP(1,1))))
  nf = size(freqArray)
  freqArray = obsCont%obs(1)%realCompP(1,1)%freq
  !If there was an FSCInputFileName specified in the namelist,
  !  the fsc observers were created from this file and this is 
  !  where we need to extract some data from: 
  if(trim(obsCont%FSCInputFileName).ne.'')then
    unit = GetStreamNumber()    
    !First we have to get some information from the input file that we did not store:
    open(unit=unit,file=trim(globalFolderName)//trim(obsCont%FSCInputFileName), &
         access='sequential', form='formatted', status='old')
    read(unit,*)nsurf
    allocate(titles(nsurf),ncolls(nsurf))
    do i=1,nsurf
      read(unit,*)titles(i)
      read(unit,*)ncolls(i)
      do j=1,ncolls(i)
        read(unit,*)
      end do
    end do
    close(unit)
  else !This is an observer field defined by parameters in the namelist
    allocate(titles(1),ncolls(1))
    nsurf     = 1
    titles(1) = 'Observer'
    ncolls(1) = obsCont%nbObs
    !We have to write out a observer point location file so the people
    !  with FSC know where the observer grid positions were. It is a 
    !  plot3D type file.
    dummyArray(1) = 1.
    call WriteMovingObsContSurface(obsCont,dummyArray,trim(globalFolderName)// &
                                   trim(obsCont%FSCOutputFileName)//'.x')
  end if    
  unit   = GetStreamNumber()
  obsNum = 0
  open(unit=unit,file=trim(globalFolderName)//trim(obsCont%FSCOutputFileName), &
       status='replace')

  if( freqArray(1)/=0.0 ) then
    firstFreqIndex = 1
  else
    firstFreqIndex = 2
  end if
  
  write(unit,'(A)')trim(IntegerToString(nf-firstFreqIndex+1))
  !This is necessary for correct formatting:
  do i=firstFreqIndex,nf-1
    write(unit,fmt='(2A)',advance='no')trim(RealToString(freqArray(i))),' '
  end do
  write(unit,'(A)')trim(RealToString(freqArray(nf-firstFreqIndex+1)))
  do i=1,nsurf
    write(unit,'(A)')trim(titles(i))
    write(unit,'(A)')trim(IntegerToString(ncolls(i)))
    do j=1,ncolls(i)
      obsNum=obsNum+1
      write(unit,'(A,A)')"'COLLOCATION POINT ',",trim(IntegerToString(j))      
      do k=firstFreqIndex,nf
        tempDenom = rho*2.0*pi*freqArray(k)
        if(FSCThicknessFlag)then
          vxr = -obsCont%obs(obsNum)%imagCompP(THICK_PGX,1)%f(k) / tempDenom
          vyr = -obsCont%obs(obsNum)%imagCompP(THICK_PGY,1)%f(k) / tempDenom
          vzr = -obsCont%obs(obsNum)%imagCompP(THICK_PGZ,1)%f(k) / tempDenom
          vxi = obsCont%obs(obsNum)%realCompP(THICK_PGX,1)%f(k)  / tempDenom
          vyi = obsCont%obs(obsNum)%realCompP(THICK_PGY,1)%f(k)  / tempDenom
          vzi = obsCont%obs(obsNum)%realCompP(THICK_PGZ,1)%f(k)  / tempDenom
          write(unit,'(8(ES14.7,1x))') obsCont%obs(obsNum)%realCompP(THICK_APTH,1)%f(k), &
                                 obsCont%obs(obsNum)%imagCompP(THICK_APTH,1)%f(k), &
                                 vxr, vxi, vyr, vyi, vzr, vzi
        else if(FSCLoadingFlag)then
          vxr = -obsCont%obs(obsNum)%imagCompP(LOAD_PGX,1)%f(k) / tempDenom
          vyr = -obsCont%obs(obsNum)%imagCompP(LOAD_PGY,1)%f(k) / tempDenom
          vzr = -obsCont%obs(obsNum)%imagCompP(LOAD_PGZ,1)%f(k) / tempDenom
          vxi = obsCont%obs(obsNum)%realCompP(LOAD_PGX,1)%f(k)  / tempDenom
          vyi = obsCont%obs(obsNum)%realCompP(LOAD_PGY,1)%f(k)  / tempDenom
          vzi = obsCont%obs(obsNum)%realCompP(LOAD_PGZ,1)%f(k)  / tempDenom
          write(unit,'(8(ES14.7,1x))') obsCont%obs(obsNum)%realCompP(LOAD_APTH,1)%f(k), &
                                 obsCont%obs(obsNum)%imagCompP(LOAD_APTH,1)%f(k), &
                                 vxr, vxi, vyr, vyi, vzr, vzi
        else !This is a total noise calculation
          vxr = -obsCont%obs(obsNum)%imagCompP(TOTAL_PGX,1)%f(k) / tempDenom
          vyr = -obsCont%obs(obsNum)%imagCompP(TOTAL_PGY,1)%f(k) / tempDenom
          vzr = -obsCont%obs(obsNum)%imagCompP(TOTAL_PGZ,1)%f(k) / tempDenom
          vxi = obsCont%obs(obsNum)%realCompP(TOTAL_PGX,1)%f(k)  / tempDenom
          vyi = obsCont%obs(obsNum)%realCompP(TOTAL_PGY,1)%f(k)  / tempDenom
          vzi = obsCont%obs(obsNum)%realCompP(TOTAL_PGZ,1)%f(k)  / tempDenom
          write(unit,'(8(ES14.7,1x))') obsCont%obs(obsNum)%realCompP(TOTAL_APTH,1)%f(k), &
                                 obsCont%obs(obsNum)%imagCompP(TOTAL_APTH,1)%f(k), &
                                 vxr, vxi, vyr, vyi, vzr, vzi
        end if
      end do                                          
    end do
  end do
  close(unit)
end subroutine WriteFSCOutput

!**
!SUBROUTINE WriteObsContOASPLData
!  Binary support: no
!**
subroutine WriteObsContOASPLData(obsCont)
  implicit none
  type(observerContainer),intent(in)::obsCont
  real,dimension(:,:,:),allocatable::OASPLdB,OASPLdBA
  character(len=4096),dimension(:),allocatable::titles
  integer::i,nbObs,nbSources,j
  nbObs=obsCont%nbObs
  if(nbObs==1)then
    if(OASPLdBFlag)then
      call WriteSingleObsTimeHistoryData(obsCont%obs(1)%OASPLdB,trim(obsCont%OASPLpath)//'dB', &
                                         ' OASPL (dB)', 'SPL ')
    end if
    if(OASPLdBAFlag)then
      call WriteSingleObsTimeHistoryData(obsCont%obs(1)%OASPLdBA,trim(obsCont%OASPLpath)//'dBA', &
                                         'OASPL (dBA)', 'SPL ')
    end if
  else
    !First write the grid file:
    if(getSize(obsCont%CBList).gt.0)then      
      call WriteMovingObsContSurface(obsCont,obsCont%segTimeArray,trim(obsCont%OASPLpath)//'.x')
    else    
      call WriteStationaryObsSurface(obsCont,size(obsCont%segTimeArray),&
                                       trim(obsCont%OASPLpath)//'.x')
    end if
    if(OASPLdBFlag)then
      nbSources = size(obsCont%obs(1)%OASPLdB)
      allocate(OASPLdB(nbObs,nbSources,size(obsCont%segTimeArray)), &
               titles(nbSources+1))
      titles(1)='SPL Time (s)'
      do i=1,nbObs
        do j=1,nbSources
          OASPLdB(i,j,:)= obsCont%obs(i)%OASPLdB(j)%f
        end do
      end do
      do j=1,nbSources
        titles(j+1) = GettHTitle(obsCont%obs(1)%OASPLdB(j))
      end do
      call WriteTimeHistoryFunctionFile(OASPLdB,titles,obsCont%segTimeArray, &
                                        obsCont%nbObsDim1,obsCont%nbObsDim2, &
                                        trim(obsCont%OASPLpath)//'dB')
      deallocate(OASPLdB,titles)
    end if
    if(OASPLdBAFlag)then
      nbSources = size(obsCont%obs(1)%OASPLdBA)
      allocate(OASPLdBA(nbObs,nbSources,size(obsCont%segTimeArray)), &
               titles(nbSources+1))
      titles(1)='SPL Time (s)'
      do i=1,nbObs
        do j=1,nbSources
          OASPLdBA(i,j,:)= obsCont%obs(i)%OASPLdBA(j)%f
        end do
      end do
      do j=1,nbSources
        titles(j+1) = GettHTitle(obsCont%obs(1)%OASPLdBA(j))
      end do
      call WriteTimeHistoryFunctionFile(OASPLdBA,titles,obsCont%segTimeArray, &
                                        obsCont%nbObsDim1,obsCont%nbObsDim2,  &
                                        trim(obsCont%OASPLpath)//'dBA')
      deallocate(OASPLdBA,titles)
    end if  
  end if
end subroutine WriteObsContOASPLData

!**
!SUBROUTINE WriteSpectralFunctionFile
!  Binary support: yes
!**
subroutine WriteSpectralFunctionFile(obsCont,funcData,freqArray,path)
  implicit none
  type(observerContainer),intent(in)::obsCont
  !funcData(nbObs,nbFunctions,nf)
  real,dimension(:,:,:),intent(in)::funcData
  character(len=*),intent(in)::path  
  real,dimension(:),intent(in)::freqArray
  integer::unitNum,nbFunctions,obsNum,it,k,nf,stat
  real(kind=4):: temp
  nbFunctions = size(funcData,2)
  nf          = size(funcData,3)  
  unitNum     = GetStreamNumber()
  if(ASCIIOutputFlag)then
    open(unitNum, file=trim(path), form='FORMATTED')  
      write(unitNum,'(4I12)')obsCont%NbObsDim1, &
                             obsCont%NbObsDim2, &
                             nf,                &
                             nbFunctions+1
      !The following nested do loops writes plot-3d structured
      !  format.  It increments obsNum, then "it", then k.
      write(unitNum,'(G15.6)') ((freqArray(it),&
                                obsNum =1, obsCont%nbObs), &
                                it=1, nf)        
      write(unitNum,'(G15.6)') (((funcData(obsNum,k,it),&
                                 obsNum =1, obsCont%nbObs), &
                                 it=1, nf), &
                                 k =1, nbFunctions)
    close(unitNum) 
  else  !It is a binary file
    call CreateBinaryFile (unitNum, trim(path), .false., stat)
    if (stat/=0) then
      call Error('Could not open the file '//trim(path))
    end if    
    call WriteBinaryInteger(unitNum,obsCont%nbObsDim1,stat)
    call WriteBinaryInteger(unitNum,obsCont%nbObsDim2,stat)
    call WriteBinaryInteger(unitNum,nf,stat)
    call WriteBinaryInteger(unitNum,nbFunctions + 1,stat)
    if (stat/=0) then
      call Error('Could not write to the file '//trim(path))
    end if
    do it = 1, nf
      temp = freqArray(it)
      do obsNum = 1,obsCont%nbObs
        call WriteBinaryReal(unitNum,temp,stat)
      end do
    end do
    if (stat/=0) then
      call Error('Could not write to the file '//trim(path))
    end if
    do k = 1,nbFunctions
      do it = 1, nf
        do obsNum = 1,obsCont%nbObs
          temp = funcData(obsNum,k,it)
          call WriteBinaryReal(unitNum,temp,stat)
        end do
      end do
    end do
    if (stat/=0) then
      call Error('Could not write to the file '//trim(path))
    end if 
    call CloseBinaryFile(unitNum)     
  end if      
end subroutine WriteSpectralFunctionFile

!**
!SUBROUTINE WriteTimeHistoryFunctionFile
!  Binary support: yes
!**
subroutine WriteTimeHistoryFunctionFile(func,fTitles,time, &
                                        nbObsDim1,nbObsDim2,path)
  implicit none
  character(len=*),intent(in)::path
  integer::nt,unitNum,obsNum,it,k,nbSources,nbObs,nbObsDim1,nbObsDim2,stat
  !func(obsNum,noiseSource,nt)
  real(kind=4):: temp
  real,dimension(:,:,:),intent(in)::func
  !time(nt)
  real,dimension(:),intent(in)::time
  character(len=4096),dimension(:),intent(in)::fTitles
  unitNum   = GetStreamNumber()
  nbSources = size(func,2)
  nt        = size(func,3)
  nbObs     = size(func,1)
  
  if (ASCIIOutputFlag) then
    open(unitNum, file=trim(path)//'.fn', form='FORMATTED')
    write(unitNum,'(4I12)')nbObsDim1,             &
                           nbObsDim2,             &
                           nt,                    &
                           nbSources+1
    !The following nested do loops write plot-3d structured
    !  format.  It increments obsNum, then 'it', then k.    6
    write(unitNum,'(G15.6)')((time(it),           &
                              obsNum = 1, nbObs), &
                              it=1, nt)    
    write(unitNum,'(G15.6)')(((func(obsNum,k,it), &
                               obsNum = 1, nbObs),&
                               it=1, nt),         &
                               k =1, nbSources)
    close(unitNum)                               
    else  !It is a binary file
      call CreateBinaryFile (unitNum, trim(path)//'.fn', .false., stat)
      if (stat/=0) then
        call Error('Could not open the file '//trim(path))
      end if    
      call WriteBinaryInteger(unitNum,nbObsDim1,stat)
      call WriteBinaryInteger(unitNum,nbObsDim2,stat)
      call WriteBinaryInteger(unitNum,nt,stat)
      call WriteBinaryInteger(unitNum,nbSources + 1,stat)
      if (stat/=0) then
        call Error('Could not write to the file '//trim(path))
      end if
      do it = 1, nt
        temp = time(it)
        do obsNum = 1,nbObs
          call WriteBinaryReal(unitNum,temp,stat)
        end do
      end do
      if (stat/=0) then
        call Error('Could not write to the file '//trim(path))
      end if
      do k = 1,nbSources
        do it = 1, nt
          do obsNum = 1,nbObs
            temp = func(obsNum,k,it)          
            call WriteBinaryReal(unitNum,temp,stat)
          end do
        end do
      end do
      if (stat/=0) then
        call Error('Could not write to the file '//trim(path))
      end if 
      call CloseBinaryFile(unitNum)     
  end if            
  call WriteNameFile(fTitles,trim(path)//'.nam')                                                              
end subroutine WriteTimeHistoryFunctionFile

!**
!SUBROUTINE WriteMovingObsContSurface
!  Binary support: yes
!**
subroutine WriteMovingObsContSurface(obsCont,array,path)
  implicit none
  type(observerContainer),intent(in)::obsCont
  integer::nPoints,it,i,j,obsNum,unit,k,stat
  character(len=*),intent(in)::path
  !We need a temporary coordinate array that records the position in time
  type(vector), dimension(:,:,:), allocatable::tempCoordinates
  !This is a transformation matrix that will give us the position of any vector in time.
  type(matrix4)::transformationMatrix
  real(kind=4):: temp
  real,dimension(:),intent(in)::array
  nPoints=size(array)
  allocate(tempCoordinates(obsCont%NbObsDim1,obsCont%NbObsDim2,nPoints))
  !This do loop calculates the transformation matrix for each time step and records
  ! where each point on the oberver is at that time.
  do it=1, nPoints
    transformationMatrix=observerTransformationMatrix(obsCont%CBList,array(it))
    obsNum=0
    do j=1, obsCont%NbObsDim2
      do i=1, obsCont%NbObsDim1
        obsNum=ObsNum+1
        tempCoordinates(i,j,it)=transformationMatrix * &
                                obsCont%obs(obsNum)%coordinates
      end do
    end do
  end do
  !Open the grid file and write out the position information for the observer surface
  unit=GetStreamNumber()
  if (ASCIIOutputFlag) then
    open(unit, file=trim(path), form='FORMATTED')
!    write(unit,'(1x,"1")')  !ksb - added to make this a "multigrid" file, needed by FSC.
    write(unit,'(3I12)')   size(tempCoordinates, 1),  &
                           size(tempCoordinates, 2),  &
                           size(tempCoordinates, 3)
    !This increments i, then j, then it, then k
    write(unit,'(G15.6)') ((((tempCoordinates(i,j,it)%A(k),          &
                       i =1, size(tempCoordinates, 1)), &
                       j =1, size(tempCoordinates, 2)), &
                       it=1, size(tempCoordinates, 3)), &
                       k =1, 3)
      close(unit)
    else
      call CreateBinaryFile (unit, trim(path), .false., stat)
      if (stat/=0) then
        call Error('Could not open the file '//trim(path))
      end if 
      call WriteBinaryInteger (unit, size(tempCoordinates, 1), stat)
      call WriteBinaryInteger (unit, size(tempCoordinates, 2), stat)
      call WriteBinaryInteger (unit, size(tempCoordinates, 3), stat)
      if (stat/=0) then
        call Error('Could not write to the file '//trim(path))
      end if
      do k=1, 3
        do it=1, size(tempCoordinates, 3)
          do j=1, size(tempCoordinates, 2)
            do i=1, size(tempCoordinates, 1)
              temp = tempCoordinates(i,j,it)%A(k)
              call WriteBinaryReal (unit, temp, stat)
              if (stat/=0) then
                call Error('Could not write to the file '//trim(path))
              end if
            end do
          end do
        end do
      end do
      call CloseBinaryFile(unit)
  end if
  !deallocate the temporary array
  if (allocated(tempCoordinates)) deallocate(tempCoordinates)  
end subroutine WriteMovingObsContSurface

!**
!SUBROUTINE WriteStationaryObsSurface
!  Binary support: yes
!**
subroutine WriteStationaryObsSurface(obsCont,nPoints,path)
  implicit none 
  integer,intent(in)::nPoints
  type(observerContainer),intent(in)::obsCont
  character(len=*),intent(in)::path
  integer::nDim1,nDim2,it,k,unitNum,obsNum,stat
  real(kind=4):: temp
  nDim1 = obsCont%NbObsDim1
  nDim2 = obsCont%NbObsDim2
  unitNum = GetStreamNumber()     
  if (ASCIIOutputFlag) then
    open(unitNum, file=trim(path), form='FORMATTED')
    write(unitNum,'(3I12)')    nDim1,nDim2,nPoints
    !The following nested do loops writes plot-3d structured
    !  format.  It increments obsNum, then "it", then k.
    do k=1,3
      do it=1,nPoints
        do obsNum=1,obsCont%nbObs
          write(unitNum,'(G15.6)')obsCont%obs(obsNum)%coordinates%A(k)
        end do
      end do
    end do
    close(unitNum)
    else
      call CreateBinaryFile (unitNum, trim(path), .false., stat)
      if (stat/=0) then
        call Error('Could not open the file '//trim(path))
      end if
      call WriteBinaryInteger (unitNum, nDim1,   stat)
      call WriteBinaryInteger (unitNum, nDim2,   stat)
      call WriteBinaryInteger (unitNum, nPoints, stat)    
      if (stat/=0) then
        call Error('Could not write to the file '//trim(path))
      end if
      do k = 1,3
        do it = 1,nPoints
          do obsNum = 1,obsCont%nbObs
            temp = obsCont%obs(obsNum)%coordinates%A(k)
            call WriteBinaryReal(unitNum,temp,stat)
            if (stat/=0) then
              call Error('Could not write to the file '//trim(path))
            end if  
          end do
        end do
      end do
      call CloseBinaryFile(unitNum)      
  end if
end subroutine WriteStationaryObsSurface

!**
!SUBROUTINE WriteObsContSigmaSurface
!  Binary support: yes
!**
!  This routine writes out a sigma surface for an observer container.
!CALLS:
! - WriteBinaryReal3DArray (linux/windowsIO.f90)
!**
subroutine WriteObsContSigmaSurface(obsCont, nFunc, gridUnit, dataUnit)
  type(observerContainer), intent(in)::obsCont
  integer, intent(in):: gridUnit, dataUnit, nFunc
  integer::i, j, k, iTime, stat, iMax,jMax,nt
  real, dimension(:,:,:,:), allocatable::obsSigma
  real, dimension(:), allocatable::time
  type(matrix4)::transformMatrix
  type(vector)::tempVector
  real(kind=4), dimension(:,:,:), allocatable:: temp3DArray
  iMax = obsCont%NbObsDim1
  jMax = obsCont%NbObsDim2
  nt   = obsCont%obs(1)%nt
  allocate(time(nt))
  time = CreateEvenlySpacedArray(obsCont%obs(1)%tMin,obsCont%obs(1)%tMax,nt)
  if (iMax < 1 .or. jMax < 1 .or. nt < 1 ) then 
    if (debugLevel >=1) write(*,*) "WARNING: No data calculated on this observer."
    return
  end if
  allocate(obsSigma(nFunc, obsCont%NbObsDim1, obsCont%NbObsDim2, obsCont%obs(1)%nt))
  allocate(temp3DArray(obsCont%nbObsDim1, obsCont%NbObsDim2, obsCont%obs(1)%nt))
  !nFunc - 1,2,3 - x,y,z location
  !nFunc - 4 - source time
  !nFunc - 5 - observer time
  !nFunc - 6,7,8,etc. - the actual functions
  obsSigma=0.0
  do iTime=1, nt
    transformMatrix=observerTransformationMatrix(obsCont%CBList, time(iTime))
    k=1
    do j=1, jMax
      do i=1, iMax
        tempVector=transformMatrix*obsCont%obs(k)%coordinates
        obsSigma(1:3, i, j, iTime) = tempVector%A
        obsSigma(5,   i, j, iTime) = time(iTime)
        !At the observer location, source time and observer time are equal,
        !becuase if a source was placed exactly at the observer location,
        !the travel time of the wave would be equal to zero.
        obsSigma(4,   i, j, iTime) = time(iTime)
        k=k+1 
      end do
    end do
  end do
  call Message ("Writing sigma data for "//trim(obsCont%title))
  do k=1,3
    temp3DArray = obsSigma(k,:,:,:)
    call WriteBinaryReal3DArray (gridUnit, temp3DArray, stat)
    if (stat /= 0) then
      call Error("Error writing observer sigma: out of space.")
    end if
  end do
  do k=4,nFunc
    temp3DArray = obsSigma(k,:,:,:)  
    call WriteBinaryReal3DArray (dataUnit, temp3DArray, stat)
    if (stat /= 0) then
      call Error("Error writing observer sigma: out of space.")
    end if
  end do
  deallocate(obsSigma, temp3DArray)  
end subroutine WriteObsContSigmaSurface

subroutine WriteMultipleObsFreqDomainData(obsCont,fD1,fD2,segInc,path)
  implicit none
  type(observerContainer),intent(in)::obsCont
  !fD1(obsNum,noiseSource,segNum)
  type(freqDomain),dimension(:,:,:),intent(in),optional::fD1,fD2
  integer,intent(in)::segInc
  character(len=*),intent(in)::path
  integer::i,k,segNum,nbSegs,start,finish,nbFuncs,nf
  real, dimension(:,:,:),allocatable::f  
  character(len=4096),dimension(:),allocatable::fTitles
  real,dimension(:),allocatable::freqArray    
  if(present(fD1))then
    nf        = GetNf(fD1(1,1,1))
    allocate(freqArray(nf))
    freqArray = fD1(1,1,1)%freq
    nbSegs    = size(fD1,3)    
  else if(present(fD2))then
    nf        = GetNf(fD2(1,1,1))
    allocate(freqArray(nf))
    freqArray = fD2(1,1,1)%freq
    nbSegs    = size(fD2,3)
  else
    call Warning('WriteSingleObsFreqDomainData called with no data.')
  end if  
  if(present(fD1).and.present(fD2))then
    nbFuncs = size(fD1,2) + size(fD2,2)
    allocate(fTitles(nbFuncs+1))
    fTitles(1)='Frequency'
    do i=1,nbFuncs/2
      fTitles(i+1) = trim(GetFreqDomainTitle(fD1(1,i,1)))
    end do
    do i=1,size(fD2,2)
      fTitles(i+size(fD1,2)+1) = trim(GetFreqDomainTitle(fD2(1,i,1)))
    end do
    else if(present(fD1))then
      nbFuncs = size(fD1,2)
      allocate(fTitles(nbFuncs+1))
      fTitles(1) = 'Frequency'
      do i=1,size(fD1,2)
        fTitles(i+1) = GetFreqDomainTitle(fD1(1,i,1))
      end do 
      else
        nbFuncs = size(fD2,2)
        allocate(fTitles(nbFuncs+1))
        fTitles(1) = 'Frequency'
        do i=1,size(fD2,2)
          fTitles(i+1) = GetFreqDomainTitle(fD2(1,i,1))
        end do
  end if   
  if(segInc.ne.0)then
    start = segInc
    finish = segInc
    else
      start  = 1
      finish = nbSegs
  end if  
  allocate(f(obsCont%nbObs,nbFuncs,nf))
  do segNum = start,finish
    do i=1,obsCont%nbObs
      if(present(fD1))then
        do k=1,size(fD1,2)
          f(i,k,:) = fD1(i,k,segNum)%f
        end do
        k = size(fD1,2)
      end if
      if(present(fD2))then
        do k = k+1,nbFuncs
          f(i,k,:) = fD2(i,k-size(fD1,2),segNum)%f
        end do
      end if
    end do
    if(obsCont%nbSegments==1.or.obsCont%nbSegments==0)then
      call WriteSpectralFunctionFile(obsCont,f,freqArray,trim(path)//'.fn')
    else
      call WriteSpectralFunctionFile(obsCont,f,freqArray,trim(path)//'_segment'// &
                                   trim(IntegerToString(segNum))//'.fn')
    end if
  end do  
  call WriteNameFile(fTitles,trim(path)//'.nam')  
  deallocate(f,fTitles,freqArray)    
end subroutine WriteMultipleObsFreqDomainData

!**
!SUBROUTINE WriteMultipleObsPressureData
!  This routine writes out pressure time histories
!  for an observer container with multiple observers.
!**
subroutine WriteMultipleObsPressureData(obsCont)
  implicit none
  type(observerContainer),intent(in)::obsCont
  real,dimension(:),allocatable::time
  real,dimension(:,:,:),allocatable::func
  character(len=4096)::path
  character(len=4096),dimension(:),allocatable::funcTitles
  integer::i,j    
  path = trim(obsCont%pressurePath)
  allocate(time(obsCont%nt))   
  time = CreateEvenlySpacedArray(obsCont%tMin,obsCont%tMax,obsCont%nt)
  !First write the grid file:
  if(getSize(obsCont%CBList).gt.0)then      
    call WriteMovingObsContSurface(obsCont,time,trim(path)//'.x')
  else    
    call WriteStationaryObsSurface(obsCont,obsCont%nt,trim(path)//'.x')
  end if
  !Now write the function file:
  allocate(func(obsCont%nbObs,GetNewSizeOfPPrime(),obsCont%nt), &
           funcTitles(GetNewSizeOfPPrime()+1))           
  do i=1,obsCont%nbObs
    do j=1,GetNewSizeOfPPrime()
      func(i,j,:) = obsCont%obs(i)%pPrime(j)%f
    end do
  end do
  do j=1,GetNewSizeOfPPrime()
    funcTitles(j+1) = trim(GettHTitle(obsCont%obs(1)%pPrime(j)))//'Acoustic Pressure'
  end do
  funcTitles(1) = 'Observer Time (s)'
  call WriteTimeHistoryFunctionFile(func,funcTitles,time,obsCont%nbObsDim1, &
                                    obsCont%nbObsDim2,trim(path))
  deallocate(func,time,funcTitles)
end subroutine WriteMultipleObsPressureData

subroutine WriteObsContFreqDomainData(obsCont,fD1,fD2,path)
  implicit none
  type(observerContainer),intent(in)::obsCont
  !fD1(obsNum,noiseSource,segNum)
  !fD2(obsNum,noiseSource,segNum)  
  type(freqDomain),dimension(:,:,:),intent(in),optional::fD1,fD2
  character(len=*),intent(in)::path
  integer::segInc  
  segInc=obsCont%segmentIncrement
  if(obsCont%nbObs==1)then
    if(present(fD1).and.present(fD2))then
      call WriteSingleObsFreqDomainData(fD1(1,:,:),fD2(1,:,:),path=path)
    else
      call WriteSingleObsFreqDomainData(fD1(1,:,:),path=path)
    end if
  else
    if(present(fD1).and.present(fD2))then
      call WriteMultipleObsFreqDomainData(obsCont,fD1,fD2,segInc,path)
    else
      call WriteMultipleObsFreqDomainData(obsCont,fD1,segInc=segInc,path=path)
    end if
  end if  
end subroutine WriteObsContFreqDomainData

subroutine WriteObsContPressureData(obsCont)
  implicit none
  type(observerContainer),intent(in)::obsCont  
  if(obsCont%nbObs==1) then
    call WriteSingleObsTimeHistoryData(obsCont%obs(1)%pPrime(1:GetNewSizeOfPPrime()),&
                                       trim(obsCont%pressurePath),                   &
                                       ' Acoustic Pressure (Pa)', 'Observer ')
  else  
    call WriteMultipleObsPressureData(obsCont)
  end if
end subroutine WriteObsContPressureData

!**
!SUBROUTINE WriteObsContPNLData
!  This routine writes out PNL and tonal corrected PNL (PNLT)
!**
subroutine WriteObsContPNLData(obsCont)
  implicit none
  type(observerContainer),intent(in)::obsCont
  real,dimension(:,:,:),allocatable::PNL,PNLT
  character(len=4096),dimension(:),allocatable::titles
  character(len=1024):: BBtitle
  integer::i,j,k,nbSources,unitNum,obsNum

  if (globalPeggNoiseFlag) then
    BBtitle = 'Broadband-Pegg'
  else if (globalBPMNoiseFlag) then
    BBtitle = 'Broadband-BPM'
  end if 
  
  if(PNLFlag)then
    if(obsCont%nbObs==1)then
      call WriteSingleObsTimeHistoryData(obsCont%obs(1)%PNL,trim(obsCont%SPLpath)//'_PNL', &
                                        ' PNLdB', 'SPL ')
    else
      !First write the grid file:
      if(getSize(obsCont%CBList).gt.0)then      
        call WriteMovingObsContSurface(obsCont,obsCont%segTimeArray, &
                                       trim(obsCont%SPLpath)//'_PNL.x')
        else    
          call WriteStationaryObsSurface(obsCont,size(obsCont%segTimeArray), &
                                         trim(obsCont%SPLpath)//'_PNL.x')
      end if
      nbSources = size(obsCont%obs(1)%PNL)
      allocate(PNL(obsCont%nbObs,nbSources,size(obsCont%segTimeArray)), &
               titles(nbSources+1))
      titles(1)='Time'
      do i=1, obsCont%nbObs
        do j=1,nbSources
          PNL(i,j,:)= obsCont%obs(i)%PNL(j)%f
        end do
      end do
      do j=1,nbSources
        titles(j+1) = GettHTitle(obsCont%obs(1)%PNL(j))
      end do
      call WriteTimeHistoryFunctionFile(PNL,titles,obsCont%segTimeArray, &
                                        obsCont%nbObsDim1,obsCont%nbObsDim2,&
                                        trim(obsCont%SPLpath)//'_PNLdB')
      deallocate(PNL,titles)
    end if
  end if
  !This is for writing out tonal corrected PNL:
  if(PNLTFlag)then
    if(obsCont%nbObs==1)then
      call WriteSingleObsTimeHistoryData(obsCont%obs(1)%PNLT,trim(obsCont%SPLpath)//'_PNLT', &
                                         ' PNLTdB', 'SPL ')
    else
      !First write the grid file:
      if(getSize(obsCont%CBList).gt.0)then      
        call WriteMovingObsContSurface(obsCont,obsCont%segTimeArray, &
                                       trim(obsCont%SPLpath)//'_PNLT.x')
        else    
          call WriteStationaryObsSurface(obsCont,size(obsCont%segTimeArray), &
                                         trim(obsCont%SPLpath)//'_PNLT.x')
      end if
      nbSources = size(obsCont%obs(1)%PNLT)
      allocate(PNLT(obsCont%nbObs,nbSources,size(obsCont%segTimeArray)), &
               titles(nbSources+1))
      titles(1)='Time'

      do i=1, obsCont%nbObs
        do j=1,nbSources
          PNLT(i,j,:)= obsCont%obs(i)%PNLT(j)%f
        end do
      end do
      do j=1,nbSources
        titles(j+1) = GettHTitle(obsCont%obs(1)%PNLT(j))
      end do
   
      call WriteTimeHistoryFunctionFile(PNLT,titles,obsCont%segTimeArray, &
                                        obsCont%nbObsDim1,obsCont%nbObsDim2, &
                                        trim(obsCont%SPLpath)//'_PNLTdB')
      deallocate(PNLT,titles)
    end if
  end if
  if(EPNLFlag .and. (all(PNLTdBDrop).or.forceEPNL))then
    nbSources = size(obsCont%obs(1)%EPNL)
    if(obsCont%nbObs==1)then
      unitNum = GetStreamNumber()
      open(unit=unitNum,file=trim(obsCont%SPLpath)//'_EPNL.txt',status='replace')
        do j=1,nbSources
          if (.not.associated(obsCont%obs(1)%mspBB)) then
            write(unitNum,*)trim(PPRIMETITLEARRAY(j))//' EPNLdB =',obsCont%obs(1)%EPNL(j)
          else 
            write(unitNum,*)trim(InsertBBTitle(obsCont%obs(1),BBtitle,j,nbSources))//' EPNLdB =', obsCont%obs(1)%EPNL(j)
          end if
        end do                      
      close(unitNum)
    else
      allocate(titles(nbSources))
      call WriteStationaryObsSurface(obsCont,1,trim(obsCont%SPLpath)//'_EPNL.x')
      unitNum = GetStreamNumber()
      open(unit=unitNum,file=trim(obsCont%SPLpath)//'_EPNLdB.fn',status='replace')
      write(unitNum,'(4I12)')obsCont%nbObsDim1,obsCont%nbObsDim2,1,nbSources
      write(unitNum,'(G15.6)')((obsCont%obs(obsNum)%EPNL(k), &
                                 obsNum = 1, obsCont%nbObs),&
                                 k = 1, nbSources)
      close(unitNum)
      !do j=1,nbSources
      !  if (.not. associated(obsCont%obs(1)%mspBB)) then
      !    titles(j) = trim(PPRIMETITLEARRAY(j))//' EPNLdB'
      !    !ksb debug:
      !    print*,'sourceNb=',j,', line A'
      !  else 
      !    titles(j) = trim(InsertBBTitle(obsCont%obs(1),BBtitle,j,nbSources))//' EPNLdB'
      !    print*,'sourceNb=',j,', line B'
      !  end if
      !end do
      !ksb debug: 5/29/2015 - problem with the above namefile generation in parallel.  The loop below only works if the
      !PNLT Flag is also true.  !!ksb debug: Check this for BB inclusion
      do j=1,nbSources
        titles(j) = GettHTitle(obsCont%obs(1)%PNLT(j))
        titles(j) = trim(titles(j))//', EPNLdB'
      end do
      call WriteNameFile(titles,trim(obsCont%SPLpath)//'_EPNLdB.nam')
      deallocate(titles)
    end if                                  
  end if
end subroutine WriteObsContPNLData


subroutine WriteObsContSELData(obsCont)
  implicit none
  type(observerContainer),intent(in)::obsCont
  real,dimension(:,:,:),allocatable::PNL,PNLT
  character(len=4096),dimension(:),allocatable::titles
  character(len=1024):: BBTitle
  integer::i,j,k,nbSources,unitNum,obsNum
 
  if (globalPeggNoiseFlag) then
    BBtitle = 'Broadband-Pegg'
  else if (globalBPMNoiseFlag) then
    BBtitle = 'Broadband-BPM'
  end if  

  nbSources = size(obsCont%obs(1)%SEL)  
  allocate(titles(nbSources))
  !do j=1,nbSources
  !  if (.not.associated(obsCont%obs(1)%mspBB)) then
  !    titles(j) = trim(PPRIMETITLEARRAY(j))//' SELdBA'
  !  else 
  !    titles(j) = trim(InsertBBTitle(obsCont%obs(1),BBtitle,j,nbSources))//' SELdBA'
  !  end if
  !end do
  !ksb debug: 5/29/2015 - problem with the above namefile generation in parallel.  The loop below only works if the
  !PNLT Flag is also true.  !!ksb debug: Check this out for BB inclusion.
  do j=1,nbSources
    titles(j) = GettHTitle(obsCont%obs(1)%PNLT(j))
    titles(j) = trim(titles(j))//', SELdBA'
  end do

  if(obsCont%nbObs==1)then
    unitNum = GetStreamNumber()
    open(unit=unitNum,file=trim(obsCont%SPLpath)//'_SEL.txt',status='replace')
      do j=1,nbSources
        write(unitNum,*)trim(titles(j))//' =',obsCont%obs(1)%SEL(j)
      end do
    close(unitNum)
  else
    call WriteStationaryObsSurface(obsCont,1,trim(obsCont%SPLpath)//'_SEL.x')
    unitNum = GetStreamNumber()
    open(unit=unitNum,file=trim(obsCont%SPLpath)//'_SELdBA.fn',status='replace')
    write(unitNum,'(4I12)')obsCont%nbObsDim1,obsCont%nbObsDim2,1,nbSources
    write(unitNum,'(G15.6)')((obsCont%obs(obsNum)%SEL(k), &
                               obsNum = 1, obsCont%nbObs),&
                               k = 1, nbSources)
    close(unitNum)
    call WriteNameFile(titles,trim(obsCont%SPLpath)//'_SELdBA.nam')
    deallocate(titles)
  end if
end subroutine WriteObsContSELData


!**
!SUBROUTINE WriteObsContFreqRangeData
!  This routine writes out SPL/SPLdBA data versus segment time
!  for each frequency range specified in the namelist.
!**
subroutine WriteObsContFreqRangeData(obsCont)
  implicit none
  type(observerContainer),intent(in)::obsCont
  real,dimension(:,:,:),allocatable::SPLdB,SPLdBA
  character(len=4096),dimension(:),allocatable::titles
  integer::nbSources,i,k,j  
  if(obsCont%nbObs==1)then
    call WriteSingleObsFreqRangeData(obsCont%obs(1),obsCont%fRangeTitle,trim(obsCont%SPLpath))
  else    
    !First write the grid file:
    if(getSize(obsCont%CBList).gt.0)then      
      call WriteMovingObsContSurface(obsCont,obsCont%segTimeArray, &
                                     trim(obsCont%SPLpath)//'_freqRanges.x')
    else    
      call WriteStationaryObsSurface(obsCont,size(obsCont%segTimeArray),&
                                     trim(obsCont%SPLpath)//'_freqRanges.x')
    end if
    do k=1,obsCont%nbFreqRanges
      if(SPLdBFlag)then
        nbSources = size(obsCont%obs(1)%dBRanges,1)
        allocate(SPLdB(obsCont%nbObs,nbSources,size(obsCont%segTimeArray)), &
                 titles(nbSources+1))
        titles(1)='SPL Time (s)'
        do i=1,obsCont%nbObs
          do j=1,nbSources
            SPLdB(i,j,:)= obsCont%obs(i)%dBRanges(j,k)%f
          end do
        end do
        do j=1,nbSources
          titles(j+1) = GettHTitle(obsCont%obs(1)%dBRanges(j,k))
        end do
        call WriteTimeHistoryFunctionFile(SPLdB,titles,obsCont%segTimeArray, &
                                        obsCont%nbObsDim1,obsCont%nbObsDim2, &
                                        trim(obsCont%SPLpath)//'_'//         &
                                        trim(obsCont%fRangeTitle(k))//'_dB')
        deallocate(SPLdB,titles)
      end if
      if(SPLdBAFlag)then
        nbSources = size(obsCont%obs(1)%dBARanges,1)
        allocate(SPLdBA(obsCont%nbObs,nbSources,size(obsCont%segTimeArray)), &
                 titles(nbSources+1))
        titles(1)='SPL Time (s)'
        do i=1,obsCont%nbObs
          do j=1,nbSources
            SPLdBA(i,j,:)= obsCont%obs(i)%dBARanges(j,k)%f
          end do
        end do
        do j=1,nbSources
          titles(j+1) = GettHTitle(obsCont%obs(1)%dBARanges(j,k))
        end do
        call WriteTimeHistoryFunctionFile(SPLdBA,titles,obsCont%segTimeArray,  &
                                          obsCont%nbObsDim1,obsCont%nbObsDim2, &
                                          trim(obsCont%SPLpath)//'_'//         &
                                          trim(obsCont%fRangeTitle(k))//'_dBA')
        deallocate(SPLdBA,titles)
      end if      
    end do
  end if
end subroutine WriteObsContFreqRangeData

subroutine WriteObsContOctaveData(obsCont)
  implicit none
  type(observerContainer),intent(in)::obsCont
  type(freqDomain),dimension(:,:,:),allocatable::SPLdB,SPLdBA
  integer::nbObs,nbSources,nbSegs,i,j,k,nf
  nbObs     = obsCont%nbObs
  if(SPLdBFlag)then
    nf = GetNf(obsCont%obs(1)%octaveFiltdB(1,1))
    else if(SPLFlag)then
    nf = GetNf(obsCont%obs(1)%octaveFiltdBA(1,1))
  end if   
  if(obsCont%nbObs.gt.1)then
    call WriteStationaryObsSurface(obsCont,nf, &
                 trim(obsCont%SPLpath)//'_octFilt_spectrum.x')    
   end if
  if(SPLdBFlag.and.SPLdBAFlag)then
    nbSources = size(obsCont%obs(1)%octaveFiltdB,1)  
    nbSegs    = size(obsCont%obs(1)%octaveFiltdB,2)
    allocate(SPLdB(nbObs,nbSources,nbSegs), &
             SPLdBA(nbObs,nbSources,nbSegs))
    do i=1,nbObs
      do j=1,nbSources
        do k=1,nbSegs
          nullify(SPLdB(i,j,k)%f,SPLdB(i,j,k)%freq,  &
                  SPLdBA(i,j,k)%f,SPLdBA(i,j,k)%freq) 
          call CopyFreqDomain(obsCont%obs(i)%octaveFiltdB(j,k), SPLdB(i,j,k))
          call CopyFreqDomain(obsCont%obs(i)%octaveFiltdBA(j,k), SPLdBA(i,j,k))
        end do
      end do
    end do
    call WriteObsContFreqDomainData(obsCont,SPLdB,SPLdBA, &
                                    trim(obsCont%SPLpath)//'_octFilt_spectrum')
    deallocate(SPLdB,SPLdBA)
  else if(SPLdBFlag)then
    nbSources = size(obsCont%obs(1)%octaveFiltdB,1)    
    nbSegs    = size(obsCont%obs(1)%octaveFiltdB,2)
    allocate(SPLdB(nbObs,nbSources,nbSegs))
    do i=1,nbObs
      do j=1,nbSources
        do k=1,nbSegs
          nullify(SPLdB(i,j,k)%f,SPLdB(i,j,k)%freq)
          call CopyFreqDomain(obsCont%obs(i)%octaveFiltdB(j,k), SPLdB(i,j,k))
        end do
      end do
    end do
    call WriteObsContFreqDomainData(obsCont,SPLdB, &
                                    path=trim(obsCont%SPLpath)//'_octFilt_spectrum')
    deallocate(SPLdB)
  else if(SPLdBAFlag)then
    nbSources = size(obsCont%obs(1)%octaveFiltdBA,1)
    nbSegs    = size(obsCont%obs(1)%octaveFiltdBA,2)
    allocate(SPLdBA(nbObs,nbSources,nbSegs))
    do i=1,nbObs
      do j=1,nbSources
        do k=1,nbSegs
          nullify(SPLdBA(i,j,k)%f,SPLdBA(i,j,k)%freq)
          call CopyFreqDomain(obsCont%obs(i)%octaveFiltdBA(j,k), SPLdBA(i,j,k))
        end do
      end do
    end do
    call WriteObsContFreqDomainData(obsCont,SPLdBA, &
                                    path=trim(obsCont%SPLpath)//'_octFilt_spectrum')
    deallocate(SPLdBA)
  end if  
end subroutine WriteObsContOctaveData


subroutine WriteObsContSpectralData(obsCont)
  implicit none
  type(observerContainer),intent(in)::obsCont
  type(freqDomain),dimension(:,:,:),allocatable::SPLdB,SPLdBA,phase, &
                                                 imagCompP,realCompP
  integer::nbObs,nbSources,nbSegs,i,j,k,nf
  nbObs     = obsCont%nbObs
  nbSegs    = obsCont%nbSegments  
  if(SPLdBFlag)then
    nf = GetNf(obsCont%obs(1)%dB(1,1))
  else if(SPLdBAFlag)then
    nf = GetNf(obsCont%obs(1)%dBA(1,1))
  else if(phaseFlag)then
    nf = GetNf(obsCont%obs(1)%phase(1,1))
  else if(complexPressureFlag)then
    nf = GetNf(obsCont%obs(1)%realCompP(1,1))
  else
    call Error('WriteObsContSpectralData called with unknown data', &
               'to be written out. Stopping.')
  end if
  if(obsCont%nbObs.gt.1)then
    call WriteStationaryObsSurface(obsCont,nf, &
                 trim(obsCont%SPLpath)//'_spectrum.x')
  end if
  !Write the SPLdB and/or SPLdBA data:
  if(SPLdBFlag.and.SPLdBAFlag)then
    nbSources = size(obsCont%obs(1)%dB,1)
    nbSegs    = size(obsCont%obs(1)%dB,2)
    allocate(SPLdB(nbObs,nbSources,nbSegs), &
             SPLdBA(nbObs,nbSources,nbSegs))
    do i=1,nbObs
      do j=1,nbSources
        do k=1,nbSegs
          nullify(SPLdB(i,j,k)%freq)
          nullify(SPLdB(i,j,k)%f)
          nullify(SPLdBA(i,j,k)%freq)
          nullify(SPLdBA(i,j,k)%f)
          call CopyFreqDomain(obsCont%obs(i)%dB(j,k), SPLdB(i,j,k))
          call CopyFreqDomain(obsCont%obs(i)%dBA(j,k), SPLdBA(i,j,k))
        end do
      end do
    end do
    call WriteObsContFreqDomainData(obsCont,SPLdB,SPLdBA, &
                                     trim(obsCont%SPLpath)//'_spectrum')
    deallocate(SPLdB,SPLdBA)
  else if(SPLdBFlag)then
    nbSources = size(obsCont%obs(1)%dB,1)
    nbSegs    = size(obsCont%obs(1)%dB,2)
    allocate(SPLdB(nbObs,nbSources,nbSegs))
    do i=1,nbObs
      do j=1,nbSources
        do k=1,nbSegs
          nullify(SPLdB(i,j,k)%f)
          nullify(SPLdB(i,j,k)%freq)
          call CopyFreqDomain(obsCont%obs(i)%dB(j,k), SPLdB(i,j,k))
        end do
      end do
    end do
    call WriteObsContFreqDomainData(obsCont,SPLdB, &
                                    path=trim(obsCont%SPLpath)//'_spectrum')
    deallocate(SPLdB)
  else if(SPLdBAFlag)then
    nbSources = size(obsCont%obs(1)%dBA,1)
    nbSegs    = size(obsCont%obs(1)%dBA,2)
    allocate(SPLdBA(nbObs,nbSources,nbSegs))
    do i=1,nbObs
      do j=1,nbSources
        do k=1,nbSegs
          nullify(SPLdBA(i,j,k)%f)
          nullify(SPLdBA(i,j,k)%freq)
          call CopyFreqDomain(obsCont%obs(i)%dBA(j,k), SPLdBA(i,j,k))
        end do
      end do
    end do
    call WriteObsContFreqDomainData(obsCont,SPLdBA, &
                                   path=trim(obsCont%SPLpath)//'_spectrum')
    deallocate(SPLdBA)
  end if
  !Write the phase data:
  if(phaseFlag)then
    nbSources = size(obsCont%obs(1)%phase,1)
    nbSegs    = size(obsCont%obs(1)%phase,2)
    allocate(phase(nbObs,nbSources,nbSegs))
    do i=1,nbObs
      do j=1,nbSources
        do k=1,nbSegs
          nullify(phase(i,j,k)%freq)
          nullify(phase(i,j,k)%f)
          call CopyFreqDomain(obsCont%obs(i)%phase(j,k), phase(i,j,k))
        end do
      end do
    end do
    call WriteObsContFreqDomainData(obsCont,phase, &
                                    path=trim(obsCont%phasePath)//'_spectrum')
    deallocate(phase) 
  end if
  !Write the complex pressure data:
  if(complexPressureFlag)then
    nbSources = size(obsCont%obs(1)%realCompP,1)
    nbSegs    = size(obsCont%obs(1)%realCompP,2)
    allocate(imagCompP(nbObs,nbSources,nbSegs), &
             realCompP(nbObs,nbSources,nbSegs))
    do i=1,nbObs
      do j=1,nbSources
        do k=1,nbSegs
          nullify(imagCompP(i,j,k)%freq)
          nullify(imagCompP(i,j,k)%f)
          nullify(realCompP(i,j,k)%freq)
          nullify(realCompP(i,j,k)%f)
          call CopyFreqDomain(obsCont%obs(i)%imagCompP(j,k), imagCompP(i,j,k))
          call CopyFreqDomain(obsCont%obs(i)%realCompP(j,k), realCompP(i,j,k))
        end do
      end do
    end do
    call WriteObsContFreqDomainData(obsCont,imagCompP,realCompP, &
                                    trim(obsCont%complexPressurePath)//'_spectrum')
    deallocate(imagCompP,realCompP) 
  end if   
end subroutine WriteObsContSpectralData

subroutine WriteObsContRadDist(obsCont)
  implicit none
  type(observerContainer), intent(in)::obsCont
  type(timeHistory), dimension(1):: radDist
  real, dimension(:,:,:), allocatable:: radDistance
  character(len=4096),dimension(:),allocatable::titles
  integer:: i,j, k, unitNum, nbSources

  if(obsCont%nbObs==1)then
    radDist(1)%nt = size(obsCont%obs(1)%radDistance)
    call SettHTitle(radDist(1),obsCont%obs(1)%radDistance(1)%title)
    allocate(radDist(1)%t(radDist(1)%nt))
    allocate(radDist(1)%f(radDist(1)%nt))
    do i=1,size(obsCont%obs(1)%radDistance)
      radDist(1)%t(i) = obsCont%obs(1)%radDistance(i)%t(1)
      radDist(1)%f(i) = obsCont%obs(1)%radDistance(i)%f(1)
    end do
    call WriteSingleObsTimeHistoryData(radDist, trim(obsCont%SPLpath)//'_radDistance',  &
                                       ' ', 'Observer ')
    deallocate(radDist(1)%t, radDist(1)%f)
  else
    !First write the grid file:
    if(getSize(obsCont%CBList).gt.0)then      
      call WriteMovingObsContSurface(obsCont,obsCont%segTimeArray, &
                                     trim(obsCont%SPLpath)//'_radDistance.x')
    else    
      call WriteStationaryObsSurface(obsCont,size(obsCont%segTimeArray), &
                                       trim(obsCont%SPLpath)//'_radDistance.x')
    end if

    nbSources = 1
    allocate(radDistance(obsCont%nbObs,nbSources,size(obsCont%segTimeArray)), &
             titles(nbSources+1))
    titles(1)='Time'
    do i=1, obsCont%nbObs
      do j=1,nbSources
         do k=1,size(obsCont%obs(1)%radDistance)
           !ksb debug: radDistance(i,j,k)= obsCont%obs(i)%radDistance(j)%f(k) !ksb debug - I think this is where the problem is 5/27/2015
           radDistance(i,j,k)= obsCont%obs(i)%radDistance(k)%f(1) !ksb debug:  I don't think we are using the radDistance
                                                                  !array of timeHistories correctly. Need to check for better way..
         end do
      end do
    end do
    do j=1,nbSources
      titles(j+1) = GettHTitle(obsCont%obs(1)%radDistance(j))
    end do
    call WriteTimeHistoryFunctionFile(radDistance,titles,obsCont%segTimeArray, &
                                      obsCont%nbObsDim1,obsCont%nbObsDim2,&
                                      trim(obsCont%SPLpath)//'_radDistance')
    deallocate(radDistance, titles)
  end if
end subroutine WriteObsContRadDist
  
subroutine WriteObsContBBmsp(obsCont)
  implicit none
  type(observerContainer), intent(in):: obsCont

  type(freqDomain), dimension(:,:,:), allocatable:: mspBB
  integer:: i, j, k, m, nf, nbObs, nbSegs, nbSources
  real:: intensity

  nbObs     = obsCont%nbObs
  nbSegs    = obsCont%nbSegments  
  nf = GetNf(obsCont%obs(1)%mspBB(1,1))
  
  if(obsCont%nbObs.gt.1)then
    call WriteStationaryObsSurface(obsCont,nf, &
                  trim(obsCont%SPLpath)//'_mspBB.x')
  end if

  nbSources = size(obsCont%obs(1)%mspBB,1)
  nbSegs    = size(obsCont%obs(1)%mspBB,2)
  allocate(mspBB(nbObs,nbSources,nbSegs))
  do i=1,nbObs
    do j=1,nbSources
      do k=1,nbSegs
        nullify(mspBB(i,j,k)%f)
        nullify(mspBB(i,j,k)%freq)
        call CopyFreqDomain(obsCont%obs(i)%mspBB(j,k), mspBB(i,j,k)) 
      end do
    end do
  end do
  call WriteObsContFreqDomainData(obsCont,mspBB, &
                                  path=trim(obsCont%SPLpath)//'_mspBB')

  if (AtmAtten) then
    do i=1,nbObs
      do j=1,nbSources
        do k=1,nbSegs
          nullify(mspBB(i,j,k)%f)
          nullify(mspBB(i,j,k)%freq)
          call CopyFreqDomain(obsCont%obs(i)%mspBBAtmAtten(j,k), mspBB(i,j,k)) 
        end do
      end do
    end do
    call WriteObsContFreqDomainData(obsCont,mspBB, &
                                    path=trim(obsCont%SPLpath)//'_mspBB_AtmAtten')
  end if

  deallocate(mspBB)

end subroutine WriteObsContBBmsp

end module ObserverContainerObject
