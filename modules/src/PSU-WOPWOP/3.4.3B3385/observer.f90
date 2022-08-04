module ObserverObject
  
  use TimeHistoryObject
  use FrequencyDomainObject
  use FrequencyAnalysisModule
  use ObserverUtility
  use ConstantsModule
  use MPIModule
  use strings
  
  implicit none  
  
  type observer
    character(len=4096)::title
    integer::obsIndex
    !The following has 1 dimension to store the noise source. 
    !  It can be allocated to store either a pressure time history or 
    !  pressure gradient information.
    type(timeHistory),dimension(:),pointer::pPrime=>null()
    type(timeHistory),dimension(:,:),pointer::pPrimeSegs=>null()
    !Time parameters
    integer::nt, nbFreqRanges
    real::tMin,tMax, dt
    !iBlank array for visualization, single dimension for time
    integer,dimension(:),pointer::iBlankArray=>NULL()
    !Actual coordinates of observers
    type(vector)::coordinates
    !nbSources - keep this information to use for allocation in the following arrays
    integer:: nbSources
    !realCompP(noiseSource,segNum)...
    type(freqDomain),dimension(:,:),pointer::realCompP=>null(),imagCompP=>null(),phase=>null(), &
      msp=>null(),dB=>null(),dBA=>null(), octaveFiltMSP=>null(),octaveFiltdB=>null(),octaveFiltdBA=>null(), &
      mspBB=>null(), mspBBAtmAtten=>null()
    !OASPLdB(noiseSource)...
    type(timeHistory),dimension(:),pointer::OASPLdB=>null(),OASPLdBA=>null(),PNL=>null(),PNLT=>null(),radDistance=>null()
    !dBRanges(noiseSource,rangeNumber)...
    type(timeHistory),dimension(:,:),pointer::dBRanges=>null(),dBARanges=>null()
    !EPNL(noiseSource)
    real,dimension(:),pointer::EPNL=>null(), SEL=>null()
  end type observer
  
  !ksb debug: 10/10/2017
  logical:: firstObsSetup
  
  interface operator (+)
    module procedure AddTimeHistories
  end interface

  
  private
  
  public::observer,DestroyObserver,                     &
          GetObserverNt,GetObserverTMax,GetObserverTMin,&
          GetObserverCoordinates,                       &
          SetObserverTMin,SetObserverTMax,SetObserverNt,&
          SetObserverIndex,WriteSingleObsFreqDomainData,&
          WriteNameFile, 			        &
          CreateObsWavFiles,CalcObsPNLData,             &  
          CalcObsSELData,FilterPPrime,windowTimeHistory,&
          CalcFiltComplexDFTFromPPrime,                 &
          CalcObsOASPLdBData,CalcObsOASPLdBAData,       &
          CalcTotalNoise, addThirdOctMSP, AtmosAtten,	&
          CalcObsComplexPressure,CalcObsPhaseFromComplexPressure,&
          CalcObsMSPFromComplexPressure,CalcObsDbFromMSP,        &
          CalcObsdBAFromMSP,CalcOctaveFiltMSP,                   &
          ResetPPrimeAndIndices,SetObserverTitle, InsertBBTitle, &
          CalcObsFreqRangedBAData,CalcObsFreqRangedBData,        &
          WriteSingleObsTimeHistoryData,WriteSingleObsFreqRangeData,    &
          SendObserverResultsToMaster, ReceiveObserverResultsFromSlave, &
          ResetObserver, Destroy2DTHArray, Destroy2DFDArray, NullifyObserver, firstObsSetup

contains
  
subroutine ResetObserver(obs, tMin, tMax, nt, dt)
  type(observer), intent(inout)::obs
  real, intent(in)::tMin, tMax, dt
  integer, intent(in)::nt
  integer::i
  call SetObserverTMin(obs, tMin)
  call SetObserverTMax(obs, tMax)
  call SetObserverNt(obs, nt)
  call SetObserverTitle(obs, 'Single Observer')
  if (.not.associated(obs%pPrime)) then
    allocate(obs%pPrime(GetSizeOfPPrime())) 
    do i=1, GetSizeOfPPrime()
      nullify(obs%pPrime(i)%f)
      nullify(obs%pPrime(i)%t)
    end do
  end if
  do i=1, GetSizeOfPPrime()
    if (.not.associated(obs%pPrime(i)%f)) then
      allocate(obs%pPrime(i)%f(nt))
    end if
    obs%pPrime(i)%f = 0.0
    call SetupTimeHistory(obs%pPrime(i), tMin,tMax,nt,'Slave Observer')
  end do
  if (.not.associated(obs%iBlankArray)) then
    allocate(obs%iBlankArray(nt))
  end if
  obs%iBlankArray = 0
  !ksb debug:
  call SetpPrimeIndices()
end subroutine ResetObserver

subroutine NullifyObserver(obs)
  implicit none
  type(observer),intent(inout)::obs
  nullify(obs%pPrime, obs%pPrimeSegs,&
          obs%iBlankArray,&
          obs%realCompP, obs%imagCompP,&
          obs%phase,&
          obs%msp,&
          obs%dB, obs%dBA,&
          obs%octaveFiltMSP,&
          obs%octaveFiltdB,obs%octaveFiltdBA,&
          obs%OASPLdB,obs%OASPLdBA,&
          obs%PNL,obs%PNLT,obs%EPNL,obs%SEL,&
          obs%dBRanges,obs%dBARanges)
end subroutine NullifyObserver

subroutine SendObserverResultsToMaster(obs, master, tag)
  type(observer), intent(in)::obs
  integer, intent(in)::master, tag
  integer:: nbSources
  !ksb debug:
  call SendString(obs%title, master, tag)
  call SendInteger(obs%obsIndex, master, tag)
  call SendPressureTHToMaster(obs, master, tag)
  
  nbSources = obs%nbSources
  call SendInteger(nbSources, master, tag)
  !end ksb debug:

  call SendOASPLAndPNLToMaster(obs, master, tag)
  call SenddBRangesToMaster(obs, master, tag)
  call SendOctaveSpectrumsToMaster(obs, master, tag) 
  call SendFDPressureToMaster(obs, master, tag)

end subroutine SendObserverResultsToMaster

!type observer
  !  character(len=4096)::title
  !  integer::obsIndex
  !  !The following has 1 dimension to store the noise source. 
  !  !  It can be allocated to store either a pressure time history or 
  !  !  pressure gradient information.
  !  type(timeHistory),dimension(:),pointer::pPrime=>null()
  !  type(timeHistory),dimension(:,:),pointer::pPrimeSegs=>null()
  !  !Time parameters
  !  integer::nt, nbFreqRanges
  !  real::tMin,tMax, dt
  !  !iBlank array for visualization, single dimension for time
  !  integer,dimension(:),pointer::iBlankArray=>NULL()
  !  !Actual coordinates of observers
  !  type(vector)::coordinates
  !  !nbSources - keep this information to use for allocation in the following arrays
  !  integer:: nbSources
  !  !realCompP(noiseSource,segNum)...
  !  type(freqDomain),dimension(:,:),pointer::realCompP=>null(),imagCompP=>null(),phase=>null(), &
  !    msp=>null(),dB=>null(),dBA=>null(), octaveFiltMSP=>null(),octaveFiltdB=>null(),octaveFiltdBA=>null(), &
  !    mspBB=>null(), mspBBAtmAtten=>null()
  !  !OASPLdB(noiseSource)...
  !  type(timeHistory),dimension(:),pointer::OASPLdB=>null(),OASPLdBA=>null(),PNL=>null(),PNLT=>null(),radDistance=>null()
  !  !dBRanges(noiseSource,rangeNumber)...
  !  type(timeHistory),dimension(:,:),pointer::dBRanges=>null(),dBARanges=>null()
  !  !EPNL(noiseSource)
  !  real,dimension(:),pointer::EPNL=>null(), SEL=>null()
  !end type observer

subroutine ReceiveObserverResultsFromSlave(obs, master, sender)
  type(observer), intent(inout)::obs
  integer, intent(in)::master, sender
  !ksb debug:
  integer:: nbSources
  call ReceiveString(obs%title, master, sender)
  call ReceiveInteger(obs%obsIndex, master, sender)
  call ReceivePressureTHFromSlave(obs, master, sender)  
  
  call ReceiveInteger(nbSources, master, sender)
  obs%nbSources = nbSources
  !end ksb debug:

  call ReceiveOASPLAndPNLFromSlave(obs, master, sender)
  call ReceivedBRangesFromSlave(obs, master, sender)
  call ReceiveOctaveSpectrumsFromSlave(obs, master, sender)
  call ReceiveFDPressureFromSlave(obs, master, sender)

end subroutine ReceiveObserverResultsFromSlave

subroutine Send2DimTimeHistory(THArray, master, tag)
  type(timeHistory), dimension(:,:)::THArray
  integer, intent(in)::master, tag
  integer::arraySize1, arraySize2, i, j
  arraySize1 = size(THArray, 1)
  arraySize2 = size(THArray, 2)
  call SendInteger(arraySize1, master, tag)
  call SendInteger(arraySize2, master, tag)
  do i=1, arraySize1
    do j=1, arraySize2
      call SendTimeHistoryToMaster(THArray(i,j), master, tag)
    end do
  end do
end subroutine Send2DimTimeHistory

subroutine Send2DimFreqDomain(FDArray, master, tag)
  type(freqDomain), dimension(:,:)::FDArray
  integer, intent(in)::master, tag
  integer::arraySize1, arraySize2, i, j
  arraySize1 = size(FDArray, 1)
  arraySize2 = size(FDArray, 2)
  call SendInteger(arraySize1, master, tag)
  call SendInteger(arraySize2, master, tag)
  do i=1, arraySize1
    do j=1, arraySize2
      call SendFreqDomainToMaster(FDArray(i,j), master, tag)
    end do
  end do
end subroutine Send2DimFreqDomain

subroutine SendFDPressureToMaster(obs, master, tag)
  type(observer), intent(in)::obs
  integer, intent(in)::master, tag
  if ((spectrumFlag.and.complexPressureFlag).or.FSCOutputFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending Frequency Domain Pressure (imag)'
    end if
    call Send2DimFreqDomain(obs%imagCompP, master, tag)
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending Frequency Domain Pressure (real)'
    end if
    call Send2DimFreqDomain(obs%realCompP, master, tag)
  end if
  if (spectrumFlag.and.phaseFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending Frequency Domain Pressure (phase)'
    end if
    call Send2DimFreqDomain(obs%phase, master, tag)
  end if
  if (spectrumFlag.and.SPLdBFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending Frequency Domain Pressure (dB)'
    end if
    call Send2DimFreqDomain(obs%dB, master, tag)
  end if
  if (spectrumFlag.and.SPLdBAFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending Frequency Domain Pressure (dBA)'
    end if
    call Send2DimFreqDomain(obs%dBA, master, tag)
  end if
end subroutine SendFDPressureToMaster

subroutine SendOctaveSpectrumsToMaster(obs, master, tag)
  type(observer), intent(in)::obs
  integer, intent(in)::master, tag
  if (SPLdBFlag.and.octaveFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending octave filtered (dB)'
    end if
    call Send2DimFreqDomain(obs%octaveFiltdB, master, tag)
  end if
  if (SPLdBAFlag.and.octaveFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending octave filtered (dBA)'
    end if
    call Send2DimFreqDomain(obs%octaveFiltdBA, master, tag)
  end if
  if (globalPeggNoiseFlag.or.globalBPMNoiseFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending mspBB'
    end if
    call Send2DimFreqDomain(obs%mspBB, master, tag)
  end if
  if ((globalPeggNoiseFlag.or.globalBPMNoiseFlag) .and. AtmAtten) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending mspBBAtmAtten'
    end if
    call Send2DimFreqDomain(obs%mspBBAtmAtten, master, tag)
  end if 
end subroutine SendOctaveSpectrumsToMaster

subroutine SenddBRangesToMaster(obs, master, tag)
  type(observer), intent(in)::obs
  integer, intent(in)::master, tag
  if (SPLdBFlag.and.obs%nbFreqRanges.ne.0) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending Frequency Ranges (dB)'
    end if
    call Send2DimTimeHistory(obs%dBRanges, master, tag)
  end if
  if (SPLdBAFlag.and.obs%nbFreqRanges.ne.0) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending Frequency Ranges (dBA)'
    end if
    call Send2DimTimeHistory(obs%dBARanges, master, tag)
  end if
end subroutine SenddBRangesToMaster
    
subroutine SendOASPLAndPNLToMaster(obs, master, tag)
  type(observer), intent(in)::obs
  integer, intent(in)::master, tag
  integer::arraySize, i
  if(OASPLdBFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending OASPL (dB)'
    end if
    arraySize = size(obs%OASPLdB)
    call SendInteger(arraySize, master, tag)
    do i=1, arraySize
      call SendTimeHistoryToMaster(obs%OASPLdB(i), master, tag)
    end do
  end if
  if (OASPLdBAFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending OASPL (dBA)'
    end if
    arraySize = size(obs%OASPLdBA)
    call SendInteger(arraySize, master, tag)
    do i=1, arraySize
      call SendTimeHistoryToMaster(obs%OASPLdBA(i), master, tag)
    end do
  end if
  if (PNLFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending PNL (dB)'
    end if
    arraySize = size(obs%PNL)
    call SendInteger(arraySize, master, tag)
    do i=1, arraySize
      call SendTimeHistoryToMaster(obs%PNL(i), master, tag)
    end do
  end if
  if (PNLTFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending PNLT (dB)'
    end if
    arraySize = size(obs%PNLT)
    call SendInteger(arraySize, master, tag)
    do i=1, arraySize
      call SendTimeHistoryToMaster(obs%PNLT(i), master, tag)
    end do
    call SendLogical(PNLTdBDrop(1), master, tag)
    call SendLogical(PNLTdBDrop(2), master, tag)
  end if
  ! this probably doesn't belong here (but it really doesn't matter) KSB 5/28/2015
  if (AtmAtten) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending radDistance'
    end if
    arraySize = size(obs%radDistance)
    ! probably this could be simplified - there is no reason to send the distance for every source
    ! Distance should be the same for all of them, hence, this should not be an array. KSB 5/28/2015
    call SendInteger(arraySize, master, tag)
    do i=1, arraySize
      call SendTimeHistoryToMaster(obs%radDistance(i), master, tag)
    end do
  end if
  if (EPNLFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending EPNL (EPNLdB)'
    end if
    arraySize = size(obs%EPNL)
    call SendInteger(arraySize, master, tag)
    call SendReals(obs%EPNL, master, tag)
  end if
 if (SELFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending SEL (dBA)'
    end if
    arraySize = size(obs%SEL)
    call SendInteger(arraySize, master, tag)
    call SendReals(obs%SEL, master, tag)
    call SendLogical(SELdBDrop(1), master, tag)
    call SendLogical(SELdBDrop(2), master, tag)
  end if
end subroutine SendOASPLAndPNLToMaster
    
subroutine ReceiveOASPLAndPNLFromSlave(obs, master, sender)
  type(observer), intent(inout)::obs
  integer, intent(in)::master, sender
  integer::arraySize, i
  if (OASPLdBFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving OASPL (dB)'
    end if
    call ReceiveInteger(arraySize, master, sender)
    allocate(obs%OASPLdB(arraySize))
    do i=1, arraySize
      nullify(obs%OASPLdB(i)%t)
      nullify(obs%OASPLdB(i)%f)
      call ReceiveTimeHistoryFromSlave(obs%OASPLdB(i), master, sender)
    end do
  end if
  if (OASPLdBAFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving OASPL (dBA)'
    end if
    call ReceiveInteger(arraySize, master, sender)
    allocate(obs%OASPLdBA(arraySize))
    do i=1, arraySize
      nullify(obs%OASPLdBA(i)%t)
      nullify(obs%OASPLdBA(i)%f)
      call ReceiveTimeHistoryFromSlave(obs%OASPLdBA(i), master, sender)
    end do
  end if
  if (PNLFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving PNL (dB)'
    end if
    call ReceiveInteger(arraySize, master, sender)
    allocate(obs%PNL(arraySize))
    do i=1, arraySize
      nullify(obs%PNL(i)%t)
      nullify(obs%PNL(i)%f)
      call ReceiveTimeHistoryFromSlave(obs%PNL(i), master, sender)
    end do
  end if
  if (PNLTFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving PNLT (dB) '
    end if
    call ReceiveInteger(arraySize, master, sender)
    allocate(obs%PNLT(arraySize))
    do i=1, arraySize
      nullify(obs%PNLT(i)%t)
      nullify(obs%PNLT(i)%f)
      call ReceiveTimeHistoryFromSlave(obs%PNLT(i), master, sender)
    end do
    call ReceiveLogical(PNLTdBDrop(1), master, sender)
    call ReceiveLogical(PNLTdBDrop(2), master, sender)
  end if
  if (AtmAtten) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving radDistance'
    end if
    call ReceiveInteger(arraySize, master, sender)
    allocate(obs%radDistance(arraySize))
    do i=1, arraySize
      nullify(obs%radDistance(i)%t)
      nullify(obs%radDistance(i)%f)
      call ReceiveTimeHistoryFromSlave(obs%radDistance(i), master, sender)
    end do
  end if
  if (EPNLFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving EPNL'
    end if
    call ReceiveInteger(arraySize, master, sender)
    allocate(obs%EPNL(arraySize))
    call ReceiveReals(obs%EPNL, master, sender)
  end if
  if (SELFLag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving SEL'
    end if
    call ReceiveInteger(arraySize, master, sender)
    allocate(obs%SEL(arraySize))
    call ReceiveReals(obs%SEL, master, sender)
    call ReceiveLogical(SELdBDrop(1), master, sender)
    call ReceiveLogical(SELdBDrop(2), master, sender)
  end if
end subroutine ReceiveOASPLAndPNLFromSlave
    
subroutine SendPressureTHToMaster(obs, master, tag)
  type(observer), intent(in)::obs
  integer, intent(in)::master, tag
  integer::arraySize, i
  if (acousticPressureFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Sending Pressure'
    end if
    arraySize = size(obs%pPrime)
    call SendInteger(arraySize, master, tag)
    do i=1, arraySize
      call SendTimeHistoryToMaster(obs%pPrime(i), master, tag)
    end do
  !ksb debug:  I want it to do the pprime title array stuff next.
  !else
  !  return
  end if
  ! I don't know were to put this, so I'll pass it to the master here
  arraySize = size(pprimeTitleArray)
  call SendInteger(arraySize, master, tag)
  do i=1, arraySize
    call SendString(pprimeTitleArray(i), master, tag)
  end do
  !ksb debug: character(len=4096),dimension(:),pointer::fRangeTitle=>null(),pprimeTitleArray=>null()
end subroutine SendPressureTHToMaster
    
subroutine ReceivePressureTHFromSlave(obs, master, sender)
  type(observer), intent(inout)::obs
  integer, intent(in)::master, sender
  integer::arraySize,i 
  if (acousticPressureFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving Pressure'
    end if
    call ReceiveInteger(arraySize, master, sender)
    allocate(obs%pPrime(arraySize))
    do i=1, arraySize
      nullify(obs%pPrime(i)%t)
      nullify(obs%pPrime(i)%f)
      call ReceiveTimeHistoryFromSlave(obs%pPrime(i), master, sender)
    end do
  end if
  call ReceiveInteger(arraySize, master, sender)
  do i=1, arraySize
    call ReceiveString(pprimeTitleArray(i), master, sender)
  end do
end subroutine ReceivePressureTHFromSlave

subroutine Receive2DimTimeHistory(THArray, master, sender)
  type(TimeHistory), dimension(:,:), pointer::THArray
  integer, intent(in)::master, sender
  integer::arraySize1, arraySize2, i, j
  call ReceiveInteger(arraySize1, master, sender)
  call ReceiveInteger(arraySize2, master, sender)
  allocate(THArray(arraySize1, arraySize2))
  do i=1, arraySize1
    do j=1, arraySize2
      nullify(THArray(i,j)%t)
      nullify(THArray(i,j)%f)
      call ReceiveTimeHistoryFromSlave(THArray(i,j), master, sender)
    end do
  end do
end subroutine Receive2DimTimeHistory

subroutine Receive2DimFreqDomain(FDArray, master, sender)
  type(FreqDomain), dimension(:,:), pointer::FDArray
  integer, intent(in)::master, sender
  integer::arraySize1, arraySize2, i, j
  call ReceiveInteger(arraySize1, master, sender)
  call ReceiveInteger(arraySize2, master, sender)
  allocate(FDArray(arraySize1, arraySize2))
  do i=1, arraySize1
    do j=1, arraySize2
      nullify(FDArray(i,j)%freq)
      nullify(FDArray(i,j)%f)
      call ReceiveFreqDomainFromSlave(FDArray(i,j), master, sender)
    end do
  end do
end subroutine Receive2DimFreqDomain

subroutine ReceiveFDPressureFromSlave(obs, master, sender)
  type(observer), intent(in)::obs
  integer, intent(in)::master, sender
  if ((spectrumFlag.and.complexPressureFlag).or.FSCOutputFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving Frequency Domain Pressure (imag)'
    end if
    call Receive2DimFreqDomain(obs%imagCompP, master, sender)
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving Frequency Domain Pressure (real)'
    end if
    call Receive2DimFreqDomain(obs%realCompP, master, sender)
  end if
  if (spectrumFlag.and.phaseFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving Frequency Domain Pressure (phase)'
    end if
    call Receive2DimFreqDomain(obs%phase, master, sender)
  end if
  if (spectrumFlag.and.SPLdBFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving Frequency Domain Pressure (dB)'
    end if
    call Receive2DimFreqDomain(obs%dB, master, sender)
  end if
  if (spectrumFlag.and.SPLdBAFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving Frequency Domain Pressure (dBA)'
    end if
    call Receive2DimFreqDomain(obs%dBA, master, sender)
  end if
end subroutine ReceiveFDPressureFromSlave

subroutine ReceiveOctaveSpectrumsFromSlave(obs, master, sender)
  type(observer), intent(in)::obs
  integer, intent(in)::master, sender
  if (SPLdBFlag.and.octaveFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving Octave Filtered SPL (dB)'
    end if
    call Receive2DimFreqDomain(obs%octaveFiltdB, master, sender)
  end if
  if (SPLdBAFlag.and.octaveFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving Octave Filtered SPL (dBA)'
    end if
    call Receive2DimFreqDomain(obs%octaveFiltdBA, master, sender)
  end if
  if (globalPeggNoiseFlag.or.globalBPMNoiseFlag) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving mspBB'
    end if
    call Receive2DimFreqDomain(obs%mspBB, master, sender)
  end if
  if ((globalPeggNoiseFlag.or.globalBPMNoiseFlag) .and. AtmAtten) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving mspBBAtmAtten'
    end if
    call Receive2DimFreqDomain(obs%mspBBAtmAtten, master, sender)
  end if
end subroutine ReceiveOctaveSpectrumsFromSlave

subroutine ReceivedBRangesFromSlave(obs, master, sender)
  type(observer), intent(inout)::obs
  integer, intent(in)::master, sender
  if (SPLdBFlag.and.obs%nbFreqRanges.ne.0) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving Frequency Ranges (dB)'
    end if
    call Receive2DimTimeHistory(obs%dBRanges, master, sender)
  end if
  if (SPLdBAFlag.and.obs%nbFreqRanges.ne.0) then
    if (debugLevel.gt.11) then
      write(*,*) ProcID(), 'Receiving Frequency Ranges (dBA)'
    end if
    call Receive2DimTimeHistory(obs%dBARanges, master, sender)
  end if
end subroutine ReceivedBRangesFromSlave

!**
!SUBROUTINE CalcObsPNLData
!  This routine calculates PNL (perceived noise level) and 
!  tonal corrected PNL (PNLT) data for an observer
!**2

subroutine CalcObsPNLData(obs,timeArray)
  implicit none
  type(observer),intent(inout)::obs
  real,dimension(:),intent(in)::timeArray
  
  integer::i,j,k,nbSources,nbSegs,count,maxK,nt, index, endIndex, nbEmptyBins
  real,dimension(:),allocatable::PNLArray,PNLTArray,SPLArray,cMax
  real::PNLTM !Maximum PNLT from all segments
  real::startTime,endTime,d,summation,durationCorrection,dt,tMin,tMax
  real::SPLdBA
  type(freqDomain)::tempMSP
  character(len=1024):: BBtitle

  nbSources = size(obs%msp,1)
  nbSegs    = size(obs%pPrimeSegs,2)
  if (globalPeggNoiseFlag.or.globalBPMnoiseFlag) then
    nbSources = nbSources+1
    if (globalPeggNoiseFlag) then
      BBtitle = 'Broadband-Pegg'
    else
      BBtitle = 'Broadband-BPM'
    end if
  else
    BBtitle   = ''
  end if
  nt = size(timeArray)
  tMin = timeArray(1)
  tMax = timeArray(nt)
  if (EPNLFlag) PNLTdBDrop = .false.
  if(PNLFlag.or.EPNLFlag.or.PNLTFlag)then
    allocate(PNLArray(nbSegs))
  end if
  if(PNLFlag)then
    allocate(obs%PNL(nbSources))
    do i=1, nbSources
      nullify(obs%PNL(i)%t,obs%PNL(i)%f)
      call SetNT(obs%PNL(i),nt)
      call SettMin(obs%PNL(i),tMin)
      call SettMax(obs%PNL(i),tMax)
      if (len_trim(BBtitle).eq.0) then
        !ksb debug: call SettHTitle(obs%PNL(i),trim(obs%PPrime(i)%title))
        call SettHTitle(obs%PNL(i),trim(pprimeTitleArray(i)))
      else
        call SettHTitle(obs%PNL(i),trim(InsertBBTitle(obs, BBtitle, i, nbSources)))
      end if
      allocate(obs%PNL(i)%t(nt),obs%PNL(i)%f(nt)) 
      obs%PNL(i)%t = timeArray
    end do
  end if
  if(PNLTFlag)then
    allocate(obs%PNLT(nbSources))
    do i=1, nbSources
      nullify(obs%PNLT(i)%t,obs%PNLT(i)%f)
      call SetNT(obs%PNLT(i),nt)
      call SettMin(obs%PNLT(i),tMin)
      call SettMax(obs%PNLT(i),tMax)
      if (len_trim(BBtitle).eq.0) then
        !ksb debug: call SettHTitle(obs%PNLT(i),trim(obs%PPrime(i)%title))
        call SettHTitle(obs%PNLT(i),trim(pprimeTitleArray(i)))
      else
        call SettHTitle(obs%PNLT(i),trim(InsertBBTitle(obs, BBtitle, i, nbSources)))
      end if
      allocate(obs%PNLT(i)%t(nt),obs%PNLT(i)%f(nt)) 
      obs%PNLT(i)%t = timeArray
    end do
  end if
  if(PNLTFlag.or.EPNLFlag)then
    allocate(PNLTArray(nbSegs),cMax(nbSegs))
  end if
  if(EPNLFlag)then
    if(.not.associated(obs%EPNL)) then
      allocate(obs%EPNL(nbSources))
    end if
  end if
  do i=1,nbSources   
    do j=1,nbSegs
      ! We have to filter the octave filtered spectrum it so that 
      ! there are 24 bins, ranging from 50Hz to 10000Hz in 1/3 octaves:
      nullify(tempMSP%freq,tempMSP%f)
      if (len_trim(BBtitle).eq.0) then
        call FilterFreqDomain(obs%octaveFiltMSP(i,j), 44.7,11220.0, tempMSP)
      else
        ! We have to insert the broadband noise inbetween 
        ! the discrete frequency results and the total.
        call InsertBBArray(obs, 44.7, 11220.0, tempMSP, i, j)
      end if
      nbEmptyBins = 0
      do k=size(tempMSP%freq),1,-1
        if (tempMSP%f(k).eq.0.0) then
          nbEmptyBins = nbEmptyBins + 1
        else
          exit
        end if
      end do
      !ksb debug:  This is a temporary information line.  I need to figure out what should be
      ! done if the frequency bins are below the starting frequency, hence no data to follow.
      if( size(tempMSP%freq)-nbEmptyBins == 0 ) then  
          print*,' Error in observer.f90: CalcObsPNLData: '
          print*,' The number of empty bins is the same as the number of total bins'
          print*,' You can probably correct this by specifying more time steps to ensure'
          print*,' The maximum frequency is above 50 Hz.'
          stop
      end if
      if (allocated(SPLArray)) deallocate(SPLArray)
      allocate(SPLArray(size(tempMSP%freq)-nbEmptyBins))
      do k=1,size(SPLArray)
        SPLArray(k) = ConvertMSPTodB(tempMSP%f(k))
      end do

!      if (AtmAtten) call AtmosAtten(tempMSP%freq, SPLArray, obs%radDistance(j)%f(1))

      if(PNLFlag.or.PNLTFlag.or.EPNLFlag)then
        PNLArray(j) = CalcPNL(SPLArray)
      end if
      !This calculates the tone corrected PNL:
      if(PNLTFlag.or.EPNLFlag)then
        if (EPNLPrime) then
          PNLTArray(j) = PNLArray(j)
	    else
          call CalcPNLT(SPLArray,PNLArray(j),cMax(j),PNLTArray(j))
        end if
      end if
    end do
    if (size(SPLArray) .lt. 24 .and. trim(obs%PNLT(i)%title).eq.'Total') then
      call Warning('SPL_in does not have 24 bands - populating missing bands with spl(i)=spl(i-1)-1.5dB in PNLT calculation')
    end if
    if(PNLTFlag)then
      obs%PNLT(i)%f = PNLTArray
    end if
    if(PNLFlag)then
      obs%PNL(i)%f = PNLArray
    end if
    if(EPNLFlag)then
      PNLTM = maxval(PNLTArray)
      index = maxloc(PNLTArray,1)
      if((index.gt.2).and.(index.lt.(nbSegs-1)))then
        if(cMax(index).lt.(sum(cMax(index-2:index+2))/5.))then
          PNLTM = PNLTM - cMax(index) + sum(cMax(index-2:index+2))/5.
        end if
      end if
      ! Check whether the case achieves PNLTM-10dB before AND after PNLTM.  
      ! However, we need to be careful about overwriting the values in PNLTdBDrop.
      ! For the first source we need to check if PNLTM-10dB passes, for 
      ! the other sources we only care if they DON'T pass.
      if (i.eq.1 .or. trim(obs%PNLT(i)%title).eq.'Total') then !ksb debug: PPrime(i)%title).eq.'Total') then
        if (minval(PNLTArray(1:index-1))     .lt.(PNLTM-10.)) PNLTdBDrop(1) = .true.
        if (minval(PNLTArray(index+1:nbSegs)).lt.(PNLTM-10.)) PNLTdBDrop(2) = .true.
      else
        if (minval(PNLTArray(1:index-1))     .ge.(PNLTM-10.)) PNLTdBDrop(1) = .false.
        if (minval(PNLTArray(index+1:nbSegs)).ge.(PNLTM-10.)) PNLTdBDrop(2) = .false.
      end if
!!!   dt is never called upon   BAG 11/15/11
!      dt = timeArray(2)-timeArray(1)
	 
      if (all(PNLTdBDrop) .or. forceEPNL) then
        !First calculate durationCorrection, the duration correction factor:
        count = 0
        do j=1,nbSegs  !ksb debug: size(PNLTArray)
          if(PNLTArray(j).ge.(PNLTM-10.))then
            if(count==0)then
!ksb debug:                if(abs(PNLTArray(j-1)-(PNLTM-10.)).lt.abs(PNLTArray(j)-(PNLTM-10.)))then
!!                startTime = timeArray(j-1)
!                index = j-1
!                else
!!                startTime = timeArray(j)
!                index = j
!                end if
              index = j
              if( j>1 ) then
                if(abs(PNLTArray(j-1)-(PNLTM-10.)).lt.abs(PNLTArray(j)-(PNLTM-10.)))then
                  index = j-1
                end if
              end if  !ksb debug; end of change
            end if
            count = count + 1
          end if
        end do
        !ksb debug: if(count==0 ) then
        if(count==0 .or. (index+count-1)<1 .or. (index+count)>nbSegs )then
!          startTime = timeArray(1)
          index = 1
!          endTime   = startTime
          endIndex   = index
        else
          if(abs(PNLTArray(index+count)-(PNLTM-10.)).lt. &
             abs(PNLTArray(index+count-1)  -(PNLTM-10.)))then
!            endTime = timeArray(index+count+1)
            endIndex = index+count
          else
!            endTime = timeArray(index+count) 
	        endIndex = index+count-1
          end if
        end if
!        d = endTime - startTime
!        maxk = nint(2*d)
        summation = 0.
!        do j=1,maxk
!ksb debug: 2/25/16 - this modification is for the case with observer time expansion turned off
if( forceEPNL ) then
  if( index == 0 ) then
    index = 1
    print*,'CalcObsPNLData: index adjusted to 1'
  end if
  if( endIndex == 1 ) then
    endIndex = nbSegs
    print*,'CalcObsPNLData: endIndex adjusted to nbSegs'
  end if
end if
	    do j=index, endIndex
!          summation = summation + 10.**(PNLTArray(j+index-1)/10.)
          summation = summation + 10.**(PNLTArray(j)/10.)
        end do
        if (summation > 0 ) then  !ksb debug: added check that sumation > 0
          durationCorrection = 10.*log10(summation)-PNLTM-13.
        else
          durationCorrection = 0.
          call Warning('CalcObsPNLData: summation <= 0 in durationCorrection  computation')
        end if
        obs%EPNL(i) = PNLTM + durationCorrection
      end if
    end if    
  end do
  if (.not. all(PNLTdBDrop)) then
    ! We will only rerun the case if the FAACertFlag is on and the data is NOT aperiodic.
    ! Otherwise the user will simply be given a warning.      
    if (.not.PNLTdBDrop(1)) then
      if (.not.dataLimitLHS) then
!ksb debug: Turn off the observer time expansion - it's not working in parallel (or scalar?)
if(.not.forceEPNL) then
  print*,'The PNLTM-10dB critera for FAA Certification was not met BEFORE the'  !ksb debug
  print*,'maximum tone corrected Perceived Noise Level (PNLTM) value.' !ksb debug:
  print*,'Adjust observer time TMIN to an earlier time.' !ksb debug:
  print*,'obsIndex = ',obs%obsIndex
end if
      !  if (.not. readInPressureFlag) newSegLHS = .true.  
        call Notice('The PNLTM-10dB critera for FAA Certification was not met BEFORE the',&
                     'Maximum Perceived Noise Level, Tone (PNLTM).')
      else
        if (HitMaxExpansions) then
          call Notice(' The maximum number of observer time range expansions has been',&
                      ' reached and PNLTM-10dB could not be obtained BEFORE the peak.')
        else
          call Notice('The source time range in your aperiodic data files',&
                      'is not sufficient to capture PNLTM-10dB BEFORE the peak.')
        end if
      end if
    end if
    if (.not.PNLTdBDrop(2)) then
      if(.not.dataLimitRHS) then 
!ksb debug:  Turn of observer time expansin - it is not working right now. 2/6/2016
if(.not.forceEPNL) then
  print*,'The PNLTM-10dB critera for FAA Certification was not met AFTER the' !ksb debug:
  print*,'maximum tone corrected Percieved Noise (PNLTM) value.'  !ksb debug:
  print*,'Adjust observer time TMAX to a later time'  !ksb debug:
  print*,'obsIndex = ',obs%obsIndex
end if 
       ! if (.not. readInPressureFlag) newSegRHS = .true.
        call Notice('The PNLTM-10dB critera for FAA Certification was not met AFTER the',&
                     'Maximum Perceived Noise Level, Tone (PNLTM).')
      else
        if (HitMaxExpansions) then
          call Notice(' The maximum number of observer time range expansions has been',&
                      ' reached and PNLTM-10dB could not be obtained AFTER the peak.')
        else
          call Notice('The source time range in your aperiodic data files',&
                      'is not sufficient to capture PNLTM-10dB AFTER the peak.')
        end if
      end if
    end if
  end if 
  call Notice('Using '//trim(integertostring(tempMSP%nf))//' 1/3 octave filtered', &
              'frequency bins for EPNL/PNLT calculations')
  
  if(allocated(SPLArray)) then
    deallocate(SPLArray)
  end if
  if(allocated(PNLArray)) then
    deallocate(PNLArray)
  end if
  if(allocated(PNLTArray)) then
    deallocate(PNLTArray)
  end if
  if(allocated(cMax)) then  
    deallocate(cMax)
  end if
!ksb debug:  This is part of turning off the observer time expansion 2/6/2016
  !if (.not.(all(PNLTdBDrop) .or. forceEPNL)) then
    if(associated(obs%iBlankArray)) then
      deallocate(obs%iBlankArray)
    end if  
    !if(associated(obs%EPNL)) then
    !  deallocate(obs%EPNL)
    !end if  
 ! end if  !ksb debug: this is part of turning off the observer time expansion 2/6/2016
end subroutine CalcObsPNLData


subroutine CalcObsSELData(obs)
  implicit none
  type(observer),intent(inout)::obs
  
  integer::i,j,k,nbSources,nbSegs,maxIndex,startIndex,endIndex
  real,dimension(:),allocatable::SPLArray, LAarray
  real::SELtemp, LAmax
  type(freqDomain)::tempMSP
  character(len=1024):: BBtitle

  nbSources = size(obs%msp,1)
  nbSegs    = size(obs%pPrimeSegs,2)

  if (globalPeggNoiseFlag.or.globalBPMnoiseFlag) then
    nbSources = nbSources+1
    if (globalPeggNoiseFlag) then
      BBtitle = 'Broadband-Pegg'
    else
      BBtitle = 'Broadband-BPM'
    end if
  else
    BBtitle   = ''
  end if  
  
  if(.not.associated(obs%SEL)) then
    allocate(obs%SEL(nbSources))
  end if
  if (.not. allocated(LAarray)) then
    allocate(LAarray(nbSegs))
  end if
  
  SELdbDrop = .false.

  do i=1,nbSources
    SELtemp = 0.0
    LAarray = 0.0
    do j=1,nbSegs
      startIndex = 0
      endIndex = 1
      ! We have to filter the octave filtered spectrum it so that 
      ! there are 28 bins, ranging from 17.8Hz to 11220Hz in 1/3 octaves:
      nullify(tempMSP%freq,tempMSP%f)
      if (len_trim(BBtitle).eq.0) then
        call FilterFreqDomain(obs%octaveFiltMSP(i,j), 17.8,11220.0, tempMSP)
      else 
        ! We have to insert the broadband noise inbetween 
        ! the discrete frequency results and the total.
        call InsertBBArray(obs, 17.8, 11220.0, tempMSP, i, j)
      end if
         
      if(allocated(SPLArray)) deallocate(SPLArray)
      allocate(SPLArray(size(tempMSP%freq)))
      do k=1,size(SPLArray)
        ! Fill SPLArray but not before converting the MSP to A-weighted
        SPLArray(k) = ConvertMSPTodBA(tempMSP%f(k), tempMSP%freq(k))
      end do
!      if (AtmAtten .and. associated(obs%radDistance)) call AtmosAtten(tempMSP%freq, SPLArray, obs%radDistance(j)%f(1))

      do k=1,size(SPLArray)
        LAarray(j) = LAarray(j) + 10.0**(SPLArray(k)/10.0)
      end do
      LAarray(j) = 10.0*log10(LAarray(j)+epsilon**3.0)
    end do

    LAmax = maxval(LAarray)
    maxIndex = maxloc(LAarray,1)
    ! Find where the SPL first goes above LAmax before and
    ! after the peak.
    do k=1,maxIndex-1
      if (LAarray(k) .gt. LAmax-10.) then
        startIndex = k
        exit
      end if  
    end do
    do k=size(LAarray),maxIndex+1,-1
      if (LAarray(k) .gt. LAmax-10.) then
        endIndex = k
        exit
      end if  
    end do
    !ksb debug: 5/24/2015
    !do k=maxIndex-1,1,-1
    !  if (LAarray(k) .le. LAmax-10.) then
    !    startIndex = k
    !    exit
    !  end if  
    !end do
    !do k=maxIndex+1,size(LAarray), 1
    !  if (LAarray(k) .le. LAmax-10.) then
    !    endIndex = k
    !    exit
    !  end if  
    !end do
!ksb debug: 2/23/16 - this modification is for the case with observer time expansion turned off
if( forceSEL ) then
  if( startIndex == 0 ) startIndex = 1
  if( endIndex == 1 ) endIndex = maxIndex
end if

    if (startIndex /= 0 .and. endIndex /= 1) then 
      do k=startIndex,endIndex
        ! Summing SPLdBA contribution from each time segment whose value
        ! is greater than 10dB down from the maximum.
        SELtemp = SELtemp + 10.**(LAarray(k)/10.)
      end do
      obs%SEL(i) = 10.*log10(SELtemp+epsilon**3.0)
      SELdBDrop(1) = .true.
      SELdbDRop(2) = .true.
    else
      if( startIndex /= 0 ) then
        SELdBDrop(1) = .true.
        SELdBDrop(2) = .false.
      else
        SELdBDrop(1) = .false.
        SELdBDrop(2) = .true.
      end if
    end if
    
  end do
  call Notice('Using '//trim(integertostring(tempMSP%nf))//' 1/3 octave filtered', &
              'frequency bins for SEL calculations')

  if (.not. all(SELdBDrop)) then
    ! We will only rerun the case if the data is NOT aperiodic.
    ! Otherwise the user will simply be given a warning.      
    if (.not.SELdBDrop(1)) then
      if (.not.dataLimitLHS) then
!ksb debug:  Turn of observer time expansin - it is not working right now. 2/6/2016
if( .not.forceSEL ) then !ksb debug:
  print*,'The LAmax-10dB critera for FAA Certification was not met BEFORE the' !ksb debug:
  print*,'maximum A-weighted sound pressure Level (LAmax).'  !ksb debug
  print*,'Adjust observer time TMIN to an earlier time.' !ksb debug: 
  print*,'obsIndex = ',obs%obsIndex
end if
       ! if (.not. readInPressureFlag) newSegLHS = .true.
        call Notice('The LAmax-10dB critera for FAA Certification was not met BEFORE the',&
                     'maximum A-weighted sound pressure Level (LAmax).')
      else
        if (HitMaxExpansions) then
          call Notice(' The maximum number of observer time range expansions has been',&
                      ' reached and LAmax-10dB could not be obtained BEFORE the peak.')
        else
          call Notice('The source time range in your aperiodic data files',&
                      'is not sufficient to capture LAmax-10dB BEFORE the peak.')
        end if
      end if
    end if
    if (.not.SELdBDrop(2)) then
      if(.not.dataLimitRHS) then  
!ksb debug:  Turn of observer time expansin - it is not working right now. 2/6/2016
if( .not.forceSEL ) then !ksb debug:
  print*,'The LAmax-10dB critera for FAA Certification was not met AFTER the' !ksb debug:
  print*,'maximum A-weighted sound pressure Level (LAmax).'  !ksb debug:
  print*,'Adjust observer time TMAX to a later time.' !ksb debug:
  print*,'obsIndex = ',obs%obsIndex
end if !ksb debug:
        !if (.not. readInPressureFlag) newSegRHS = .true.
        call Notice('The LAmax-10dB critera for FAA Certification was not met AFTER the',&
                     'maximum A-weighted sound pressure Level (LAmax).')
      else
        if (HitMaxExpansions) then
          call Notice(' The maximum number of observer time range expansions has been',&
                      ' reached and LAmax-10dB could not be obtained AFTER the peak.')
        else
          call Notice('The source time range in your aperiodic data files',&
                      'is not sufficient to capture LAmax-10dB AFTER the peak.')
        end if
      end if
    end if
  end if 
  
  if(allocated(SPLArray)) then
    deallocate(SPLArray)
  end if
  if(allocated(LAarray)) then
    deallocate(LAarray)
  end if
end subroutine CalcObsSELData

function InsertBBTitle(obs, BBtitle, i, nbSources) result(title)
  implicit none
  type(observer):: obs
  character(len=1024):: BBtitle, title
  integer:: i, switch, nbSources

  if (thicknessNoiseFlag.and.loadingNoiseFlag) then
    switch = 3
  elseif (thicknessNoiseFlag.or.loadingNoiseFlag) then
    switch = 2
  else
    switch = 1
  end if

  if (i.lt.switch) then
    !ksb debug: title = obs%PPrime(i)%title
    title = PPRIMETITLEARRAY(i)
  elseif (i.eq.switch) then
    title = BBtitle
  else
    !ksb debug: title = obs%PPrime(i-1)%title
    title = PPRIMETITLEARRAY(i-1)
  end if
end function InsertBBTitle

subroutine InsertBBArray(obs, minFreq, maxFreq, tempMSP, i,j)
  implicit none
  type(observer):: obs
  type(freqDomain):: tempMSP
  type(freqDomain), pointer:: mspBB, mspArray
  real:: minFreq, maxFreq
  integer:: i, j, totalPos, switch

  ! We know that the order of the sources goes
  ! thickness noise, loading noise, total noise
  ! (with any of them being absent for a given case)
  if (thicknessNoiseFlag.and.loadingNoiseFlag) then
    switch = 3
  elseif (thicknessNoiseFlag.or.loadingNoiseFlag) then
    switch = 2
  else
    switch = 1
  end if

  totalPos = size(obs%mspBB,1)
  if (AtmAtten) then
    mspBB => obs%mspBBAtmAtten(totalPos,j)
  else
    mspBB => obs%mspBB(totalPos,j)
  end if

  if (i.lt.switch) then
    mspArray => obs%octaveFiltMSP(i,j)
  elseif (i.eq.switch) then
    mspArray => mspBB
  else
    mspArray => obs%octaveFiltMSP(i-1,j)
  end if
  call FilterFreqDomain(mspArray, minFreq,maxFreq, tempMSP)
end subroutine InsertBBArray


subroutine CalcObsOASPLdBData(obs,timeArray)
  implicit none
  type(observer),intent(inout)::obs
  real,dimension(:), intent(in)::timeArray
  real:: MSP
  integer::i,j,nbSources,nbSegs, totalPos
  nbSources = size(obs%msp,1)
  nbSegs    = size(obs%msp,2)
  if (.not.associated(obs%OASPLdB)) then
    allocate(obs%OASPLdB(nbSources))
    do i=1, nbSources
      nullify(obs%OASPLdB(i)%t)
      nullify(obs%OASPLdB(i)%f)
    end do
  end if
  if (associated(obs%mspBB)) totalPos = size(obs%mspBB,1)

  do i=1,nbSources
    call SetupTimeHistory(obs%OASPLdB(i), timeArray(1), timeArray(nbSegs), &
                          nbSegs, trim(GetFreqDomainTitle(obs%msp(i,1)))// &
                          'OASPL (dB, re20<greek>m</greek>Pa)')
    if (.not.associated(obs%OASPLdB(i)%f)) then
      allocate(obs%OASPLdB(i)%f(nbSegs))
    end if
    do j=1, nbSegs
      MSP = sum(obs%msp(i,j)%f)
      if (associated(obs%mspBB) .and. trim(GetFreqDomainTitle(obs%msp(i,j))).eq.'Total') then
        MSP = MSP + sum(obs%mspBB(totalPos,j)%f)
      end if
      obs%OASPLdB(i)%f(j) = ConvertMSPTodB(MSP)
    end do
  end do
end subroutine CalcObsOASPLdBData

subroutine CalcObsOASPLdBAData(obs,timeArray,octaveNumber)
  implicit none
  type(observer),intent(inout)::obs
  real,dimension(:), intent(in)::timeArray
  real, intent(in):: octaveNumber
  real:: AweightMSP
  integer::i,j,nbSources,nbSegs,totalPos
  
  type(freqDomain):: MSP
  
  nbSources = size(obs%msp,1)
  nbSegs    = size(obs%msp,2)
  if (.not.associated(obs%OASPLdBA)) then
    allocate(obs%OASPLdBA(nbSources))
    do i=1, nbSources
      nullify(obs%OASPLdBA(i)%t)
      nullify(obs%OASPLdBA(i)%f)
    end do    
  end if
  if (associated(obs%mspBB)) totalPos = size(obs%mspBB,1)

  do i=1,nbSources   
    call SetupTimeHistory(obs%OASPLdBA(i), timeArray(1),timeArray(nbSegs), &
                          nbSegs,trim(GetFreqDomainTitle(obs%msp(i,1)))// &
                          'OASPL (dBA, re20<greek>m</greek>Pa)')
    if (.not.associated(obs%OASPLdBA(i)%f)) then
      allocate(obs%OASPLdBA(i)%f(nbSegs))
    end if
    do j=1,nbSegs
      AweightMSP=AweightIntegrateOverFreqDomain(obs%msp(i,j))
      if (associated(obs%mspBB) .and. trim(GetFreqDomainTitle(obs%msp(i,j))).eq.'Total') then
        AweightMSP = AweightMSP + AweightIntegrateOverFreqDomain(obs%mspBB(totalPos,j))
      end if
      obs%OASPLdBA(i)%f(j) = ConvertMSPTodB(AweightMSP)
    end do
  end do
end subroutine CalcObsOASPLdBAData 

!**
!SUBROUTINE CalcObsFreqRangedBData
!  This routine calculates observer SPLdB versus segment time
!  for each frequency range specified in the namelist.
!**
subroutine CalcObsFreqRangedBData(obs,MSPRanges,timeArray)
  implicit none
  type(observer),intent(inout)::obs
  type(freqDomain),dimension(:,:,:),intent(in)::MSPRanges
  real,dimension(:),intent(in)::timeArray
  integer::nbSources,nbSegs,nbRanges,i,j,k
  nbSources = size(MSPRanges,1)
  nbSegs    = size(MSPRanges,2)
  nbRanges  = size(MSPRanges,3)
  if(.not.associated(obs%dBRanges)) then
    allocate(obs%dBRanges(nbSources,nbRanges)) 
    do i=1, nbSources
      do j=1, nbRanges
        nullify(obs%dBRanges(i,j)%t)
        nullify(obs%dBRanges(i,j)%f)
      end do
    end do
  end if
  do i=1,nbSources
    do k=1,nbRanges          
      call SetupTimeHistory(obs%dBRanges(i,k), timeArray(1),timeArray(nbSegs), &
                            nbSegs,trim(GetFreqDomainTitle(MSPRanges(i,1,1)))// &
                            'SPL (dB, re20<greek>m</greek>Pa)')                                        
      if (.not.associated(obs%dBRanges(i,k)%f)) then
        allocate(obs%dBRanges(i,k)%f(nbSegs))
      end if
      do j=1, nbSegs
        obs%dBRanges(i,k)%f(j) = ConvertMSPTodB(sum(MSPRanges(i,j,k)%f))
      end do
    end do  
  end do     
end subroutine CalcObsFreqRangedBData  

!**
!SUBROUTINE CalcObsFreqRangedBAData
!  This routine calculates observer SPLdBA versus segment time
!  for each frequency range specified in the namelist.
!**
subroutine CalcObsFreqRangedBAData(obs,MSPRanges,timeArray)
  implicit none
  type(observer),intent(inout)::obs
  !MSPRanges(noiseSource,segNum,rangeNum)
  type(freqDomain),dimension(:,:,:),intent(in)::MSPRanges
  real,dimension(:),intent(in)::timeArray
  integer::nbSources,nbSegs,nbRanges,i,j,k
  nbSources = size(MSPRanges,1)
  nbSegs    = size(MSPRanges,2)
  nbRanges  = size(MSPRanges,3)
  if(.not.associated(obs%dBARanges)) then
    allocate(obs%dBARanges(nbSources,nbRanges))
    do i=1, nbSources
      do j=1, nbRanges
        nullify(obs%dBARanges(i,j)%t)
        nullify(obs%dBARanges(i,j)%f)
      end do
    end do    
  end if
  do i=1,nbSources
    do k=1,nbRanges   
      call SetupTimeHistory(obs%dBARanges(i,k), timeArray(1),timeArray(nbSegs), &
                            nbSegs,trim(GetFreqDomainTitle(MSPRanges(i,1,1)))// &
                           'SPL (dBA, re20<greek>m</greek>Pa)')
      if (.not.associated(obs%dBARanges(i,k)%f)) then
        allocate(obs%dBARanges(i,k)%f(nbSegs))
      end if
      do j=1, nbSegs
        obs%dBARanges(i,k)%f(j) =  CalculateOASPLdBA(MSPRanges(i,j,k))
      end do
    end do  
  end do   
end subroutine CalcObsFreqRangedBAData

! **********************************
!SUBROUTINE CalcPNLT (modified version by K.S.Brentner)
!  calculates PNLT per 
!  FAR Part 36 , Appendix A2 Section A36.4 - Calculation of Effective Perceived
!  Noise Levels from Measured Data
!
!  If the logical PartH36Flag is .true., then the 50 Hz centerband frequency is used as the
!  starting 1/3 Octave bin.  Otherwise, the A36.4.1 definition, which starts with the 3rd band
!  (80 Hz centerband frequency) is used.
!
! **********************************
! Change history:
!  KSB 11/8/2012 - changed to match Bell Helicopter's procedure.  
!                  Added the logical variable PartH36Flag to start with either 50Hz or 80Hz 
!                  (1st of 3rd) 1/3 octave band.
!
! **********************************
!
subroutine CalcPNLT(splin,pnl,cMax,pnlt)
  use constantsModule, only: PartH36flag
  implicit none
  real, dimension(:), intent(in):: splin
  real, intent(in):: pnl
  real, intent(out):: pnlt
  real, intent(out):: cMax
!
! Local variables:
!
  integer:: i, nf, istart
  real, dimension(1:24):: s, spl, splp, splpp, sbar, f, c
  real, dimension(1:25):: sp
!
! check incomming SPL array to make sure it is the size we expect (24 1/3 octave
! bands from 50Hz to 10kHz - must assume that the first band is at 50Hz. 
! (or 80Hz - depends on PartH36Flag)
!
  nf = size(splin)
  spl = 0.
  if( nf >= 24 )then
    spl = splin(1:24)    ! if splin has more than 24 bands, just use the first 24.
  else
    spl(1:nf) = splin
    !ksb debug: (this is old, but sometime what to check usijng this) spl(nf+1:24) = 0.0
    !
    ! This padding ensures that if all 24 bands are not present there is not a tone
    ! correction for the last band.  KSB 8/25/2013
    !
    do i=nf,24 !ksb debug: nf+1,24  
      spl(i)=spl(i-1)-1.5
    end do
  end if
  if( PartH36Flag ) then
    istart = 1
  else
    istart = 3
  end if
! 
! A36.4.3  Correction for spectral irregularities
!
! Step 1:
!
  s(1:istart) = 0.
  do i=istart+1,24
    s(i) = spl(i)-spl(i-1)
  end do
!
! Steps 2-4:
!
  splp = spl
  do i=istart+2,23
    if ( abs(s(i)-s(i-1))>5. ) then
      if( s(i)>0. .and. s(i)>s(i-1) ) then
        splp(i) = 0.5*(spl(i-1)+spl(i+1)) 
      else if( s(i)<=0. .and. s(i-1)>0. ) then
        if (i.gt.2) splp(i-1) = 0.5*(spl(i-2)+spl(i))
      end if
    end if
  end do
  if( abs(s(24)-s(23))>5. ) then
    if( s(24)>0. .and. s(24)>s(23) ) then
      splp(24) = spl(23)+s(23) 
    end if
  end if      
!
! Step 5
!
  do i=istart+1,24
    sp(i)=splp(i)-splp(i-1)
  end do
  sp(istart) = sp(istart+1)
  sp(25)=sp(24)

!
! Step 6
!
  do i=istart,23
    sbar(i) = (sp(i)+sp(i+1)+sp(i+2))/3.
  end do
!
! Step 7
!
  splpp(istart) = spl(istart)
  do i=istart+1,24
    splpp(i) = splpp(i-1)+sbar(i-1)
  end do
!
! Step 8
!
  do i=istart,24    ! check indices
    f(i) = spl(i)-splpp(i) 
  end do
!
! Step 9
! 
  do i=istart,24
    if( f(i)>=20 ) then
      c(i) = 3.33333333
    else if( f(i)>=3 ) then
      c(i) = f(i)/6.
    else if( f(i)>=1.5 ) then
      c(i) = f(i)/3. - 0.5
    else
      c(i)=0.
    end if
    if( i>=11 .and. i<=21 ) then  ! 500<= freq <= 5000
      c(i) = 2.*c(i)
    end if
  end do 
!
! Step 10
!      
  cMax = maxval(c)
  pnlt = pnl + cMax
end subroutine CalcPNLT


!*************************************************************************************************************
!Abstract - Calculate atmsopheric absorption coefficient as a function of temperature 'Temp' in
!           degrees F and relative humidity 'Relhum' in percent for third octave bands from 'Band1'
!           (min 10, 10Hz) to 'Band2' (max 43, 20kHz).  'Alpha' in dB/1000 feet.
!
!Reference - SAE ARP 866A, modified to calculate from 10 to 20000 Hz.
!
!Limits - Temp: 1 to 100 degrees F, inclusive.
!         Relative Humidity: 1 to 100 percent, inclusive.
!
!Author - J. T. Brieger                                                                        Date - 06/19/86
!                                                                                               Rev - 11/13/08
!Rev - Benjamin A. Goldman								        Rev - 10/23/12
!********1*********2*********3*********4*********5*********6*********7*********8*********9*********0*********1
!
subroutine AtmosAtten_old(msp, radDistance)
  implicit none
  type(freqDomain), dimension(:,:):: msp
  real, dimension(:):: radDistance

  integer:: i, j, k, kMax, startIndex, endIndex, count, Icount, jMax, nbSources, nbSegs, nbFreq
  real, dimension(59):: Table
  real, dimension(BB_THIRD_OCTAVE_BANDS):: thirdOctArray
  real:: radDistanceFt, tempF, B, absHum, Humrat, freq, dBshift
  real:: Xa0, Xa1, Xa2, X01, X02, X12
  real:: Alprat, Amolmx, Alpmol, Alpcla, alpha, alpha40
  logical:: trip1, trip2, getAlprat

  data Table / 29.0  ,  0.0  ,  0.0  ,   .25 ,  0.315,  0.5  ,  0.7  ,  0.6  ,  0.84 ,  0.7  ,  0.93 ,  0.8  ,  0.975,  &
                0.9  ,  0.996,  1.0  ,  1.0  ,  1.1  ,   .97 ,  1.2  ,  0.9  ,  1.3  ,  0.84 ,  1.5  ,  0.75 ,  1.7  ,  &
                0.67 ,  2.0  ,  0.57 ,  2.3  ,  0.495,  2.5  ,  0.45 ,  2.8  ,  0.4  ,  3.0  ,  0.37 ,  3.3  ,  0.33 ,  &
                3.6  ,  0.3  ,  4.15 ,  0.26 ,  4.45 ,  0.245,  4.8  ,  0.23 ,  5.25 ,  0.22 ,  5.7  ,  0.21 ,  6.05 ,  &
                0.205,  6.5  ,  0.2  ,  7.0  ,  0.2  , 10.0  ,  0.2  /

  thirdOctArray = GetThirdOctFreqArray()
  trip1=.false.
  trip2=.false.

  nbSources = size(msp,1)
  nbSegs    = size(msp,2)
  nbFreq    = size(msp(1,1)%freq)

  ! Convert Temp from Kelvin to Fahrenheit
  tempF = 1.8*(Temperature-273.15)+32.

  ! Calculate Absolute Humidity
  !
  B=1.97274664+tempF*(-.02288074+tempF*(.00009589+tempF*(-3.E-7)))
  Abshum=10.0**(log10(RelHumidity)-B)
  do i=1,nbSources
    do j=1,nbSegs
      ! Convert radDistance from meters to feet, note that we divide by 1000.
      ! because the atmospheric attenutation cofficient (alpha) is in units of dB/1000ft
      radDistanceFt = 3.28084*radDistance(j)/1000.

      startIndex = 1
      endIndex = nbFreq
      do k=1,nbFreq
        if (msp(i,j)%freq(k).lt.44.7) startIndex = k+1
        if (msp(i,j)%freq(nbFreq-k+1).gt.22390.0) endIndex = nbFreq-k
      end do

      if (startIndex.ne.1 .or. endIndex.ne.nbFreq) then
        if ((startIndex.gt.(nbFreq) .or. endIndex.lt.1).and. .not.(trip1)) then
          call Notice (' The requested frequency range falls outside of what can be predicted by this routine. ',&
                        ' No atmospheric absorption correction will be applied.')
          trip1=.true.
        elseif (.not.(trip2)) then
          if (startIndex.gt.1) call Notice (' The lowest band frequency must be at least 44.7 Hz.  Atmospheric absorption',&
            ' correction will only be applied to bands with a lower frequency of at least 50.0 Hz.')

          if (endIndex.lt.nbFreq) call Notice (' The highest band frequency must be at most 22390.0 Hz.  Atmospheric absorption',&
            ' correction will only be applied to bands with an upper frequency of at most 22390.0 Hz.')
          trip2=.true.
        end if
      end if

      ! Calculate molecular atmospheric absorption coefficient
      !
      kMax = 1000000  ! I don't think we'll exceed 1 million frequency bins
      do k=startIndex,endIndex
        freq = msp(i,j)%freq(k)
        getAlprat = .false.
        Humrat=Abshum/sqrt(freq/1010.0)
        Icount=2
        if (Humrat.le.Table(2)) getAlprat = .true.

        if (.not. getAlprat) then
          do count=4,58,2
            Icount=count
            if (Table(count).gt.Humrat) then
              exit
            elseif (Table(count).eq.Humrat) then
              getAlprat=.true.
              exit
            end if
    	    if (count.eq.58) getAlprat = .true.
          end do	  
        end if
    
        if (getAlprat) then
          Alprat=Table(Icount+1)
        else
          count=Icount-2
          if (count.lt.4) count=count+2
          Xa1=Humrat-Table(count)
          Xa0=Humrat-Table(count-2)
          Xa2=Humrat-Table(count+2)
          X01=Table(count-2)-Table(count)
          X02=Table(count-2)-Table(count+2)
          X12=Table(count)-Table(count+2)
          Alprat=Table(count-1)*(Xa1/X01)*(Xa2/X02)
          Alprat=Alprat-Table(count+1)*(Xa0/X01)*(Xa2/X12)
          Alprat=Alprat+Table(count+3)*(Xa0/X02)*(Xa1/X12)
        end if
        Amolmx=10.0**(log10(freq)+.00468333*tempF-2.4215)
        Alpmol=Alprat*Amolmx

        ! Calculate classical atmospheric absorption coefficient
        !
        Alpcla=10.0**(2.05*log10(freq/1000.0)+(.000633*tempF)-1.45325)
    
        ! Calculate total sound attenuation coefficient
        !
        alpha=Alpcla+Alpmol
        if (freq/(2.**(1./6.)).ge.17780.0 .and. freq*(2.**(1./6.)).le.22390.0) then
          kMax = k
          alpha40 = alpha
        elseif (k.gt.kMax .and. alpha.gt.30.) then
          alpha=alpha40 !Added to limit divergence (OK as long as signal is above ambient)
        end if

        ! Since we are adding a 'delta dB' we need to first convert MSP to dB, add
        ! the atmospheric attenuation in dB, and then convert back to MSP.
        dBshift = ConvertMSPtodB(msp(i,j)%f(k)) - alpha*radDistanceFt
        msp(i,j)%f(k) = ConvertdBtoMSP(dBshift)
      end do
    end do
  end do
end subroutine AtmosAtten_old


!function [a] = atmAtten(Tin,Psin,hrin,dist,f)
! A function to return the atmospheric attenuation of sound due to the vibrational relaxation times of oxygen and nitrogen.
! NOTE:  This function does not account for spherical spreading!
!
! Usage: [a] = atmAtten(T,P,RH,d,f)
!               a - attenuation of sound for input parameters in dB
!               T - temperature in deg C
!               P - static pressure in Pa
!               RelHumidity - relative humidity in %
!               radDistance - distance of sound propagation in m
!               f - frequency of sound (may be a vector)
!
! Nathan Burnside 10/5/04
! AerospaceComputing Inc.
! nburnside@mail.arc.nasa.gov
!
! References:   Bass, et al., "Journal of Acoustical Society of America", (97) pg 680, January 1995.
!               Bass, et al., "Journal of Acoustical Society of America", (99) pg 1259, February 1996.
!              Kinsler, et al., "Fundamentals of Acoustics", 4th ed., pg 214, John Wiley & Sons, 2000.
!
! Rev - Benjamin A. Goldman				       Rev - 11/28/12
subroutine AtmosAtten(msp, radDistance)
  implicit none
  type(freqDomain), dimension(:,:):: msp
  type(timehistory), dimension(:), pointer:: radDistance

  real:: To1,To,T,Psat,h,F,FrN,FrO,alpha,dBatten,dbShift,nonZero
  integer:: i,j,k,nbSources,nbSegs,nbFreq

  nonZero   = tiny(radDistance(1)%f(1))
  nbSources = size(msp,1)
  nbSegs    = size(msp,2)
  nbFreq    = size(msp(1,1)%freq)

  To1 = 273.15 ! triple point in K
  To = 293.15 ! ref temp in K

  ! calculate saturation pressure
  Psat = 10.0**(10.79586*(1.0-(To1/Temperature))-5.02808*log10(Temperature/To1) &
        +1.50474e-4*(1.0-10.0**(-8.29692*((Temperature/To1)-1.0)))-4.2873e-4*   &
        (1.0-10.0**(-4.76955*((To1/Temperature)-1.0)))-2.2195983)

  ! calculate the absolute humidity
  h = RelHumidity*Psat

  ! Scaled relaxation frequency for Nitrogen
  FrN = (To/Temperature)**(0.5)*(9.0+280.0*h*exp(-4.17*((To/Temperature)**(1.0/3.0)-1.0)))

  ! scaled relaxation frequency for Oxygen
  FrO = (24.0+4.04e4*h*(0.02+h)/(0.391+h))

  do i=1,nbSources
    do j=1,nbSegs
      do k=1,nbFreq
        F = msp(i,j)%freq(k) ! frequency per atm

        ! attenuation coefficient in nepers/m
        alpha = F**2.0*(1.84e-11*((Temperature/To)**(0.5)) + &
                (Temperature/To)**(-2.5)*(1.275e-2*exp(-2239.1/Temperature)/(FrO+ &
                (F**2.0)/FrO) + 1.068e-1*exp(-3352.0/Temperature)/(FrN+(F**2.0)/FrN)))
	    dBatten = 10*log10(exp(2.0*alpha))*radDistance(j)%f(1)
        
        dBshift = ConvertMSPtodB(msp(i,j)%f(k)) - dBatten
        msp(i,j)%f(k) = ConvertdBtoMSP(dBshift)
        if (msp(i,j)%f(k).eq.0) msp(i,j)%f(k)=nonZero
      end do
    end do
  end do
end subroutine AtmosAtten


subroutine ResetPPrimeAndIndices(obs)
  implicit none
  type(observer),dimension(:),intent(inout)::obs
  integer::i  
  if((.not.thicknessNoiseFlag).and.(.not.loadingNoiseFlag))then
    do i=1,size(obs)
      call CopyTimeHistory(obs(i)%pPrime(TOTAL_APTH), obs(i)%pPrime(1))
    end do
    PPRIMETITLEARRAY(1)='Total'
    if(pressureGradientFlag.or.pressureGradient1AFlag)then
      do i=1,size(obs)
        call CopyTimeHistory(obs(i)%pPrime(TOTAL_PGX), obs(i)%pPrime(2))
        call CopyTimeHistory(obs(i)%pPrime(TOTAL_PGY), obs(i)%pPrime(3))
        call CopyTimeHistory(obs(i)%pPrime(TOTAL_PGZ), obs(i)%pPrime(4))
        call DestroyTimeHistories(obs(i)%pPrime(5:12))       
      end do
      PPRIMETITLEARRAY(2)='TotalPGX'
      PPRIMETITLEARRAY(3)='TotalPGY'
      PPRIMETITLEARRAY(4)='TotalPGZ'            
    else
      do i=1,size(obs)
        call DestroyTimeHistories(obs(i)%pPrime(2:3))
      end do
    end if      
  else if((.not.thicknessNoiseFlag).and.loadingNoiseFlag)then
    do i=1,size(obs)
      call CopyTimeHistory(obs(i)%pPrime(LOAD_APTH), obs(i)%pPrime(1))
      call CopyTimeHistory(obs(i)%pPrime(TOTAL_APTH), obs(i)%pPrime(2))
    end do
    PPRIMETITLEARRAY(1)='Loading'
    PPRIMETITLEARRAY(2)='Total'      
    if(pressureGradientFlag.or.pressureGradient1AFlag)then
      do i=1,size(obs)
        call CopyTimeHistory(obs(i)%pPrime(LOAD_PGX), obs(i)%pPrime(3))
        call CopyTimeHistory(obs(i)%pPrime(LOAD_PGY), obs(i)%pPrime(4))
        call CopyTimeHistory(obs(i)%pPrime(LOAD_PGZ), obs(i)%pPrime(5))
        call CopyTimeHistory(obs(i)%pPrime(TOTAL_PGX), obs(i)%pPrime(6))
        call CopyTimeHistory(obs(i)%pPrime(TOTAL_PGY), obs(i)%pPrime(7))
        call CopyTimeHistory(obs(i)%pPrime(TOTAL_PGZ), obs(i)%pPrime(8))
        call DestroyTimeHistories(obs(i)%pPrime(9:12))
      end do
      PPRIMETITLEARRAY(3)='LoadingPGX'
      PPRIMETITLEARRAY(4)='LoadingPGY'
      PPRIMETITLEARRAY(5)='LoadingPGZ' 
      PPRIMETITLEARRAY(6)='TotalPGX'
      PPRIMETITLEARRAY(7)='TotalPGY'
      PPRIMETITLEARRAY(8)='TotalPGZ'               
    else
      do i=1,size(obs)
        call DestroyTimeHistory(obs(i)%pPrime(3))
      end do
    end if      
  else
    do i=1,size(obs)
      call CopyTimeHistory(obs(i)%pPrime(THICK_APTH), obs(i)%pPrime(1))
      call CopyTimeHistory(obs(i)%pPrime(TOTAL_APTH), obs(i)%pPrime(2))
    end do
    PPRIMETITLEARRAY(1)='Thickness'        
    PPRIMETITLEARRAY(2)='Total'        
    if(pressureGradientFlag.or.pressureGradient1AFlag)then
      do i=1,size(obs)
        call CopyTimeHistory(obs(i)%pPrime(THICK_PGX), obs(i)%pPrime(3))
        call CopyTimeHistory(obs(i)%pPrime(THICK_PGY), obs(i)%pPrime(4))
        call CopyTimeHistory(obs(i)%pPrime(THICK_PGZ), obs(i)%pPrime(5))
        call CopyTimeHistory(obs(i)%pPrime(TOTAL_PGX), obs(i)%pPrime(6))
        call CopyTimeHistory(obs(i)%pPrime(TOTAL_PGY), obs(i)%pPrime(7))
        call CopyTimeHistory(obs(i)%pPrime(TOTAL_PGZ), obs(i)%pPrime(8))
        call DestroyTimeHistories(obs(i)%pPrime(9:12))
      end do
      PPRIMETITLEARRAY(3)='ThicknessPGX'
      PPRIMETITLEARRAY(4)='ThicknessPGY'
      PPRIMETITLEARRAY(5)='ThicknessPGZ' 
      PPRIMETITLEARRAY(6)='TotalPGX'
      PPRIMETITLEARRAY(7)='TotalPGY'
      PPRIMETITLEARRAY(8)='TotalPGZ'           
    else
      do i=1,size(obs)
        call DestroyTimeHistory(obs(i)%pPrime(3))
      end do
    end if       
  end if
  call ResetPPrimeIndices()
end subroutine ResetPPrimeAndIndices

subroutine FilterPPrime(pPrime,hpf,lpf)
  implicit none
  type(timeHistory),intent(inout)::pPrime
  real,intent(in)::hpf,lpf
  complex,dimension(:),allocatable:: complexDFTResult
  real::length
  allocate(complexDFTResult(GetNt(pPrime)/2+1))
  call CalcFiltComplexDFTFromPPrime(pPrime,complexDFTResult,hpf,lpf)
  call ConvertFrequencyToTime(complexDFTResult, pPrime%f)
  !Running FFTW and then iFFTW scales the pressure by nt
  !  so we have to correct it:
  length = size(pPrime%f)
  call SettHFArray(pPrime,pPrime%f/length)
  deallocate(complexDFTResult)
end subroutine FilterPPrime

subroutine CalcObsComplexPressure(obs,freqArray,startBin,df,hpf,lpf)
  implicit none
  type(observer),intent(inout)::obs
  real,intent(in)::hpf,lpf,df
  real::dt
  real,dimension(:),intent(in)::freqArray
  integer,intent(in)::startBin
  complex,dimension(:),allocatable::complexDFTResult,complexPressure,FSCCorrection
  integer::nbSources,nbSegs,nf,nt,i,j
  nbSources = GetNewSizeOfPPrime()
  nbSegs    = size(obs%pPrimeSegs,2)
  nf        = size(freqArray)
  dt        = GetDt(obs%pPrimeSegs(1,1))   
  nt        = GetNt(obs%pPrimeSegs(1,1))
  if (.not.associated(obs%realCompP)) then
    allocate(obs%realCompP(nbSources,nbSegs))
    do i=1, nbSources
      do j=1, nbSegs
        nullify(obs%realCompP(i,j)%f)
        nullify(obs%realCompP(i,j)%freq)
      end do
    end do
  end if
  if (.not.associated(obs%imagCompP)) then
    allocate(obs%imagCompP(nbSources,nbSegs))
    do i=1, nbSources
      do j=1, nbSegs
        nullify(obs%imagCompP(i,j)%f)
        nullify(obs%imagCompP(i,j)%freq)
      end do
    end do
  end if
  allocate(complexDFTResult(nt/2+1), complexPressure(nf))
  if(FSCOutputFlag)then
    !Calculate the FSCCorrection Array for FSC cases.
    allocate(FSCCorrection(nf))
    FSCCorrection(:) = exp(cmplx(0,1)*2*pi*freqArray(:)*obs%tMin)
  end if
  do i=1,nbSources
    do j=1,nbSegs   
      call CalcFiltComplexDFTFromPPrime(obs%pPrimeSegs(i,j), &
                                        complexDFTResult, &
                                        hpf,lpf)    
      complexPressure = complexDFTResult(startBin:startBin+nf-1)*dt*2.*df
      if(FSCOutputFlag)then
        complexPressure(:) = complexPressure(:)/FSCCorrection(:)
      end if
      !Handle the DC offset differently, as defined in the user manual:
      if(freqArray(1).eq.0.)then        
        complexPressure(1) = complexDFTResult(1)*dt*df
      end if
      call SetMinFreq(obs%realCompP(i,j),freqArray(1))
      call SetNf(obs%realCompP(i,j),nf)
      call SetMaxFreq(obs%realCompP(i,j),freqArray(size(freqArray)))
      call SetFreqArray(obs%realCompP(i,j),freqArray)
      call SetFreqDomainTitle(obs%realCompP(i,j),'Re_' // &
                              GettHTitle(obs%pPrimeSegs(i,j)))
      call CopyFreqDomain(obs%realCompP(i,j), obs%imagCompP(i,j))
      call SetFreqDomainFArray(obs%realCompP(i,j),real(complexPressure(:)))
      call SetFreqDomainFArray(obs%imagCompP(i,j),imag(complexPressure(:)))
      call SetFreqDomainTitle(obs%imagCompP(i,j),'Im_' // &
                              GettHTitle(obs%pPrimeSegs(i,j)))
    end do
  end do
  deallocate(complexDFTResult,complexPressure)
  if(FSCOutputFlag)then
    deallocate(FSCCorrection)
  end if  
end subroutine CalcObsComplexPressure

!**
!SUBROUTINE CalcFiltComplexDFTFromPPrime
!  This routine calculates a filtered complex pressure FFT
!  from a time history, pPRime, according to the high and low pass
!  frequency hpf and lpf.
!**
subroutine CalcFiltComplexDFTFromPPrime(pPrime,complexDFTResult,hpf,lpf)
  type(timeHistory),intent(in)::pPrime  
  real,intent(in)::hpf,lpf
  complex,dimension(:),intent(inout)::complexDFTResult
  integer::i
  real::df  
  call ConvertTimeToFrequency(pPrime%f,complexDFTResult=complexDFTResult)
  !Right now the units of complexDFTResult are "pascals"
  !The following filters the result based on high and low pass frequency values
  !  specified in the namelist:
  df=1./GetPeriod(pPrime)  
  do i=0,size(complexDFTResult)-1
    if((real(i)*df).lt.(hpf-epsilon))then
      complexDFTResult(i+1)=0.
    else if ((real(i)*df).gt.(lpf+epsilon))then
      complexDFTResult(i+1)=0.
    end if
  end do
end subroutine CalcFiltComplexDFTFromPPrime  

subroutine CalcObsdBFromMSP(dB,msp)
  implicit none  
  type(freqDomain),dimension(:,:),pointer::dB
  type(freqDomain),dimension(:,:),pointer::msp
  integer::i,j,k,nbSources,nbSegs
  nbSources = size(msp,1)
  nbSegs    = size(msp,2)  
  if(.not.associated(dB)) then
    allocate(dB(nbSources,nbSegs))
    do i=1, nbSources
      do j=1, nbSegs
        nullify(dB(i,j)%freq)
        nullify(dB(i,j)%f)
      end do
    end do
  end if
  do i=1,nbSources
    do j=1,nbSegs
      call CopyFreqDomain(msp(i,j), dB(i,j))
      !Do not convert the DC offset to dB:
      if(GetFreq(msp(i,j),1).lt.epsilon)then
        do k=2,size(dB(i,j)%f)
          dB(i,j)%f(k)=ConvertMSPTodB(msp(i,j)%f(k))
        end do  
      else
        do k=1,size(dB(i,j)%f)
          dB(i,j)%f(k)=ConvertMSPTodB(msp(i,j)%f(k))
        end do
      end if
      call SetFreqDomainTitle(dB(i,j),trim(PPRIMETITLEARRAY(i))//'_dB')
    end do
  end do
end subroutine CalcObsdBFromMSP

subroutine CalcObsdBAFromMSP(dBA,msp)
  implicit none  
  type(freqDomain),dimension(:,:),pointer::dBA
  type(freqDomain),dimension(:,:),pointer::msp
  integer::i,j,k,nbSources,nbSegs,start
  nbSources = size(msp,1)
  nbSegs    = size(msp,2)  
  if(.not.associated(dBA)) then
    allocate(dBA(nbSources,nbSegs))
    do i=1, nbSources
      do j=1, nbSegs
        nullify(dBA(i,j)%freq)
        nullify(dBA(i,j)%f)
      end do
    end do
  end if
  if(GetFreq(msp(1,1),1).lt.epsilon)then
    start=2
  else
    start=1
  end if
  do i=1,nbSources
    do j=1,nbSegs
      call CopyFreqDomain(msp(i,j), dBA(i,j))
      do k=start,size(dBA(i,j)%f)
        dBA(i,j)%f(k)=ConvertMSPTodBA(msp(i,j)%f(k), msp(i,j)%freq(k))
      end do
      call SetFreqDomainTitle(dBA(i,j),trim(PPRIMETITLEARRAY(i))//'_dBA')
    end do
  end do
end subroutine CalcObsdBAFromMSP

!**
!SUBROUTINE CalcOctaveFiltMSP
!  This routine calculates octave filtered MSP for an observer
!  based on the octave number.
!**
!ksb debug: subroutine CalcOctaveFiltMSP(obs,octaveNumber)
subroutine CalcOctaveFiltMSP(obs,octaveNumber,highPassFreq,lowPassFreq)
  implicit none
  type(observer),intent(inout)::obs
  type(freqDomain):: tempMSP
  type(freqDomain), dimension(:), allocatable:: BBtempMSP
  type(freqDomain), dimension(:,:), pointer:: mspBBarray
  real,intent(in)::octaveNumber, highPassFreq, lowPassFreq !ksb debug: added highPassFreq and lowPassFreq
  real:: maxDF, maxBB, minF, maxF !ksb debug: added minF and maxF
  integer::i,j,nbSources,nbSegs,totalPos

  nbSources = size(obs%msp,1)
  nbSegs = size(obs%msp,2)
  if (.not.associated(obs%octaveFiltMSP)) allocate(obs%octaveFiltMSP(nbSources,nbSegs))

  if (associated(obs%mspBB)) then
    if (AtmAtten) then
      mspBBarray => obs%mspBBAtmAtten
    else
      mspBBarray => obs%mspBB
    end if

    maxDF = maxval(obs%msp(1,1)%freq)
    maxBB = maxval(mspBBarray(1,1)%freq)
    totalPos = size(mspBBarray,1)
  end if

  do i=1,nbSources
    do j=1,nbSegs
      call CreateMSPPerBand(obs%msp(i,j),octaveNumber, tempMSP)
      call SetFreqDomainTitle(tempMSP, &
                              trim(GetFreqDomainTitle(obs%msp(i,j)))//'_octFilt')
      if (.not.associated(obs%mspBB)) then
        call copyFreqDomain(tempMSP, obs%octaveFiltMSP(i,j))
      elseif (associated(obs%mspBB) .and. octaveNumber.eq.3) then
        if(maxBB.gt.maxDF) then
          !ksb debug: call expandMSPArray(obs%octaveFiltMSP(i,j), tempMSP, mspBBarray(totalPos,j), maxDF, maxBB)
          ! minF = max( min(minDF,minBB), highPassFreq)  ! there seems to be an implicit assumption that minDF will
          ! always be less than minBB.  Also minDF should reflect the highPassFrequency - so maybe we only need
          ! to change the upper limit to either maxBB or the lowPassFrequency.
          !maxF = min( max(maxDF,maxBB), lowPassFreq) ! this is more general, but we already know maxBB > maxDF here.
          maxF = min( maxBB, lowPassFreq )
          call expandMSPArray(obs%octaveFiltMSP(i,j), tempMSP, mspBBarray(totalPos,j), maxDF, maxF)
        else
          call copyFreqDomain(tempMSP, obs%octaveFiltMSP(i,j))
        end if
        if (trim(obs%msp(i,j)%title).eq.'Total') obs%octaveFiltMSP(i,j) = addThirdOctMSP(obs%octaveFiltMSP(i,j), &
            mspBBarray(totalPos,j), .true.)
      elseif (octaveNumber.ne.3) then
        call Notice(' Broadband noise can only be added to the total noise when using', &
                    ' 1/3rd octave bands (octaveNumber = 3).  The results from this case', &
                    ' will not include broadband noise.')
      end if
    end do
  end do
end subroutine CalcOctaveFiltMSP

!**
!SUBROUTINE CalcTotalNoise
!  This routine calculates the total pressure time histories for an observer.
!**
subroutine CalcTotalNoise(obs)
  implicit none
  type(observer),intent(inout)::obs

  if(.not.associated(obs%pPrime(TOTAL_APTH)%f)) then
    allocate(obs%pPrime(TOTAL_APTH)%f(obs%pPrime(TOTAL_APTH)%nt))
  end if
  obs%pPrime(TOTAL_APTH)%f=obs%pPrime(THICK_APTH)%f+obs%pPrime(LOAD_APTH)%f
  call SettHTitle(obs%pPrime(TOTAL_APTH),'Total')
  if(pressureGradientFlag.or.pressureGradient1AFlag)then
    if(.not.associated(obs%pPrime(TOTAL_PGX)%f)) then
      allocate(obs%pPrime(TOTAL_PGX)%f(obs%pPrime(TOTAL_PGX)%nt))
    end if
    obs%pPrime(TOTAL_PGX)%f=obs%pPrime(THICK_PGX)%f+obs%pPrime(LOAD_PGX)%f
    call SettHTitle(obs%pPrime(TOTAL_PGX),'TotalPGX')
    if(.not.associated(obs%pPrime(TOTAL_PGY)%f)) then
      allocate(obs%pPrime(TOTAL_PGY)%f(obs%pPrime(TOTAL_PGY)%nt))
    end if
    obs%pPrime(TOTAL_PGY)%f=obs%pPrime(THICK_PGY)%f+obs%pPrime(LOAD_PGY)%f
    call SettHTitle(obs%pPrime(TOTAL_PGY),'TotalPGY')
    if(.not.associated(obs%pPrime(TOTAL_PGZ)%f)) then
      allocate(obs%pPrime(TOTAL_PGZ)%f(obs%pPrime(TOTAL_PGZ)%nt))
    end if
    obs%pPrime(TOTAL_PGZ)%f=obs%pPrime(THICK_PGZ)%f+obs%pPrime(LOAD_PGZ)%f
    call SettHTitle(obs%pPrime(TOTAL_PGZ),'TotalPGZ')
  end if
end subroutine CalcTotalNoise

subroutine CalcObsMSPFromComplexPressure(obs)
  implicit none
  type(observer),intent(inout)::obs
  integer::i,j,nf,nbSegs,nbSources
  character(len=4096):: BBTitle
  
  !We only calculate MSP for pressure data, not pressure gradient data:
  if(totalNoiseFlag)then
    nbSources=TOTAL_APTH
  else if(loadingNoiseFlag)then
    nbSources=LOAD_APTH
  else if(thicknessNoiseFlag)then
    nbSources=THICK_APTH
  end if
  nbSegs=size(obs%realCompP,2)
  if (associated(obs%msp)) then
    deallocate(obs%msp)
  end if
  allocate(obs%msp(nbSources,nbSegs))
!ksb debug:
  obs%nbSources = nbSources
  if (globalPeggNoiseFlag.or.globalBPMnoiseFlag) then
    obs%nbSources=obs%nbSources+1
  end if
!end ksb debug:
  do i=1, nbSources
    do j=1, nbSegs
      nullify(obs%msp(i,j)%freq)
      nullify(obs%msp(i,j)%f)
    end do
  end do
  nf=GetNf(obs%realCompP(1,1))  
  do i=1,nbSources
    do j=1,nbSegs
      call CopyFreqDomain(obs%realCompP(i,j), obs%msp(i,j))       
      call SetFreqDomainFArray(obs%msp(i,j), &  
           0.5*(obs%realCompP(i,j)%f**2.0 + obs%imagCompP(i,j)%f**2.0))
      call SetFreqDomainTitle(obs%msp(i,j),trim(PPRIMETITLEARRAY(i)))                
    end do
  end do
end subroutine CalcObsMSPFromComplexPressure

subroutine CalcObsPhaseFromComplexPressure(obs)
  implicit none
  type(observer),intent(inout)::obs
  integer::i,j,nbSources,nbSegs           
  !We only calculate phase for pressure data, not pressure gradient data:
  if(totalNoiseFlag)then
    nbSources=TOTAL_APTH
  else if(loadingNoiseFlag)then
    nbSources=LOAD_APTH
  else if(thicknessNoiseFlag)then
    nbSources=THICK_APTH
  end if
  nbSegs=size(obs%realCompP,2)
  if (.not.associated(obs%phase)) then
    allocate(obs%phase(nbSources,nbSegs))  
    do i=1, nbSources
      do j=1, nbSegs
        nullify(obs%phase(i,j)%freq)
        nullify(obs%phase(i,j)%f)
      end do
    end do
  end if     
  do i=1,nbSources
    do j=1,nbSegs
      call CopyFreqDomain(obs%realCompP(i,j), obs%phase(i,j))
      call SetFreqDomainFArray(obs%phase(i,j), -atan2(obs%imagCompP(i,j)%f, &
                       obs%realCompP(i,j)%f))
      call SetFreqDomainTitle(obs%phase(i,j),trim(PPRIMETITLEARRAY(i))//'_phase')             
    end do
  end do  
end subroutine CalcObsPhaseFromComplexPressure
               
subroutine CreateObsWavFiles(obs,path)
  implicit none
  type(observer),intent(in)::obs
  character(len=*),intent(in)::path
  if(thicknessNoiseFlag)then
    call CreateWavFile(obs%pPrime(THICK_APTH)%f,&
                       obs%pPrime(THICK_APTH)%t,  &
                       trim(path)//'_thickness')
  end if
  if(loadingNoiseFlag)then
    call CreateWavFile(obs%pPrime(LOAD_APTH)%f,&
                       obs%pPrime(LOAD_APTH)%t,  &
                       trim(path)//'_loading')
  end if
  if(totalNoiseFlag)then
    call CreateWavFile(obs%pPrime(TOTAL_APTH)%f,&
                       obs%pPrime(TOTAL_APTH)%t,  &
                       trim(path)//'_total')
  end if  
end subroutine CreateObsWavFiles               
                                     
subroutine DestroyObserver(obs)
  implicit none
  type(observer),intent(inout)::obs
  call DestroyTHArray(obs%pPrime)
  nullify(obs%pPrime)
  call Destroy2DTHArray(obs%pPrimeSegs)
  nullify(obs%pPrimeSegs)
  call Destroy2DFDArray(obs%realCompP) 
  nullify(obs%realCompP)
  call Destroy2DFDArray(obs%imagCompP)   
  nullify(obs%imagCompP)
  call DestroyTHArray(obs%OASPLdB)
  nullify(obs%OASPLdB)
  call DestroyTHArray(obs%OASPLdBA)
  nullify(obs%OASPLdBA)
  call DestroyTHArray(obs%PNL)
  nullify(obs%PNL)
  call DestroyTHArray(obs%PNLT)
  nullify(obs%PNLT)
  call DestroyTHArray(obs%radDistance)
  nullify(obs%radDistance)
  call Destroy2DTHArray(obs%dBRanges)
  nullify(obs%dBRanges)
  call Destroy2DTHArray(obs%dBARanges)
  nullify(obs%dBARanges)
  call Destroy2DFDArray(obs%phase) 
  nullify(obs%phase)
  call Destroy2DFDArray(obs%msp) 
  nullify(obs%msp)
  call Destroy2DFDArray(obs%mspBB) 
  nullify(obs%mspBB)
  call Destroy2DFDArray(obs%mspBBAtmAtten)
  nullify(obs%mspBBAtmAtten)
  call Destroy2DFDArray(obs%dB) 
  nullify(obs%dB)
  call Destroy2DFDArray(obs%dBA) 
  nullify(obs%dBA)
  call Destroy2DFDArray(obs%octaveFiltMSP) 
  nullify(obs%octaveFiltMSP)
  call Destroy2DFDArray(obs%octaveFiltdB) 
  nullify(obs%octaveFiltdB)
  call Destroy2DFDArray(obs%octaveFiltdBA) 
  nullify(obs%octaveFiltdBA)
  if(associated(obs%iBlankArray))then
    deallocate(obs%iBlankArray)
    nullify(obs%iBlankArray)
  end if
  if(associated(obs%EPNL))then
    deallocate(obs%EPNL)
    nullify(obs%EPNL)
  end if
  if(associated(obs%SEL))then
    deallocate(obs%SEL)
    nullify(obs%SEL)
  end if
end subroutine DestroyObserver

subroutine DestroyTHArray(array)
  type(timeHistory), dimension(:), pointer::array
  integer::i
  if(associated(array)) then
    do i=1, size(array)
      call DestroyTimeHistory(array(i))
    end do
  end if
end subroutine DestroyTHArray

subroutine Destroy2DTHArray(array)
  type(timeHistory), dimension(:,:), pointer::array
  integer::i, j
  if(associated(array)) then
    do i=1, size(array, 1)
      do j=1, size(Array, 2)
        call DestroyTimeHistory(array(i, j))
      end do
    end do
  end if
end subroutine Destroy2DTHArray

subroutine Destroy2DFDArray(array)
  type(freqDomain), dimension(:,:), pointer::array
  integer::i, j
  if(associated(array)) then
    do i=1, size(array, 1)
      do j=1, size(Array, 2)
        call DestroyFreqDomain(array(i, j))
      end do
    end do
  end if
end subroutine Destroy2DFDArray

function GetObserverCoordinates(obs) result(coordinates)
  implicit none
  type(observer),intent(in)::obs
  type(vector)::coordinates
  coordinates = obs%coordinates
end function GetObserverCoordinates

function GetObserverNt(obs)result(nt)
  implicit none
  type(observer),intent(in)::obs
  integer::nt
  nt=obs%nt
end function GetObserverNt

function GetObserverTitle(obs)result(title)
  implicit none
  type(observer),intent(in)::obs
  character(len=4096)::title
  title = trim(obs%title)
end function GetObserverTitle

function GetObserverTMax(obs)result(tMax)
  implicit none
  type(observer),intent(in)::obs
  real::tMax
  tMax = obs%tMax
end function GetObserverTMax

function GetObserverTMin(obs)result(tMin)
  implicit none
  type(observer),intent(in)::obs
  real::tMin
  tMin = obs%tMin
end function GetObserverTMin

function GetObserverDt(obs)result(Dt)
  implicit none
  type(observer),intent(in)::obs
  real::Dt
  Dt = obs%Dt
end function GetObserverDt

subroutine SetObserverTitle(obs,title)
  implicit none
  type(observer),intent(inout)::obs
  character(len=*),intent(in)::title
  obs%title = trim(title)
end subroutine SetObserverTitle

subroutine SetObserverTMin(obs,tMin)
  implicit none
  type(observer),intent(inout)::obs
  real,intent(in)::tMin
  obs%tMin = tMin
end subroutine SetObserverTMin

subroutine SetObserverTMax(obs,tMax)
  implicit none
  type(observer),intent(inout)::obs
  real,intent(in)::tMax
  obs%tMax = tMax
end subroutine SetObserverTMax

subroutine SetObserverNt(obs,nt)
  implicit none
  type(observer),intent(inout)::obs
  integer,intent(in)::nt
  obs%Nt = nt
end subroutine SetObserverNt

subroutine SetObserverDt(obs,dt)
  implicit none
  type(observer),intent(inout)::obs
  real,intent(in)::dt
  obs%Dt = dt
end subroutine SetObserverDt

subroutine SetObserverIndex(obs,index)
  implicit none
  type(observer),intent(inout)::obs
  integer,intent(in)::index
  obs%ObsIndex = index
end subroutine SetObserverIndex

!**
!SUBROUTINE WindowTimeHistory
!  This function determines the window function specified by the user,
!  and creates a new time history which is identical to the original
!  time history except for the range array. which is windowed.
!**
subroutine WindowTimeHistory(tH,windowFunction)
  implicit none
  type(timeHistory),intent(inout)::tH
  character(len=4096),intent(in)::windowFunction
  character(len=4096)::string
  real::correctionFactor
  real,dimension(:),allocatable::array
  integer:: i
  allocate(array(tH%nt))
  string = trim(windowFunction)
  call strToUpper(string,len_trim(string))
  selectcase(trim(string))
    case('BLACKMAN WINDOW')
      call CreateBlackmanWindow(array)
    case('FLAT TOP WINDOW')
      call CreateFlatTopWindow(array)
    case('HANNING WINDOW')
      call CreateHanningWindow(array)
    case('HAMMING WINDOW')
      call CreateHammingWindow(array)                        
    case default
      call Error('Unrecognized windowFunction in namelist.')
  end select
  correctionFactor = sqrt(tH%nt/sum(array**2.0))
  call SettHFArray(tH,correctionFactor*tH%f*array)
  deallocate(array)
end subroutine WindowTimeHistory

subroutine WriteNameFile(titleArray,path)
  implicit none
  character(len=*),dimension(:),intent(in)::titleArray
  character(len=*),intent(in)::path
  integer::unit,i
  unit = GetStreamNumber()
  open(unit=unit,file=trim(path),status='replace')
    do i=1,size(titleArray)
      write(unit,*)trim(titleArray(i))
    end do
  close(unit)    
 end subroutine WriteNameFile

!**
!SUBROUTINE WriteSingleObsFreqDomainData
!This routine creates a Tecplot (*.tec) file or similar file
!  that stores spectrum data for a single observer  point.  
!  The file that is created is named "path".tec. This file
!  contains thickness, loading, and total spectrums, depending 
!  on the noise flags declared in the namelist.  
!**
subroutine WriteSingleObsFreqDomainData(fD1,fD2,path)
  implicit none
  type(freqDomain),dimension(:,:),intent(in),optional::fD1,fD2
  character(len=*),intent(in)::path
  integer::i,j,segNum,nbSegs,dim1,nf
  real,dimension(:),allocatable::freqArray
  real, dimension(:,:),allocatable::f  
  character(len=4096),dimension(:),allocatable::fTitles  
  !First we have to determine what the frequency array is for output:
  if(present(fD1))then
    nf=GetNf(fD1(1,1))
    allocate(freqArray(nf))
    freqArray=fD1(1,1)%freq
    nbSegs=size(fD1,2)    
    else if(present(fD2))then
      nf=GetNf(fD2(1,1))
      allocate(freqArray(nf))
      freqArray=fD2(1,1)%freq
      nbSegs=size(fD2,2)
      else
        call Warning('WriteSingleObsFreqDomainData called with no data.')
  end if
  dim1=0 
  if(present(fD1))then
    dim1=dim1+size(fD1,1)
  end if
  if(present(fD2))then
    dim1=dim1+size(fD2,1)
  end if  
  do segNum=1,nbSegs   
    allocate(f(dim1,nf),fTitles(dim1))
    if(present(fD1))then
      do i=0,size(fD1,1)-1
        f(i+1,:)=fD1(i+1,segNum)%f
        fTitles(i+1)=GetFreqDomainTitle(fD1(i+1,segNum))
      end do
    end if    
    if(present(fD2))then
      do j=1,size(fD2,1)
        f(i+j,:)=fD2(j,segNum)%f
        fTitles(i+j)=GetFreqDomainTitle(fD2(j,segNum))
      end do
    end if
    if(nbSegs==1)then 
      call WriteOutTecplotData(trim(path)//'.tec',freqArray,'Frequency', &
                                                  f,         fTitles)
    else
      call WriteOutTecplotData(trim(path)//'_segment'//trim(IntegerToString(segNum))// &
                               '.tec',freqArray,'Frequency',f,fTitles)
    end if
    if(allocated(f))deallocate(f)
    if(allocated(fTitles))deallocate(fTitles)    
  end do
  if(allocated(freqArray))deallocate(freqArray)                                                   
end subroutine WriteSingleObsFreqDomainData  

subroutine WriteSingleObsFreqRangeData(obs,rangeTitles,path)
  implicit none
  type(observer),intent(in)::obs
  character(len=*),intent(in)::path
  integer::nbRanges,i,nbSources,j
  type(timeHistory),dimension(:),pointer::f
  character(len=*),dimension(:),pointer::rangeTitles
  nullify(f)
  if(SPLdBFlag.and.SPLdBAFlag)then
    nbSources = size(obs%dBRanges,1)
    allocate(f(2*nbSources))
    do i=1,2*nbSources
      call NullifyTimeHistory(f(i))
    end do
    nbRanges = size(obs%dBRanges,2)
  else if(SPLdBFlag)then
    nbSources = size(obs%dBRanges,1)
    allocate(f(nbSources))  
    do i=1,nbSources
      call NullifyTimeHistory(f(i))
    end do
    nbRanges = size(obs%dBRanges,2)
  else
    nbSources = size(obs%dBARanges,1)
    allocate(f(nbSources))
    do i=1,nbSources
      call NullifyTimeHistory(f(i))
    end do
    nbRanges = size(obs%dBARanges,2)
  end if
  do i=1,nbRanges
    if(SPLdBFlag.and.SPLdBAFlag)then
      do j=1,nbSources
        call CopyTimeHistory(obs%dBRanges(j,i), f(j))
        call CopyTimeHistory(obs%dBARanges(j,i), f(j+nbSources))
      end do
      else if(SPLdBFlag)then
        do j=1,nbSources
          call CopyTimeHistory(obs%dBRanges(j,i), f(j))
        end do
      else
        do j=1,nbSources
          call CopyTimeHistory(obs%dBARanges(j,i), f(j))
        end do
    end if
    call WriteSingleObsTimeHistoryData(f,trim(path)//'_'//rangeTitles(i), &
                                       'SPL (dB)', 'SPL ')
  end do  
end subroutine WriteSingleObsFreqRangeData

subroutine WriteSingleObsTimeHistoryData(tH, path, ytag, xtag)
  implicit none
  type(timeHistory),dimension(:),intent(in)::tH
  character(len=*),intent(in)::path, ytag, xtag
  character(len=4096),dimension(:),allocatable::fTitles
  real, dimension(:,:), allocatable::f
  integer::i,nbSources
  nbSources = size(tH)
  allocate(fTitles(nbSources), f(nbSources, TH(1)%nt))  
  do i=1,nbSources
    f(i,:) = TH(i)%f
    fTitles(i) = trim(GetTHTitle(tH(i)))//trim(ytag)
  end do

  call WriteOutTecplotData(trim(path)//'.tec',TH(1)%t, trim(xtag)//' Time (s)',f,fTitles)  
  deallocate(fTitles,f)                  
end subroutine WriteSingleObsTimeHistoryData

end module observerObject
