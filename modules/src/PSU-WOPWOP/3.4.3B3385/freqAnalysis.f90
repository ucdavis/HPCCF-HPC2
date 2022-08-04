! PSU-WOPWOP
! $Id: freqAnalysis.f90 3270 2013-04-19 21:23:29Z brentner $
! $LastChangedDate: 2013-04-19 17:23:29 -0400 (Fri, 19 Apr 2013) $
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


module FrequencyAnalysisModule
  use constantsModule
  use interpolateModule
  use debugModule
  implicit none

  real (kind=8), dimension(:), pointer :: mg_fftTempReal
  complex (kind=8), dimension(:), pointer :: mg_fftTempComplex
  integer (kind=8), save :: mg_fftPlan, mg_ifftPlan

contains

subroutine InitializeFreqAnalysisData()
  nullify(mg_fftTempReal)
  nullify(mg_fftTempComplex)
  mg_fftPlan=0
  mg_ifftPlan=0
end subroutine InitializeFreqAnalysisData


subroutine DestroyFreqAnalysisData ()
  call dfftw_destroy_plan (mg_fftPlan)
  call dfftw_destroy_plan (mg_ifftPlan)
  if (associated(mg_fftTempReal)) then
    deallocate (mg_fftTempReal)
    nullify (mg_fftTempReal)
  end if
  if (associated(mg_fftTempComplex)) then
    deallocate(mg_fftTempComplex)
    nullify (mg_fftTempComplex)
  end if
end subroutine DestroyFreqAnalysisData

subroutine ConvertToSpectrum(array1, array2)
  real, dimension(:,:), pointer::array1, array2
  if (associated(array1)) call ConvertToSpectrumHelper(array1)
  if (associated(array2)) call ConvertToSpectrumHelper(array2)
end subroutine ConvertToSpectrum

subroutine ConvertToSpectrumHelper(array)
  real, dimension(:,:)::array
  real, dimension(:,:), allocatable::temp
  integer::i, iMax, nt, nfreq
  iMax  = size(array,1)
  nt    = size(array,2) 
  nfreq = int(nt/2.0)
  allocate(temp(iMax,nt))
  temp = 0.0
  do i=1, iMax
    call ConvertTimeToFrequency (array(i,:), temp(i,1:nFreq), &
                                               temp(i,nFreq+1:nt))
  end do
  array = temp
  deallocate(temp)
end subroutine ConvertToSpectrumHelper

subroutine AddSignal(signal1L, signal1C, time1L, time1C, signal2L, signal2C, time2L, time2C)
  real(kind=8), dimension(:):: signal1L, signal1C, time1L, time1C, signal2L, signal2C, time2L, time2C
  real(kind=8):: dtKind8
  real, dimension(2):: spectrum1, spectrum2 
  real:: timeOffset, dt
  integer::j
  if (integrationType.eq.TIME_DOMAIN) then
    dtKind8 = time2C(1)-time2L(1)
    call precisetinterpolate(time2L, time2C, signal2L, signal2C, time1C, signal1C, dtKind8)
  else
    do j=1, size(signal1L)  
      timeOffset = abs(time2C(j)-time1C(j))
      dt = time1C(j)-time1L(j)
      spectrum1 = (/signal1L(j), signal1C(j)/)
      spectrum2 = (/signal2L(j), signal2C(j)/)
      call AddSpectrumWithPhase(spectrum1, spectrum2, timeOffset, dt)
                                !time1(j, size(time1(j,:)) )-time1L(j))
    end do
  end if
end subroutine AddSignal


subroutine AddSpectrumWithPhase(spectrum1, spectrum2, timeOffset, dt)
  real, dimension(:)::spectrum1, spectrum2
  real::timeOffset, dt
  integer::nfreq, i
  real::f
  nfreq  = int(size(spectrum1)/2.0)
  do i=1, nfreq
    f = real(i)/dt
    spectrum2(i+nFreq)=spectrum2(i+nFreq)+2.0*pi*f*timeOffset
    call AddAmplitudeAndPhase(spectrum1(i), spectrum2(i), spectrum1(i+nfreq), spectrum2(i+nfreq))
  end do
end subroutine AddSpectrumWithPhase


subroutine AddAmplitudeAndPhase(amp1, amp2, angle1, angle2)
  real::amp1, amp2, angle1, angle2
  real::x1, x2, y1, y2
  x1=amp1*cos(angle1)
  x2=amp2*cos(angle2)
  y1=amp1*sin(angle1)
  y2=amp2*sin(angle2)
  x1=x1+x2
  y1=y1+y2
  amp1  =sqrt(x1**2+y1**2)
  angle1=atan(y1/x1)  
end subroutine AddAmplitudeAndPhase

subroutine SumSignalTogether(inputTime, inputData, outputTime, outputData, iBlankArray)
  real, dimension(:,:)::inputTime, inputData
  real(kind=8), dimension(:)::outputTime, outputData
  integer, dimension(:) ::iBlankArray
  
  integer::i, nFreq
  real:: dt 
  real, dimension(:), allocatable::outputFrequency

  nfreq = int(size(inputData,2)/2.0)
  allocate(outputFrequency(size(outputData)))
  dt = outputTime(size(outputTime))-outputTime(1)
  do i=1, size(outputData)
    outputFrequency(i)=i/dt
  end do       

  do i=1, size(inputData,1)
    call AddSignalByFrequency(inputTime(i,:), inputData(i,1:nfreq), &
                              outputFrequency, outputData)
  end do
  deallocate(outputFrequency)     

end subroutine SumSignalTogether

!subroutine SumSignalTogetherVector(inputTime, inputData, outputTime, outputData, &
!                                   iBlankArray)
!  real(kind=8),    dimension(:),intent(in)  :: inputTime
!  type(vector), dimension(:),intent(in) :: inputData
!  real(kind=8),    dimension(:)   ::outputTime
!  type(vector), dimension(:), intent(inout) :: outputData
!  integer, dimension(:),intent(inout)    ::iBlankArray
!  integer::i
!  real(kind=8),dimension(:),allocatable::preciseTime, preciseDatax,preciseDatay,preciseDataz
!  real(kind=8),dimension(:),allocatable::preciseOutputTime, preciseOutputDatax,preciseOutputDatay,preciseOutputDataz
!  allocate (preciseTime(size(inputTime,1),size(inputTime,2)), &
!            preciseDatax(size(inputData,1),size(inputData,2)), &
!            preciseDatay(size(inputData,1),size(inputData,2)), &
!            preciseDataz(size(inputData,1),size(inputData,2)), &
!            preciseOutputTime(size(outputTime)), &
!            preciseOutputDatax(size(outputData)), &
!            preciseOutputDatay(size(outputData)), &
!            preciseOutputDataz(size(outputData)))
!
!  preciseTime = inputTime
!  preciseDatax = inputData%A(1)
!  preciseDatay = inputData%A(2)
!  preciseDataz = inputData%A(3)
!  preciseOutputTime = outputTime
!  preciseOutputDatax(:) = outputData(:)%A(1)
!  preciseOutputDatay(:) = outputData(:)%A(2)
!  preciseOutputDataz(:) = outputData(:)%A(3)
!  call precisetinterpolate (preciseTime, preciseDatax, &
!                            preciseOutputTime, preciseOutputDatax, &
!                            iBlankArray)
!  call precisetinterpolate (preciseTime, preciseDatay, &
!                            preciseOutputTime, preciseOutputDatay, &
!                            iBlankArray)
!  call precisetinterpolate (preciseTime, preciseDataz, &
!                            preciseOutputTime, preciseOutputDataz, &
!                            iBlankArray)
!  outputTime = preciseOutputTime
!  outputData(:)%A(1) = preciseOutputDatax
!  outputData(:)%A(2) = preciseOutputDatay
!  outputData(:)%A(3) = preciseOutputDataz
!  deallocate (preciseTime, preciseDatax,preciseDatay,preciseDataz, preciseOutputTime)
!  deallocate (preciseOutputDatax,preciseOutputDatay,preciseOutputDataz)
! end subroutine SumSignalTogetherVector


subroutine AddSignalByFrequency(sourceTime, sourceAmplitude, obsFrequency, obsSpectrum)
  real, dimension(:), intent(inout)::sourceTime, sourceAmplitude, obsFrequency
  real(kind=8), dimension(:), intent(inout):: obsSpectrum

  integer::i, nbSourceFreq, nbObsFreq
  real:: dt
  real, dimension(:), allocatable::sourceFrequency, tempAmplitude

  nbObsFreq    = size(obsSpectrum)
  nbSourceFreq = size(sourceAmplitude)
  ! Need to allocate the temporary arrays
  allocate(sourceFrequency(nbSourceFreq))
  allocate(tempAmplitude  (nbObsFreq))
  tempAmplitude   = 0.0
  sourceFrequency = 0.0
  dt = sourceTime(2)-sourceTime(1)
  do i=1, nbSourceFreq
    sourceFrequency(i)=i/(sourceTime(size(sourceTime))-sourceTime(1))
  end do
  call tinterpolate (sourceFrequency, sourceAmplitude, obsFrequency, tempAmplitude)
  ! Add this elements contribution to the total history then proceed to the next one
  obsSpectrum = obsSpectrum + tempAmplitude**2
  ! Deallocate the allocatables
  deallocate(sourceFrequency)
  deallocate(tempAmplitude)
end subroutine AddSignalByFrequency

subroutine ResetFFT (nbTime)    
  integer :: nbTime
  call dfftw_destroy_plan (mg_fftPlan)
  call dfftw_destroy_plan (mg_ifftPlan)
  if (associated (mg_fftTempReal)) deallocate (mg_fftTempReal)
  if (associated (mg_fftTempComplex)) deallocate (mg_fftTempComplex)
  allocate (mg_fftTempReal(nbTime))
  allocate (mg_fftTempComplex(0:nbTime/2)) ! See fftw docs for size info
  call dfftw_plan_dft_r2c_1d (mg_fftPlan, nbTime, mg_fftTempReal, &
                              mg_fftTempComplex, FFTW_MEASURE)
  call dfftw_plan_dft_c2r_1d (mg_ifftPlan, nbTime, mg_fftTempComplex, &
                              mg_fftTempReal, FFTW_MEASURE+FFTW_PRESERVE_INPUT)
end subroutine ResetFFT

subroutine ConvertTimeToFrequency (pressure, amplitude, phase, complexDFTResult)
  real, dimension(:), intent(in)::pressure
  real, dimension(:), intent(out), optional :: amplitude, phase
  complex, dimension(:), intent(out), optional :: complexDFTResult
  integer::nbTime
  real,parameter::roottwo=1.41421356237
  
  ! Allocate p if we haven't already: this code optimizes for lots of calls
  ! using the same size p, a common case.
  nbTime=size(pressure)
  if (nbTime < 2) then
    call Error ("Not enough pressure time history to compute Fourier transform.")
  end if
  
  if (associated(mg_fftTempReal)) then
    if (size(mg_fftTempReal) /= nbTime) then
      call ResetFFT (nbTime)
    end if
  else
    call ResetFFT (nbTime)
  end if
  ! Copy over the data to our temp storage:
  mg_fftTempReal = pressure(1:nbTime)
  
  ! Calculate spectra using FFTW
  call dfftw_execute (mg_fftPlan)

  ! Note that ubound and lbound are used to account for the arrays that may
  ! start at zero, rather than one. This makes sense in this context, since the
  ! zeroth harmonic is the DC offset component.
  
  ! Convert to amplitude and phase
  if (present(amplitude)) then
    amplitude = 0
    ! Calculate the normalized RMS amplitude:
    amplitude(1:) = abs(mg_fftTempComplex(0:)) / (real(nbTime/2+1)*roottwo)
  end if
  
  if (present(phase)) then
    phase = 0
    ! Now the phase: negate it to conform to the engineering sign convention
    phase(1:) = -atan2(aimag(mg_fftTempComplex(0:)), &
              real(mg_fftTempComplex(0:)))
  end if
  
  if (present(complexDFTResult)) then
    complexDFTResult = 0
    complexDFTResult(1:) = mg_fftTempComplex(0:)
  end if
    
end subroutine ConvertTimeToFrequency


subroutine ConvertFrequencyToTime (complexDFTResult, pressure)
  complex, dimension(:), intent(in)::complexDFTResult
  real, dimension(:), intent(out)::pressure
  integer::npts
  real,parameter::rootwo=1.41421356237
 
  ! Copy over the data to our temp storage 
  mg_fftTempComplex = complexDFTResult
  mg_fftTempComplex(0) = 0.0
  ! Calculate reals using FFTW (note that we are using the inverse plan):
  call dfftw_execute (mg_ifftPlan)

  ! Copy over to pressure, making certain that our array sizes are OK:
  npts = min(size(mg_fftTempReal), size(pressure))
  pressure = 0
  pressure(1:npts) = mg_fftTempReal(1:npts)
  
end subroutine ConvertFrequencyToTime

subroutine CreateBlackmanWindow(window)
  implicit none
  real,dimension(:)::window
  integer::i, nt
  nt=size(window)
  do i=1,nt
    window(i)=0.42-.5*cos((2*pi)*real(i)/real(nt)) &
              +0.08*cos((4*pi)*real(i)/real(nt))
  end do
end subroutine CreateBlackmanWindow

subroutine CreateFlatTopWindow(window)
  implicit none
  real,dimension(:)::window
  integer::i, nt
  nt=size(window)
  do i=1,nt
    window(i)=1.-1.93*cos((2*pi)*real(i)/real(nt)) &
                  +1.29*cos((4*pi)*real(i)/real(nt)) &
                  -0.388*cos((6*pi)*real(i)/real(nt))&
                  +0.0322*cos((8*pi)*real(i)/real(nt))
  end do
end subroutine CreateFlatTopWindow

subroutine CreateHammingWindow(window)
  implicit none
  real,dimension(:)::window
  integer::i, nt
  nt = size(window)
  do i=1,nt
    window(i)=0.54-0.46*cos((2*pi)*real(i)/real(nt))
  end do
end subroutine CreateHammingWindow

subroutine CreateHanningWindow(window)
  implicit none
  real,dimension(:)::window
  integer::i, nt
  nt = size(window)
  do i=1,nt
    window(i)=0.5-0.5*cos((2.0*pi)*real(i)/real(nt))
  end do
end subroutine CreateHanningWindow


!**
!FUNCTION Aweight(f) result(weightA)
!  This function calculates the A weighting value used in the calculation of OASPLdBA.
!  This method is based off of ACS516 notes in chapter VI.  Note that A weighting 
!  includes the C weighting factor.
!
!  Large frequencies cause this routine process (Aweight and Cweight) to crash.  To overcome
!  this I made the computations use doubl precision.  KSB 2/6/2013
!
!**
function Aweight(fs) result(weightAs)
  implicit none
  real, intent(in):: fs  
  real:: weightAs
  real (kind=8)::weightA,k3,f2,f3,f
  f = fs  ! work in double precision
  k3=1.562339
  f2=107.65265
  f3=737.86223
  weightA=k3*f**4/((f**2+f2**2)*(f**2+f3**2))*Cweight(f)
  weightAs=weightA  ! send back a single precision answer
end function Aweight

!**
!FUNCTION Cweight(f) result(weightC)
!  This function calculates the C weighting value used in the calculation of OASPLdBA.
!  weightC is used in the calculation of the weightA value.
!**
function Cweight(f) result(weightC) 
  implicit none
  real (kind=8), intent(in):: f
  real (kind=8)::weightC,k1,G1,f4
  k1=2.242881E16
  G1=20.598997
  f4=12194.22
  weightC=k1*f**4/(((f**2+G1**2)**2)*((f**2+f4**2)**2))
end function Cweight


!**
!FUNCTION ConvertMSPTodB(meanSquareP) result(dB)
!  This function converts a mean-square pressure value to 
!  decibels.  If the mean square pressure is zero, it sets
!  the SPL equal to zero instead of taking the log of zero.    
!ARGUMENTS
! - meanSquareP: the mean square pressure value
!**
function ConvertMSPTodB(MSP) result(dB)
  implicit none
  real,intent(in)::MSP
  real::dB
  dB=10.*log10((MSP+epsilon**3)/(P_ref)**2)
end function ConvertMSPTodB

function ConvertMSPTodBArray(MSP) result(dB)
  implicit none
  real, dimension(:), intent(in)::MSP
  real, dimension(size(MSP))::dB
  dB=10.*log10((MSP+epsilon**3)/(P_ref)**2)
end function ConvertMSPTodBArray

function ConvertMSPTodBA(MSP,freq) result(dBA)
  implicit none
  real,intent(in)::MSP,freq
  real::dBA
  dBA=10.0*log10((MSP*Aweight(freq)+epsilon**3)/P_ref**2.0)
end function ConvertMSPTodBA

function ConvertdBtoMSP(dB) result(MSP)
  implicit none
  real, intent(in):: dB
  real:: MSP
  MSP = (P_ref**2.)*10.**(dB/10.)
end function ConvertdBtoMSP

function sinc(x) result(y)
  real, intent(in)::x
  real::y
  if(x==0.0) then
    y=1.0
  else
    y=sin(x)/x
  end if
end function sinc


!**
!The theory for this subroutine is from FARA36.4.7.2
!This subroutine was written by Leonard Lopes
!**
function CalcPNL(SPL) result(PNL)
  implicit none
  real,dimension(:),intent(in)::SPL
  integer::i, maxi
  real::maxnoy, noisiness, PNL  
  real, dimension(24)::SPLa, SPLb, SPLc, SPLd, SPLe, Mb, Mc, Md, Me
  real, dimension(24)::noy  
  !**Comment the following lines out. Dr. Brentner said it's okay to use
  !  less than 24 1/3 octave filtered bands for PNL calculations
  !  because higher harmonics of BPF output by WOPWOP are not high fideltiy
  !  anyway, and are also less than 10dB attenuated. -JPE 2/5/2008
  !if(size(SPL).ne.24)then
  !  call Error('In CalcPNL, the size of the SPL array does not equal 12')
  !end if
                   
  data SPLa / 91.00, 85.9, 87.3, 79.9, 79.8, 76.0, 74.0, 74.9, 94.6,  &
              1000.00, 1000.00, 1000.00, 1000.00, 1000.00, 1000.00, 1000.00, 1000.00, &
              1000.00, 1000.00, 1000.00, 1000.00, 1000.00, 44.3, 50.7 /
                   
  data SPLb / 64.00, 60.0, 56.0, 53.0, 51.0, 48.0, 46.0, 44.0, 42.0,  &
              40.0, 40.00, 40.00, 40.00, 40.00, 38.00, 34.00, 32.00,  &
              30.00, 29.00, 29.00, 30.00, 31.00, 37.0, 41.0 /
                   
  data SPLc / 52.00, 51.0, 49.0, 47.0, 46.0, 45.0, 43.0, 42.0, 41.0,  &
              40.00, 40.00, 40.00, 40.00, 40.00, 38.00, 34.00, 32.00, &
              30.00, 29.00, 29.00, 30.00, 31.00, 34.0, 37.0 /
                   
  data SPLd / 49.00, 44.0, 39.0, 34.0, 30.0, 27.0, 24.0, 21.0, 18.0,  &
              16.00, 16.00, 16.00, 16.00, 16.0, 15.00, 12.00, 9.00, 5.00, &
              4.00, 5.00, 6.00, 10.0, 17.0, 21.0 /
                   
  data SPLe / 55.0, 51.0, 46.0, 42.0, 39.0, 36.0, 33.0, 30.0, 27.0,   &
              25.0, 25.0, 25.0, 25.0, 25.0, 23.0, 21.0, 18.0, &
              15.0, 14.0, 14.0, 15.0, 17.0, 23.0, 29.0 /
 
  data Mb / 0.043478, 0.040570, 0.036831, 0.036831, 0.035336, 0.033333, &
            0.033333, 0.032051, 0.030675, 0.030103, 0.030103, 0.030103, &
            0.030103, 0.030103, 0.030103, 0.029960, 0.029960, 0.029960, &
            0.029960, 0.029960, 0.029960, 0.029960, 0.042285, 0.042285 /
            
  data Mc / 0.030103, 0.030103, 0.030103, 0.030103, 0.030103, 0.030103, &
            0.030103, 0.030103, 0.030103, 0.0,  0.0, &
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, 0.0, 0.029960, 0.029960 /
            
  data Md / 0.079520, 0.068160, 0.068160, 0.059640, 0.053013, 0.053013, &
            0.053013, 0.053013, 0.053013, 0.053013, 0.053013, 0.053013, &
            0.053013, 0.053013, 0.059640, 0.053013, 0.053013, 0.047712, &
            0.047712, 0.053013, 0.053013, 0.068160, 0.079520, 0.059640 /
            
  data Me / 0.058098, 0.058098, 0.052288, 0.047534, 0.043573, 0.043573, &
            0.040221, 0.037349, 0.034859, 0.034859, 0.034859, 0.034859, &
            0.034859, 0.034859, 0.034859, 0.040221, 0.037349, 0.034859, &
            0.034859, 0.034859, 0.034859, 0.037349, 0.037349, 0.043573 /     
            
  do i=1, size(SPL)
    if(SPL(i).ge.SPLa(i)) then
      noy(i)=10.**(Mc(i)*(SPL(i)-SPLc(i)))
    else if(SPL(i).ge.SPLb(i).and.SPL(i).lt.SPLa(i)) then
      noy(i)=10.**(Mb(i)*(SPL(i)-SPLb(i)))
    else if(SPL(i).ge.SPLe(i).and.SPL(i).lt.SPLb(i)) then
      noy(i)=0.3*10.**(Me(i)*(SPL(i)-SPLe(i)))
    else if(SPL(i).ge.SPLd(i).and.SPL(i).lt.SPLe(i)) then
      noy(i)=0.1*10.**(Md(i)*(SPL(i)-SPLd(i)))
    else
      noy(i)=0.0
    end if
  end do
  maxi=0
  maxnoy=-1000.00
  do i=1, size(SPL)
    if(maxnoy.lt.noy(i)) then
      maxi=i
      maxnoy=noy(i)
    end if
  end do  
  noisiness = 0.85*maxnoy + 0.15*sum(noy)    
  PNL = 40.0 + (10.0/log10(2.0))*log10(noisiness+epsilon)

end function CalcPNL

  
end module FrequencyAnalysisModule

