module FrequencyDomainObject

  use FrequencyAnalysisModule
  use TimeHistoryObject
 
  implicit none

  type freqDomain
    !title: the title of this frequency domain type
    !freq: the frequency array (domain)
    !f: the function array (range)
    !nf: the number of frequency bins, which is also the size
    !    of the two arrays
    !minFreq,maxFreq: the minimum and maximum frequency in the spectrum 
    character(len=4096)::title
    real::minFreq,maxFreq,df
    integer::nf
    real,dimension(:),pointer::freq=>null(),f=>null()
  end type freqDomain   
  
  private
  public::freqDomain,FilterfreqDomain,GetNf,          &
          DestroyfreqDomain,Destroy3DFreqDomain,      &
          SetFreqArray, addThirdOctMSP,               &
          SetMaxFreq,SetMinFreq,SetNf,                &
          SetfreqDomainTitle,GetfreqDomainTitle,      &
          GetMinFreq,GetMaxFreq,GetDf,expandMSPArray, &
          SetfreqDomainFArray,IntegrateOverFreqDomain,&
          AweightIntegrateOverFreqDomain,             &
          CalculateOASPLdBA, GetThirdOctFreqArray,    &
          CreateMSPPerBand,GetFreq,CopyFreqDomain,    &
          ReceiveFreqDomainFromSlave, SendFreqDomainToMaster
          
contains

subroutine ReceiveFreqDomainFromSlave(FD, master, sender)
  type(FreqDomain), intent(inout)::FD
  integer, intent(in)::master, sender
  logical::flag
  integer::n
  call ReceiveString(FD%title, master, sender)
  call ReceiveInteger(FD%nf, master, sender)
  call ReceiveReal(FD%minFreq, master, sender)
  call ReceiveReal(FD%maxFreq, master, sender)
  nullify(FD%f)
  call ReceiveLogical(flag, master, sender)
  if (flag) then
    call ReceiveInteger(n, master, sender)
    allocate(FD%f(n))
    call ReceiveReals(FD%f, master, sender)
  end if
  nullify(FD%freq)
  call ReceiveLogical(flag, master, sender)
  if (flag) then
    call ReceiveInteger(n, master, sender)
    allocate(FD%freq(n))
    call ReceiveReals(FD%freq, master, sender)
  end if
end subroutine ReceiveFreqDomainFromSlave

subroutine SendFreqDomainToMaster(FD, master, tag)
  type(freqDomain), intent(in)::FD
  integer, intent(in)::master, tag
  logical::flag
  call SendString(FD%title, master, tag)
  call SendInteger(FD%nf, master, tag)
  call SendReal(FD%minFreq, master, tag)
  call SendReal(FD%maxFreq, master, tag)
  flag = associated(FD%f)
  call SendLogical(flag, master, tag)
  if (flag) then
    call SendInteger(FD%nf, master, tag)
    call SendReals(FD%f, master, tag)
  end if
  flag = associated(FD%freq)
  call SendLogical(flag, master, tag)
  if (flag) then
    call SendInteger(FD%nf, master, tag)
    call SendReals(FD%freq, master, tag)
  end if
end subroutine SendFreqDomainToMaster

!**
!SUBROUTINE Destroy3DFreqDomain
!  This routine destroys a 3-D input frequency
!  domain array, fD.  It calls DestroyFreqDomain
!  for each element in the array.
!**
subroutine Destroy3DFreqDomain(fD)
  implicit none
  type(freqDomain),dimension(:,:,:),intent(inout)::fD
  integer::i,j,k
  do i=1,size(fD,1)
    do j=1,size(fD,2)
      do k=1,size(fD,3)      
        call DestroyFreqDomain(fD(i,j,k))
      end do
    end do
  end do
end subroutine Destroy3DFreqDomain

!**
!FUNCTION IntegrateOverFreqDomain
!  This function adds all the values in fD%f to produce
!  the integral result.  It assumes a unit frequency bin
!  spacing of 1.
!**
function IntegrateOverFreqDomain(fD) result(integral)
  implicit none
  type(freqDomain),intent(in)::fD
  real::integral
  if(fD%nf.lt.1)then
    integral=0.
  else
    integral=sum(fD%f)
  end if
end function IntegrateOverFreqDomain

!**
!FUNCTION CalculateOASPLdBA
!  This function calculates the sound pressure level, in 
!  A-weighted decibels, for the entire fD named MSP.
!**
function CalculateOASPLdBA(MSP) result(dBA)
  implicit none
  type(freqDomain),intent(in)::MSP
  real::dBA,AmeanSquareP 
  AmeanSquareP = AweightIntegrateOverFreqDomain(MSP)
  dBA = ConvertMSPTodB(AmeanSquareP)
end function CalculateOASPLdBA

!**
!FUNCTION AweightIntegrateOverFreqDomain
!  This function integrates over the input frequency domain
!  fD using the A-weighted method.
!**
function AweightIntegrateOverFreqDomain(fD) result(sum)
  implicit none
  type(freqDomain),intent(in)::fD
  real::sum
  integer::i
  sum = 0.0
  do i=1,FD%nf
    sum=sum+FD%f(i)*Aweight(FD%freq(i))
  end do
end function AweightIntegrateOverFreqDomain

!**
!FUNCTION CalcOctaveBandCenterFreq
!  This function calculates the center frequency of any ith band.
!  i=0 is the fundamental frequency, 1000Hz.  +i are frequencies above
!  1000Hz and -i are frequencies below 1000Hz.
!ARGUMENTS
! - i: the ith band (see above and manual)
! - octaveNumber: (1/octaveNumber) octave filtering
! - r: value depends on octave approximation flag (see user manual)
!**
function CalcOctaveBandCenterFreq(i,r,octaveNumber) result(centerFreq)
  implicit none
  real,intent(in)::r,octaveNumber
  integer,intent(in)::i
  real::centerFreq
  centerFreq = 1000.*r**(real(i)/octaveNumber)
end function CalcOctaveBandCenterFreq

!**
!SUBROUTINE CalcOctaveBandCenterFreqArray
!  This subroutine calculates and returns an octave band center frequency array
!  according to the maximum and minimum frequencies of the input fD.
!  If the max and min frequencies of fD are in the middle of the octave 
!  band, the band is still included.
!ARGUMENTS
! - fD: the fD which the octave band array is being filled for
! - centerFreqArray: a real, 1-D array of center frequencies
! - octaveNumber: (1/octaveNumber) octave filtering
!**                                
subroutine CalcOctaveBandCenterFreqArray(fD,centerFreqArray,octaveNumber)
  implicit none
  type(freqDomain),intent(in)::fD
  real,intent(in)::octaveNumber
  real,dimension(:),pointer::centerFreqArray
  real::center,lowLimit,highLimit,r,factor
  integer::i,nbBands,iBand 
  !The function GetR() returns the value of r depending on whether octaveApproxFlag
  !  is turned on or off.  See the post-processing documentation for more information.
  r      = GetR(octaveApproxFlag) 
  factor = r**(1/(2*octaveNumber))
  !This is where we determine what band number the frequency range
  !starts in.  i is the band number.  When i=0 this is the fundamental
  !frequency, 1000Hz.  +i are bands above 1000Hz, and -i are bands below 1000Hz.
  i = 0  !Start at the fundamental frequency
  center    = 1000.
  lowLimit  = center/factor
  highLimit = center*factor 
  if(fD%minFreq.lt.lowLimit)then           
    !if the minimum frequency is less than 1000Hz, start moving
    !backwards until the starting band is found:
    do while(fD%minFreq.lt.lowLimit/(factor**2)) ! lowLimit/(factor**2) is the next lower band low limit     
      i = i-1
      center   = CalcOctaveBandCenterFreq(i,r,octaveNumber)      
      lowLimit = center/factor
      if(lowLimit.lt.GetDf(fD))then
        exit
      end if
    end do       
    !else, move forward until the first band is found
    else                               
      do while(fD%minFreq.gt.highLimit)
        i=i+1
        center    = CalcOctaveBandCenterFreq(i,r,octaveNumber)
        highLimit = center*factor
      end do
  end if  
  iBand = i
  !The value of i is now the first band in the input frequency range.
  !Now we have to count the bands.
  nbBands = 1
  center = CalcOctaveBandCenterFreq(iBand,r,octaveNumber)
  do while(fD%maxFreq.gt.(center*factor))
    i = i+1
    center  = CalcOctaveBandCenterFreq(i,r,octaveNumber)
    nbBands = nbBands+1
  end do
  !Now we can allocate the arrays according to how many bands 
  !there are
  if(.not.associated(centerFreqArray))then
    allocate(centerFreqArray(nbBands))
  end if
  !Now we have to fill the array:
  do i=1,nbBands
    centerFreqArray(i) = CalcOctaveBandCenterFreq(iBand,r,octaveNumber)
    iBand = iBand+1
  end do  
end subroutine CalcOctaveBandCenterFreqArray

!**
!FUNCTION CreateMSPPerBand
!  This function calculates the mean-square pressure for each
!  octave band in MSP.
!ARGUMENTS
! - MSP: the narrow band mean-square pressure spectrum.
! - octaveNumber: the filtering number used to create MSPPerBand
!**
subroutine CreateMSPPerBand(MSP, octaveNumber, MSPPerBand)
  implicit none
  type(freqDomain),intent(in)::MSP
  real,intent(in)::octaveNumber
  type(freqDomain)::MSPPerBand
  real,dimension(:), pointer::centerFreqArray,MSPArray
  type(freqDomain)::bandSpec
  integer::i,nf,nf1 !ksb debug startIndex
  real::factor,lowLimit,highLimit, minFreq, maxFreq
  nullify(centerFreqArray,MSPArray,bandSpec%f,bandSpec%freq) 
  !First, calculate the center frequency array
  call CalcOctaveBandCenterFreqArray(MSP, centerFreqArray, octaveNumber)
  nf = size(centerFreqArray)
  allocate(MSPArray(nf))
  factor = GetR(octaveApproxFlag)**(1/(2*octaveNumber))
  !Create a temporary spectrum (bandSpec) for each octave band and calculate
  !the SPL for each:  
  do i=1,size(centerFreqArray)
    lowLimit    = centerFreqArray(i)/factor
    highLimit   = centerFreqArray(i)*factor
    call FilterFreqDomain(MSP,lowLimit,highLimit, bandSpec)        
    MSPArray(i) = IntegrateOverFreqDomain(bandSpec)
  end do
  !startIndex = 1  !ksb debug - I don't think this logic is correct and I'm not sure why it is needed.!
  !do i=nf,1,-1
  !  if(MSPArray(i).eq.0.0)then !ksb debug: The problem is that MSPArray(i) can be 0.0 but i+1 and i-1 can have values
  !    startIndex = i + 1
  !    exit
  !  end if
  !end do
  !nf1 = nf-startIndex+1
  !if (startIndex .eq. size(centerFreqArray)+1) then
  !  minFreq = 0.
  !  maxFreq = 0.
  !else
  !  minFreq = centerFreqArray(startIndex)
  !  maxFreq = centerFreqArray(size(centerFreqArray))
  !end if
  !nf = size(centerFreqArray)-startIndex+1
  nf1 = nf
  if(.not.associated(MSPPerBand%f)) then
    allocate(MSPPerBand%f(nf))
  end if
  if(.not.associated(MSPPerBand%freq)) then
    allocate(MSPPerBand%freq(nf))
  end if
  call SetFreqDomainTitle(MSPPerBand, 'MSPPerOctaveBand')
  call SetFreqDomainFArray(MSPPerBand, MSPArray)  !ksb debug (startIndex:))
  call SetFreqArray(MSPPerBand, centerFreqArray)  !ksb debug (startIndex:))
  call SetNf(MSPPerBand, nf)
  call SetMinFreq(MSPPerBand, centerFreqArray(1)) !ksb debug minFreq)
  call SetMaxFreq(MSPPerBand, centerFreqArray(nf))!ksb debug maxFreq)
  deallocate(MSPArray,centerFreqArray)
end subroutine CreateMSPPerBand

!SUBROUTINE expandMSPArray
! Expands MSP array to accommodate broadband data at frequencies
! higher than is present in the discrete frequency data.
subroutine expandMSPArray(dfMSP, origMSP, bbMSP, maxDF, maxBB)
  implicit none
  type(freqDomain), intent(inout):: origMSP
  type(freqDomain):: dfMSP, bbMSP, BBtempMSP
  real:: maxDF, maxBB

  ! Filter out any BB content that lies in the octaveFiltMSP freq domain
  ! by setting the high pass frequency to the lower limit of the third octave above 
  ! the last third octave in octaveFiltMSP. (Add a little bit to the low-pass freq
  ! to make sure nothing is missed.)
  
  !ksb debug:  I would like to change this to be more general and not assume
  ! that minDF < minBB.  Currently it only expands the upper bound. 1/16/16
  call FilterFreqDomain(bbMSP, maxDF*(2.**(1./6.)), maxBB+5., BBtempMSP)

  ! create a new MSP array whose length is the sum of
  ! of the octaveFiltMSP and mspBB arrays
  call SetNf(dfMSP, origMSP%nf+BBtempMSP%nf)
  call SetMaxFreq(dfMSP, BBtempMSP%maxFreq)
  allocate(dfMSP%freq(dfMSP%nf), dfMSP%f(dfMSP%nf))
  dfMSP%freq(1:origMSP%nf) = origMSP%freq
  dfMSP%freq(origMSP%nf+1:origMSP%nf+BBtempMSP%nf) = BBtempMSP%freq
  dfMSP%f(:) = 0.0
  dfMSP%f(1:origMSP%nf) = origMSP%f

  call DestroyFreqDomain(BBtempMSP)
end subroutine expandMSPArray

!FUNCTION addThirdOctMSP
!  Add MSP by third-octave bands.
function addThirdOctMSP(msp, mspBB, useMSP) result(totalMSP)
  implicit none 
  type(freqDomain), target:: msp, mspBB, totalMSP
  type(freqDomain), pointer:: band=>null()
  integer:: i,j, jMin
  logical:: useMSP

  ! Use either the msp or mspBB array to set the frequency range
  ! as the foundation of the copy routine
  if (useMSP) then
    call CopyFreqDomain(msp,totalMSP)
    band => mspBB
  else
    call CopyFreqDomain(mspBB,totalMSP)
    band => msp
  end if
 
  jMin = 1
  do i=1,size(totalMSP%freq)
    do j=jMin,size(band%freq)
      if (totalMSP%freq(i)/(2.**(1./6.0)).lt.band%freq(j) .and. &
          totalMSP%freq(i)*(2.**(1./6.0)).gt.band%freq(j)) then
        totalMSP%f(i) = totalMSP%f(i) + band%f(j)
        jMin = j
        exit
      end if
    end do
  end do
end function addThirdOctMSP

!**
!FUNCTION CopyFreqDomain
!  This routine makes a complete copy of the input freq domain
!  fD1 into the time history fDcopy.  This routine is necessary because
!  the fD type contains pointers which need their own unique memory space.
!**
subroutine CopyFreqDomain(fD1, fDcopy)
  implicit none
  type(freqDomain),intent(in)::fD1
  type(freqDomain), intent(inout)::fDcopy
  
  if (associated(fDcopy%freq)) deallocate(fDcopy%freq)
  if (associated(fD1%freq)) then
    allocate(fDcopy%freq(size(fD1%freq)))  
    fDcopy%freq   =fD1%freq
  end if

  if (associated(fDcopy%f)) deallocate(fDcopy%f) 
  if (associated(fD1%f)) then
    allocate(fDcopy%f(size(fD1%f)))    
    fDcopy%f	    =fD1%f
  end if

  fDcopy%nf     =fD1%nf
  fDcopy%minFreq=fD1%minFreq
  fDcopy%maxFreq=fD1%maxFreq
  fDcopy%title  =fD1%title
end subroutine CopyFreqDomain

!**
!SUBROUTINE DestroyFreqDomain
!  This routine nullifies the pointers contained
!  in the fD type.
!**
subroutine DestroyFreqDomain(fD)
  implicit none
  type(freqDomain),intent(inout)::fD  
  if(associated(fD%freq))then
    deallocate(FD%freq)
    nullify(fD%freq)
  end if  
  if(associated(fD%f))then
    deallocate(FD%f)
    nullify(fD%f)
  end if  
end subroutine DestroyFreqDomain

!**
!FUNCTION FilterFreqDomain
!  Checks the input freq domain, fD, to see if any frequencies are below the 
!  high pass frequency or above the high pass frequency specified by
!  the user.  If they are, these frequency points are removed from the
!  freq domain.
!ARGUMENTS
!  spec - the spectrum type that we are filtering
!  hpf - the high pass frequency
!  lpf - the low pass frequency
!RESULT
!  fDfiltered - the new filtered spectrum type
!**
subroutine FilterFreqDomain(fD, hpf, lpf, fDfiltered)
  implicit none
  type(freqDomain),intent(in)::fD
  real,intent(in)::hpf,lpf
  type(freqDomain),intent(inout)::fDfiltered  
  integer::sum,i, minIndex
  call SetFreqDomainTitle(fDfiltered, GetFreqDomainTitle(fD))
  sum = 0

   !Here we are counting to see how many frequency points
   !are within the hpf and lpf boundaries:
   ! (Start from the beginning of the original spectrum and find
   ! the first point that is higher than the high pass frequency)
  do i=1, fD%nf  !ksb debug size(fD%freq)
    if((fD%freq(i).ge.(hpf)).and.(fD%freq(i).le.(lpf)))then
      sum = sum+1
      if (sum.eq.1) minIndex = i
    end if
  end do
  if(sum==0)then
    fDfiltered%nf      = 0
    fDfiltered%minFreq = hpf
    fDfiltered%maxFreq = lpf
    !This is simply to eliminate access violations:
    allocate(fDFiltered%f(1),fDfiltered%freq(1))
    fDFiltered%f(1)    = 0.
    fdFiltered%freq(1) = 0.
  else
    !Sum is the amount of frequency bins in the new spectrum  
    call SetNf(fDfiltered,sum)
    allocate(fDfiltered%f(sum),fDfiltered%freq(sum))
    !then set the values for the filtered spectrum:
    do i=minIndex,size(fD%freq)
      if(fD%freq(i).ge.(hpf-epsilon))then
        call SetFreqDomainFArray(fDfiltered,fD%f(i:(i+sum-1)))
        call SetFreqArray(fDfiltered,fD%freq(i:(i+sum-1)))
        call SetMinFreq(fDfiltered,fD%freq(i))
        call SetMaxFreq(fDfiltered,fD%freq(i+sum-1))
        exit
      end if
    end do
  end if  
end subroutine FilterFreqDomain

!**
!FUNCTION GetMaxFreq
!  This function returns the maximum frequency
!  contained in the input frequency domain range, fR%maxFreq.
!**
function GetMaxFreq(fR)result(max)
  implicit none
  type(freqDomain),intent(in)::fR
  real::max
  max=fR%maxFreq
end function GetMaxFreq

!**
!FUNCTION GetMaxFreq
!  This function returns the minimum frequency
!  contained in the input frequency domain range, fR%minFreq.
!**
function GetMinFreq(fR)result(min)
  implicit none
  type(freqDomain),intent(in)::fR
  real::min
  min=fR%minFreq
end function GetMinFreq

!**
!FUNCTION GetMaxFreq
!  This function returns the number of frequency bins
!  contained in the input frequency domain, fD%nf.
!**
function GetNf(fD)result(nf)
  implicit none
  type(freqDomain),intent(in)::fD
  integer::nf
  nf=fD%nf
end function GetNf  

!**
!FUNCTION GetR()
!  R is used in the calculation of the octave band center frequencies.
!  This function calculates R based on an approximation flag.
!  See the user manual for a detailed explanation of the 
!  calculation of R.
!**
function GetR(approxFlag) result(r)
  implicit none
  logical,intent(in)::approxFlag
  real::r
  !calculate the value of R based on octaveApproxFlag (see Manual)
  if(approxFlag)then  
    r = 2.                    
  else
    r = 10**(0.3)
  end if
end function GetR

!**
!FUNCTION GetMaxFreq
!  This function returns the frequency bin spacing
!  of the input frequency domain, fD.  It is mandatory
!  that each frequency bin is spaced equally or this routine
!  does not return anything useful.  For example, do not use
!  this routine for an octave-filtered frequency domain.
!**
function GetDf(fD) result(df)
  implicit none
  type(freqDomain),intent(in)::fD
  real::df
  df = (fD%maxFreq - fD%minFreq)/(real(fD%nf) - 1.)
end function GetDf

!**
!FUNCTION GetFreq
!  This function returns the ith value of the frequency
!  domain frequency array, fD%freq(i).
!**
function GetFreq(fD,i) result(freq)
  implicit none
  type(freqDomain),intent(in)::fD
  integer,intent(in)::i
  real::freq
  freq=fD%freq(i)
end function GetFreq

!**
!FUNCTION GetFreqDomainTitle
!  This function returns the title of the frequency
!  domain, fD%title.
!**
function GetFreqDomainTitle(fD) result(title)
  implicit none
  type(freqDomain),intent(in)::fD
  character(len=4096)::title
  title=fD%title
end function GetFreqDomainTitle

function GetThirdOctFreqArray() result(freqArray)
  implicit none
  real, dimension(BB_THIRD_OCTAVE_BANDS):: frCen, freqArray
  data frCen         /  50.  ,    63.   ,    80.   ,    100.  ,   125.   , &
                       160.  ,   200.   ,   250.   ,    315.  ,   400.   , &
                       500.  ,   630.   ,   800.   ,   1000.  ,  1250.   , &
                      1600.  ,  2000.   ,  2500.   ,   3150.  ,  4000.   , &
                      5000.  ,  6300.   ,  8000.   ,  10000.  , 12500.   , &
                     16000.  , 20000.   /!,  31500.  , 40000.   /

    freqArray = frCen
    
end function GetThirdOctFreqArray

!**
!SUBROUTINE SetFreqArray
!  This routine sets the input frequency domain
!  frequency array equal to the input array, freq.
!**
subroutine SetFreqArray(fD,freq)
  implicit none
  type(freqDomain),intent(inout)::fD
  real,dimension(:),intent(in)::freq
  fD%nf = size(freq)
  if(.not.associated(fD%freq)) then
    allocate(fD%freq(fD%nf))
  end if
  fD%freq=freq
end subroutine SetFreqArray

!**
!SUBROUTINE SetMaxFreq
!  This routine sets the maximum frequency of a freq
!  domain, fD%maxFreq.  Make sure this value is equal
!  to the last entry of the frequency array, fD%freq(nf).
!**
subroutine SetMaxFreq(fD,maxFreq)
  implicit none
  type(freqDomain),intent(inout)::fD
  real,intent(in)::maxFreq
  fD%maxFreq=maxFreq
end subroutine SetMaxFreq

!**
!SUBROUTINE SetMinFreq
!  This routine sets the minimum frequency of a freq
!  domain, fD%minFreq.  Make sure this value is equal
!  to the first entry of the frequency array, fD%freq(1).
!**
subroutine SetMinFreq(fD,minFreq)
  implicit none
  type(freqDomain),intent(inout)::fD
  real,intent(in)::minFreq
  fD%minFreq=minFreq
end subroutine SetMinFreq

!**
!SUBROUTINE SetNf
!  This routine sets the number of frequency bins
!  in a frequency domain, fD%nf.
!**
subroutine SetNf(fD,nf)
  implicit none 
  type(freqDomain),intent(inout)::fD
  integer,intent(in)::nf
  fD%nf=nf
end subroutine SetNf

!**
!SUBROUTINE SetFreqDomainFArray
!  This routine sets the input frequency domain
!  function array, fD%f, equal to the input array f.
!**
subroutine SetFreqDomainFArray(fD,f)
  implicit none
  type(freqDomain),intent(inout)::fD
  real,dimension(:),intent(in)::f
  if (associated(FD%f)) then
    deallocate(fD%f)
  end if    
  allocate(fD%f(size(f)))
  fD%f = f
end subroutine SetFreqDomainFArray

!**
!SUBROUTINE SetFreqDomainTitle
!  This routine sets the freq domain title, fD%title
!** 
subroutine SetFreqDomainTitle(fR,title)
  implicit none 
  type(freqDomain),intent(inout)::fR
  character(len=*),intent(in)::title
  fR%title=trim(title)
end subroutine SetFreqDomainTitle

end module FrequencyDomainObject 
