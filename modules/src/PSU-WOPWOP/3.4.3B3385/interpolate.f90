! $LastChangedDate: 2013-10-17 22:12:03 -0400 (Thu, 17 Oct 2013) $
! $Id: interpolate.f90 3298 2013-10-18 02:12:03Z brentner $
!********************************************************
!TYPE: integrate module
!  This module incorporates all the functions used for integration
! and extrapolation of data. It is mainly used by patches via
! integratePatch routine.
! - most of this material was elsewhere in the code but moved here
!   by Leonard Lopes summer '02
!********************************************************
module interpolateModule
use mathModule
use constantsModule, only:iBlankFlag
implicit none

integer::interpolatePosition

contains

!******************************************************************************
!SUBROUTINE differentiateArrayReal(array,arraytua,darray)
!  this subroutine calculate the differental of an array of reals and returns
!  the result darray
!ARGUEMENTS:
!  - array: the array that we need the derivative of
!  - arraytua: the time array that we need to take the derivative with respect to
!  - darray: the result we are calculating
!******************************************************************************
subroutine differentiateArrayReal(array,arraytau,darray)
  real, dimension(:), intent(inout)::array,arraytau,darray

  integer::n,i

  n=size(array)
  darray(1)=(array(2)-array(1))/(arraytau(2)-arraytau(1))
  do i=2, n-1
     darray(i)=(array(i+1)-array(i-1))/(arraytau(i+1)-arraytau(i-1))
  end do
  darray(n)=(array(n)-array(n-1))/(arraytau(n)-arraytau(n-1))

end subroutine DifferentiateArrayReal

!******************************************************************************
!SUBROUTINE differentiateArrayVector(array,arraytua,darray)
!  this subroutine calculate the differental of an array of vectors and returns
!  the result darray
!ARGUEMENTS:
!  - array: the array that we need the derivative of
!  - arraytua: the time array that we need to take the derivative with respect to
!  - darray: the result we are calculating
!******************************************************************************
subroutine differentiateArrayVector(array,arraytau,darray)
  type(vector), dimension(:), intent(inout)::array, darray
  real, dimension(:), intent(inout)::arraytau
  integer::i,n
  
  n=size(array)
  do i=1, n-1
    darray(i)%A(1)=(array(i+1)%A(1)-array(i)%A(1))/(arraytau(i+1)-arraytau(i)) 
    darray(i)%A(2)=(array(i+1)%A(2)-array(i)%A(2))/(arraytau(i+1)-arraytau(i))
    darray(i)%A(3)=(array(i+1)%A(3)-array(i)%A(3))/(arraytau(i+1)-arraytau(i))
  end do
  darray(n)%A(1)=(array(n)%A(1)-array(n-1)%A(1))/(arraytau(n)-arraytau(n-1))
  darray(n)%A(2)=(array(n)%A(2)-array(n-1)%A(2))/(arraytau(n)-arraytau(n-1))
  darray(n)%A(3)=(array(n)%A(3)-array(n-1)%A(3))/(arraytau(n)-arraytau(n-1))  

end subroutine DifferentiateArrayVector

!*******************************************************************************
!SUBROUTINE interpolateVector(tau, arrayTau,arrayValues,res)
!  this routine takes time, array of time history, array of values and interpolates
!  to find the result resulting vector
!ARGUEMENTS:
!  - tau: the time we need to interpolate at
!  - arrayTau: the complete time history, we need this to interpolate accurately
!  - res: the vector result we need
!********************************************************************************
function interpolateVector(time,arrayTime,arrayValues,index) result(res)
  real, intent(in)::time
  real,dimension(:), intent(in)::arrayTime
  type(vector), dimension(:), intent(in)::arrayValues
  integer, intent(inout)::index

  type(vector)::res
  integer::n, i

  n=size(arrayTime)
  i=index
  call hunt(arrayTime,n,time,i)! i is the index of the value in arrayTau just below tau.
  if(i<=0) then
    res = arrayValues(1)
  else if (i>n) then
    res = arrayValues(n)
  else if(i==n) then
    res=(arrayValues(i))
  else
    res=(arrayTime(i+1)-time)/(arrayTime(i+1)-arrayTime(i))*arrayValues(i)+&
       &(time-arrayTime(i))/(arrayTime(i+1)-arrayTime(i))*arrayValues(i+1)
    index=i  
  end if

end function interpolateVector

!*******************************************************************************
!SUBROUTINE interpolateReal(tau, arrayTau,arrayValues,res)
!  this routine takes time, array of time history, array of values and interpolates
!  to find the result resulting real value
!ARGUEMENTS:
!  - tau: the time we need to interpolate at
!  - arrayTau: the complete time history, we need this to interpolate accurately
!  - res: the real result we need
!********************************************************************************
function interpolateReal(time,arrayTime,arrayValues,index) result(res)
  real, intent(in)::time
  real, dimension(:), intent(inout)::arrayTime, arrayValues 
  integer, intent(inout)::index

  integer::n,i 
  real::res

  n=size(arrayTime)
  i=index
  call hunt(arrayTime,n,time,i)! i is the index of the value in arrayTau just below tau.
  if(i<=0) then
    res = arrayValues(1)
  else if (i>n) then
    res = arrayValues(n)
  else if(i==n) then
    res=(arrayValues(i))
  else
    res=(arrayTime(i+1)-time)/(arrayTime(i+1)-arrayTime(i))*arrayValues(i)+&
       &(time-arrayTime(i))/(arrayTime(i+1)-arrayTime(i))*arrayValues(i+1)
    index=i  
  end if

end function interpolateReal

!**********************************************************************************
!SUBROUTINE tinterpolate(tLocal,qdSLocal,tGlobal,qdSPatch)
!  Given the integral contribution qdSLocal (or a function f) at the observer
!time tLocal (or at point x), the subroutine interpolates the value at 
!observer time tGlobal (or the value f(x0) at x0) and sums it to qdS-patch
!ARGUMENTS:
!  - tLocal:array of real, independant variable x(i)
!  - qdSLocal:array of real, dependant variable f(x(i))
!  - tGlobal:array of real: points x0(i) at which the interpolation will be 
!performed
!  - qdSPatch:array of real, sum of the original value qdSPatch(i) and the 
!interpolated value f(x0(i))
!REMARKS:
!  - new subroutine for interpolation (more efficient)
!  - extrapolation is not allowed: the subroutine stops if extrapolation is needed
!**********************************************************************************
subroutine precisetinterpolate(tLocalL, tLocalR, qdSLocalL, qdSLocalR, &
                               tGlobal,qdSPatch, dtGlobal, iBlankArray)
  real(kind=8), intent(in), dimension(:)::tLocalL, tLocalR, qdSLocalL, qdSLocalR, tGlobal
  real(kind=8), intent(inout), dimension(:)::qdSPatch
  integer, dimension(:), optional::iBlankArray
  real(kind=8)::dtGlobal
  
  integer:: i, it, ntGlobal, nbNodes, minGlobal, maxGlobal, nMinGlobal, nMaxGlobal
  real(kind=8)::t, Atemp, tMinLocal, tMaxLocal
  
  nbNodes = size(tLocalL)
  ntGlobal = size(tGlobal)
  
  tMinLocal = minval(tLocalL)
  tMaxLocal = maxval(tLocalR)
  nMinGlobal = floor((tMinLocal-tGlobal(1))/dtGlobal)+1
  nMaxGlobal = ceiling((tMaxLocal-tGlobal(1))/dtGlobal)+1
  
 
  do i=1, nbNodes  
    nMinGlobal = floor((tLocalL(i)-tGlobal(1))/dtGlobal)+1
    nMaxGlobal = ceiling((tLocalR(i)-tGlobal(1))/dtGlobal)+1

    if (nMinGlobal <= ntGlobal .and. nMaxGlobal >= 1) then  
      ! We need the loop not to crash.  These bounds will help speed things up
      ! but we can still rely on the checks within the loop to skip unnecessary 
      ! interpolation times.
      if (nMinGlobal < 1) nMinGlobal = 1
      if (nMaxGlobal > ntGlobal) nMaxGlobal = ntGlobal
      
      Atemp = (qdSLocalR(i)-qdSLocalL(i))/(tLocalR(i)-tLocalL(i))    
      do it=nMinGlobal,nMaxGlobal  
        t = tGlobal(it)      
        if( tLocalL(i) <= t) then          
          if( tLocalR(i) > t) then
            qdSPatch(it) = qdSPatch(it) + qdSLocalL(i)+Atemp*(t-tLocalL(i))
          else
           ! If the Global time is ahead of the Local time (no way to interpolate)
             if (( t > tLocalR(i)) .and. iBlankFlag) then
               if(present(iBlankArray)) iBlankArray(it)=0
               exit
             end if  
          end if
        else ! If the Global time is behind the Local time (no way to interpolate)
           if(present(iBlankArray)) iBlankArray(it)=0
        end if
      end do   
    end if        
  end do
  return

end subroutine precisetinterpolate


function SourceTimeInRange(tMinLocal, tMaxLocal, tMinGlobal, ntGlobal, dtGlobal) result(inRange)
  implicit none
  real(kind=8):: tMinLocal, tMaxLocal, tMinGlobal, dtGlobal
  real(kind=8):: temp
  integer, intent(in):: ntGlobal
  integer::minGlobal, maxGlobal
  logical:: inRange
  
  inRange = .false.
  if (tMinLocal.ne.huge(temp) .and. tMaxLocal.ne.-huge(temp)) then
    ! This minGlobal can be (+1) when only dealing with the interpolation
    ! but because we're using it to check if the aperiodic data is in
    ! range of the observer time array we need it to be (-1) to ensure
    ! the pressure gradients still work the same as before.
    ! maxGlobal must be +8 to ensure that, if we expand the observer
    ! time range, we have a 2nd-order-accurate solution at the first new
    ! observer time step.
    minGlobal = floor((tMinLocal-tMinGlobal)/dtGlobal)-1
    if (newSegLHS.or.newSegRHS) then
      maxGlobal = ceiling((tMaxLocal-tMinGlobal)/dtGlobal)+8
    else
      maxGlobal = ceiling((tMaxLocal-tMinGlobal)/dtGlobal)+3
    end if
    if (minGlobal.le.ntGlobal .and. maxGlobal.ge.1) then
      inRange = .true.      
    end if
  else
    inRange = .true.
  end if
end function SourceTimeInRange


subroutine preciseFreqTinterpolate(tLocal,qdSLocal,tGlobal,qdSPatch,iBlankArray)
  real(kind=8), intent(in), dimension(:)::tLocal, qdSLocal, tGlobal
  real(kind=8), intent(inout), dimension(:)::qdSPatch
  integer, dimension(:), optional::iBlankArray
  
  integer:: i, i1, it, ntLocal, ntGlobal
  real(kind=8)::t
  real(kind=8), allocatable, dimension(:),save::Atemp, Btemp

  ntLocal = size(tLocal)
  ntGlobal= size(tGlobal)
 !Setup up temporary arrays - only compute once
  if (.not. allocated(Atemp)) then
    allocate( Atemp(ntLocal-1), Btemp(ntLocal-1) )
  else if (size(Atemp) /= ntLocal-1) then
    deallocate (Atemp, Btemp)
    allocate( Atemp(ntLocal-1), Btemp(ntLocal-1) )
  end if
  do it=1,ntLocal-1
    Btemp(it)= (qdSLocal(it+1)-qdSLocal(it))/(tLocal(it+1)-tLocal(it))
!ksb    Atemp(it)= qdSLocal(it) - Btemp(it)*tLocal(it)
  end do
!monotonicity test and interpolation:
  i1 = 1
  do it=1,ntGlobal
    t = tGlobal(it)
    do i=i1, ntLocal-1    
      if( tLocal(i+1)>=t ) then  
        if( tLocal(i)>t ) then
          if(iBlankFlag.and.present(iBlankArray)) iBlankArray(it)=0
          exit
        end if
!ksb        qdSPatch(it) = qdSPatch(it) + Atemp(i) + Btemp(i)*t
        qdSPatch(it) = qdSPatch(it) + qdsLocal(i)+Btemp(i)*(t-tLocal(i))
        i1 = i
        exit
      end if
    end do
    if (( t > tLocal(ntLocal-1)) .and. iBlankFlag) then
      if(present(iBlankArray)) iBlankArray(it)=0
    end if
  end do
  deallocate( Atemp, Btemp )
  return

end subroutine preciseFreqTinterpolate


subroutine tinterpolate(tLocal,qdSLocal,tGlobal,qdSPatch,iBlankArray)
  real, intent(in), dimension(:)::tLocal, qdSLocal, tGlobal
  real, intent(inout), dimension(:)::qdSPatch
  integer, dimension(:), optional::iBlankArray
  
  integer:: i, i1, it, ntLocal, ntGlobal
  real::t
  real, allocatable, dimension(:)::Atemp, Btemp

  ntLocal = size(tLocal)
  ntGlobal= size(tGlobal)
!Setup up temporary arrays - only compute once
  allocate( Atemp(ntLocal-1), Btemp(ntLocal-1) )
  do it=1,ntLocal-1
    Btemp(it)= (qdSLocal(it+1)-qdSLocal(it))/(tLocal(it+1)-tLocal(it))
    Atemp(it)= qdSLocal(it) - Btemp(it)*tLocal(it)
  end do
!monotonicity test and interpolation:
  i1 = 1
  do it=1,ntGlobal
    t = tGlobal(it)
    do i=i1, ntLocal-1    
      if( tLocal(i+1)>t ) then  
        if( tLocal(i)>t ) then
          if(iBlankFlag.and.present(iBlankArray)) iBlankArray(it)=0
          exit
        end if
        qdSPatch(it) = qdSPatch(it) + Atemp(i) + Btemp(i)*t
        i1 = i
        exit
      end if
    end do
    if (( t > tLocal(ntLocal-1)) .and. iBlankFlag) then
      if(present(iBlankArray)) iBlankArray(it)=0
    end if
  end do
  deallocate( Atemp, Btemp )
  return

end subroutine tinterpolate


end module interpolateModule
