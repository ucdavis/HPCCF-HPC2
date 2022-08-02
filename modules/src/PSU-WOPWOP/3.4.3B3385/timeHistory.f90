module TimeHistoryObject
  use DebugModule
  use MPIModule
  use constantsModule, only: PNLTdBDrop
  
  implicit none

  type timeHistory
    !title: the title of this time history
    !t: the time array (domain)
    !f: the function array (range)
    !nt: the number of time points, which is also the size
    !    of the two arrays
    !tMin,tMax: the minimum and maximum time in the time history
    character(len=4096)::title
    real,dimension(:), pointer::t=>null(),f=>null()
    integer::nt
    real::tMin,tMax
  end type timeHistory
  
  interface operator (+)
    module procedure AddTimeHistories
  end interface 
  
  private
  public::timeHistory,SetupTimeHistory,NullifyTimeHistory,                  &
          CreateTimeHistorySegment,CheckAndFinalizeSegParameters,           &
          AddTimeHistories,CopyTimeHistory,                                 & 
          GetPeriod,GetNt,GettHTitle,GetTimeValue,                          &
          GetDt,GetTMin,GetTMax,GetFunctionValue,                           &
          SettHFArray,SettHTitle,SetTMin,SetTMax,SetNt,SetTimeArray,        &
          DestroyTimeHistory,DestroyTimeHistories,                          &
          ReceiveTimeHistoryFromSlave, SendTimeHistoryToMaster
          
contains

subroutine ReceiveTimeHistoryFromSlave(th, master, sender)
  type(timeHistory), intent(inout)::th
  integer, intent(in)::master, sender
  logical::flag
  integer::n
  call ReceiveString(th%title, master, sender)
  call ReceiveInteger(th%nt, master, sender)
  call ReceiveReal(th%tMin, master, sender)
  call ReceiveReal(th%tMax, master, sender)
  nullify(th%t)
  call ReceiveLogical(flag, master, sender)
  if (flag) then
    call ReceiveInteger(n, master, sender)
    allocate(th%t(n))
    call ReceiveReals(th%t, master, sender)
  end if
  nullify(th%f)
  call ReceiveLogical(flag, master, sender)
  if (flag) then
    call ReceiveInteger(n, master, sender)
    allocate(th%f(n))
    call ReceiveReals(th%f, master, sender)
  end if
end subroutine ReceiveTimeHistoryFromSlave

subroutine SendTimeHistoryToMaster(th, master, tag)
  type(timeHistory), intent(in)::th
  integer, intent(in)::master, tag
  logical::flag
  call SendString(th%title, master, tag)
  call SendInteger(th%nt, master, tag)
  call SendReal(th%tMin, master, tag)
  call SendReal(th%tMax, master, tag)
  flag = associated(th%t)
  call SendLogical(flag, master, tag)
  if (flag) then
    call SendInteger(size(th%t), master, tag)
    call SendReals(th%t, master, tag)
  end if
  flag = associated(th%f)
  call SendLogical(flag, master, tag)
  if (flag) then
    call SendInteger(size(th%f), master, tag)
    call SendReals(th%f, master, tag)
  end if
end subroutine SendTimeHistoryToMaster

!**
!FUNCTION AddTimeHistories
!  This function adds two time histories.  It checks to make sure
!  that they are of the same length, and also checks to make sure 
!  that they have the same domain arrays.  It returns a new time history.
!  This function operates on the assumption that all the variables
!  in the resulting tHsum type are the same as those in the tH1 and tH2 
!  types except for the range array and the title.
!**
function AddTimeHistories(tH1,tH2) result(tHsum)
  implicit none
  type(timeHistory),intent(in)::tH1
  type(timeHistory),intent(in)::tH2
  type(timeHistory)::tHsum
  !Check to see if the time history arrays are the same length:
  if(size(tH1%f).ne.size(tH2%f))then
    call Error('In AddTimeHistories, cannot add time histories of unequal length.', &
               'The length of the first history is '//trim(IntegerToString(GetNt(tH1)))//&
               ' and the length of the second is '//trim(IntegerToString(GetNt(tH2)))//'.')
  end if
  !Check to see if the time histories have unequal domain arrays:
  if(GetTMin(tH1).ne.GetTMin(tH2) .or. &
     GetTMax(tH1).ne.GetTMax(tH2) .or. &
     GetNt(tH1)  .ne.GetNt(tH2) )then
    call Error('In AddTimeHistories, cannot add time histories with unequal time arrays.')
  end if
  nullify(tHsum%f,tHsum%t)
  allocate(tHsum%f(tH1%nt),tHsum%t(tH1%nt))
  !All the variables of tHsum are the same as the variables of 
  !tH1 and tH2 except for the title and the range array, f.
  tHsum%t    = tH1%t  
  tHsum%f    = tH1%f+tH2%f
  tHsum%tMin = tH1%tMin
  tHsum%tMax = tH1%tMax
  tHsum%nt   = tH1%nt  
end function AddTimeHistories
 
 
!**
!SUBROUTINE CheckAndFinalizeSegParameters
!  This routine checks time segmenting parameters input as arguments, 
!  and if they're invalid it recalculates the next closest valid parameters 
!  and warns the user.  It checks the segment size, and segment step size.
!NOTES
! - A "segment" is a section out of an entire time history.  Any time history
!   can be split into segments.  Each segment is also a time history type.
!   See the PSU-WOPWOP post processing manual for more detailed information
!   regarding segments.
!ARGUMENTS:
!  tMin, tMax - the maximum and minimum time of the entire time period
!  nt - the total number of points in the entire time period
!  segmentSize - the length (in seconds) of each desired segment
!  segmentStepSize - the length (in seconds) to increment each segment
!  segmentIncrement - 0 if all segments are to be analyzed, 1 if only one
!                     segment is to be analyzed (for example, if the input segment
!                     size is equal to (tMax-tMin)
!  ntPerSegment - the number of time points per segment  
!**
subroutine CheckAndFinalizeSegParameters(tMin,tMax,nt,segmentSize,     &
                                         segmentStepSize,nbSegments,   &
                                         segmentIncrement,ntPerSegment)
  implicit none
  real,intent(in)::tMin,tMax
  integer,intent(inout)::nt
  real,intent(inout)::segmentSize,segmentStepSize
  integer,intent(out)::nbSegments,ntPerSegment
  integer,intent(inout)::segmentIncrement
  real::period,dt
  period = tMax-tMin
  dt = period/(real(nt)-1.)
  !If the specified segment size is "very close" to the total time,
  !the total time history is analyzed:
  if(abs(segmentSize-period).lt.1e-5)then
    nbSegments = 0
    ntPerSegment = nt
    call Warning('In CheckAndFinalizeSegmentParameters,',                        &
                 'The requested segment size is equal to the period of the',     &
                 'time history.  The whole time history will be analyzed.')
  !If the specified segment size is greater than the total time, 
  !  the total time history will by analyzed:
  else if(segmentSize.gt.period)then
    nbSegments = 0
    ntPerSegment = nt
    call Warning('In CheckAndFinalizeSegmentParameters,',                        &
                 'The requested segment size is larger than the period',         &
                 'of the time history. The whole time history will be analyzed.') 
  !If the specified segment step size is "very close" to the total time,
  !only the first segment is analyzed.             
  else if(abs(segmentStepSize-period).lt.1e-5)then
    nbSegments       = 1
    segmentIncrement = 1
    ntPerSegment = (segmentSize/dt)+1
    call Warning('In CheckAndFinalizeSegmentParameters,',                        &
                 'The requested segment step size is equal to the period of the',&
                 'time history.  The first segment will be analyzed.')
  !If the specified segment step size is greater than the total time,
  !only analyze the first segment:                
  else if(segmentStepSize.gt.period)then
    nbSegments       = 1
    segmentIncrement = 1
    ntPerSegment = (segmentSize/dt)+1    
    call Warning('In CheckAndFinalizeSegmentParameters,',                        &
                 'The requested segment step size is greater than the period of',&
                 'the time history.  The first segment will be analyzed.')   
  !If both the segment size and segment step size less than the total time, we have 
  !to make sure that both the segment size and the segment step size will work
  !with the discrete time history data points, and the result is equally sized 
  !and spaced segments:
  else
    call CheckSegmentSize(tMax,tMin,segmentSize, nt, dt, ntPerSegment)
    call CheckSegmentStepSize(tMin,tMax,nt,segmentSize, segmentStepSize, nbSegments)

    !nbSegments = nint((period-segmentSize)/segmentStepSize)+1    
    !ntPerSegment = (segmentSize/dt)+1
  end if
                                         
end subroutine CheckAndFinalizeSegParameters

!**
!SUBROUTINE CheckSegmentSize
!  This subroutine checks the segment size specified by the user.
!  The segment size must be an integer multiple of the observer time dt.
!  If it is invalid, nt is adjusted (and hence dt) so that this requirement
!  is met.
!ARGUMENTS
! - tMax,tMin: the maximum and minimum time of the entire time history
! - segmentSize: the length of each segment, in seconds
! - nt, dt: the number of time steps in the time history, and the 
!   observer time history step step
! - ntPerSeg: number of time steps per segment
! - nSeg: number of segments in the time history
!**
subroutine CheckSegmentSize(tMax,tMin,segmentSize, nt, dt, ntPerSeg) 
  implicit none
  real,intent(in)::tMin,tMax,segmentSize
  integer,intent(inout)::nt
  real,intent(inout):: dt
  
  integer:: ntPerSeg, nSeg

  ! Check to see if an integer multiple of dt fits in segmentSize.  If not, change
  ! nt (and hence dt) to make sure there is an integer multiple of dt in segmentSize.
  If( segmentSize/dt /= nint(segmentSize/dt) ) then
    ntPerSeg = nint(segmentSize/dt)
    nSeg = int((tMax-tMin)/segmentSize)
    nt = nSeg*ntPerSeg + 1
    dt = (tMax-tMin)/(nt-1)
    call Warning('In CheckSegmentSize, the segment size was not an integer multiple of',&
                 'dt, so nt (and thus dt) has been adjusted to ensure segmentSize is',&
                 'now an integer multiple of dt.')
    end if
    
end subroutine CheckSegmentSize

!**
!SUBROUTINE CheckSegmentStepSize
!  This routine checks the segment step size specified by the user.
!  The step size must be an integer multiple of the time step of the time history
!  because we are working with discrete time points.  If the step size
!  specified by the user is invalid, it warns the user and sets the step 
!  size to the the next closest (smaler) step size that is valid
!ARGUMENTS
! - tMax,tMin: the maximum and minimum time of the entire time history
! - nt: the number of time points in the time history defined by tMax and tMin
! - segmentSize: the lenght of a segment (which can be different from the segmentStepSize
! - segmentStepSize: the incremental distance between the starting point of each segment,
!                    in seconds
! - nbSegments:  The number of segments that fit in the time hisory.
!**
subroutine CheckSegmentStepSize(tMin,tMax,nt,segmentSize,segmentStepSize,nbSegments)
  implicit none
  real,intent(in)::tMin,tMax,segmentSize
  integer,intent(in)::nt
  real,intent(inout)::segmentStepSize 
  integer, intent(out):: nbSegments
  real:: dt,period
  period = tMax-tMin
  dt = period/(nt-1)

  if( segmentStepSize/dt /= int(segmentStepSize/dt) ) then
    ! need to adjust the segmentStepSize to be an integer multiple of dt.  
    ! Round down the segmentStepSize/dt to nearest integer and then multiply by dt.
    segmentStepSize = int(segmentStepSize/dt)*dt
    ! Send a warning to let users know SegmentStepSize was changed.
    call Warning('In CheckSegmentStepSize, the segment step size is not valid. ',  &
                 'The SegmentStepSize must be an integer multiple of the time step', &
                 'size, dt.  The SegmentStepSize was changed (reduced) to meet this', &
                 'requirement.')
  end if
  nbSegments = int(period/segmentStepSize)
  ! It is conceivable that a segmentStepSize could be bigger than the segmentSize, thus
  ! we need to check to see if the fraction of a period left is >= a segmentSize.  If it
  ! is then increase the number of segments by 1.
  ! ksb debug: I don't think this really works correctly for and needs some further work.
  ! (the code runs, but I don't think it does what I wanted.  I would like to update it
  ! do do things like overlap, etc. - but that really hasn't been incorporated yet.)
  if( nbSegments*segmentStepSize + segmentSize <= period ) then
    nbSegments = nbSegments+1
  end if

end subroutine CheckSegmentStepSize

!**
!FUNCTION CreateTimeHistorySegment
!  This function creates a time history segment out of the total time history tH.
!  It creates the segment whose number is segNum, where the first segment is 
!  segNum=1, the second is segNum=2, and so on.  It creates the segment according
!  to the segment parameters segmentSize and segmentStepSize.
!**
subroutine CreateTimeHistorySegment(tH,segmentSize,segmentStepSize,     &
                                    segNum, tHsegment)
  implicit none
  real,intent(in)::segmentSize,segmentStepSize
  integer,intent(in)::segNum
  type(timeHistory),intent(in)::tH
  type(timeHistory), intent(inout)::tHsegment
  integer::ntPerSegmentStep,nt
  real::dt,tMin,tMax
  character(len=4096)::title 
  dt    = GetDt(tH)  
  ntPerSegmentStep = nint(segmentStepSize/dt)  
  tMin  = tH%tMin+real((segNum-1)*segmentStepSize)
  tMax  = tMin + segmentSize
  nt = segmentSize/dt + 1
  title = trim(tH%title)//'_Segment#'//IntegerToString(segNum)
  call SetupTimeHistory(tHsegment,tMin,tMax,nt,title)
  !Fill the range array:
  if (associated(THsegment%f)) then
    if(size(THsegment%f).ne.nt) then
      nullify(THsegment%f)
    end if
  end if
  allocate(THsegment%f(THSegment%nt))
  tHsegment%f = tH%f((ntPerSegmentStep*(segNum-1)+1): &
                   (ntPerSegmentStep*(segNum-1)   &
                    +tHsegment%nt))
end subroutine CreateTimeHistorySegment

!**
!SUBROUTINE DestroyTimeHistory
!  This routine nullifies the pointers contained
!  in the tH type.
!**
subroutine DestroyTimeHistory(tH)
  implicit none
  type(timeHistory),intent(inout)::tH  
  if(associated(tH%f))then
    deallocate(tH%f)
    nullify(tH%f)
  end if    
  if(associated(tH%t))then
    deallocate(tH%t)
    nullify(tH%t)  
  end if  
end subroutine DestroyTimeHistory

subroutine CopyTimeHistory(tH1, tHcopy)
  implicit none
  type(timeHistory),intent(in)::tH1
  type(timeHistory),intent(inout)::tHcopy
  
  if( associated(tHcopy%t) ) deallocate(tHcopy%t)
  if( associated(tH1%t) ) then
    allocate(tHcopy%t(th1%nt))
    tHcopy%t = tH1%t
  end if
  
  if( associated(tHcopy%f) ) deallocate(tHcopy%f)
  if( associated(tH1%f) ) then
    allocate(tHcopy%f(th1%nt))
    tHcopy%f = tH1%f
  end if
 
  tHcopy%nt    = tH1%nt
  tHcopy%tMin  = tH1%tMin
  tHcopy%tMax  = tH1%tMax
  tHcopy%title = tH1%title
end subroutine CopyTimeHistory

!**
!SUBROUTINE DestroyTimeHistories
!  This routine destroys an array of time histories
!  tH by calling the routine DestroyTimeHistory
!  for each time history in the tH array.
!**
subroutine DestroyTimeHistories(tH)
  implicit none
  type(timeHistory),dimension(:),intent(inout)::tH
  integer::i
  do i=1,size(tH)
    call DestroyTimeHistory(tH(i))
  end do
end subroutine DestroyTimeHistories

!**
!FUNCTION GetDt
!  This routine calculates and returns the time step
!  of an input time history according to the length
!  of the time history (tMax - tMin) and the number of time points 
!  it contains (nt).
!**
function GetDt(tH) result(dt)
  implicit none
  type(timeHistory)::tH
  real::dt
  if(tH%nt==1)then
    dt = 0.
  else
    dt = (tH%tMax-tH%tMin)/(tH%nt-1)
  end if  
end function GetDt

!**
!FUNCTION GetFunctionValue
!  This function returns the function value corresponding
!  to the nt entry of the tH%f array.  It returns tH%f(nt).
!**
function GetFunctionValue(tH,nt) result(f)
  implicit none
  type(timeHistory)::tH
  integer::nt
  real::f
  f = tH%f(nt)
end function GetFunctionValue

!**
!FUNCTION GetNt
!  This function returns the number of time points
!  in the input time history, tH%nt.
!**
function GetNt(tH)result(nt)
  implicit none
  type(timeHistory),intent(in)::tH
  integer::nt
  nt = tH%nt
end function GetNt

!**
!FUNCTION GetPeriod
!  This function returns the length of the time history,
!  tH%tMax - tH%tMin
!**
function GetPeriod(tH) result(period)
  implicit none
  type(timeHistory)::tH
  real::period
  period = GetTMax(tH)-GetTMin(tH)
end function GetPeriod

!**
!FUNCTION GettHTitle
!  This function returns the character string
!  tH%title which is the title of the time history.
!**
function GettHTitle(tH) result(title)
  implicit none
  type(timeHistory)::tH
  character(len=4096)::title
  title = tH%title
end function GettHTitle

!**
!FUNCTION GetTimeValue
!  This function returns the time value corresponding
!  to the nt entry of the tH%t array.  It returns tH%t(nt).
!**
function GetTimeValue(tH,nt) result(time)
  implicit none
  type(timeHistory)::tH
  integer::nt
  real::time
  time = tH%t(nt)
end function GetTimeValue

!**
!FUNCTION GetTMax
!  This function returns the maximum time of the
!  time history, tH%tMax.
!**
function GetTMax(tH) result(tMax)
  implicit none
  type(timeHistory)::tH
  real::tMax
  tMax = tH%tMax
end function GetTMax

!**
!FUNCTION GetTMin
!  This function returns the minimum time of the
!  time history, tH%tMin.
!**
function GetTMin(tH) result(tMin)
  implicit none
  type(timeHistory)::tH
  real::tMin
  tMin = tH%tMin
end function GetTMin

!**
!SUBROUTINE SettHFArrayTMin
!  This routine sets the function array
!  of the input time history, tH%f.  It first checks
!  to ensure that the length of the input function array
!  is equal to the number of time points in the time history,
!  tH%nt.
!**
subroutine SettHFArray(tH,fArray)
  implicit none
  type(timeHistory),intent(inout)::tH
  real,dimension(:),intent(in)::fArray
  if(size(fArray).ne.th%nt)then
    call Error('In SettHFArray, the length of the input array (' &
                //trim(IntegerToString(size(fArray)))// &
               ') does not equal nt ('//trim(IntegerToString(GetNt(tH)))//').')
  end if
  if (.not.associated(TH%f)) then
    allocate(tH%f(th%nt))
  end if
  tH%f = fArray
end subroutine SettHFArray

!**
!SUBROUTINE SetNt
!  This routine sets the number of time history points
!  in a time history, tH%nt.
!**
subroutine SetNt(tH,nt)
  implicit none
  type(timeHistory),intent(inout)::tH
  integer,intent(in)::nt
  tH%nt = nt
end subroutine SetNt   

!**
!SUBROUTINE SettHTitle
!  This routine sets the time history title, tH%title
!** 
subroutine SettHTitle(tH,title)
  implicit none
  type(timeHistory),intent(inout)::tH
  character(len=*),intent(in)::title
  tH%title = trim(title)
end subroutine SettHTitle

!**
!SUBROUTINE SetTimeArray
!  This routine sets the time time array
!  of the time history, tH%t
!** 
subroutine SetTimeArray(tH,timeArray)
  implicit none
  type(timeHistory),intent(inout)::tH
  real,dimension(:),intent(in)::timeArray
  if (.not.associated(TH%t)) then
    allocate(tH%t(size(timeArray)))
  end if
  tH%t = timeArray
end subroutine SetTimeArray

!**
!SUBROUTINE SetTMax
!  This routine sets the maximum time of a time
!  history, tH%tMax.  Make sure this value is equal
!  to the last entry of the time array, tH%t(nt).
!** 
subroutine SetTMax(tH,tMax)
  implicit none
  type(timeHistory),intent(inout)::tH
  real,intent(in)::tMax
  tH%tMax = tMax
end subroutine SetTMax

!**
!SUBROUTINE SetTMax
!  This routine sets the minimum time of a time
!  history, tH%tMin.  Make sure this value is equal
!  to the first entry of the time array, tH%t(1).
!** 
subroutine SetTMin(tH,tMin)
  implicit none
  type(timeHistory),intent(inout)::tH
  real,intent(in)::tMin
  tH%tMin = tMin
end subroutine SetTMin

subroutine SetupTimeHistory(tH,tMin,tMax,nt,title) 
  implicit none
  real,intent(in)::tMin,tMax
  integer,intent(in)::nt
  character(len=*),intent(in)::title
  integer::i
  type(timeHistory)::tH
  call SetNt(tH,nt)
  call SetTMin(tH,tMin)
  call SetTMax(tH,tMax)
  if (associated(TH%t)) then
    if(size(TH%t).ne.nt) then
      nullify(TH%t)
    end if
  end if
  if (.not.associated(TH%t)) then
    allocate(TH%t(nt))
  end if
  do i=1,nt
    tH%t(i) = tMin+(i-1)*GetDt(tH)
  end do
  call SettHTitle(tH,trim(title))
end subroutine SetupTimeHistory
 
subroutine NullifyTimeHistory(tH)
  implicit none
  type(timeHistory),intent(inout)::tH
  nullify(tH%t,tH%f)
end subroutine NullifyTimeHistory

end module TimeHistoryObject


