! NAME
! $Id: serialMPIModule.f90 3307 2014-03-03 04:03:50Z brentner $
! $LastChangedDate: 2014-03-02 23:03:50 -0500 (Sun, 02 Mar 2014) $
!   mpiModule -- A module containing the data and functions relative
!                to MPI routine.

module MPIModule
  implicit none
 
  integer:: ierr, myid, numprocs, master
  integer:: status(10)
  
contains
    
  subroutine NotParallel()
  end subroutine NotParallel
  
  subroutine BothParallel()
  end subroutine BothParallel

  subroutine ObserverParallel()
  end subroutine ObserverParallel

  subroutine SourceParallel()
  end subroutine SourceParallel
  
  subroutine InitializeParallelParameters()
  end subroutine InitializeParallelParameters
  
  subroutine BroadcastParallelParameters()
  end subroutine BroadcastParallelParameters
  
  function SourceParallelized() result(res)
    logical::res
    res=.false.    
  end function SourceParallelized
  
  function ObserverParallelized() result(res)
    logical::res
    res=.false.
  end function ObserverParallelized
  
  function BothParallelized() result(res)
    logical::res
    res=.false.
  end function BothParallelized
  
  subroutine KillAllProcesses()
    implicit none
    stop
  end subroutine KillAllProcesses
  
  function NotParallelized() result(res)
    logical::res
    res=.true.
  end function NotParallelized
  
  function Parallelized() result(res)
    logical::res
    res=.false.
  end function Parallelized
    
  subroutine InitializeProgram()
  end subroutine InitializeProgram

  function GetTime() result(res)
    real(kind=8)::res
    integer::count,rate,max
    call system_clock(count,rate,max)
    res=dble(count)/dble(rate)
  end function GetTime
   
  subroutine FinalizeProgram()
  end subroutine FinalizeProgram

  subroutine Barrier()
  end subroutine Barrier

  function NotMaster() result(flag)
    logical::flag
    flag=.false.
  end function NotMaster
  
  function IsMaster() result(flag)
    logical::flag
    flag=.true.
  end function IsMaster

  function MasterDoesComputation() result(flag)
    logical::flag
    flag=.true.
  end function MasterDoesComputation
  
  function GetMaster() result(value)
    integer::value 
    value=0
  end function GetMaster

  function ProcID() result(value)
    integer::value
    value=0
  end function ProcID
  
  subroutine SetMasterValue(value)
    integer::value
    value=0
  end subroutine SetMasterValue
  
  function GetNbProcs() result(nbProcs)
    integer::nbProcs
    nbProcs=1
  end function GetNbProcs
  
  function GetNbSlaves() result(nbSlaves)
    integer::nbSlaves
    nbSlaves = 0
  end function GetNbSlaves

  function GetAnySource() result(source)
    integer::source
    source=0
  end function GetAnySource
  
  function GetAnyTag() result(tag)
    integer::tag
    tag=0
  end function GetAnyTag

  subroutine SetTagValue(tag)
    integer::tag
    tag=0
  end subroutine SetTagValue
 
  subroutine SetSenderValue(source)
    integer::source
    source=0
  end subroutine SetSenderValue
  
  subroutine BroadcastLogical (logicalIn)
    logical, intent(inout)::logicalIn
  end subroutine BroadcastLogical
  
  subroutine BroadcastLogicals (logicalIn)
    logical, dimension(:), intent(inout)::logicalIn
  end subroutine BroadcastLogicals
  
  subroutine BroadcastInteger (integerIn)
    integer, intent(inout)::integerIn
  end subroutine BroadcastInteger
  
  subroutine BroadcastIntegers(integerArray)
    integer, dimension(:), intent(inout)::integerArray
  end subroutine BroadcastIntegers
  
  subroutine BroadcastReal    (realIn)
    real, intent(inout)::realIn
  end subroutine BroadcastReal
  
  subroutine BroadcastReals   (realArray)
    real, dimension(:), intent(inout)::realArray
  end subroutine BroadcastReals
  
  subroutine BroadcastIntegerss  (intArray)
    integer, dimension(:,:), intent(inout)::intArray
    integer::arraySize
    arraySize=size(intArray,1)*size(intArray,2)
  end subroutine BroadcastIntegerss
  
  subroutine BroadcastRealss  (realArray)
    real, dimension(:,:), intent(inout)::realArray
    integer::arraySize
    arraySize=size(realArray,1)*size(realArray,2)
  end subroutine BroadcastRealss
  
  subroutine BroadcastRealsss (realArray)
    real, dimension(:,:,:), intent(inout)::realArray
    integer::arraySize
    arraySize=size(realArray,1)*size(realArray,2)*size(realArray,3)
  end subroutine BroadcastRealsss

  subroutine BroadcastRealssss (realArray)
    real, dimension(:,:,:,:), intent(inout)::realArray
    integer::arraySize
    arraySize=size(realArray,1)*size(realArray,2)*size(realArray,3)*size(realArray,4)
  end subroutine BroadcastRealssss
  
  subroutine BroadcastString  (stringIn)
    character(len=*), intent(inout)::stringIn
  end subroutine BroadcastString
  
  subroutine ReceiveLogical   (logicalIn,    tag, source)
    logical, intent(inout)::logicalIn
    integer, intent(in)   ::source, tag
  end subroutine ReceiveLogical

  subroutine SendLogical      (logicalIn,    tag, source)
    logical, intent(in)::logicalIn
    integer, intent(in)::source, tag
  end subroutine SendLogical
    
  subroutine ReceiveReals     (RealArray,    tag, source)
    real, dimension(:)::realArray
    integer, intent(in)::source, tag
  end subroutine ReceiveReals
    
  subroutine ReceiveRealss     (RealArray,    tag, source)
    real, dimension(:,:)::realArray
    integer, intent(in)::source, tag
  end subroutine ReceiveRealss

  subroutine ReceiveRealsss   (RealArray,    tag, source)
    real, dimension(:,:,:)::realArray
    integer, intent(in)::source, tag
  end subroutine ReceiveRealsss

  subroutine ReceiveRealssss  (RealArray,    tag, source)
    real, dimension(:,:,:,:)::realArray
    integer, intent(in)::source, tag
  end subroutine ReceiveRealssss

  subroutine SendReals        (realArray,    tag, source)
    real, dimension(:), intent(in)::realArray
    integer, intent(in)::tag, source
  end subroutine SendReals

  subroutine SendRealss       (realArray,    tag, source)
    real, dimension(:,:), intent(in)::realArray
    integer, intent(in)::tag, source
  end subroutine SendRealss

  subroutine SendRealsss      (realArray,    tag, source)
    real, dimension(:,:,:), intent(in)::realArray
    integer, intent(in)::tag, source
  end subroutine SendRealsss

  subroutine SendRealssss     (realArray,    tag, source)
    real, dimension(:,:,:,:), intent(in)::realArray
    integer, intent(in)::tag, source
  end subroutine SendRealssss

  subroutine ReceiveString      (stringIn,   tag, source)
    character(len=*)::stringIn
    integer, intent(in)::source, tag
  end subroutine ReceiveString

  subroutine SendString         (StringIn,   tag, source)
    character(len=*), intent(in)::StringIn
    integer, intent(in)::source, tag
  end subroutine SendString
  
  subroutine ReceiveReal      (realIn,       tag, source)
    real::realIn
    integer, intent(in)::source, tag
  end subroutine ReceiveReal

  subroutine SendReal         (realIn,       tag, source)
    real, intent(in)::realIn
    integer, intent(in)::source, tag
  end subroutine SendReal
  
  subroutine ReceiveIntegers  (integerArray, tag, source)
    integer, dimension(:)::integerArray
    integer, intent(in)::source, tag
  end subroutine ReceiveIntegers

  subroutine ReceiveIntegersss(integerArray, tag, source)
    integer, dimension(:,:,:)::integerArray
    integer, intent(in)::source, tag
  end subroutine ReceiveIntegersss

  subroutine SendIntegers     (integerArray, tag, source)
    integer, dimension(:), intent(in)::integerArray
    integer, intent(in)::tag, source
  end subroutine SendIntegers

  subroutine SendIntegersss     (integerArray, tag, source)
    integer, dimension(:,:,:), intent(in)::integerArray
    integer, intent(in)::tag, source
  end subroutine SendIntegersss

  subroutine ReceiveInteger   (integerIn,    tag, source)
    integer::integerIn
    integer, intent(in)::source, tag
  end subroutine ReceiveInteger

  subroutine SendInteger      (integerIn,    tag, source)
    integer, intent(in)::integerIn
    integer, intent(in)::source, tag
  end subroutine SendInteger
    
  subroutine AllocateAndBroadcastReals(array)
    real, dimension(:), pointer::array
  end subroutine AllocateAndBroadcastReals

  subroutine AllocateAndBroadcastRealss(array)
    real, dimension(:), pointer::array
  end subroutine AllocateAndBroadcastRealss

  subroutine AllocateAndBroadcastRealsss(array)
    real, dimension(:), pointer::array
  end subroutine AllocateAndBroadcastRealsss
  
  subroutine AllocateAndBroadcastRealssss(array)
    real, dimension(:), pointer::array
  end subroutine AllocateAndBroadcastRealssss
  
  subroutine AllocateAndBroadcastStrings(array)
    character(len=*), dimension(:), pointer::array
  end subroutine AllocateAndBroadcastStrings
  
end module MPIModule
