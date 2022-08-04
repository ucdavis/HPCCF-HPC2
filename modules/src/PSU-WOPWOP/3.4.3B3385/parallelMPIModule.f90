! NAME
! $LastChangedDate: 2018-07-31 13:46:39 -0400 (Tue, 31 Jul 2018) $
! $Id: parallelMPIModule.f90 3385 2018-07-31 17:46:39Z brentner $
!   mpiModule -- A module containing the data and functions relative
!                to MPI routine.

module MPIModule
  use mpi
  implicit none
  !include 'mpif.h'
 
  integer:: ierr, myid, numprocs, master
  integer:: status(MPI_STATUS_SIZE)
  
  logical:: NOT_PARALLEL, &
            OBSERVER_PARALLEL, &
            SOURCE_PARALLEL, &
            BOTH_PARALLEL
contains

  

  subroutine InitializeProgram()
    call MPI_INIT(ierr)
    call MPI_COMM_RANK (MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, numprocs, ierr)
    master            = 0
  end subroutine InitializeProgram
  
  subroutine KillAllProcesses()
    integer:: errorcode, ierror
    call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
  end subroutine KillAllProcesses

  subroutine InitializeParallelParameters()
    NOT_PARALLEL      = .false.
    OBSERVER_PARALLEL = .false.
    SOURCE_PARALLEL   = .false.
    BOTH_PARALLEL     = .false.
  end subroutine InitializeParallelParameters
  
  subroutine BroadcastParallelParameters()
    call BroadcastLogical(NOT_PARALLEL)
    call BroadcastLogical(OBSERVER_PARALLEL)
    call BroadcastLogical(SOURCE_PARALLEL)
    call BroadcastLogical(BOTH_PARALLEL)
  end subroutine BroadcastParallelParameters
  
  subroutine NotParallel()
    NOT_PARALLEL=.true.
  end subroutine NotParallel
  
  subroutine BothParallel()
    BOTH_PARALLEL=.true.
  end subroutine BothParallel

  subroutine ObserverParallel()
    OBSERVER_PARALLEL=.true.
  end subroutine ObserverParallel

  subroutine SourceParallel()
    SOURCE_PARALLEL=.true.
  end subroutine SourceParallel
  
  function SourceParallelized() result(res)
    logical::res
    res=SOURCE_PARALLEL    
  end function SourceParallelized
  
  function ObserverParallelized() result(res)
    logical::res
    res=OBSERVER_PARALLEL    
  end function ObserverParallelized
  
  function BothParallelized() result(res)
    logical::res
    res=BOTH_PARALLEL    
  end function BothParallelized
  
  function NotParallelized() result(res)
    logical::res
    res=NOT_PARALLEL    
  end function NotParallelized
  
  function Parallelized() result(res)
    logical::res
    res=.not.NotParallelized()    
  end function Parallelized
  
  function GetTime() result(res)
    real(kind=8)::res
    res=MPI_WTIME()
  end function GetTime
  
  subroutine FinalizeProgram()
   call MPI_FINALIZE(ierr)
  end subroutine FinalizeProgram

  subroutine Barrier()
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine Barrier

  function NotMaster() result(flag)
    logical::flag
    flag=.not.IsMaster()
  end function NotMaster
  
  function IsMaster() result(flag)
    logical::flag
    if(myid==0) then
      flag=.true.
    else
      flag=.false.
    end if
  end function IsMaster
  
  function ProcId() result(id)
    integer:: id
    id = myid
  end function
  
  function GetMaster() result(value)
    integer::value  
    value=master
  end function GetMaster

  subroutine SetMasterValue(value)
    integer::value
    value=master
  end subroutine SetMasterValue
  
  function MasterDoesComputation() result(flag)
    logical::flag
    if(GetNbProcs().eq.1) then
      flag=.true.
    else
      flag=.false.
    end if
  end function MasterDoesComputation
  
  function GetNbProcs() result(nbProcs)
    integer::nbProcs
    nbProcs = numprocs
  end function GetNbProcs
  
  function GetNbSlaves() result(nbSlaves)
    integer::nbSlaves
    nbSlaves = numprocs-1
  end function GetNbSlaves

  function GetAnySource() result(source)
    integer::source
    source = MPI_ANY_SOURCE
  end function GetAnySource
  
  function GetAnyTag() result(tag)
    integer::tag
    tag = MPI_ANY_TAG
  end function GetAnyTag

  subroutine SetTagValue(tag)
    integer::tag
    tag    = status(MPI_TAG) 
  end subroutine SetTagValue
 
  subroutine SetSenderValue(source)
    integer::source
    source = status(MPI_SOURCE)
  end subroutine SetSenderValue
   
  subroutine CheckIErr(ierr)
    integer::ierr, stringLen,ierror
    character(len=1024)::errorString
    if (ierr.ne.0)  then
      call MPI_ERROR_STRING(ierr, errorString, stringLen,ierror)
      write(*,*) trim(errorString) 
    end if
  end subroutine CheckIErr
  
  subroutine BroadcastLogical (logicalIn)
    logical, intent(inout)::logicalIn
    call MPI_BCAST(logicalIn,    1,                   MPI_LOGICAL,   0, MPI_COMM_WORLD,ierr)
  end subroutine BroadcastLogical
  
  subroutine BroadcastLogicals (logicalIn)
    logical, dimension(:), intent(inout)::logicalIn
    call MPI_BCAST(logicalIn,    size(logicalIn),     MPI_LOGICAL,   0, MPI_COMM_WORLD,ierr)
  end subroutine BroadcastLogicals
  
  subroutine BroadcastInteger (integerIn)
    integer, intent(inout)::integerIn
    call MPI_BCAST(integerIn,    1,                   MPI_INTEGER,   0, MPI_COMM_WORLD,ierr)
    call CheckIErr(ierr)
  end subroutine BroadcastInteger
  
  subroutine BroadcastIntegers(integerArray)
     integer, dimension(:), intent(inout)::integerArray
    call MPI_BCAST(integerArray, size(integerArray),  MPI_INTEGER,   0, MPI_COMM_WORLD,ierr)
    call CheckIErr(ierr)
  end subroutine BroadcastIntegers
  
  subroutine BroadcastReal    (realIn)
    real, intent(inout)::realIn
    call MPI_BCAST(realIn,       1,                   MPI_REAL,      0, MPI_COMM_WORLD,ierr)
    call CheckIErr(ierr)
  end subroutine BroadcastReal
  
  subroutine BroadcastReals   (realArray)
    real, dimension(:), intent(inout)::realArray
    call MPI_BCAST(realArray,    size(realArray),     MPI_REAL,      0, MPI_COMM_WORLD,ierr)
    call CheckIErr(ierr)
  end subroutine BroadcastReals
  
  subroutine BroadcastIntegerss  (intArray)
    integer, dimension(:,:), intent(inout)::intArray
    integer::arraySize
    arraySize=size(intArray,1)*size(intArray,2)
    call MPI_BCAST(intArray,    arraySize,           MPI_INTEGER,    0, MPI_COMM_WORLD,ierr)
    call CheckIErr(ierr)
  end subroutine BroadcastIntegerss
  
  subroutine BroadcastRealss  (realArray)
    real, dimension(:,:), intent(inout)::realArray
    integer::arraySize
    arraySize=size(realArray,1)*size(realArray,2)
    call MPI_BCAST(realArray,    arraySize,           MPI_REAL,      0, MPI_COMM_WORLD,ierr)
    call CheckIErr(ierr)
  end subroutine BroadcastRealss
  
  subroutine BroadcastRealsss (realArray)
    real, dimension(:,:,:), intent(inout)::realArray
    integer::arraySize
    arraySize=size(realArray,1)*size(realArray,2)*size(realArray,3)
    call MPI_BCAST(realArray,    arraySize,           MPI_REAL,      0, MPI_COMM_WORLD,ierr)
    call CheckIErr(ierr)
  end subroutine BroadcastRealsss

  subroutine BroadcastRealssss (realArray)
    real, dimension(:,:,:,:), intent(inout)::realArray
    integer::arraySize
    arraySize=size(realArray,1)*size(realArray,2)*size(realArray,3)*size(realArray,4)
    call MPI_BCAST(realArray,    arraySize,           MPI_REAL,      0, MPI_COMM_WORLD,ierr)
    call CheckIErr(ierr)
  end subroutine BroadcastRealssss
  
  subroutine BroadcastString  (stringIn)
    character(len=*), intent(inout)::stringIn
    call MPI_BCAST(stringIn,     len(stringIn),      MPI_CHARACTER, 0, MPI_COMM_WORLD,ierr)
    call CheckIErr(ierr)
  end subroutine BroadcastString
  
  subroutine BroadcastStrings  (stringArray)
    character(len=*), dimension(:), intent(inout)::stringArray
    integer::arraySize
    arraySize=len(stringArray(1))*size(stringArray)
    call MPI_BCAST(stringArray, arraySize, MPI_CHARACTER, 0, MPI_COMM_WORLD,ierr)
    call CheckIErr(ierr)
  end subroutine BroadcastStrings
  
  subroutine ReceiveLogical   (logicalIn,    tag, source)
    logical, intent(inout)::logicalIn
    integer, intent(in)   ::source, tag
    call MPI_RECV(logicalIn,    1, MPI_LOGICAL, source, tag, MPI_COMM_WORLD, status, ierr)
    call CheckIErr(ierr)
  end subroutine ReceiveLogical

  subroutine SendLogical      (logicalIn,    tag, source)
    logical, intent(in)::logicalIn
    integer, intent(in)::source, tag
    call MPI_SSEND(logicalIn,    1, MPI_LOGICAL,  tag, source, MPI_COMM_WORLD, ierr)
    call CheckIErr(ierr)
  end subroutine SendLogical
  
  subroutine ReceiveString       (stringIn,     tag, source)
    character(len=4096), intent(inout)  ::stringIn
    integer, intent(in)::source, tag
    call MPI_RECV(stringIn, 4096, MPI_CHARACTER, source, tag, MPI_COMM_WORLD, status, ierr)
    call CheckIErr(ierr)
  end subroutine ReceiveString
 
  subroutine SendString       (stringIn,     tag, source)
    character(len=4096), intent(in)  ::stringIn
    integer, intent(in)::source, tag
    call MPI_SSEND(stringIn, 4096, MPI_CHARACTER, tag, source, MPI_COMM_WORLD, ierr)
    call CheckIErr(ierr)
  end subroutine SendString
    
  subroutine ReceiveReals     (RealArray,    tag, source)
    real, dimension(:)::realArray
    integer, intent(in)::source, tag
    call MPI_RECV(realArray, size(realArray), MPI_REAL,&
                  source, tag, MPI_COMM_WORLD, status, ierr)
  end subroutine ReceiveReals
  
  subroutine ReceiveRealss   (RealArray,    tag, source)
    real, dimension(:,:)::realArray
    integer, intent(in)::source, tag
    integer::arraySize
    arraySize=size(realArray,1)*size(realArray,2)
    call MPI_RECV(realArray,    arraySize,          MPI_REAL,    &
                  source, tag, MPI_COMM_WORLD, status, ierr)
  end subroutine ReceiveRealss
    
  subroutine ReceiveRealsss   (RealArray,    tag, source)
    real, dimension(:,:,:)::realArray
    integer, intent(in)::source, tag
    integer::arraySize
    arraySize=size(realArray,1)*size(realArray,2)*size(realArray,3)
    call MPI_RECV(realArray,    arraySize,          MPI_REAL,    &
                  source, tag, MPI_COMM_WORLD, status, ierr)
  end subroutine ReceiveRealsss
    
  subroutine ReceiveRealssss  (RealArray,    tag, source)
    real, dimension(:,:,:,:)::realArray
    integer, intent(in)::source, tag
    integer::arraySize
    arraySize=size(realArray,1)*size(realArray,2)*size(realArray,3)*size(realArray,4)
    call MPI_RECV(realArray,    arraySize,          MPI_REAL,    &
                  source, tag, MPI_COMM_WORLD, status, ierr)
  end subroutine ReceiveRealssss

  subroutine SendReals        (realArray,    tag, source)
    real, dimension(:), intent(in)::realArray
    integer, intent(in)::tag, source
    call MPI_SSEND(realArray,    size(realArray),    MPI_REAL,    &
                  tag, source, MPI_COMM_WORLD, ierr)
  end subroutine SendReals

  subroutine SendRealss       (realArray,    tag, source)
    real, dimension(:,:), intent(in)::realArray
    integer, intent(in)::tag, source
    integer::arraySize
    arraySize=size(realArray,1)*size(realArray,2)
    call MPI_SSEND(realArray,    arraySize,          MPI_REAL,    &
                  tag, source, MPI_COMM_WORLD, ierr)
  end subroutine SendRealss

  subroutine SendRealsss      (realArray,    tag, source)
    real, dimension(:,:,:), intent(in)::realArray
    integer, intent(in)::tag, source
    integer::arraySize
    arraySize=size(realArray,1)*size(realArray,2)*size(realArray,3)
    call MPI_SSEND(realArray,    arraySize,          MPI_REAL,    &
                  tag, source, MPI_COMM_WORLD, ierr)
  end subroutine SendRealsss

  subroutine SendRealssss     (realArray,    tag, source)
    real, dimension(:,:,:,:), intent(in)::realArray
    integer, intent(in)::tag, source
    integer::arraySize
    arraySize=size(realArray,1)*size(realArray,2)*size(realArray,3)*size(realArray,4)
    call MPI_SSEND(realArray,    arraySize,          MPI_REAL,    &
                  tag, source, MPI_COMM_WORLD, ierr)
  end subroutine SendRealssss

  subroutine ReceiveReal      (realIn,       tag, source)
    real::realIn
    integer, intent(in)::source, tag
    call MPI_RECV(realIn,       1,                  MPI_REAL,    &
                  source, tag, MPI_COMM_WORLD, status, ierr)
  end subroutine ReceiveReal

  subroutine SendReal         (realIn,       tag, source)
    real, intent(in)::realIn
    integer, intent(in)::source, tag
    call MPI_SSEND(realIn,       1,                  MPI_REAL,    &
                  tag, source, MPI_COMM_WORLD, ierr)
  end subroutine SendReal
  
  subroutine ReceiveIntegers  (integerArray, tag, source)
    integer, dimension(:)::integerArray
    integer, intent(in)::source, tag
    call MPI_RECV(integerArray, size(integerArray), MPI_INTEGER, &
                  source, tag, MPI_COMM_WORLD, status, ierr)
  end subroutine ReceiveIntegers

  subroutine ReceiveIntegersss(integerArray, tag, source)
    integer, dimension(:,:,:)::integerArray
    integer, intent(in)::source, tag
    integer::arraySize
    arraySize=size(integerArray,1)*size(integerArray,2)*size(integerArray,3)
    call MPI_RECV(integerArray, arraySize,          MPI_INTEGER, &
                  source, tag, MPI_COMM_WORLD, status, ierr)
  end subroutine ReceiveIntegersss

  subroutine SendIntegers     (integerArray, tag, source)
    integer, dimension(:), intent(in)::integerArray
    integer, intent(in)::tag, source
    call MPI_SSEND(integerArray, size(integerArray), MPI_INTEGER, &
                  tag, source, MPI_COMM_WORLD, ierr)
  end subroutine SendIntegers

  subroutine SendIntegersss   (integerArray, tag, source)
    integer, dimension(:,:,:), intent(in)::integerArray
    integer, intent(in)::tag, source
    integer::arraySize
    arraySize=size(integerArray,1)*size(integerArray,2)*size(integerArray,3)
    call MPI_SSEND(integerArray, arraySize, MPI_INTEGER, &
                  tag, source, MPI_COMM_WORLD, ierr)
  end subroutine SendIntegersss

  subroutine ReceiveInteger   (integerIn,    tag, source)
    integer::integerIn
    integer, intent(in)::source, tag
    call MPI_RECV(integerIn,    1,                  MPI_INTEGER, &
                  source, tag, MPI_COMM_WORLD, status, ierr)
  end subroutine ReceiveInteger

  subroutine SendInteger      (integerIn,    tag, source)
    integer, intent(in)::integerIn
    integer, intent(in)::source, tag
    call MPI_SSEND(integerIn,    1,                  MPI_INTEGER, &
                  tag, source, MPI_COMM_WORLD, ierr)
  end subroutine SendInteger
  
  subroutine SmartBroadcastLogicals(array)
    character(len=*), dimension(:), pointer::array
    logical::allocatedFlag
    integer::arraySize
    
    allocatedFlag = .false.
    arraySize = 0
    if(isMaster()) then
      allocatedFlag = associated(array)
      if (allocatedFlag) then
        arraySize = size(array)
      end if
    end if
    call BroadcastLogical(allocatedFlag)
    call BroadcastInteger(arraySize)
    if (NotMaster()) then
      if (allocatedFlag.and.arraySize.ne.0) then
        allocate(array(arraySize))
      else
        nullify(array)
      end if
    end if
    if (allocatedFlag.and.arraySize.ne.0) then
      call BroadcastStrings(array)
    end if
  end subroutine SmartBroadcastLogicals
  
  subroutine SmartBroadcastReals(array)
    real, dimension(:), pointer::array
    logical::allocatedFlag
    integer::arraySize
    
    allocatedFlag = .false.
    arraySize = 0
    if(isMaster()) then
      allocatedFlag = associated(array)
      if (allocatedFlag) then
        arraySize = size(array)
      end if
    end if
    call BroadcastLogical(allocatedFlag)
    call BroadcastInteger(arraySize)
    if (NotMaster()) then
      if (allocatedFlag.and.arraySize.ne.0) then
        allocate(array(arraySize))
      else
        nullify(array)
      end if
    end if
    if (allocatedFlag) then
      call BroadcastReals(array)
    end if
  end subroutine SmartBroadcastReals
    
end module mpiModule
