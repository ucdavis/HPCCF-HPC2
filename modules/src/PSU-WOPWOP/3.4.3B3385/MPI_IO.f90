! PSU-WOPWOP
! $Id: WindowsIO.f90 3298 2013-10-18 02:12:03Z brentner $
! $LastChangedDate: 2013-10-17 22:12:03 -0400 (Thu, 17 Oct 2013) $
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

module IOModule
  use MPIModule
  use debugModule
  use DFLIB ! For MAKEDIRQQ
  implicit none

  public :: OpenBinaryFile, CreateBinaryFile, CloseBinaryFile, RewindBinaryFile
  public :: ReadBinaryInteger, ReadBinaryReal, ReadBinaryReal1DArray, &
            ReadBinaryString, ReadBinaryReal2DArray, ReadBinaryReal3DArray, &
            ReadBinaryInteger2DArray, MakeDirectory
  public :: RestartBinaryFileAtByte, GetBinaryFilePosition
  public :: WriteBinaryString, WriteBinaryReal, WriteBinaryInteger, &
            WriteBinaryInteger1DArray,WriteBinaryInteger2DArray, &
            WriteBinaryShortInteger1DArray, &
            WriteBinaryReal1DArray, WriteBinaryReal2DArray, &
            WriteBinaryReal3DArray

  ! We'll go ahead and hard-code a 5000 file-pointer limit to reduce the
  ! complexity of the code:
  integer, dimension(5000) :: mg_filePosition, fh
  
  public :: InitializeStreamNumberCounter, GetStreamNumber
  
  integer::streamNumberCounter
            
  contains

  subroutine InitializeStreamNumberCounter()
    streamNumberCounter=100
  end subroutine InitializeStreamNumberCounter

  function GetStreamNumber() result(streamNumber)
    integer::streamNumber
    streamNumber=streamNumberCounter
    streamNumberCounter=streamNumberCounter+1
  end function GetStreamNumber

  
  subroutine MakeDirectory (directory)
    
    character (len=*), intent(in):: directory
    logical :: success
    integer :: err

    success = MAKEDIRQQ (trim(directory))
    if (.not. success) then
      err = GETLASTERRORQQ ()
      if (err /= ERR$EXIST) then
        call Error ("Could not create '"//trim(directory)//"'")
        stop
      end if
    end if

  end subroutine MakeDirectory
  
  
  subroutine OpenBinaryFile (unit, filename, convertFlag, ierr)
    integer, intent(in) :: unit
    integer, intent(inout):: ierr
    character (len=*), intent(in):: filename
    logical, intent(in)::convertFlag
    character(len=20)::convert

    ! Determine if conversion is necessary (to big_endian).
    if (convertFlag) then
      convert = "big_endian"
      print*, 'convert is not supported in MPI_IO (yet?)' !stop
      stop
    else
      convert = "native"  !default conversion
    end if

    call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_RDONLY, &
      MPI_INFO_NULL, fh(unit), ierr)

    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = 0
    end if

  end subroutine OpenBinaryFile
  
  
  subroutine CreateBinaryFile (unit, filename, convertFlag, ierr)
    integer, intent(in) :: unit
    integer, intent(inout):: ierr
    character (len=*), intent(in):: filename
    logical, intent(in)::convertFlag
    character(len=20)::convert

    ! Determine if conversion is necessary (to big_endian).
    if (convertFlag) then
      convert = "big_endian"
      print*, 'convert is not supported in MPI_IO (yet?)' !stop
    else
      convert = "native"  !default conversion
    end if
    
    call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
      MPI_INFO_NULL, fh(unit), ierr)
  
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = 0
    end if

  end subroutine CreateBinaryFile

  subroutine CloseBinaryFile (unit)
    integer, intent(in) :: unit
    integer:: ierr
    call MPI_FILE_CLOSE(fh(unit),ierr)
  end subroutine CloseBinaryFile
  
  subroutine RewindBinaryFile (unit)
    integer, intent(in) :: unit
    integer(kind=MPI_OFFSET_KIND):: offset
    integer:: ierr
    
    offset=0
    call MPI_FILE_SEEK(fh(unit),offset,MPI_SEEK_SET,ierr)

    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = 0
    end if
    
  end subroutine RewindBinaryFile

  subroutine ReadBinaryInteger (unit, data, ierr)
    integer, intent(in) :: unit
    integer, intent(inout) :: ierr
    integer, intent(out) :: data
    integer, dimension(MPI_STATUS_SIZE):: status
    
    call MPI_FILE_READ(fh(unit),data,1,MPI_INTEGER, status, ierr)
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + 4
    end if

  end subroutine ReadBinaryInteger


  subroutine ReadBinaryReal (unit, data, ierr)
    integer, intent(in) :: unit
    integer, intent(inout) :: ierr
    real(kind=4), intent(out) :: data
    integer, dimension(MPI_STATUS_SIZE):: status

    call MPI_FILE_READ(fh(unit),data,1,MPI_REAL, status, ierr)
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + 4
    end if

  end subroutine ReadBinaryReal

  subroutine ReadBinaryReal1DArray (unit, data, ierr)
    integer, intent(in) :: unit
    integer, intent(inout) :: ierr
    real(kind=4), dimension(:), intent(out) :: data
    integer, dimension(MPI_STATUS_SIZE):: status
    
    call MPI_FILE_READ(fh(unit),data,size(data),MPI_REAL, status, ierr)
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + 4*size(data)
    end if

  end subroutine ReadBinaryReal1DArray

  subroutine ReadBinaryReal2DArray (unit, data, ierr)
    integer, intent(in) :: unit
    integer, intent(inout) :: ierr
    real(kind=4), dimension(:,:), intent(out) :: data
    integer, dimension(MPI_STATUS_SIZE):: status

    call MPI_FILE_READ(fh(unit),data,size(data),MPI_REAL, status, ierr)
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + &
                              4*size(data,1)*size(data,2)
    end if

  end subroutine ReadBinaryReal2DArray

  subroutine ReadBinaryReal3DArray (unit, data, ierr)
    integer, intent(in) :: unit
    integer, intent(inout) :: ierr
    real(kind=4), dimension(:,:,:), intent(out) :: data
    integer, dimension(MPI_STATUS_SIZE):: status

    call MPI_FILE_READ(fh(unit),data,size(data),MPI_REAL, status, ierr)
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + &
                              4*size(data,1)*size(data,2)*size(data,3)
    end if
      
  end subroutine ReadBinaryReal3DArray

  subroutine ReadBinaryString (unit, data, ierr)
    integer, intent(in) :: unit
    integer, intent(inout) :: ierr
    character(len=*), intent(inout)::data
    integer, dimension(MPI_STATUS_SIZE):: status

    call MPI_FILE_READ(fh(unit),data,len(data),MPI_CHARACTER, status, ierr)
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + len(data)
    end if

  end subroutine ReadBinaryString

  subroutine ReadBinaryInteger2DArray (unit, data, ierr)
    integer, intent(in) :: unit
    integer, intent(inout) :: ierr
    integer, dimension(:,:), intent(out) :: data
    integer, dimension(MPI_STATUS_SIZE):: status

    call MPI_FILE_READ(fh(unit),data,size(data),MPI_INTEGER, status, ierr) 
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + &
                              4*size(data,1)*size(data,2)
    end if

  end subroutine ReadBinaryInteger2DArray

  subroutine RestartBinaryFileAtByte (unit, byte, status)
    integer, intent(in) :: unit
    integer, intent(in) :: byte
    integer, intent(out) :: status
    integer(kind=MPI_OFFSET_KIND):: offset
    character :: oneByte
    integer :: i, ierror
    
    ! Set the file position to byte
    offset = byte
    call MPI_FILE_SEEK(fh(unit),offset,MPI_SEEK_SET,status)
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = byte
    end if
    
  end subroutine RestartBinaryFileAtByte

  function GetBinaryFilePosition (unit) result (byte)
    integer, intent(in) :: unit
    integer :: byte

    if (unit <= size(mg_filePosition)) then
      byte = mg_filePosition(unit)
    else
      byte = 0
    end if
    
  end function GetBinaryFilePosition


  subroutine WriteBinaryString (unit, outputData, ierr)
    integer, intent(in) :: unit
    character (len=*), intent(in) :: outputData
    integer, intent(out) :: ierr 
    integer, dimension(MPI_STATUS_SIZE):: status
    
    call MPI_FILE_WRITE(fh(unit),outputData,len(outputData),MPI_CHARACTER, status, ierr)
    !ksb debug: write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryString

  subroutine WriteBinaryInteger (unit, outputData, ierr)
    integer, intent(in)  :: unit
    integer, intent(in)  :: outputData
    integer, intent(out) :: ierr 
    integer, dimension(MPI_STATUS_SIZE):: status
    
    call MPI_FILE_WRITE(fh(unit), outputData, 1, MPI_INTEGER, status, ierr)
    !ksb debug: write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryInteger

  subroutine WriteBinaryInteger1DArray (unit, outputData, ierr)
    integer, intent(in) :: unit
    integer, dimension(:), intent(in) :: outputData
    integer, intent(out) :: ierr 
    integer, dimension(MPI_STATUS_SIZE):: status
    
    call MPI_FILE_WRITE(fh(unit), outputData, size(outputData), MPI_INTEGER, status, ierr)
    !ksb debug: write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryInteger1DArray

  subroutine WriteBinaryShortInteger1DArray (unit, outputData, ierr)
  integer, intent(in) :: unit
    integer(kind=selected_int_kind(4)), dimension(:), intent(in) :: outputData
    integer, intent(out) :: ierr 
    integer, dimension(MPI_STATUS_SIZE):: status

    call MPI_FILE_WRITE(fh(unit), outputData, size(outputData), MPI_INTEGER, status, ierr)
    !ksb debug: write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryShortInteger1DArray

  subroutine WriteBinaryReal1DArray (unit, outputData, ierr)
    integer, intent(in) :: unit
    real(kind=4), dimension(:), intent(in) :: outputData
    integer, intent(out) :: ierr 
    integer, dimension(MPI_STATUS_SIZE):: status

    call MPI_FILE_WRITE(fh(unit), outputData, size(outputData), MPI_INTEGER, status, ierr)
    !ksb debug: write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryReal1DArray
  
  subroutine WriteBinaryInteger2DArray (unit, outputData, ierr)
    integer, intent(in) :: unit
    integer, dimension(:,:), intent(in) :: outputData
    integer, intent(out) :: ierr 
    integer, dimension(MPI_STATUS_SIZE):: status
    
    call MPI_FILE_WRITE(fh(unit), outputData, size(outputData), MPI_INTEGER, status, ierr)
    !ksb debug: write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryInteger2DArray

  subroutine WriteBinaryReal (unit, outputData, ierr)
    integer, intent(in) :: unit
    real(kind=4), intent(in) :: outputData
    integer, intent(out) :: ierr 
    integer, dimension(MPI_STATUS_SIZE):: status
    
    call MPI_FILE_WRITE(fh(unit), outputData, 1, MPI_REAL, status, ierr)
    !ksb debug: write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryReal

  subroutine WriteBinaryReal2DArray (unit, outputData, ierr)
    integer, intent(in) :: unit
    real(kind=4), dimension(:,:), intent(in) :: outputData
    real(kind=4), dimension(:,:), allocatable :: outputdata2 !Temp array is used to avoid stack overflow
    integer, intent(out) :: ierr 
    integer, dimension(MPI_STATUS_SIZE):: status

    call MPI_FILE_WRITE(fh(unit), outputData, size(outputData), MPI_REAL, status, ierr)
    !ksb debug: allocate(outputdata2(size(outputData,1),size(outputData,2)))
    !outputdata2=outputData
    !write (unit=unit,iostat=status) outputData2
  end subroutine WriteBinaryReal2DArray


  subroutine WriteBinaryReal3DArray (unit, outputData, ierr)
    integer, intent(in) :: unit
    real(kind=4), dimension(:,:,:), intent(in) :: outputData
    real(kind=4), dimension(:,:,:), allocatable :: outputdata2 !Temp array is used to avoid stack overflow
    integer, intent(out) :: ierr 
    integer, dimension(MPI_STATUS_SIZE):: status

    call MPI_FILE_WRITE(fh(unit), outputData, size(outputData), MPI_REAL, status, ierr)
    !ksb debug: allocate(outputdata2(size(outputData,1),size(outputData,2),size(outputData,3)))
    !outputdata2=outputData
    !write (unit=unit,iostat=status) outputData2
  end subroutine WriteBinaryReal3DArray


  function IsGood (num) result (res)
    real, intent(in) :: num
    logical :: res
    integer :: test
     res = .true.
     return
!    res = .not. isnan(num)
    
  end function IsGood
  !function IsGood (num) result (res)
  !  real, intent(in) :: num
  !  logical :: res
  !  !
  !  ! The IEEE standard says that NaN /= NaN - so this a check for NaN
  !  if (num /= num) then
  !    res = .false.
  !  else
  !    res = .true.
  !  endif
  !  return
  !end function IsGood


end module IOModule
