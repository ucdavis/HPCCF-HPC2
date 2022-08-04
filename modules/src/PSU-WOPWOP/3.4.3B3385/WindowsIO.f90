! PSU-WOPWOP
! $Id: WindowsIO.f90 3340 2015-03-27 21:08:44Z brentner $
! $LastChangedDate: 2015-03-27 17:08:44 -0400 (Fri, 27 Mar 2015) $
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
  integer, dimension(5000) :: mg_filePosition
  
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
  
  
  subroutine OpenBinaryFile (unit, filename, convertFlag, status)
    integer, intent(in) :: unit
    integer, intent(inout):: status
    character (len=*), intent(in):: filename
    logical, intent(in)::convertFlag
    character(len=20)::convert

    ! Determine if conversion is necessary (to big_endian).
  if (convertFlag) then
    convert = "big_endian"
  else
    convert = "native"  !default conversion
  end if
  open(unit=unit, file=trim(filename), form='unformatted',access='stream', &
       status="old",convert=trim(convert), iostat=status)

  if (unit <= size(mg_filePosition)) then
    mg_filePosition(unit) = 0
  end if

  end subroutine OpenBinaryFile
  
  
  subroutine CreateBinaryFile (unit, filename, convertFlag, status)
    integer, intent(in) :: unit
    integer, intent(inout):: status
    character (len=*), intent(in):: filename
    logical, intent(in)::convertFlag
    character(len=20)::convert

    ! Determine if conversion is necessary (to big_endian).
    if (convertFlag) then
      convert = "big_endian"
    else
      convert = "native"  !default conversion
    end if

    open(unit=unit, file=trim(filename), form='unformatted', access='stream', &
         status="replace",convert=trim(convert), iostat=status)
  
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = 0
    end if

  end subroutine CreateBinaryFile

  subroutine CloseBinaryFile (unit)
    integer, intent(in) :: unit
    close(unit)
  end subroutine CloseBinaryFile

  subroutine RewindBinaryFile (unit)
    integer, intent(in) :: unit
    rewind(unit)
  end subroutine RewindBinaryFile

  subroutine ReadBinaryInteger (unit, data, status)
    integer, intent(in) :: unit
    integer, intent(inout) :: status
    integer, intent(out) :: data
    
    read(unit,iostat=status) data
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + 4
    end if

  end subroutine ReadBinaryInteger


  subroutine ReadBinaryReal (unit, data, status)
    integer, intent(in) :: unit
    integer, intent(inout) :: status
    real(kind=4), intent(out) :: data

    read(unit, iostat=status) data
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + 4
    end if

  end subroutine ReadBinaryReal

  subroutine ReadBinaryReal1DArray (unit, data, status)
    integer, intent(in) :: unit
    integer, intent(inout) :: status
    real(kind=4), dimension(:), intent(out) :: data

    read(unit, iostat=status) data
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + 4*size(data)
    end if

  end subroutine ReadBinaryReal1DArray

  subroutine ReadBinaryReal2DArray (unit, data, status)
    integer, intent(in) :: unit
    integer, intent(inout) :: status
    real(kind=4), dimension(:,:), intent(out) :: data
  
    read(unit, iostat=status) data
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + &
                              4*size(data,1)*size(data,2)
    end if

  end subroutine ReadBinaryReal2DArray

  subroutine ReadBinaryReal3DArray (unit, data, status)
    integer, intent(in) :: unit
    integer, intent(inout) :: status
    real(kind=4), dimension(:,:,:), intent(out) :: data

    read(unit, iostat=status) data
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + &
                              4*size(data,1)*size(data,2)*size(data,3)
    end if
      
  end subroutine ReadBinaryReal3DArray

  subroutine ReadBinaryString (unit, data, status)
    integer, intent(in) :: unit
    integer, intent(inout) :: status
    character(len=*), intent(inout)::data

    read(unit, iostat=status) data
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + len(data)
    end if

  end subroutine ReadBinaryString

  subroutine ReadBinaryInteger2DArray (unit, data, status)
    integer, intent(in) :: unit
    integer, intent(inout) :: status
    integer, dimension(:,:), intent(out) :: data
  
    read(unit, iostat=status) data
      
    if (unit <= size(mg_filePosition)) then
      mg_filePosition(unit) = mg_filePosition(unit) + &
                              4*size(data,1)*size(data,2)
    end if

  end subroutine ReadBinaryInteger2DArray

  subroutine RestartBinaryFileAtByte (unit, byte, status)
    integer, intent(in) :: unit
    integer, intent(in) :: byte
    integer, intent(out) :: status
    character :: oneByte
    integer :: i

    read(unit, pos=byte, iostat=status) oneByte
      
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


  subroutine WriteBinaryString (unit, outputData, status)
    integer, intent(in) :: unit
    character (len=*), intent(in) :: outputData
    integer, intent(out) :: status
    write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryString

  subroutine WriteBinaryInteger (unit, outputData, status)
    integer, intent(in)  :: unit
    integer, intent(in)  :: outputData
    integer, intent(out) :: status
    write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryInteger

  subroutine WriteBinaryInteger1DArray (unit, outputData, status)
    integer, intent(in) :: unit
    integer, dimension(:), intent(in) :: outputData
    integer, intent(out) :: status
    write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryInteger1DArray

  subroutine WriteBinaryShortInteger1DArray (unit, outputData, status)
  integer, intent(in) :: unit
    integer(kind=selected_int_kind(4)), dimension(:), intent(in) :: outputData
    integer, intent(out) :: status
    write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryShortInteger1DArray

  subroutine WriteBinaryReal1DArray (unit, outputData, status)
    integer, intent(in) :: unit
    real(kind=4), dimension(:), intent(in) :: outputData
    integer, intent(out) :: status
    write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryReal1DArray
  
  subroutine WriteBinaryInteger2DArray (unit, outputData, status)
    integer, intent(in) :: unit
    integer, dimension(:,:), intent(in) :: outputData
    integer, intent(out) :: status
    write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryInteger2DArray

  subroutine WriteBinaryReal (unit, outputData, status)
    integer, intent(in) :: unit
    real(kind=4), intent(in) :: outputData
    integer, intent(out) :: status
    write (unit=unit,iostat=status) outputData
  end subroutine WriteBinaryReal

  subroutine WriteBinaryReal2DArray (unit, outputData, status)
    integer, intent(in) :: unit
    real(kind=4), dimension(:,:), intent(in) :: outputData
    real(kind=4), dimension(:,:), allocatable :: outputdata2 !Temp array is used to avoid stack overflow
    integer, intent(out) :: status

    allocate(outputdata2(size(outputData,1),size(outputData,2)))
    outputdata2=outputData
    write (unit=unit,iostat=status) outputData2
  end subroutine WriteBinaryReal2DArray


  subroutine WriteBinaryReal3DArray (unit, outputData, status)
    integer, intent(in) :: unit
    real(kind=4), dimension(:,:,:), intent(in) :: outputData
    real(kind=4), dimension(:,:,:), allocatable :: outputdata2 !Temp array is used to avoid stack overflow
    integer, intent(out) :: status

    allocate(outputdata2(size(outputData,1),size(outputData,2),size(outputData,3)))
    outputdata2=outputData
    write (unit=unit,iostat=status) outputData2
  end subroutine WriteBinaryReal3DArray


  function IsGood (num) result (res)
    real, intent(in) :: num
    logical :: res
    integer :: test
     res = .true.
     return
!    res = .not. isnan(num)
    
  end function IsGood

end module IOModule
