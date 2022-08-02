module observerUtility
! PSU-WOPWOP
! $LastChangedDate: 2014-08-07 14:36:12 -0400 (Thu, 07 Aug 2014) $
! $Id: obsUtility.f90 3321 2014-08-07 18:36:12Z brentner $
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



  use COBObject
  use mathModule
  use constantsModule
  use MPIModule
  use debugModule
  use IOModule
  
  implicit none

  ! An enumeration of possible function values:
  integer, parameter :: THICKNESS_OASPL=1,&
                        LOADING_OASPL=2,&
                        TOTAL_OASPL=3,&
                        THICKNESS_OASPLDBA=4,&
                        LOADING_OASPLDBA=5,&
                        TOTAL_OASPLDBA=6,&
                        NUM_FUNCTIONS=6

  ! The names of the functions enumerated above. Note that the order of this
  ! list must match the order of the list above.
  character (len=50), dimension (NUM_FUNCTIONS), parameter :: functionNames = &
                    (/ "Thickness OASPL   ",&
                       "Loading OASPL     ",&
                       "Total OASPL       ",&
                       "Thickness OASPLdBA",&
                       "Loading OASPLdBA  ",&
                       "Total OASPLdBA    "/)
  
  ! A type to store structured function data
  type StructuredStorage
    real, dimension(:,:,:),   pointer::A   ! For single-zoned data
    real, dimension(:,:,:,:), pointer::B   ! For multi-zoned data
  end type
  
contains

!**
!SUBROUTINE WriteFunctionFile(path, time, function1, function2, function3,
!                             function4, function5, function6)
!   Needs commenting.
!ARGUMENTS:
! - path:the path the file will be written to 
! - time:the observer time array
! - function(1-6): Thickness, loading, total noise,
!                  Thickness, loading, total phase, respectively
!                  (Phase functions are optional)
!PROCESS:
!CALLS:
! - createBinaryFile       (linux/windowsIO.f90)
! - writeBinaryInteger     (linux/windowsIO.f90)
! - writeBinaryReal        (linux/windowsIO.f90)
! - writeBinaryReal1DArray (linux/windowsIO.f90)
!**
subroutine WriteFunctionFile (path, time, function1, function2, function3, &
                                          function4, function5, function6)
  implicit none                                          
  character(len=*), intent(in)::path
  real(kind=4):: temp
  real, dimension(:)::time
  real, dimension(:,:,:), pointer, optional:: function1, function2, function3, &
                                              function4, function5, function6
  logical :: f1 = .false., f2 = .false., f3 = .false., f4 = .false., &
             f5 = .false., f6 = .false.
  type(StructuredStorage), dimension (:), pointer :: validFunctions
  
  integer :: stat, iMax, jMax, kMax, nFun, f, i, j, k
  nFun = 0 !this is counting how many functions are present,
           !from Function 1 through Function 6
  if (present (function1)) then
    if (associated (function1)) then
      nFun = nFun+1
      f1 = .true.
      iMax = size(function1,1)
      jMax = size(function1,2)
      kMax = size(function1,3)
    end if
  end if
  if (present (function2)) then
    if (associated (function2)) then
      nFun = nFun+1
      f2 = .true.
      iMax = size(function2,1)
      jMax = size(function2,2)
      kMax = size(function2,3)
    end if
  end if
  if (present (function3)) then
    if (associated (function3)) then
      nFun = nFun+1
      f3 = .true.
      iMax = size(function3,1)
      jMax = size(function3,2)
      kMax = size(function3,3)
    end if
  end if
  if (present (function4)) then
    if (associated (function4)) then
      nFun = nFun+1
      f4 = .true.
      iMax = size(function4,1)
      jMax = size(function4,2)
      kMax = size(function4,3)
    end if
  end if
  if (present (function5)) then
    if (associated (function5)) then
      nFun = nFun+1
      f5 = .true.
      iMax = size(function5,1)
      jMax = size(function5,2)
      kMax = size(function5,3)
    end if
  end if
  if (present (function6)) then
    if (associated (function6)) then
      nFun = nFun+1
      f6 = .true.
      iMax = size(function6,1)
      jMax = size(function6,2)
      kMax = size(function6,3)
    end if
  end if

  if (nFun == 0) return

  ! Now we jump through some hoops to make unstructured output work.
  ! All functions must be output in the same write statement, which is difficult
  ! in their present format. Set it up so that a loop can be used:
  allocate (validFunctions(nFun))
  f= 0
  if (f1) then
    f=f+1
    validFunctions(f)%A => function1
  end if
  if (f2) then
    f=f+1
    validFunctions(f)%A => function2
  end if
  if (f3) then
    f= f+1
    validFunctions(f)%A => function3
  end if
  if (f4) then
    f = f+1
    validFunctions(f)%A => function4
  end if
  if (f5) then
    f =f+1
    validFunctions(f)%A => function5
  end if
  if (f6) then
    f = f+1
    validFunctions(f)%A => function6
  end if
  ! So, we now have an array called validFunctions which can be looped over to
  ! provide the functions to be output

 if (.not.ASCIIOutputFlag) then
    call CreateBinaryFile (116, trim(path), .false., stat)
    if (stat/=0) then
      call Error('Could not open the file '//trim(path))
      stop
    end if
    call WriteBinaryInteger (116, iMax,   stat)
    call WriteBinaryInteger (116, jMax,   stat)
    call WriteBinaryInteger (116, kMax,   stat)
    call WriteBinaryInteger (116, nFun+1, stat)
      
    !The following nested do-loop was created to write out the time array
    !into this binary file.  This was creating the problem opening the
    !binary pressure.fn file in Fieldview.
    do k=1, kMax
      do j=1, jMax
        do i=1, iMax
          call WriteBinaryReal (116, temp, stat)
          time(k) = temp
        end do
      end do
    end do 
    
    !This is what was removed.  This was not writing the time array correctly
    !to the binary file.
    !call WriteBinaryReal1DArray (116, time, stat)
    
    do f=1, nFun
      do k=1, kMax
        do j=1, jMax
          do i=1, iMax
            call WriteBinaryReal(116, temp, stat)
            validFunctions(f)%A(i,j,k) = temp
            if (stat/=0) then
              call Error('Could not write to the file '//trim(path))
              stop
            end if
          end do
        end do
      end do
    end do
    close(116)
    
  else
    open (116, file=trim(path), form="formatted", iostat=stat)
    if (stat /= 0) then
      write (*,*) "ERROR: The file ", trim(path), " could not be written to."
      return
    end if
    write(116,'(4I12)') iMax, jMax, kMax, nFun+1
    write(116,'(e15.8)') ((( time(k), i=1, iMax), j=1, jMax), k=1, kMax),  &
               ((((validFunctions(f)%A(i,j,k),i=1, iMax), j=1, jMax), k=1, kMax),f=1,nFun)
    close(116)
  end if
    
end subroutine WriteFunctionFile



!**
!SUBROUTINE writeOutSingleFile()
!   This routine takes up to seven arrays, and writes an ASCII file
!   with their values.  If any of the y arrays are not the same
!   length as the x array, it warns the user.  If there is no data
!   available when this routine is called, it warns the user also.
!ARGUMENTS:
! - xArray: usually the observer time array
! - yArray(1-6): thickness, loading, total,
!                thickness phase, loading phase, total phase
! - yArray(1-6)Name: The associated names
!CALLS:
! - warnings and errors
!**
subroutine writeOutSingleFile(path, xArray,  xAxisTitle, &
            yArray1, yArray1Name,yArray2, yArray2Name,&
            yArray3, yArray3Name,yArray4, yArray4Name,&
            yArray5, yArray5Name,yArray6, yArray6Name)
  implicit none                                    
  character(len=*), intent(in)::path
  character(len=*), intent(in), optional::xAxisTitle
  real, dimension(:), intent(in)::xArray
  real, dimension(:), optional,intent(in):: yArray1, yArray2, yArray3, &
                                          yArray4, yarray5, yArray6
  
  character (len=*), optional ::  yArray1Name, yArray2Name, yArray3Name, &
                                  yArray4Name, yArray5Name, yArray6Name 
                                                           
  integer::i, start, unitNumber
  logical:: hasf1, hasF2, hasF3, hasF4, hasF5, hasF6

  unitNumber=GetStreamNumber()
              
  start = lbound(xArray,1)
  
  hasf1  = .false.;hasF2  = .false.;hasF3  = .false.;hasF4  = .false.
  hasF5  = .false.;hasF6  = .false.
  if (present(yArray1)) then 
    if (size(yArray1).gt.0) then
      if (size(yArray1) == size(xArray)) then
        hasf1 = .true.
      else
        call Warning ("In WriteOutSingleFile, The size of the " &
                      //trim(xAxisTitle)//" array does not match ",&
                      "the size of the "//trim(yArray1Name)//" array.")
      end if
    end if
  end if
  
  if (present(yArray2)) then
    if (size(yArray2).gt.0)  then
      if (size(yArray2) == size(xArray)) then
        hasF2 = .true.
      else
        call Warning ("The size of the"//trim(xAxisTitle)//" array does not match",&
                      "the size of the "//trim(yArray2Name)//" array.")
      end if
    end if
  end if
    
  if (present(yArray3)) then
    if (size(yArray3).gt.0)  then
      if (size(yArray3) == size(xArray))then
        hasF3 = .true.
      else
        call Warning ("The size of the"//trim(xAxisTitle)//" array does not match",&
                      "the size of the "//trim(yArray3Name)//" array.")
      end if
    end if
  end if
  if (present(yArray4)) then
    if (size(yArray4).gt.0)  then
      if (size(yArray4) == size(xArray)) then
        hasF4 = .true.
      else
        call Warning ("The size of the"//trim(xAxisTitle)//" array does not match",&
                      "the size of the "//trim(yArray4Name)//" array.")
      end if
    end if
  end if
  if (present(yArray5)) then
    if (size(yArray5).gt.0)  then
      if (size(yArray5) == size(xArray))then
        hasF5 = .true.
      else
        call Warning ("The size of the"//trim(xAxisTitle)//" array does not match",&
                      "the size of the "//trim(yArray5Name)//" array.")
      end if
    end if
  end if
  if (present(yArray6)) then
    if (size(yArray6).gt.0)  then
      if (size(yArray6) == size(xArray))then
        hasF6 = .true.
      else
        call Warning ("The size of the"//trim(xAxisTitle)//" array does not match",&
                      "the size of the "//trim(yArray6Name)//" array.")
      end if
    end if
  end if  

  if (hasf1 .or. hasF2 .or. hasF3 .or. hasF4 .or. hasF5 .or. hasF6) then
      open(unitNumber,file=path,status="REPLACE",form='FORMATTED')
      
      write(unitNumber,*)'TITLE="PSU-WOPWOP"'      
      
      write(unitNumber, '(A,$)') 'VARIABLES= "'//trim(xAxisTitle)//'"'
      if(hasf1) then
        write (unitNumber, '(A,$)') ', "'//trim(yArray1Name)//'"'
      end if
      if(hasF2) then
        write (unitNumber, '(A,$)') ', "'//trim(yArray2Name)//'"'
      end if
      if(hasF3) then
        write (unitNumber, '(A,$)') ', "'//trim(yArray3Name)//'"'
      end if
      if(hasF4) then
        write (unitNumber, '(A,$)') ', "'//trim(yArray4Name)//'"'
      end if
      if(hasF5) then
        write (unitNumber, '(A,$)') ', "'//trim(yArray5Name)//'"'
      end if
      if(hasF6) then
        write (unitNumber, '(A,$)') ', "'//trim(yArray6Name)//'"'
      end if
       
            
      write(unitNumber, '(A,$)')
      
      !This is the title
        write (unitNumber, '(A,$)') 'TEXT X=100 Y=95 AN=HEADRIGHT T="'//trim(xAxisTitle)//' vs \\n'
      if(hasf1) then
        write (unitNumber, '(A,$)') trim(yArray1Name)//'\\n'
      end if
      if(hasF2) then
        write (unitNumber, '(A,$)') trim(yArray2Name)//'\\n'
      end if
      if(hasF3) then 
        write (unitNumber, '(A,$)') trim(yArray3Name)//'\\n'
      end if
      if(hasF4) then
        write (unitNumber, '(A,$)') trim(yArray4Name)//'\\n'
      end if
      if(hasF5) then
        write (unitNumber, '(A,$)') trim(yArray5Name)//'\\n'
      end if
      if(hasF6) then
        write (unitNumber, '(A,$)') trim(yArray6Name)//'\\n'
      end if

          
      write(unitNumber,*)'"'
      write (unitNumber, '(A,$)')

      
      do i=start, size(xArray)
        write (unitNumber,'(e15.8," ",$)') xArray(i)
        if (hasf1) write (unitNumber,'(e15.8," ",$)') yArray1(i)
        if (hasF2) write (unitNumber,'(e15.8," ",$)') yArray2(i)
        if (hasF3) write (unitNumber,'(e15.8," ",$)') yArray3(i)
        if (hasF4) write (unitNumber,'(e15.8," ",$)') yArray4(i)
        if (hasF5) write (unitNumber,'(e15.8," ",$)') yArray5(i)
        if (hasF6) write (unitNumber,'(e15.8," ",$)') yArray6(i)
                       
        write (unitNumber,*)
      end do  
    close(unitNumber)
  else
    if (debugLevel >= 1) write (*,*) "WARNING: WriteOutSingleFile called with no data when writing the file"
    if (debugLevel >= 1) write (*,'(2A)') "          ", trim(path)
  end if
end subroutine writeOutSingleFile

  
subroutine writeOutTecplotData(path,xArray,xAxisTitle,f,fTitles)
  implicit none                                    
  character(len=*), intent(in)::path
  character(len=*), intent(in), optional::xAxisTitle
  real, dimension(:), intent(in)::xArray
  real, dimension(:,:),intent(in):: f  
  character(len=*),dimension(:),intent(in) ::  fTitles                                                           
  integer::i,j, unitNumber
  unitNumber=GetStreamNumber()  
    open(unitNumber,file=path,status="REPLACE",form='FORMATTED')      
      write(unitNumber,*)'TITLE="PSU-WOPWOP"'     
      write(unitNumber, '(A,$)') 'VARIABLES= "'//trim(xAxisTitle)//'"'
      do i=1,size(f,1)
        write (unitNumber, '(A,$)') ', "'//trim(fTitles(i))//'"'
      end do               
      write(unitNumber, '(A,$)')      
      !This is the title
      write (unitNumber, '(A,$)') 'TEXT X=100 Y=95 AN=HEADRIGHT T="'//trim(xAxisTitle)//' vs \\n'
      do i=1,size(f,1)
        write (unitNumber, '(A,$)') trim(fTitles(i))//'\\n'
      end do                
      write(unitNumber,*)'"'
      write (unitNumber, '(A,$)')      
      do i=1, size(xArray)
        write (unitNumber,'(e15.8," ",$)') xArray(i)
        do j=1,size(f,1)
          write (unitNumber,'(e15.8," ",$)') f(j,i)
        end do
        write (unitNumber,*)
      end do  
    close(unitNumber)
end subroutine writeOutTecplotData  
 
subroutine CreateWavFile(pressure, timeHistory, path)
  real, dimension(:), intent(in)::pressure, timeHistory
  character(len=*), intent(in)::path
  integer::nPoints,chunkSize,subChunk1Size,subChunk2Size
  real::period, vref
  real, dimension(:), allocatable::tVal, val
  integer::nout, iv, i0, i,unit, stat
  integer(kind=selected_int_kind(4)), allocatable:: ival(:)
  integer(kind=selected_int_kind(4)),dimension(2)::array1,array2
  !array1 = (wFormatTag,nChannels)
  !array2 = (nBlockAlign,nBitsPerSample)
  integer :: nSamplesPerSec=44100, nAvgBytesPerSec=88200
  real::vmin, vmax, vave, v, t, vref1  
  array1(1) = 1; array1(2) = 1
  array2(1) = 2; array2(2) = 16
  unit = GetStreamNumber()
  nPoints=size(timeHistory)
  if(.not.(nPoints==size(timeHistory)) .and. debugLevel >=1 ) then
    print*, 'WARNING: Points in time history do not match points in acoustic history.  '
    print*, '         Stopping wav file generation.                                    '
    print*, '         Cannot create ', path//'.wav '
    return
  end if
  allocate(tVal(nPoints))
  period = timeHistory(nPoints)-timeHistory(1)
  vref=0.0
  ! period is the length of actual time the signal takes
  ! vref is the scale factor, if vref=0, then scale by max values of this signal
  nout=period*nSamplesPerSec
  allocate(ival(nout))
  allocate(val(nout))
  ! Determine the min and max values of the data.  Also set up to do the 
  ! interpolation needed to get the sample rate right.
  vmin = 10e20
  vmax = -vmin
  do i=1,nPoints
    vmin = min(vmin,pressure(i))
    vmax = max(vmax,pressure(i))
    tval(i) = (timeHistory(i)-timeHistory(1))/(timeHistory(nPoints)-timeHistory(1))
  end do
  vave = .5*(vmax+vmin)
  if( vref == 0. ) then
     vref1 = .5*(vmax-vmin)
  else
    vref1 = vref
  end if
  !  Now we need to interpolate to get the data at the prescibed sample rate 
  !  and write it out using a C function.
  i0=1
  do i=1,nout
    t = float(i-1)/float(nout-1)
    ! Use linear interpolation
    call hunt(tval,nPoints,t,i0)
    if(i==nout)then
      v = pressure(i0)
      else
    v=pressure(i0)+(t-tval(i0))*(pressure(i0+1)-pressure(i0))/(tval(i0+1)-tval(i0))
    end if
    ! 16 bit data is stored as integers ranging from -32768 to +32767, 0 being no sound
    ! Normalize the data and clip if it exceeds the range -32768<val<32767.
    iv = (v-vave)*32768/vref1
    ival(i) = iv
    if(iv.gt.32767) then
      ival(i)=32767
    else if(iv.lt.-32768) then
      ival(i)=-32768
    end if
  end do
  chunkSize = nout*2+36
  subChunk1Size = 16
  subChunk2Size = 2*nout
  call CreateBinaryFile    (unit, trim(path)//'.wav', .false., stat)
    call WriteBinaryString (unit,'RIFF',            stat)
    call WriteBinaryInteger(unit,chunkSize,         stat)
    call WriteBinaryString (unit,'WAVEfmt ',        stat)
    call WriteBinaryInteger(unit,subChunk1Size,     stat)
    call WriteBinaryShortInteger1DArray(unit,array1,stat)
    call WriteBinaryInteger(unit,nSamplesPerSec,    stat)        
    call WriteBinaryInteger(unit,nAvgBytesPerSec,   stat)    
    call WriteBinaryShortInteger1DArray(unit,array2,stat)
    call WriteBinaryString (unit, 'data',           stat)
    call WriteBinaryInteger(unit,subChunk2Size,     stat)
    call WriteBinaryShortInteger1DArray(unit,ival,  stat)
  call CloseBinaryFile(unit)    
  deallocate(ival)  
end subroutine CreateWavFile

end module observerUtility
