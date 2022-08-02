! $Id: debug.f90 3360 2016-02-01 04:58:14Z brentner $
module debugModule
  use MPIModule
  use strings
  integer::debugLevel
  integer:: BufferUnitNum

contains
  
  subroutine SetDebugLevel(input)
    integer::input
    debugLevel=input
    if (debugLevel.gt.1.and.IsMaster()) then
      write(*,*) 'Setting debug level to ', debugLevel
    end if    
  end subroutine SetDebugLevel

  subroutine Message(string1, string2, string3, string4, &
                     string5, string6, string7)
    character(len=*), intent(in)::string1
    character(len=*), intent(in), optional::string3, string4, string5, &
                                            string6, string7, string2
    if (debugLevel.gt.2.and.IsMaster()) then
       write (*,'(2A)') ' ', string1
       if (present(string2)) write (*,'(2A)') ' ', string2
       if (present(string3)) write (*,'(2A)') ' ', string3
       if (present(string4)) write (*,'(2A)') ' ', string4
       if (present(string5)) write (*,'(2A)') ' ', string5
       if (present(string6)) write (*,'(2A)') ' ', string6
       if (present(string7)) write (*,'(2A)') ' ', string7
    end if
  end subroutine Message
  
  subroutine Notice(string1, string2, string3, string4, &
                     string5, string6, string7)
    character(len=*), intent(in)::string1
    character(len=*), intent(in), optional::string3, string4, string5, &
                                            string6, string7, string2
    if (debugLevel.gt.2.and.IsMaster()) then
       write (*,'(2A)') '  NOTICE: ', string1
       if (present(string2)) write (*,'(2A)') '          ', string2
       if (present(string3)) write (*,'(2A)') '          ', string3
       if (present(string4)) write (*,'(2A)') '          ', string4
       if (present(string5)) write (*,'(2A)') '          ', string5
       if (present(string6)) write (*,'(2A)') '          ', string6
       if (present(string7)) write (*,'(2A)') '          ', string7
    end if
  end subroutine Notice  
    
  subroutine Warning(string1, string2, string3, string4, &
                     string5, string6, string7)
    character(len=*), intent(in)::string1
    character(len=*), intent(in), optional::string3, string4, string5, &
                                            string6, string7, string2
    if (debugLevel.gt.0.and.IsMaster()) then
       write (*,'(2A)') ' WARNING: ', string1
       if (present(string2)) write (*,'(2A)') '          ', string2
       if (present(string3)) write (*,'(2A)') '          ', string3
       if (present(string4)) write (*,'(2A)') '          ', string4
       if (present(string5)) write (*,'(2A)') '          ', string5
       if (present(string6)) write (*,'(2A)') '          ', string6
       if (present(string7)) write (*,'(2A)') '          ', string7
    end if
  end subroutine Warning
    
  subroutine Error(string1, string2, string3, string4, &
                   string5, string6, string7)
    character(len=*), intent(in)::string1
    character(len=*), intent(in), optional::string2, string3, string4, &
                                            string5, string6, string7
    if (debugLevel.gt.0) then
       write (*,*)
       write (*,'(2A)') ' ERROR: ', string1
       if (present(string2)) write (*,'(2A)') '        ', string2
       if (present(string3)) write (*,'(2A)') '        ', string3
       if (present(string4)) write (*,'(2A)') '        ', string4
       if (present(string5)) write (*,'(2A)') '        ', string5
       if (present(string6)) write (*,'(2A)') '        ', string6
       if (present(string7)) write (*,'(2A)') '        ', string7
    end if
    call KillAllProcesses()
  end subroutine Error
   


  subroutine Write1DBufferSize(title,i)
    integer::i
    character(len=*)::title
    if (debugLevel.lt.12) then
      return
    else
      write(GetBufferUnitNum(),*) trim(title)
      write(GetBufferUnitNum(),*) 'Size :', trim(integertostring(i))
    end if 
  end subroutine Write1DBufferSize
  
  subroutine Write2DBufferSize(title,i,j)
    integer::i,j
    character(len=*)::title
    if (debugLevel.lt.12) then
      return
    else
      write(GetBufferUnitNum(),*) trim(title)
      write(GetBufferUnitNum(),*) 'Size :', trim(integertostring(i)) ,' by ',&
                              trim(integertostring(j))
    end if
  end subroutine Write2DBufferSize

  subroutine Write3DBufferSize(title,i,j,k)
    integer::i,j,k
    character(len=*)::title
    if (debugLevel.lt.12) then
      return
    else
      write(GetBufferUnitNum(),*) trim(title)
      write(GetBufferUnitNum(),*) 'Size :', trim(integertostring(i)) ,' by ',&
                              trim(integertostring(j)), ' by ', &
                              trim(integertostring(k))
    end if
  end subroutine Write3DBufferSize
  
  subroutine Write4DBufferSize(title,i,j,k,n)
    integer::i,j,k,n
    character(len=*)::title
    if (debugLevel.lt.12) then
      return
    else
      write(GetBufferUnitNum(),*) trim(title)
      write(GetBufferUnitNum(),*) 'Size :', trim(integertostring(i)) ,' by ',&
                              trim(integertostring(j)), ' by ', &
                              trim(integertostring(k)), ' by ', &
                              trim(integertostring(n))
    end if
  end subroutine Write4DBufferSize
    
  subroutine SetBufferUnitNum(num)
    integer::num
    BufferUnitNum=num
  end subroutine SetBufferUnitNum

  function GetBufferUnitNum() result(num)
    integer::num
    num=BufferUnitNum
  end function GetBufferUnitNum
  
end module debugModule