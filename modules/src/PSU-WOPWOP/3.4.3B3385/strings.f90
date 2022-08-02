module strings
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
contains

    subroutine strToUpper(string,len)
        character(len=*)::string
        integer, intent(in), optional::len
        integer::i, rLen

        if (present(len)) then
        rLen = len
        else
        rLen = len_trim(string)
        end if
  
        do i=1, rLen
        if(LGE(string(i:i),"a").and.LGE("z",string(i:i))) then
            string(i:i) = char(iand(ichar(string(i:i)),223))
        end if
        end do

    end subroutine strToUpper
    
    function RealToString(input) result(output)
        real, intent(in)::input
        character(len=15)::output
        write(output, '(G15.7)') input
        output = adjustl(output)
    end function RealToString
  
    function IntegerToString(input) result(output)
        integer, intent(in)::input
        character(len=30)::output
        write(output,*) input
        output = adjustl(output)
    end function IntegerToString

end module strings