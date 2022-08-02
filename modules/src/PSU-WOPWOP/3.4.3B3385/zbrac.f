C PSU-WOPWOP
C $Id: zbrac.f 3298 2013-10-18 02:12:03Z brentner $
! $LastChangedDate: 2013-10-17 22:12:03 -0400 (Thu, 17 Oct 2013) $
C
C This code has been developed through support from The Pennsylvania State 
C University and the following U.S. Government Contracts:
C 
C 9/1/00 - 8/31/01  NASA Cooperative Agreement NCC-1-406  
C (NASA Langley Research Center; H. E Jones, Technical Officer) 
C
C 2/1/01 - 1/31/06  NASA Ames Training Grant No. NGT-2-52275 
C (NASA Ames Research Center; Y. H. Yung, Technical Officer) 
C
C 5/15/2003 - 5/14/2004  NASA LaRC Purchase Order L-71744D 
C (NASA Langley Research Center; C. L. Burley, Technical Officer) 
C
C 3/1/2003 - 2/28/2006 NASA Langley Grant NAG-03025 
C (NASA Langley Research Center, D. P. Lockard, Technical Officer)
C
C
C Written by Guillaume Bres, Guillaume Perez, Leonard Lopes, Hsuan-Nien Chen, 
C Christopher Hennes, Rui Cheng and Benjamin Goldman. 
C Faculty advisor Dr. Kenneth S. Brentner.
C 
      SUBROUTINE zbrac(func,x1,x2,succes,f1,f2)
      INTEGER NTRY
      REAL x1,x2,func,FACTOR
      EXTERNAL func
      PARAMETER (FACTOR=1.6,NTRY=50)
      INTEGER j
      REAL f1,f2
      LOGICAL succes
      if(x1.eq.x2) then
        succes = .false.
        return
      end if
      f1=func(x1)
      f2=func(x2)
      succes=.true.
      do 11 j=1,NTRY
        if(f1*f2.lt.0.)return
        if(abs(f1).lt.abs(f2))then
          x1=x1+FACTOR*(x1-x2)
          f1=func(x1)
        else
          x2=x2+FACTOR*(x2-x1)
          f2=func(x2)
        endif
11    continue
      succes=.false.
      return
      END