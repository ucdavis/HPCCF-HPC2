C PSU-WOPWOP
C $Id: hunt.f 3298 2013-10-18 02:12:03Z brentner $
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
      SUBROUTINE HUNT(XX,N,X,JLO)
C ----- Added by Chris Hennes, 6/1/2004 -----
      REAL X
      INTEGER N, INC, JLO, JHI, JM
      REAL XX
      DIMENSION XX(N)
C -------------------------------------------
      LOGICAL ASCND
      ASCND=XX(N).GT.XX(1)
      IF(JLO.LE.0.OR.JLO.GT.N)THEN
        JLO=0
        JHI=N+1
        GO TO 3
      ENDIF
      INC=1
      IF(X.GE.XX(JLO).EQV.ASCND)THEN
1       JHI=JLO+INC
        IF(JHI.GT.N)THEN
          JHI=N+1
        ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
          JLO=JHI
          INC=INC+INC
          GO TO 1
        ENDIF
      ELSE
        JHI=JLO
2       JLO=JHI-INC
        IF(JLO.LT.1)THEN
          JLO=0
        ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
          JHI=JLO
          INC=INC+INC
          GO TO 2
        ENDIF
      ENDIF
3     IF(JHI-JLO.EQ.1)RETURN
      JM=(JHI+JLO)/2
      IF(X.GT.XX(JM).EQV.ASCND)THEN
        JLO=JM
      ELSE
        JHI=JM
      ENDIF
      GO TO 3
      END
