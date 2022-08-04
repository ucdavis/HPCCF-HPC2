C PSU-WOPWOP
C $Id: zbrent.f 3298 2013-10-18 02:12:03Z brentner $
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
      FUNCTION zbrent(func,x1,x2,tol,fa,fb)
      INTEGER ITMAX
      REAL zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      !fa=func(a)
      !fb=func(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
        print*, fa,'FA',fb,'FB'
        write (*,*) "ERROR: root must be bracketed for zbrent"
        stop
      end if
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm) 
        endif
        fb=func(b)
11    continue
      write (*,*) "ERROR: zbrent exceeding maximum iterations"
      stop
      zbrent=b
      return
      END
