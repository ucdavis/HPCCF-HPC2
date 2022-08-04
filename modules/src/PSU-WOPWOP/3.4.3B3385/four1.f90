! $LastChangedDate: 2006-07-24 11:44:56 -0400 (Mon, 24 Jul 2006) $
! $Id: four1.f90 569 2006-07-24 15:44:56Z chennes $
SUBROUTINE four1(data1,nn,isign)
  implicit none
  INTEGER isign,nn
  real data1(2*nn)
  INTEGER i,istep,j,m,mmax,n
  real tempi,tempr
  double precision theta,wi,wpi,wpr,wr,wtemp
  n=2*nn
  j=1
  do i=1,n,2
     if(j.gt.i)then
       tempr=data1(j)
       tempi=data1(j+1)
       data1(j  )=data1(i  )
       data1(j+1)=data1(i+1)
       data1(i  )=tempr
       data1(i+1)=tempi
     endif
     m=nn
     do while ((m.ge.2).and.(j.gt.m))
        j=j-m
        m=m/2
     end do
        j=j+m
  end do
  mmax=2
  do while (n.gt.mmax)
     istep=2*mmax
     theta=6.28318530717959d0/(isign*mmax)
     wpr=-2.d0*sin(0.5d0*theta)**2
     wpi=sin(theta)
     wr=1.d0
     wi=0.d0
     do m=1,mmax,2
        do i=m,n,istep
           j=i+mmax
           tempr=sngl(wr)*data1(j  )-sngl(wi)*data1(j+1)
           tempi=sngl(wr)*data1(j+1)+sngl(wi)*data1(j)
           data1(j  )=data1(i  )-tempr
           data1(j+1)=data1(i+1)-tempi
           data1(i  )=data1(i  )+tempr
           data1(i+1)=data1(i+1)+tempi
        end do
        wtemp=wr
        wr=wr*wpr-wi   *wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
     end do
     mmax=istep
  end do
  return
END SUBROUTINE four1
