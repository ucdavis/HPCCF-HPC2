! $Id: realfft.f90 569 2006-07-24 15:44:56Z chennes $
! $LastChangedDate: 2006-07-24 11:44:56 -0400 (Mon, 24 Jul 2006) $
subroutine realfft(data1, n, isig)
  implicit none
   integer isig, n
   real data1(n)
   integer i,i1,i2,i3,i4,n2p3
   real c1,c2,h1i,h1r,h2i,h2r,wis,wrs
   double precision theta, wi, wpi, wpr, wr, wtemp
   theta=3.141592653589793d0/dble(n/2)
   c1=0.5
   if(isig.eq.1) then 
      c2=-0.5
      call four1(data1, n/2, +1)
   else
      c2=0.5
      theta=-theta
   endif
   wpr=-2.0d0*sin(0.5d0*theta)**2
   wpi=sin(theta)
   wr=1.0d0+wpr
   wi=wpi
   n2p3=n+3
   do i=2, n/4
        i1=2*i-1
      i2=i1+1
      i3=n2p3-i2
      i4=i3+1
      wrs=sngl(wr)
      wis=sngl(wi)
      h1r= c1*(data1(i1)+data1(i3))
      h1i= c1*(data1(i2)-data1(i4))
      h2r=-c2*(data1(i2)+data1(i4))
      h2i= c2*(data1(i1)-data1(i3))
      data1(i1)= h1r+wrs*h2r-wis*h2i
      data1(i2)= h1i+wrs*h2i+wis*h2r
      data1(i3)= h1r-wrs*h2r+wis*h2i
      data1(i4)=-h1i+wrs*h2i+wis*h2r
      wtemp=wr
      wr=wr*wpr-wi   *wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
   enddo
   if(isig.eq.1) then
      h1r=data1(1)
      data1(1)=h1r+data1(2)
      data1(2)=h1r-data1(2)
   else
      h1r=data1(1)
      data1(1)=c1*(h1r+data1(2))
      data1(2)=c1*(h1r-data1(2))
      call four1(data1,n/2,-1)
   end if
return
end subroutine realfft
