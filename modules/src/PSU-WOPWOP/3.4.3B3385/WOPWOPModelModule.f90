! $LastChangedDate: 2011-11-01 15:24:38 -0400 (Tue, 01 Nov 2011) $
! $Id: WOPWOPModelModule.f90 3169 2011-11-01 19:24:38Z bzg5066 $

module modelModule
  use debugModule


  type model
    integer::nbPoints
  end type model
  
  type(model)::mg_Cyl


contains
  
subroutine CreateTELoading(loadingArray, L, k, omega, nKey, deltaT, incVel)
  real(kind=8), dimension(:,:)::loadingArray
  real(kind=8)::L,k,omega, deltaT
  integer::nKey
  real(kind=8), dimension(3)::incVel
  if(size(loadingArray,2).gt.nKey) then
  end if
  L=L
  k=k
  omega=omega
  deltaT=deltaT
  incVel=incVel
  call Error('This routine should never be called', &
             'You need to make with LGMAP option')
  stop
end subroutine CreateTELoading

subroutine CreateCylLoading(temp2DArray, Cl, Cd, Ck, St,& 
                                p, e, nKey, deltaT, incidentVelocity) 
  real(kind=8), dimension(:,:)::temp2DArray
  real(kind=8)::Cl, Cd, Ck, St, p, e, deltaT
  integer::nKey
  real(kind=8), dimension(3)::incidentVelocity
  incidentVelocity=incidentVelocity
  temp2DArray=temp2DArray
  Cl=Cl
  Cd=Cd
  Ck=Ck
  St=St
  p=p
  e=e
  nKey=nKey
  deltaT=deltaT
  call Error('This routine should never be called', &
             'You need to make with LGMAP option')
  stop
end subroutine CreateCylLoading
  
subroutine CreateMemberGlobalCylinder(coord, diameter, correlationLength, &
                                      startPoint, endPoint, incVel)
  real(kind=8), dimension(:,:), pointer::coord
  real(kind=8), dimension(3)::startPoint, endPoint, incVel
  real(kind=8)::diameter, correlationLength
  coord=coord
  diameter=diameter
  correlationLength=correlationLength
  startPoint=startPoint
  endPoint=endPoint
  incVel=incVel
  call Error('This routine should never be called', &
             'You need to make with LGMAP option')
  stop
end subroutine CreateMemberGlobalCylinder

end module modelModule
