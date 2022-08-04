! $LastChangedDate: 2007-03-26 09:45:26 -0400 (Mon, 26 Mar 2007) $
module CylinderObject
  use frequencyAnalysisModule
  use constantsModule

  type cylinder
    real::diameter, length
    real, dimension(3)::lengthVector, normal, drag
    integer::nbPoints
  end type 

contains

subroutine CylinderLoading(cyl, loadingArray, Cl, Cd, Ck, St, p, e, &
                       nKey, deltaT, velocityVector)
  type(cylinder)::cyl
  real, dimension(:,:)::loadingArray
  real, intent(in)::Cl, Cd, Ck, St, p, e, deltaT
  integer::nKey
  real, dimension(3)::velocityVector
  
  real:: inVel, phaseOffset
  integer::nbHarmonics, i, j
  real, dimension(:), allocatable::drag, lift, phase, liftTime, dragTime
  real, dimension(3)::dragVector, normalVector
  
  loadingArray=0.0
  call CheckCylinderConstants(nKey, deltaT, Cl, Cd, Ck, St, p, e)
  inVel  = absolute(velocityVector) - &
           abs(dot_product(cyl%lengthVector/cyl%length, velocityVector))
  call random_seed()
  nbHarmonics=int(nKey/2.0)
  allocate(drag (nbHarmonics))
  allocate(lift (nbHarmonics))
  allocate(phase(nbHarmonics))
  allocate(liftTime(nKey))
  allocate(dragTime(nKey))
  dragVector = velocityVector/absolute(velocityVector)
  if (dragVector(1).ne.cyl%drag(1).and. &
      dragVector(2).ne.cyl%drag(2).and. &
      dragVector(3).ne.cyl%drag(3)) then
    call Warning('Loading cylinder model drag vector does not match ',&
                 'geometry drag vector')
  end if
  call calculateNormalVector(normalVector, dragVector, cyl%lengthVector)
  if (normalVector(1).ne.cyl%normal(1).and. &
      normalVector(2).ne.cyl%normal(2).and. &
      normalVector(3).ne.cyl%normal(3)) then
    call Warning('Loading cylinder model normal vector does not match ',&
                 'geometry normal vector')
  end if
  call createAmplitude(drag, lift, phase, nbHarmonics, St, p, e, &
                       deltaT, cyl%diameter, inVel)
  liftTime=0.0
  dragTime=0.0
  phaseOffset=0.0 !pi*sin((cyl%length*j/cyl%nbPoints)*(2.0*pi*kappa)/cyl%diameter)
!  call ifft(lift, phase+phaseOffset,        liftTime)
!  call ifft(drag, phase+phaseOffset+pi/2.0, dragTime) 
  call NormalizeForce(dragTime, liftTime)
  dragTime = dragTime*(.5*Cd*Ck*rho*cyl%diameter*(inVel**2))
  liftTime = liftTime*(.5*Cl*Ck*rho*cyl%diameter*(inVel**2))
  do i=1, nKey
    loadingArray(:,i)=liftTime(i)*normalVector + dragTime(i)*dragVector 
  end do
  deallocate(drag, lift, phase, liftTime, dragTime)
end subroutine CylinderLoading

subroutine CreateCylinder(cyl, coord, diameter, correlationLength, &
                          startPoint, endPoint, incVel)
  type(cylinder)::cyl
  real, dimension(:,:), pointer::coord
  real, dimension(3)::startPoint, endPoint, incVel
  real::diameter, correlationLength
  real::rem
  integer::nbSegments, i
  real, dimension(3)::dS
  cyl%diameter     = diameter
  cyl%drag         = incVel/absolute(incVel)
  cyl%lengthVector = endPoint - startPoint
  cyl%length       = absolute(cyl%lengthVector)
  call calculateNormalVector(cyl%normal, cyl%drag, cyl%lengthVector)
  rem = modulo(cyl%length, correlationLength*diameter)
  if(rem < (correlationLength*diameter/2)) then
    nbSegments=int(cyl%length/(correlationLength*diameter))
  else
    nbSegments=int(cyl%length/(correlationLength*diameter))+1
  end if
  if(nbSegments==0) nbSegments=1
  allocate(coord(3, nbSegments+2))
  dS=(cyl%lengthVector)/nbSegments
  do i=1, nbSegments+1
     coord(:,i) = startPoint+(i-1)*dS
  end do
  cyl%nbPoints=nbSegments+1
  coord(:,cyl%nbPoints+1) = cyl%normal  
end subroutine CreateCylinder

subroutine createAmplitude(dragAmplitude, liftAmplitude, phase,  &
                           nbHarmonics, St, p, e, deltaT, d, incVel)
  real, dimension(:), intent(inout)::dragAmplitude, liftAmplitude, phase
  integer, intent(in)::nbHarmonics
  real::St, p, e, deltaT, d, incVel
  
  integer::i
  real::S, f, A, B 
  real, dimension(:), allocatable::FS
  
  liftAmplitude = 0.0
  dragAmplitude = 0.0
  phase         = 0.0
  if(incVel.lt.(1e-3*c)) return
  B=-(St**p)*(e*(1.0-p)-1.0)/(e-1.0)
  allocate(FS(nbHarmonics))
  do i=1, nbHarmonics
    f=i/deltaT
    S=f*d/incVel
    FS(i)=((( S )**(e-1))*((B+( S )**p)**(-e)))
  end do
  call ComputeA(FS, d/(deltaT*incVel), A)
  call random_seed()
  do i=1, nbHarmonics
     f=i/deltaT
     S=f*d/incVel
     dragAmplitude(i)=(A*(( S )**(e-1))*((B+( S )**p)**(-e)))
     liftAmplitude(i)=(A*((2*S)**(e-1))*((B+(2*S)**p)**(-e)))
     call random_number(phase(i))
  end do
  phase=2.0*pi*phase-pi
  deallocate(FS)
end subroutine createAmplitude

subroutine NormalizeForce(input1, input2)
  real, dimension(:), intent(inout)::input1, input2
  integer::i
  real::sum
  sum=0.0
  do i=1, size(input1)
    sum=sum+((input1(i))**2)
  end do
  sum=sqrt(sum/size(input1))
  input1=input1/sum
  input2=input2/sum
end subroutine NormalizeForce

subroutine calculateNormalVector(normal, drag, length)
  real, dimension(3)::normal, drag, length
  normal = crossProduct(length, drag)
  if(absolute(normal).eq.0.0) then
    normal=0.0
  else
    normal = normal/absolute(normal)
  end if
end subroutine calculateNormalVector

subroutine ComputeA(F,delta,A)
  real, dimension(:)::F
  real::A,delta, sum
  integer::i
  sum=0.0
  do i=1, size(F)
    sum=sum+F(i)
  end do
  sum=sum*delta
  A=sum**(-1.0)  
end subroutine ComputeA

function absolute(V) result(res)
  real, dimension(3)::V
  real::res
  res=sqrt(V(1)**2 + V(2)**2 + V(3)**2)  
end function absolute

function crossProduct(V1,V2) result(res)
  real, dimension(3),intent(in)::V1,V2
  real, dimension(3)::res
  res(1)=V1(2)*V2(3)-V1(3)*V2(2)
  res(2)=V1(3)*V2(1)-V1(1)*V2(3)
  res(3)=V1(1)*V2(2)-V1(2)*V2(1)    
end function crossProduct

subroutine CheckCylinderConstants(nKey, deltaT, Cl, Cd, Ck, St, p, e)
  real   :: deltaT, Cl, Cd, Ck, St, p, e
  integer:: nKey
  if (nKey <= 0 ) then
    call Error('Invalid number of key locations in aperiodic loading file.',&
               'Stopping')
    stop
  end if
  if (deltaT <=0.0) then
    call Error('Invalid deltaT in cylinder loading file.', &
               'DeltaT read in is: '//RealToString(deltaT),   &
               'Stopping')
    stop
  end if
  if (Cl <=0.0) then
    call Error('Invalid lift coefficient in cylinder loading file.', &
               'Lift coefficient read in is: '//RealToString(Cl),   &
               'Stopping')
    stop
  end if
  if (Cd <=0.0) then
    call Error('Invalid drag coefficient in cylinder loading file.', &
               'Drag coefficient read in is: '//RealToString(Cd),    &
               'Stopping')
    stop
  end if    
  if (Ck <=0.0) then
    call Error('Invalid turbulent coefficient in cylinder loading file.', &
               'Turbulent coefficient read in is: '//RealToString(Ck),   &
               'Stopping')
    stop
  end if
  if (St <=0.0) then
    call Error('Invalid Strouhal number in cylinder loading file.', &
               'Strouhal number read in is: '//RealToString(St),   &
               'Stopping')
    stop
  end if
  if ((p <= 0.0) .or. (e <= 0.0)) then
    call Error('Invalid spectrum coefficients in cylinder loading file.', &
               'Spectrum coefficients read in are: '//RealToString(p)//' '&
                                                    //RealToString(e),    &
               'Stopping')
    stop
  end if
end subroutine CheckCylinderConstants
 
end module CylinderObject
