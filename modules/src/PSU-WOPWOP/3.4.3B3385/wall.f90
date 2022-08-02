! $Id: wall.f90 3298 2013-10-18 02:12:03Z brentner $
! $LastChangedDate: 2013-10-17 22:12:03 -0400 (Thu, 17 Oct 2013) $
module wallObject
 use mathModule
 use MPIModule
 
 type wall
   private
   ! This vector defines which way the wall is facing
   type(vector)::normalVector, pointOnPlane
 end type wall

  integer::mg_NbWall, mg_wallIteration, mg_currentWallNb
  type(wall), dimension(:), allocatable::mg_wallArray
 
contains

  subroutine CreateWall(unitNumber, w)
    type(Wall), intent(inout)::w
    integer,    intent(in)   ::unitNumber

    real(kind=4), dimension(3)::normalVector, pointOnPlane
    
    NameList / WallIn / normalVector, pointOnPlane
    
    read(unitNumber, nml=WallIn)
    w%normalVector = vectorSetValue(normalVector)
    !
    ! w%normalVector should be a unit normal vector - normalize here KSB 12/2/2012
    !    
    w%normalVector = w%normalVector/vectorAbsolute(w%normalVector)
    w%pointOnPlane = vectorSetValue(pointOnPlane)
    
  
  end subroutine CreateWall
  !
  ! IMPORTANT this assumes the input vector is on the correct side of the wall
  function GetImagePosition(w, inputVector) result(outputVector)
    type(wall),   intent(in) ::w
    type(vector), intent(in) ::inputVector
    type(vector)             ::outputVector
    outputVector=inputVector+2.0*((w%pointOnPlane-inputVector)*w%normalVector)*w%normalVector
  end function GetImagePosition
  !
  ! This computes the velocity of the wall observer in the case of a uniformly translating
  ! observer. (V')  
  !
  function GetImageVelocity(w,inputVelocity) result(outputVelocity)
    type(wall),  intent(in) :: w
    type(vector),intent(in) :: inputVelocity
    type(vector)            :: outputVelocity
    outputVelocity = inputVelocity-2.0*(inputVelocity*w%normalVector)*w%normalVector
  end function GetImageVelocity

  subroutine BroadcastWall(w)
    type(wall), intent(inout)::w
    call BroadcastVector(w%normalVector)
    call BroadcastVector(w%pointOnPlane)
  end subroutine BroadcastWall
  
end module wallObject
