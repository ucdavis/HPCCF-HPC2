module MathModule
! PSU-WOPWOP
! $LastChangedDate: 2014-01-28 18:13:47 -0500 (Tue, 28 Jan 2014) $
! $Id: mathModule.f90 3306 2014-01-28 23:13:47Z brentner $
!
! This code has been developed through support from The Pennsylvania State 
! University and the following U.S. Government Contracts:
! 
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


!**************************************************************************************
!    MODULE OF THE MATRIX AND VECTOR DEFINITIONS                           *
!**************************************************************************************  
!
!This module defines matrix and vector, and provides the following operations:
!  - assignement for matrix M and vector V
!  - add: M1+M2, V1+V2 
!  - product: M1*M2, real*M, real*V, M*V, dot product V1*V2, cross product
!V1.cross.V2
!  - substraction: M1-M2, V1-V2
!  - Negation: -M, -V
!  - Norm of a vector: abs(V)
!  - Transpose: Transpose(V)
!REMARKS:- Matrix and vector can be visualized with the subroutines
!PrintMatrix, PrintMatrix4, PrintVector and PrintVector

!LAST MODIFIED: G.BRES 24/08/01
!**************************************************************************************

  use MPIModule
  use constantsModule
  implicit none

  type matrix
     real, dimension(3,3)::A
  End type matrix

  type matrix4
     real, dimension(4,4)::A
  End type matrix4

  Type vector
     real, dimension(3)::A
  End type vector
  
  Type vectorKind8
     real(kind=8), dimension(3)::A
  End type vectorKind8

  Type vector4
     real, dimension(4)::A
  End type vector4

  Interface operator (+)
     module procedure matrixAdd, vectorAdd
  end interface
  
  Interface operator (*)
     module procedure matrixMultiply, matrix4Multiply, realMatrixMultiply, realVectorMultiply, &
                      matrixVectorMultiply, matrix4Vector4Multiply, matrix4VectorMultiply, vectorDotProduct
  end interface
  interface operator (/)
     module procedure vectorRealDivide
  End interface
  Interface operator (.cross.)
     module procedure vectorCrossProduct
  End interface
  Interface operator (-)
     module procedure matrixNegate, matrixSubtract, vectorNegate, vectorSubtract
  End interface
  interface assignment (=)
    module procedure vector4ToVectorAssign, vectorToVector4Assign
  end interface

  type(matrix),  parameter:: Imatrix =matrix( reshape( (/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/) ))
  type(matrix4), parameter:: Imatrix4=matrix4(reshape( (/1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1./),(/4,4/) ))

contains
  !************************************************************************************
  !           MATRIX DEFINITIONS                          *
  !************************************************************************************

  !************************************************************************************
  !FUNCTION matrixSetValue(A) result(M)
  !  The function sets a value to a matrix with a provided array
  !ARGUMENTS:
  !  - A: array of real 
  !RESULT:
  !  - M: matrix
  !REMARKS:
  !  - the function tests if the array in argument A is (3:3), and provides an 
  !error message if not
  !***********************************************************************************
  function matrixSetValue(A) result(M)
    real,dimension(:,:),intent(in)::A
    type (matrix)::M

    if ((size(A,1)/=3) .or. (size(A,2)/=3))  then
       write(*,*)"Error in array size: array must be (3:3)"
    else
       M%A=A
    end if
  end function matrixSetValue

  function matrix4SetValue(mat,vec)  result(mat4)
     type(matrix), intent(in):: mat
   type(vector), intent(in):: vec
   type(matrix4):: mat4

   mat4%A(1:3,1:3)= mat%A
   mat4%A(1:3,4)  = vec%A
   mat4%A(4,1:4) = (/0.,0.,0.,1./)

   return
  end function matrix4SetValue
     

  !************************************************************************************
  !FUNCTION matrixAdd(M1,M2) result(res)
  !  The function adds 2 Matrix
  !ARGUMENTS:
  !  - M1,M2: matrix  
  !RESULT:
  !  - res: matrix
  !************************************************************************************
  function matrixAdd(M1,M2) result(res)
    type (matrix),intent(in)::M1,M2
    type (matrix)::res

    res%A=M1%A + M2%A
  end function matrixAdd

  !************************************************************************************
  !FUNCTION matrixNegate(M) result(res)
  !  The function negates a matrix: M to -M
  !ARGUMENTS:
  !  - M: matrix
  !RESULT:
  !  - res: matrix
  !************************************************************************************
  function matrixNegate(M) result(res)
    type (matrix),intent(in)::M
    type (matrix)::res

    res%A=-M%A
  end function matrixNegate

  !************************************************************************************
  !FUNCTION matrixSubstract(M1,M2) result(res)
  !  The function substract a matrix from another
  !ARGUMENTS:
  !  - M1,M2: matrix
  !RESULT:
  !  - res: matrix
  !************************************************************************************
  function matrixSubtract(M1,M2) result(res)
    type (matrix),intent(in)::M1,M2
    type (matrix)::res

    res%A=M1%A - M2%A
  end function matrixSubtract

  !************************************************************************************
  !FUNCTION realMatrixMultiply(r,M) result(res)
  !  The function multiplies a matrix by a real
  !ARGUMENTS:
  !  - r: real 
  !  - M: matrix
  !RESULT:
  !  - res: matrix
  !*************************************************************************************
  function realMatrixMultiply(r,M) result(res)
    type (matrix),intent(in)::M
    type (matrix)::res
    real,intent(in)::r

    res%A=r*M%A
  end function realMatrixMultiply

  !************************************************************************************
  !FUNCTION matrixMultiply(A) result(res)
  !  The function realize the product of 2 matrix  
  !ARGUMENTS:
  !  - M1,M2: matrix
  !RESULT:
  !  - res: matrix
  !REMARKS:
  !  - the intrinsic function MatMul is used. The arguments of this function must 
  !be arrays of the suitable size.
  !************************************************************************************ 
  function matrixMultiply(M1,M2) result(res)
    type (matrix),intent(in)::M1,M2
    type (matrix)::res
    res%A=matmul(M1%A,M2%A)
  end function matrixMultiply
  
  function matrixInverse(M) result(res)
    type(matrix), intent(in):: M
    type(matrix):: res
    real:: det
    det = M%A(1,1)*(M%A(2,2)*M%A(3,3)-M%A(2,3)*M%A(3,2))+ &
          M%A(1,2)*(M%A(2,3)*M%A(3,1)-M%A(2,1)*M%A(3,3))+ &
          M%A(1,3)*(M%A(2,1)*M%A(3,2)-M%A(2,2)*M%A(3,1))
    
    if (det.ne. 0.) then
      res = realMatrixMultiply(1./det, matrixTranspose(M))
    else
      call Error('Attempted to invert a singular matrix.')
    end if
  
  end function matrixInverse
  
  function matrixTranspose(M) result(res)
    type(matrix), intent(in):: M
    type(matrix):: res
    integer:: i,j
    
    do i=1,3
      do j=1,3
        res%A(i,j) = M%A(j,i)
      end do
    end do
  
  end function matrixTranspose

  !************************************************************************************
  !FUNCTION matrix4Multiply(A) result(res)
  !  The function realize the product of 2 matrix4, 4X4 matrices  
  !ARGUMENTS:
  !  - M1,M2: matrix4
  !RESULT:
  !  - res: matrix4
  !REMARKS:
  !  - We use a special form of 4*4 matrices and we have optimized the
  !         function for this kind of matrices
  !***********************************************************************************
  function matrix4Multiply(M1,M2) result(res)
    type (matrix4), intent(in):: M1,M2
    type (matrix4) :: res
    res%A(1,1)=M1%A(1,1)*M2%A(1,1)+M1%A(1,2)*M2%A(2,1)+M1%A(1,3)*M2%A(3,1)
    res%A(1,2)=M1%A(1,1)*M2%A(1,2)+M1%A(1,2)*M2%A(2,2)+M1%A(1,3)*M2%A(3,2)
    res%A(1,3)=M1%A(1,1)*M2%A(1,3)+M1%A(1,2)*M2%A(2,3)+M1%A(1,3)*M2%A(3,3)
    res%A(2,1)=M1%A(2,1)*M2%A(1,1)+M1%A(2,2)*M2%A(2,1)+M1%A(2,3)*M2%A(3,1)
    res%A(2,2)=M1%A(2,1)*M2%A(1,2)+M1%A(2,2)*M2%A(2,2)+M1%A(2,3)*M2%A(3,2)
    res%A(2,3)=M1%A(2,1)*M2%A(1,3)+M1%A(2,2)*M2%A(2,3)+M1%A(2,3)*M2%A(3,3)
    res%A(3,1)=M1%A(3,1)*M2%A(1,1)+M1%A(3,2)*M2%A(2,1)+M1%A(3,3)*M2%A(3,1)
    res%A(3,2)=M1%A(3,1)*M2%A(1,2)+M1%A(3,2)*M2%A(2,2)+M1%A(3,3)*M2%A(3,2)
    res%A(3,3)=M1%A(3,1)*M2%A(1,3)+M1%A(3,2)*M2%A(2,3)+M1%A(3,3)*M2%A(3,3)

    res%A(1,4)=M1%A(1,1)*M2%A(1,4)+M1%A(1,2)*M2%A(2,4)+M1%A(1,3)*M2%A(3,4)+M1%A(1,4)
    res%A(2,4)=M1%A(2,1)*M2%A(1,4)+M1%A(2,2)*M2%A(2,4)+M1%A(2,3)*M2%A(3,4)+M1%A(2,4)
    res%A(3,4)=M1%A(3,1)*M2%A(1,4)+M1%A(3,2)*M2%A(2,4)+M1%A(3,3)*M2%A(3,4)+M1%A(3,4)

    res%A(4,1)=0.
    res%A(4,2)=0.
    res%A(4,3)=0.
    res%A(4,4)=1.
  end function matrix4Multiply

  !************************************************************************************
  !FUNCTION matrixVectorMultiply(M,V) result(res)
  !  The function left multiplies a matrix by a vector
  !ARGUMENTS:
  !  - M: matrix
  !  - V: vector
  !RESULT:
  !  - res: matrix
  !************************************************************************************
  function matrixVectorMultiply(M,V) result(res)
    type (matrix),intent(in)::M
    type (vector),intent(in)::V
    type (vector)::res
    real::sum
    integer::i,k

    do i=1,3
       sum=0
       do k=1,3
          sum=sum + (M%A(i,k) * v%A(k))
       end do
       res%A(i)=sum
    end do
  end function matrixVectorMultiply

  !************************************************************************************
  !FUNCTION matrix4Vector4Multiply(M,V) result(res)
  !  The function left multiplies a matrix4 by a vector4
  !ARGUMENTS:
  !  - M: matrix4
  !  - V: vector4
  !RESULT:
  !  - res: matrix4
  !************************************************************************************
  function matrix4Vector4Multiply(M,V) result(res)
    type (matrix4),intent(in)::M
    type (vector4),intent(in)::V
    type (vector4)::res

    res%A(1)=M%A(1,1)*V%A(1)+M%A(1,2)*V%A(2)+M%A(1,3)*V%A(3)+M%A(1,4)
    res%A(2)=M%A(2,1)*V%A(1)+M%A(2,2)*V%A(2)+M%A(2,3)*V%A(3)+M%A(2,4)
    res%A(3)=M%A(3,1)*V%A(1)+M%A(3,2)*V%A(2)+M%A(3,3)*V%A(3)+M%A(3,4)
    res%A(4)=1
  end function matrix4Vector4Multiply

  !************************************************************************************
  !FUNCTION matrix4VectorMultiply(M,V) result(res)
  !  The function left multiplies a matrix4 by a vector
  !ARGUMENTS:
  !  - M: matrix4
  !  - V: vector
  !RESULT:
  !  - res: vector
  !************************************************************************************
  function matrix4VectorMultiply(M,V) result(res)
    type (matrix4),intent(in)::M
    type (vector),intent(in)::V
    type (vector)::res

    res%A(1)=M%A(1,1)*V%A(1)+M%A(1,2)*V%A(2)+M%A(1,3)*V%A(3)+M%A(1,4)
    res%A(2)=M%A(2,1)*V%A(1)+M%A(2,2)*V%A(2)+M%A(2,3)*V%A(3)+M%A(2,4)
    res%A(3)=M%A(3,1)*V%A(1)+M%A(3,2)*V%A(2)+M%A(3,3)*V%A(3)+M%A(3,4)

  end function matrix4VectorMultiply

  !************************************************************************************
  !FUNCTION Transpose(M) result(res)
  !  The function transposes a matrix M
  !ARGUMENTS:
  !  - M: matrix
  !RESULT:
  !  - res: matrix
  !***********************************************************************************
  function transpose(M) result(res)
    type (matrix),intent(in)::M
    type (matrix)::res
    integer::i,j

    do i=1,3
       do j=1,3
          res%A(i,j) = M%A(j,i)
       end do
    end do
  end function transpose

  !*************************************************************************************
  !SUBROUTINE PrintMatrix(M) / PrintMatrix(M) 
  !  The subroutine visualizes a matrix M / a Matrix4 M
  !ARGUMENTS:
  !  - M: matrix / Matrix4
  !REMARKS:
  !  - the format for each element of the matrix is G11.4: this means that 11 characters 
  !are available to Print the real (accounting the signs - and .) with 4 characters devoted 
  !to decimals. 
  !*************************************************************************************  
  subroutine PrintMatrix(M)
    implicit none
    type (matrix),intent(in)::M
    integer::i

    do i=1,3
       Print "(G11.4,A,G11.4,A,G11.4)",M%A(i,1),"",M%A(i,2),"",M%A(i,3)
    end do
    Print*," "
  end subroutine PrintMatrix

  subroutine PrintMatrix4(M)
    implicit none
    type (matrix4),intent(in)::M
    integer::i

    do i=1,4
       Print "(G11.4,A,G11.4,A,G11.4,A,G11.4)",M%A(i,1),"",M%A(i,2),"",M%A(i,3),"", M%A(i,4)
    end do
    Print*," "
  end subroutine PrintMatrix4

  !************************************************************************************
  !           VECTOR DEFINITIONS                                     *
  !************************************************************************************

  !************************************************************************************
  !FUNCTION vectorSetValue(A) result(V)
  !  The function sets a value to a vector with a provided array
  !ARGUMENTS:
  !  - A: array of real 
  !RESULT:
  !  - V: vector
  !REMARKS:
  !  - the function tests if the array in argument A is (3), and provides an 
  !error message if not
  !************************************************************************************
  function vectorSetValue(A) result(V)
    implicit none
    real,dimension(:),intent(in)::A
    type (vector)::V
    if (size(A)/=3)  then
       Print*,'Error in array size: array must be (3)'
    else
       V%A=A
    end if
  end function vectorSetValue

  function vectorSetCoordinates(x,y,z) result(V)
    real, intent(in)::x, y, z
    type(vector)::V
    V%A(1)=x
    V%A(2)=y
    V%A(3)=z
  end function vectorSetCoordinates
  
  function vectorSetKind8Coordinates(x,y,z) result(V)
    real(kind=8), intent(in)::x, y, z
    type(vectorKind8)::V
    V%A(1)=x
    V%A(2)=y
    V%A(3)=z
  end function vectorSetKind8Coordinates


  !************************************************************************************
  !FUNCTION vectorAdd(V1,V2) result(res)
  !  The function adds 2 Vectors
  !ARGUMENTS:
  !  - V1,V2: vector  
  !RESULT:
  !  - res: vector
  !************************************************************************************
  function vectorAdd(V1,V2) result(res)
    type (vector),intent(in)::V1,V2
    type (vector)::res

    res%A=V1%A + V2%A
  end function vectorAdd

  !************************************************************************************
  !FUNCTION vectorNegate(V) result(res)
  !  The function negates a Vector: V to -V
  !ARGUMENTS:
  !  - V: vector  
  !RESULT:
  !  - res: vector
  !************************************************************************************
  function vectorNegate(V) result(res)
    type (vector),intent(in)::V
    type (vector)::res
  
    res%A=-V%A

  end function vectorNegate

  !************************************************************************************
  !FUNCTION VectSubstract(V1,V2) result(res)
  !  The function substract a Vector to another
  !ARGUMENTS:
  !  - V1,V2: vector  
  !RESULT:
  !  - res: vector
  !************************************************************************************
  function vectorSubtract(V1,V2) result(res)
    type (vector),intent(in)::V1,V2
    type (vector)::res

    res%A=V1%A - V2%A
  end function vectorSubtract

  !************************************************************************************
  !FUNCTION realVectorMultiply(r,V) result(res)
  !  The function multiplies a Vector by a real
  !ARGUMENTS:
  !  - r: real
  !  - V: vector  
  !RESULT:
  !  - res: vector
  !************************************************************************************
  function realVectorMultiply(r,V) result(res)
    type (vector),intent(in)::v
    type (vector)::res
    real,intent(in)::r

    res%A=r*V%A
  end function realVectorMultiply
  
  !************************************************************************************
  !FUNCTION vectorRealDivide(V,r) result(res)
  !  The function divides a Vector by a real
  !ARGUMENTS:
  !  - r: real
  !  - V: vector  
  !RESULT:
  !  - res: vector
  !************************************************************************************
  function vectorRealDivide(V,r) result(res)
    type (vector),intent(in)::v
    type (vector)::res
    real,intent(in)::r

    res%A=V%A / r
  end function vectorRealDivide

  !************************************************************************************
  !FUNCTION vectorDotProduct(V1,V2) result(res)
  !  The function computes the dot product of 2 Vectors
  !ARGUMENTS:
  !  - V1,V2: vector  
  !RESULT:
  !  - res: real
  !REMARK:
  !  - The intrinsic function vectorDotProduct is used: the arguments must be 2 rank-one 
  !arrays.
  !**************************************************************************************
  function vectorDotProduct(V1,V2) result(res)
    type (vector),intent(in)::V1,V2
    real::res

    res=dot_product(V1%A,V2%A)

  end function vectorDotProduct

  function vectorDotProductDouble(V1,V2) result(res)
    type (vector),intent(in)::V1,V2
    real(kind=selected_real_kind(12))::res

    res=V1%A(1)*V2%A(1)+V1%A(2)*V2%A(2)+V1%A(3)*V2%A(3)

  end function vectorDotProductDouble

  !************************************************************************************
  !FUNCTION vectorCrossProduct(V1,V2) result(res)
  !  The function computes the cross product of 2 Vectors
  !ARGUMENTS:
  !  - V1,V2: vector  
  !RESULT:
  !  - res: vector
  !**************************************************************************************
  function vectorCrossProduct(V1,V2) result(res)
    type (vector),intent(in)::V1,V2
    type (vector)::res

    res%A(1)=V1%A(2)*V2%A(3)-V1%A(3)*V2%A(2)
    res%A(2)=V1%A(3)*V2%A(1)-V1%A(1)*V2%A(3)
    res%A(3)=V1%A(1)*V2%A(2)-V1%A(2)*V2%A(1)    
  end function vectorCrossProduct

  !************************************************************************************
  !FUNCTION vectorAbsolute(V) result(res)
  !  The function computes the norm of a vector
  !ARGUMENTS:
  !  - V: vector  
  !RESULT:
  !  - res: real
  !************************************************************************************
  function vectorAbsolute(V) result(res)
    type (vector),intent(in)::V
    real::res

    res=sqrt((V%A(1))**2 + (V%A(2))**2 + (V%A(3))**2)  
  end function vectorAbsolute
  
  function inverseTangent(r1, r2) result(res)
    real,intent(in)::r1, r2
    real::res
    
    if (r1.ne.0.0) then
      res = atan(r2/r1)
    else
      ! mathematically we're simply solving for tangent at pi/2, 
      ! but we have to work around the computation limitations
      ! of dividing by zero.
      if (r2.gt.0.0) then
        res = pi/2.0
      else if (r2.lt.0.0) then
        res = -pi/2.0
      else
        res = 0.0
      end if
    end if 
  end function inverseTangent  

  !*************************************************************************************
  !SUBROUTINE PrintVector(V) / PrintVector(V)
  !  The subroutine visualizes a vector V / a vector4 V
  !ARGUMENTS:
  !  - V: vector / vector4 
  !REMARKS:
  !  - the format for each element of the vector is G11.4: this means that 11 characters 
  !are available to Print the real (accounting the signs - and .) with 4 characters devoted 
  !to decimals. 
  !*************************************************************************************
  subroutine PrintVector(V)
    implicit none
    type (vector),intent(in)::V
    integer::i

    do i=1,3
       Print "(G15.6)",V%A(i)
    end do
    Print*," "
  end subroutine PrintVector

  subroutine PrintVector4(V)
    implicit none
    type (vector4),intent(in)::V
    integer::i

    do i=1,4
       Print "(G11.4)",V%A(i)
    end do
    Print*," "
  end subroutine PrintVector4

  subroutine vector4ToVectorAssign(vec, vec4) 
  type(vector), intent(out):: vec
    type(vector4), intent(in):: vec4

  vec%A=vec4%A(1:3)

  end subroutine vector4ToVectorAssign

  subroutine vectorToVector4Assign(vec4, vec)
    type(vector4), intent(out):: vec4
  type(vector), intent(in):: vec

  vec4%A(1:3)=vec%A
  vec4%A(4) = 1.

  end subroutine vectorToVector4Assign

function createRotationMatrix(eulerAxis,angle) result(res)
  type(vector), intent(in)::eulerAxis
  type(matrix)::res
  real::angle
 
  if ((eulerAxis%A(1)==1.0).and.(eulerAxis%A(2)==0.0).and.(eulerAxis%A(3)==0.0)) then
    res%A(1,1)= 1.;             res%A(1,2)=0.;             res%A(1,3)=0.
    res%A(2,1)= 0.;             res%A(2,2)=cos(angle);     res%A(2,3)=-sin(angle)
    res%A(3,1)= 0.;             res%A(3,2)=-res%A(2,3);    res%A(3,3)=res%A(2,2)
  else if ((eulerAxis%A(2)==1.0).and.(eulerAxis%A(1)==0.0).and.(eulerAxis%A(3)==0.0)) then
    res%A(1,1)=cos(angle);     res%A(1,2)=0.;              res%A(1,3)=sin(angle)
    res%A(2,1)=0.;             res%A(2,2)=1.;              res%A(2,3)=0.
    res%A(3,1)=-res%A(1,3);     res%A(3,2)=0.;              res%A(3,3)=res%A(1,1)
  else if ((eulerAxis%A(3)==1.0).and.(eulerAxis%A(1)==0.0).and.(eulerAxis%A(2)==0.0)) then
    res%A(1,1)=cos(angle);     res%A(1,2)=-sin(angle);    res%A(1,3)=0.
    res%A(2,1)=-res%A(1,2);     res%A(2,2)=res%A(1,1);      res%A(2,3)=0.
    res%A(3,1)=0.;             res%A(3,2)=0.;             res%A(3,3)=1.
  else  
    res=MatRotation(eulerAxis,angle) 
  end if 
return

end function createRotationMatrix

!***********************************************************
!SUBROUTINE MatRotation(vect3, angle) result(rotmat)
!  To build the second element of a baseChange, we need to 
! create a subroutine which is able to calculate the matrix of a
! rotation given its axis (vect3) and its angle (angle). The 
! subroutine creates the rotation matrix and stores it in rotmat.
! To have more informations about the algorithm used to create 
! this matrix of rotation : read the corresponding documentation:
! Guillaume Bres master thesis.
!*************************************************************
function MatRotation(eulerAxis,angle) Result(rotmat)
  type(vector),intent(in)::eulerAxis
  real :: angle , normv1, normv3
  real, dimension(3,3):: temp2DArray
  type (vector) :: vect1, vect2, vect3
  type(matrix) :: P,invP, rotmat 
  
  vect3=eulerAxis
  !Step1 Calculate of a base containing vect3  
  if (vect3%A(3)/=0) then 
    vect1%A(1)=1
    vect1%A(2)=1
    vect1%A(3)=-vect3%A(1)/vect3%A(3)-vect3%A(2)/vect3%A(3)
    normv1=vectorAbsolute(vect1) 
    normv3=vectorAbsolute(vect3)
    vect1=realVectorMultiply((1/normv1),vect1)
    vect3=realVectorMultiply((1/normv3),vect3)
    vect2=vectorCrossProduct(vect3,vect1)
  end if
  if ((vect3%A(3)==0) .and. ((vect3%A(1)/=0) .or. (vect3%A(2)/=0))) then
    vect1%A(1)=0.
    vect1%A(2)=0.  
    vect1%A(3)=1.
    normv3=vectorAbsolute(vect3)
    vect3=realVectorMultiply((1./normv3),vect3)
    vect2=vectorCrossProduct(vect3,vect1) 
  end if  
  if ((vect3%A(1)==0.) .and. (vect3%A(2)==0.) .and. (vect3%A(3)==0.)) then
    write(*,*) 'error, wrong vector input'
  end if   
  !Step2 Calculate the matrix P that changes the usual base in the base containing vect3
  temp2DArray(:,1)=vect1%A
  temp2DArray(:,2)=vect2%A
  temp2DArray(:,3)=vect3%A
  P=matrixSetValue(temp2DArray)
  !Step3 Calculate the inverse of the P matrix
  invP=Transpose(P)
  !Step4 Definition of the matrix of rotation in the base containing vect3
  temp2DArray(1,1)=cos(angle);         temp2DArray(1,2)=-sin(angle);      temp2DArray(1,3)=0.
  temp2DArray(2,1)=-temp2DArray(1,2);  temp2DArray(2,2)=temp2DArray(1,1); temp2DArray(2,3)=0.
  temp2DArray(3,1)=0.;                 temp2DArray(3,2)=0.;               temp2DArray(3,3)=1.   
  !Step5 Calculate the rotation matrix in the usual base
  rotmat=matrixMultiply(P,matrixMultiply(matrixSetValue(temp2DArray),invP))
  
end function MatRotation

subroutine BroadcastVectorsss (vectorArray)
  type(vector), dimension(:,:,:), intent(inout)::vectorArray
  real, dimension(:,:,:,:), allocatable::tempArray
  integer::i, j, k, arraySize1, arraySize2, arraySize3
  
  if (IsMaster()) then
    arraySize1=size(vectorArray,1)
    arraySize2=size(vectorArray,2)
    arraySize3=size(vectorArray,3)
  end if
  call BroadcastInteger(arraySize1)
  call BroadcastInteger(arraySize2)
  call BroadcastInteger(arraySize3)
  allocate(tempArray(arraySize1, arraySize2, arraySize3, 3))
  if (IsMaster()) then
    do i=1, arraySize1
      do j=1, arraySize2
        do k=1, arraySize3
          tempArray(i,j,k,:)=vectorArray(i,j,k)%A(:)
        end do
      end do
    end do
  end if
  call BroadcastRealssss(tempArray)
  if (NotMaster()) then
    do i=1, arraySize1
      do j=1, arraySize2
        do k=1, arraySize3
          vectorArray(i,j,k)%A(:)=tempArray(i,j,k,:)
        end do
      end do
    end do
  end if
  deallocate(tempArray)
    
end subroutine BroadcastVectorsss

subroutine SendVectorsss (vectorArray, proc, master)
  type(vector), dimension(:,:,:), intent(in)::vectorArray
  real, dimension(:,:,:,:), allocatable::tempArray
  integer::i, j, k, arraySize1, arraySize2, arraySize3, proc, master
  arraySize1=size(vectorArray,1)
  arraySize2=size(vectorArray,2)
  arraySize3=size(vectorArray,3)
  call SendInteger(arraySize1,         proc, master)
  call SendInteger(arraySize2,         proc, master)
  call SendInteger(arraySize3,         proc, master)
  allocate(tempArray(arraySize1, arraySize2, arraySize3, 3))
  do i=1, arraySize1
    do j=1, arraySize2
      do k=1, arraySize3
        tempArray(i,j,k,:)=vectorArray(i,j,k)%A(:)
      end do
    end do
  end do
  call SendRealssss(tempArray,          proc, master)
  deallocate(tempArray)
end subroutine SendVectorsss

subroutine ReceiveVectorsss(vectorArray, proc, master)
  type(vector), dimension(:,:,:), intent(inout)::vectorArray
  real, dimension(:,:,:,:), allocatable::tempArray
  integer::i, j, k, arraySize1, arraySize2, arraySize3, proc, master
  call ReceiveInteger(arraySize1,       proc, master)
  call ReceiveInteger(arraySize2,       proc, master)
  call ReceiveInteger(arraySize3,       proc, master)
  allocate(tempArray(arraySize1, arraySize2, arraySize3, 3))
  call ReceiveRealssss(tempArray,       proc, master)
  do i=1, arraySize1
    do j=1, arraySize2
      do k=1, arraySize3
        vectorArray(i,j,k)%A(:)=tempArray(i,j,k,:)
      end do
    end do
  end do
  deallocate(tempArray)
end subroutine ReceiveVectorsss

subroutine BroadcastVectorss (vectorArray)
  type(vector), dimension(:,:), intent(inout)::vectorArray
  real, dimension(:,:,:), allocatable::tempArray
  integer::i, j, arraySize1, arraySize2
  
  if (IsMaster()) then
    arraySize1=size(vectorArray,1)
    arraySize2=size(vectorArray,2)
  end if
  call BroadcastInteger(arraySize1)
  call BroadcastInteger(arraySize2)
  allocate(tempArray(arraySize1, arraySize2, 3))
  if (IsMaster()) then
    do i=1, arraySize1
      do j=1, arraySize2
        tempArray(i,j,:)=vectorArray(i,j)%A(:)
      end do
    end do
  end if
  call BroadcastRealsss(tempArray)
  if (NotMaster()) then
    do i=1, arraySize1
      do j=1, arraySize2
        vectorArray(i,j)%A(:)=tempArray(i,j,:)
      end do
    end do
  end if
  deallocate(tempArray)
    
end subroutine BroadcastVectorss

subroutine SendVectorss (vectorArray, proc, master)
  type(vector), dimension(:,:), intent(in)::vectorArray
  real, dimension(:,:,:), allocatable::tempArray
  integer::i, j, arraySize1, arraySize2, proc, master
  arraySize1=size(vectorArray,1)
  arraySize2=size(vectorArray,2)
  call SendInteger(arraySize1,         proc, master)
  call SendInteger(arraySize2,         proc, master)
  allocate(tempArray(arraySize1, arraySize2, 3))
  do i=1, arraySize1
    do j=1, arraySize2
      tempArray(i,j,:)=vectorArray(i,j)%A(:)
    end do
  end do
  call SendRealsss(tempArray,          proc, master)
  deallocate(tempArray)
end subroutine SendVectorss

subroutine ReceiveVectorss (vectorArray, proc, master)
  type(vector), dimension(:,:), intent(inout)::vectorArray
  real, dimension(:,:,:), allocatable::tempArray
  integer::i, j, arraySize1, arraySize2, proc, master
  call ReceiveInteger(arraySize1,       proc, master)
  call ReceiveInteger(arraySize2,       proc, master)
  allocate(tempArray(arraySize1, arraySize2, 3))
  call ReceiveRealsss(tempArray,        proc, master)
  do i=1, arraySize1
    do j=1, arraySize2
      vectorArray(i,j)%A(:)=tempArray(i,j,:)
    end do
  end do
  deallocate(tempArray)
end subroutine ReceiveVectorss

subroutine BroadcastVectors (vectorArray)
  type(vector), dimension(:), intent(inout)::vectorArray
  real, dimension(:,:), allocatable::tempArray
  integer::i, arraySize
  
  if (IsMaster()) arraySize=size(vectorArray,1)
  call BroadcastInteger(arraySize)
  allocate(tempArray(arraySize, 3))
  if (IsMaster()) then
    do i=1, arraySize
      tempArray(i,:)=vectorArray(i)%A(:)
    end do
  end if
  call BroadcastRealss(tempArray)
  if (NotMaster()) then
    do i=1, arraySize
      vectorArray(i)%A(:)=tempArray(i,:)
    end do
  end if
  deallocate(tempArray)
    
end subroutine BroadcastVectors

subroutine SendVectors (vectorArray, proc, master)
  type(vector), dimension(:), intent(in)::vectorArray
  real, dimension(:,:), allocatable::tempArray
  integer::i, arraySize, proc, master
  arraySize=size(vectorArray,1)
  call SendInteger(arraySize,             proc, master)
  allocate(tempArray(arraySize, 3))
  do i=1, arraySize
    tempArray(i,:)=vectorArray(i)%A(:)
  end do
  call SendRealss(tempArray,              proc, master)
  deallocate(tempArray)
end subroutine SendVectors

subroutine ReceiveVectors (vectorArray, proc, master)
  type(vector), dimension(:), intent(inout)::vectorArray
  real, dimension(:,:), allocatable::tempArray
  integer::i, arraySize, proc, master
  call ReceiveInteger(arraySize,          proc, master)
  allocate(tempArray(arraySize, 3))
  call ReceiveRealss(tempArray,           proc, master)
  do i=1, arraySize
    vectorArray(i)%A(:)=tempArray(i,:)
  end do
  deallocate(tempArray)
end subroutine ReceiveVectors

subroutine BroadcastVector(vectorIn)
  type(vector), intent(inout)::vectorIn
  real, dimension(3)::temp
  if (IsMaster()) temp=vectorIn%A(:)
  call BroadcastReals(temp)
  if (NotMaster()) vectorIn%A(:)=temp
end subroutine BroadcastVector
  
  subroutine SmartBroadcastVectors(array)
    type(vector), dimension(:), pointer::array
    logical::allocatedFlag
    integer::arraySize
    
    allocatedFlag = .false.
    arraySize = 0
    if(isMaster()) then
      allocatedFlag = associated(array)
      if (allocatedFlag) then
        arraySize = size(array)
      end if
    end if
    call BroadcastLogical(allocatedFlag)
    call BroadcastInteger(arraySize)
    if (NotMaster()) then
      if (allocatedFlag.and.arraySize.ne.0) then
        allocate(array(arraySize))
      else
        nullify(array)
      end if
    end if
    if (allocatedFlag) then
      call BroadcastVectors(array)
    end if
  end subroutine SmartBroadcastVectors

subroutine SendVector(vectorIn, proc, master)
  type(vector), intent(in)::vectorIn
  real, dimension(3)::temp
  integer::proc, master
  temp=vectorIn%A(:)
  call SendReals(temp, proc, master)
end subroutine SendVector

subroutine ReceiveVector(vectorIn, proc, master)
  type(vector), intent(inout)::vectorIn
  real, dimension(3)::temp
  integer::proc, master
  call ReceiveReals(temp, proc, master)
  vectorIn%A(:)=temp
end subroutine ReceiveVector

function CreateEvenlySpacedArray(min,max,n)result(array)
  real,intent(in)::min,max
  integer,intent(in)::n
  real(kind=8),dimension(:),pointer::array
  real(kind=8)::delta, minKind8, maxKind8
  integer::i
  allocate(array(n))
  minKind8 = min
  maxKind8 = max
  array(1) = minKind8
  if (n.gt.1) then
    delta=(maxKind8-minKind8)/(n-1)
    do i=2,n
      array(i)=min+(i-1)*delta
    end do
  end if
end function CreateEvenlySpacedArray


  !FUNCTION L2norm(array1,array2,arraysize) result(res)
  !  The function computes the L2norm of two array
  !ARGUMENTS:
  !  - array1,array2: array  
  !  - arraysize: integer
  !RESULT:
  !  - res: real
  !************************************************************************************
  function L2norm(array1,array2,arraysize) result(res)
    real, dimension(:), pointer::array1,array2
    integer, intent(in) ::arraysize
    integer ::i
    real::res,total

    total=0
    do i=1, arraySize
      if (array1(i) .eq. array2(i) ) then
        cycle
      end if
      if (array1(i) .eq. 0) then
        !total=total+( (array1(i)-array2(i) )/array2(i) )**2
        total=total+1
      else
        total=total+( (array1(i)-array2(i) )/array1(i) )**2
      end if
    end do  

    res=total/arraysize
    if (res .eq. 0) then
      return
    end if
    res=sqrt(total/arraysize)
  end function L2norm
  
  

end module mathModule

