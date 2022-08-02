! PSU-WOPWOP
! $Id: COB.f90 3374 2017-10-10 01:33:08Z brentner $
! $LastChangedDate: 2017-10-09 21:33:08 -0400 (Mon, 09 Oct 2017) $
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


!***************************************************
!Type: Change of Base module
!  This module encapsulates all the routines and structures
! for the COB data structure.  A CBStructure is the COB data
! read in from the namelist, it is stored here and, when we 
! integrate, this is turned into an actual change of base matrix 
! in the file COBObject.  A CBStructure is actually a
! node in a linked list, it has a data member that points to the 
! next CBStructure, as well as a data member for nonPeriodicData that
! is for data read in from a file.
! - This file and routines created by Leonard Lopes during
!   the summer of 2002.
!*****************************************************
module COBObject
  use MPIModule
  use constantsModule
  use interpolateModule
  use debugModule
  use mathModule
  use strings
  !***************************************************************************
  !TYPE: nonPeriodicData
  !   This is a part of CBStructure, we are grouping all of the non periodic data
  !  here instead of in CBStructure.
  !****************************************************************************
  type nonPeriodicData
     integer::iB
     integer::nTau
     real, dimension(:), pointer::tau =>NULL()
     real, dimension(:), pointer::angle =>NULL()
     real, dimension(:), pointer::dangle =>NULL()
     real, dimension(:), pointer::ddangle =>NULL()
     type(vector), dimension(:), pointer::axis =>NULL()
     type(vector), dimension(:), pointer::trans =>NULL()
     type(vector), dimension(:), pointer::dtrans =>NULL()
     type(vector), dimension(:), pointer::ddtrans =>NULL()
  end type nonPeriodicData

  type functionData
     type(vector)::AH
     type(vector)::VH
     type(vector)::Y0
     real::OMEGA
     real::PSI0
     real::A0
     real,dimension(100)::A ! OK, 100 is overkill... memory is cheap
     real,dimension(100)::B
  end type functionData

  !**************************************************************************
  !TYPE: CBStructure
  !   This is the COB data before it gets changed into an actual change of base matrix
  !  when the integration happens, it contains all the periodic data plus a non
  !  periodic data structure when the COB is non periodic with time.
  !*************************************************************************** 
  type CBStructure
     character(len=4096)::Title
     integer::iB
     type(nonPeriodicData)::nonPeriodicCB
     type(functionData)::funcData, storedData
     logical::rotation,windFrame,flightPath
     integer::AxisType
     integer::TranslationType
     integer::AngleType
     integer::flightPathAccel
     integer::FlightPathProfile
     character(len=4096)::filename
     type(vector)::axisvalue
     type(vector)::translationvalue
     real::anglevalue
     type(CBStructure),pointer::next
  end type CBStructure

  type baseChange
     type(matrix)::rotation
     type(vector)::position, velocity, omega
     type(vector)::acceleration, omegaDot
  end type baseChange

  type(baseChange)::resultantBase
  !this is the CB num, that cooresponds to rotation.
  integer::CBrotation
  type(CBStructure), pointer::rotationCB
  !These are two index values  needed so that we can save time when interpolating
  integer::sourceIndex, observerIndex

  integer, parameter :: NO_VALUE = 0, &
                        TIME_INDEPENDENT = 1, &
                        KNOWN_FUNCTION = 2, &
                        PERIODIC = 3, &
                        APERIODIC = 4, &
                        
                        CONSTANT_VELOCITY = 0, &
                        CONSTANT_ACCELERATION = 1, &
                        
                        CUSTOM = 0, &                        
                        TAKEOFF = 1, &
                        OVERFLIGHT = 2, &
                        APPROACH = 3
                        
contains
  
  subroutine WriteNonPerDebugInfo(NonPData,unitnum)
    type(nonPeriodicData), intent(in)::nonPData
    integer::unitnum

    write(unitnum,*) '*** NonPeriodicDebugInfo ***'
    write(unitnum,*) 'iB= ', trim(integertostring(NonPData%iB))
    write(unitnum,*) 'nTau= ', trim(integertostring(NonPData%nTau))
    if (associated(NonPData%tau)) then
      write(unitnum,*) 'Size of tau array= ', trim(integertostring(size(NonPData%tau)))
      write(unitNum,*) NonPData%tau
    end if
    if (associated(NonPData%angle)) then
      write(unitnum,*) 'Size of angle array= ', trim(integertostring(size(NonPData%angle)))
      write(unitNum,*) NonPData%angle
    end if
    if (associated(NonPData%dangle)) then
      write(unitnum,*) 'Size of dangle array= ', trim(integertostring(size(NonPData%dangle)))
      write(unitNum,*) NonPData%dangle
    end if
    if (associated(NonPData%ddangle)) then
      write(unitnum,*) 'Size of ddangle array= ', trim(integertostring(size(NonPData%ddangle)))
      write(unitNum,*) NonPData%ddangle
    end if
    if (associated(NonPData%axis)) then
      write(unitnum,*) 'Size of axis array= ', trim(integertostring(size(NonPData%axis)))
      write(unitNum,*) NonPData%axis
    end if
    if (associated(NonPData%trans)) then
      write(unitnum,*) 'Size of trans array= ', trim(integertostring(size(NonPData%trans)))
      write(unitNum,*) NonPData%trans
    end if
    if (associated(NonPData%dtrans)) then
      write(unitnum,*) 'Size of dtrans array= ', trim(integertostring(size(NonPData%dtrans)))
      write(unitNum,*) NonPData%dtrans
    end if
    if (associated(NonPData%ddtrans)) then
      write(unitnum,*) 'Size of ddtrans array= ', trim(integertostring(size(NonPData%ddtrans)))
      write(unitNum,*) NonPData%ddtrans
    end if
    write(unitnum,*) '--- End NonPericDebug Info ---'
 
  end subroutine WriteNonPerDebugInfo
    
  subroutine WriteFuncDebugInfo(func,unitnum)
    type(functionData), intent(in)::func
    integer::unitnum

    write(unitnum,*) '*** COBFunctionDebugInfo ***'
    write(unitnum,*) 'AH=', func%AH
    write(unitnum,*) 'VH=', func%VH
    write(unitnum,*) 'Y0=', func%Y0
    write(unitnum,*) 'OMEGA=', func%OMEGA
    write(unitnum,*) 'PSI0=', func%PSI0
    write(unitnum,*) 'A0=', func%A0
    write(unitnum,*) 'A=', func%A
    write(unitnum,*) 'B=', func%B
    write(unitnum,*) '--- End COB Function Data Debug Info ---' 
  end subroutine WriteFuncDebugInfo
  
  recursive subroutine WriteCBListDebugInfo(CB,unitnum)
    type(CBStructure), pointer::CB
    integer::unitnum
    write(unitnum,*) '*** COB Debug Info ***'
    write(unitnum,*) 'Title=', trim(CB%title)
    write(unitnum,*) 'iB=', trim(integertostring(CB%iB))
    write(unitnum,*) 'Rotation? ', CB%rotation
    write(unitnum,*) 'WindFrame? ', CB%windFrame
    write(unitnum,*) 'AxisType= ', CB%AxisType
    write(unitnum,*) 'TranslationType= ', CB%TranslationType
    write(unitnum,*) 'AngleType= ', CB%AngleType
    write(unitnum,*) 'Filename=', trim(CB%filename)
    write(unitnum,*) 'AxisValue=', CB%axisValue
    write(unitnum,*) 'TranslationValue=', CB%translationValue
    write(unitnum,*) 'AngleValue=', CB%angleValue
    write(unitnum,*) 'NonPeriodicData'
    call WriteNonPerDebugInfo(CB%nonPeriodicCB,unitnum)
    write(unitnum,*) 'FunctionData'
    call WriteFuncDebugInfo(CB%funcData,unitnum)
    write(unitnum,*) 'Next CB associated?', associated(CB%next) 
    if (associated(CB%next)) then
      call WriteCBListDebugInfo(CB%next,unitnum)
    end if
    write(unitnum,*) '--- End CB Debug Info ---'
  end subroutine WriteCBListDebugInfo

  subroutine initializeInterpolateIndex()
    sourceIndex=1
    observerIndex=1
  end subroutine initializeInterpolateIndex

  subroutine CalculateResultantChangeOfBase(CBList,n,tau)
    use mathModule, only: Imatrix4
    real::tau
    integer::n
    type(CBStructure), pointer::CBList
    type(baseChange)::firstBase
    type(vector)::zero
    
    zero%A(:)=0
    resultantBase%rotation=Imatrix
    resultantBase%position=zero
    resultantBase%velocity=zero
    resultantBase%acceleration=zero
    resultantBase%omega=zero
    resultantBase%omegaDot=zero
    firstBase%rotation=Imatrix
    firstBase%position=zero
    firstBase%velocity=zero
    firstBase%acceleration=zero
    firstBase%omega=zero
    firstBase%omegaDot=zero
    if(associated(CBList)) then
         call reduceChangeOfBase(CBList, 1, n, tau, firstBase, resultantBase)
    end if

  end subroutine CalculateResultantChangeOfBase

  recursive subroutine ReduceChangeOfBase(CBList, i, n, tau, CBOfThisBaseInOrig, CBOfNBaseInOrig)
    type(CBStructure), pointer::CBList
    type(baseChange)::CBOfThisBaseInOrig, CBOfNextBaseInOrig, CBOfNBaseInOrig
    real, intent(inout)::tau
    real:: tau0
    integer::i, n

    type(vector)::rotationVector
    
    if (.not. associated(CBList)) then
      return
    end if
    if (CBList%flightPath) then
      call updateFuncData(CBList%funcData, CBList%nonPeriodicCB, tau, tau0)
    else
      tau0 = 0.0
    end if

    rotationVector=axis(CBList,tau,sourceIndex)
    CBOfNextBaseInOrig%rotation    =CBOfThisBaseInOrig%rotation*&
         &(createRotationMatrix(rotationVector,angle(CBList,tau,sourceIndex)))
    CBOfNextBaseInOrig%omega       =CBOfThisBaseInOrig%rotation*(dangle(CBList,tau,sourceIndex)*rotationVector)
    CBOfNextBaseInOrig%omegaDot    =CBOfThisBaseInOrig%rotation*(ddangle(CBList,tau,sourceIndex)*rotationVector) + &
         &(CBOfNBaseInOrig%omega.cross.CBOfNextBaseInOrig%omega)
    CBOfNextBaseInOrig%position    =CBOfThisBaseInOrig%rotation*readTranslation(CBList,tau,tau0,sourceIndex)
    CBOfNextBaseInOrig%velocity    =CBOfThisBaseInOrig%rotation*readVelocity(CBList,tau,tau0,sourceIndex)
    CBOfNextBaseInOrig%acceleration=CBOfThisBaseInOrig%rotation*readAcceleration(CBList,tau,sourceIndex) + &
         &(CBOfNBaseInOrig%omegaDot.cross.CBOfNextBaseInOrig%position) + &
         &(CBOfNBaseInOrig%omega.cross.&
         &  (CBOfNBaseInOrig%omega.cross.CBOfNextBaseInOrig%position)) + &
         &(2.*(CBOfNBaseInOrig%omega.cross.CBOfNextBaseInOrig%velocity))
    !Forward Sums
    CBOfNBaseInOrig%omega          =CBOfNBaseInOrig%omega    + CBOfNextBaseInOrig%omega
    CBOfNBaseInOrig%omegaDot       =CBOfNBaseInOrig%omegaDot + CBOfNextBaseInOrig%omegaDot
    
    
    if(i<n) then
      call reduceChangeOfBase(CBList%next, i+1, n, tau, CBOfNextBaseInOrig, CBOfNBaseInOrig)
      !Reverse Sums
      CBOfNBaseInOrig%position     =CBOfNBaseInOrig%position + CBOfNextBaseInOrig%position
      CBOfNBaseInOrig%velocity     =CBOfNBaseInOrig%velocity + CBOfNextBaseInOrig%velocity + &
           &(CBOfThisBaseInOrig%omega.cross.CBOfNBaseInOrig%position)
      CBOfNBaseInOrig%acceleration =CBOfNBaseInOrig%acceleration + CBOfNextBaseInOrig%acceleration
    else
       CBOfNBaseInOrig%rotation     =CBOfNextBaseInOrig%rotation
       CBOfNBaseInOrig%position     =CBOfNextBaseInOrig%position
       CBOfNBaseInOrig%velocity     =CBOfNextBaseInOrig%velocity
       CBOfNBaseInOrig%acceleration =CBOfNextBaseInOrig%acceleration
    end if

  end subroutine ReduceChangeOfBase

  function observerTransformationMatrix(CBList,tau) result(T)
    use mathModule, only: Imatrix4, matrix4SetValue    
    implicit none
    real, intent(in)::tau
    real:: tau0
    integer::i
    type(CBStructure), pointer::CBList
    type(CBStructure)::CBData
    type(matrix4)::T
    type(matrix)::rotation
    if (.not. associated(CBList)) then
      T = Imatrix4
      return
    end if
    T=Imatrix4
    do i=1, getSize(CBList)
      if (CBList%flightPath) then
        call updateFuncData(CBList%funcData, CBList%nonPeriodicCB, tau, tau0)
      else
        tau0 = 0.0
      end if
      
      call getCB(i, CBList, CBData)
      rotation=createRotationMatrix(axis(CBData,tau,observerIndex),angle(CBData,tau,observerIndex))
      T=T*matrix4SetValue(rotation,readTranslation(CBData,tau,tau0,observerIndex))
    end do
  end function observerTransformationMatrix
  
  subroutine updateFuncData(funcData, nonPerCB, tau, tau0)
    type(nonPeriodicData):: nonPerCB
    type(functionData):: funcData
    real, intent(in):: tau
    real, intent(out):: tau0
    integer:: i,n,index

    n = size(nonPerCB%tau)
    call hunt(nonPerCB%tau,n,tau,i)! i is the index of the value in arrayTau just below tau.
    if(i.lt.1) then
      index = 1
      funcData%AH = vectorSetCoordinates(0.0, 0.0, 0.0)
    else if (i.gt.n) then
      index = n
      funcData%AH = vectorSetCoordinates(0.0, 0.0, 0.0)
    else
      index = i
      funcData%AH = nonPerCB%ddtrans(index)
    end if
    tau0        = nonPerCB%tau(index)
    funcData%Y0 = nonPerCB%trans(index)
    funcData%VH = nonPerCB%dtrans(index)
    
    if (i.lt.i.and.vectorAbsolute(nonPerCB%ddtrans(1)).ne.0.0 .or. &
         i.gt.n.and.vectorAbsolute(nonPerCB%ddtrans(n)).ne.0.0) then
      call Warning('The source time range is outside the time range of the',&
                   'flight path.  When this occurs the velocity remains',&
                   'constant but the acceleration is set to zero.')
    end if
 
  end subroutine updateFuncData


  subroutine CalculateVelaccel(eta,vel,accel)
    type(vector), intent(in)::eta
    type(vector)::etaB0, vel, accel, omegaCrossEta
 
    etaB0=resultantBase%rotation*eta
    omegaCrossEta=resultantBase%omega.cross.etaB0
    vel  =resultantBase%velocity+omegaCrossEta
    accel=(resultantBase%acceleration) + &
          (resultantBase%omegaDot.cross.etaB0) + &
          (resultantBase%omega.cross.omegaCrossEta)

  end subroutine calculateVelaccel
  
  !************************************************
  !SUBROUTINE Position(eta) result(res)
  !  This function finds the position of a point associated
  ! with the T matrix defined in this file by multiplying the
  ! initial location to a new location.
  !ARGUEMENTS:
  ! - eta: The initial location of the point on the blade.
  ! - res: The vector that is the position of the point at the 
  !         time associated with the T matrix.
  !*****************************************************
  function Position(eta) result(res)
    type(vector)::eta,res

    res=resultantBase%position + resultantBase%rotation*eta

  end function Position

  !************************************************************************
  ! This function gives the expression of a vector belonging to the frame N 
  ! into the frame 0
  ! N.B : to use vectorBaseChange, we need to have CB(N/0), this matrix is calculated
  ! in the velocity subroutines, thus we need to have calculated the velocity
  ! at the given tau before we use this function 
  !**********************************************************************
  function vectorBaseChange(V) result(res)
    type(vector) :: V, res

    res=resultantBase%rotation*V

  end function vectorBaseChange

  !************************************************************************
  ! This function differentiate a vector of constant magnitude 
  ! belonging to the frame N.
  !************************************************************************
  function vectorDerivative(V) result(res) 
    type(vector) :: V,res

    res=(resultantBase%omega.cross.vectorBaseChange(V))

  end function vectorDerivative

  !********************************************************************************
  !SUBROUTINE CreateCBData(CBData,unitNumber)
  !   This subroutine reads in a CB Structure from the namelist, it initializes all the
  !  data in the CB and also calls createNonPeriodicCB(CBData,unitNumber) when
  !  the CB is nonperiodic.
  !ARGUEMENTS:
  !  - CBData: the CBData that we want to create.
  !  - unitNumber: cooresponds to the namelist where the COB data is contained.
  !REMARKS:
  !  3.0.4.1b: added default definitions
  !SUBROUTINES USED:
  !  createNonPeirodicCB(CBData,unitNumber) in this file
  !*********************************************************************************
  subroutine CreateCBData(CBData,unitNumber)
    type(CBStructure), intent(inout)::CBData
    integer, intent(in)::unitNumber

    integer::iB, nbWayPoints
    character(len=4096)::Title, filename
    real::OMEGA, PSI0, A0
    real,dimension(100) :: A, B
    character(len=4096)::AxisType, TranslationType, AngleType, flightPathAccel, flightPathProfile
    logical::rotation, windframe, flightPath
    real, dimension(3):: AH, VH, Y0, axisvalue, translationvalue
    real::anglevalue

    !default values
    Title               = 'Default COB'
    iB                  = 0
    windFrame           = .false.
    rotation            = .false.
    flightPath          = .false.
    axistype            = 'TimeIndependent'
    translationtype     = 'TimeIndependent'
    angletype           = 'TimeIndependent'
    flightPathAccel     = 'ConstantVelocity'
    flightPathProfile   = 'Custom'
    AH(:)               = 0.0
    VH(:)               = 0.0
    Y0(:)               = 0.0
    OMEGA               = 0.0
    PSI0                = 0.0
    A0                  = 0.0
    A                   = 0.0
    B                   = 0.0
    filename            = 'void'
    axisvalue(:)        = 0.0
    axisValue(3)        = 1.0
    translationvalue(:) = 0.0
    angleValue          = 0.0
    nbWayPoints         = 0
    nullify(CBData%nonPeriodicCB%tau)
    nullify(CBData%nonPeriodicCB%angle)
    nullify(CBData%nonPeriodicCB%dangle)
    nullify(CBData%nonPeriodicCB%ddangle)
    nullify(CBData%nonPeriodicCB%axis)
    nullify(CBData%nonPeriodicCB%trans)
    nullify(CBData%nonPeriodicCB%dtrans)
    nullify(CBData%nonPeriodicCB%ddtrans)
    call ReadInCBData(unitNumber, CBData, Title, iB, windFrame, rotation, flightPath, axisType, &
                      translationType, angleType, flightPathAccel, AH, VH, Y0, OMEGA, PSI0, A0, A,  &
                      B, fileName, axisValue, translationValue, angleValue, nbWayPoints, flightPathProfile)
    if (CBData%filename/='void') then
         call CreateNonPeriodicCB(CBData%nonPeriodicCB,CBData%filename,unitNumber)
    end if
    if(flightPath) then
      call BuildFlightPathProfile(CBData%nonPeriodicCB, unitNumber, CBData%flightPathProfile, iB, nbWayPoints)
      
      if(nbWayPoints.gt.1) then
        call BuildFlightPath(CBData%nonPeriodicCB, CBData%flightPathAccel)
      else
        call Error('A flight path must include data for at least two waypoints.  Stopping')
      end if
    end if
    if (rotation) then
       CBrotation=iB 
       allocate(rotationCB)
       call CopyOneCB(rotationCB,CBData)
    end if
   
  end subroutine createCBData

  subroutine ReadInCBData(unitNumber, CBData, Title, iB, windFrame, rotation, flightPath, axisType, &
                      translationType, angleType, flightPathAccel, AH, VH, Y0, OMEGA, PSI0, A0, A, &
                      B,fileName, axisValue, translationValue, angleValue, nbWayPoints, flightPathProfile)
    type(CBStructure), intent(inout)::CBData
    integer, intent(in)::unitNumber
    integer::iB, nbWayPoints
    character(len=4096)::Title, filename
    character(len=4096)::AxisType, TranslationType, AngleType, flightPathAccel, flightPathProfile    
    real::anglevalue, OMEGA, PSI0, A0, A1, A2, B1, B2
    real, dimension(3):: AH, VH, Y0, axisvalue, translationvalue    
    real,dimension(100) :: A, B
    logical::rotation, windframe, flightPath

    namelist / CB / Title, rotation, windframe, flightPath, iB, AxisType, TranslationType, &
                     AngleType, AH, VH, Y0, OMEGA, PSI0, A0, A1, A2, B1, B2, A, B, &
                     filename, axisvalue, translationvalue, anglevalue, &
                     flightPathAccel, flightPathProfile, nbWayPoints
    
    ! These aren't passed in, so we need to remember to set them to zero     
    A1 = 0.0
    A2 = 0.0
    B1 = 0.0
    B2 = 0.0
         
    read(unitNumber, nml=CB)

    if (A1 /= 0.0 .and. A(1) /= 0.0 .and. A1 /= A(1)) then
      call Warning ("If specifying A1, do not set A as well (or set them equal...)")
    end if
    if (A2 /= 0.0 .and. A(2) /= 0.0 .and. A2 /= A(2)) then
      call Warning ("If specifying A1, do not set A as well (or set them equal...)")
    end if
    if (B1 /= 0.0 .and. B(1) /= 0.0 .and. B1 /= B(1)) then
      call Warning ("If specifying A1, do not set A as well (or set them equal...)")
    end if
    if (B2 /= 0.0 .and. B(2) /= 0.0 .and. B2 /= B(2)) then
      call Warning ("If specifying A1, do not set A as well (or set them equal...)")
    end if
   
    CBData%Title     = trim(Title)
    CBData%iB        = iB
    CBData%rotation  = rotation
    CBData%windFrame = windframe
    CBData%flightPath = flightPath
    
    call strToUpper (AxisType,        len_trim(AxisType))
    call strToUpper (AngleType,       len_trim(AngleType))
    call strToUpper (TranslationType, len_trim(TranslationType))
    call strToUpper (flightPathAccel, len_trim(flightPathAccel))
    call strToUpper (flightPathProfile, len_trim(flightPathProfile))
 
    if (trim(AxisType) == "TIMEINDEPENDENT") then
      CBData%AxisType=TIME_INDEPENDENT
    else if (trim(AxisType) == "KNOWNFUNCTION") then
      CBData%AxisType=KNOWN_FUNCTION
    else if (trim(AxisType) == "PERIODIC") then
      CBData%AxisType=PERIODIC
    else if (trim(AxisType) == "NONPERIODIC" .or. trim(AxisType) == "APERIODIC") then
      CBData%AxisType=APERIODIC
    else
      CBData%AxisType=TIME_INDEPENDENT
      write(*,*) 'Unknown axis type defaulting to TIME INDEPENDENT'
    end if  
    if (trim(AngleType) == "TIMEINDEPENDENT") then
      CBData%AngleType=TIME_INDEPENDENT
    else if (trim(AngleType) == "KNOWNFUNCTION") then
      CBData%AngleType=KNOWN_FUNCTION
    else if (trim(AngleType) == "PERIODIC") then
      CBData%AngleType=PERIODIC
    else if (trim(AngleType) == "NONPERIODIC" .or. trim(AngleType) == "APERIODIC" ) then
      CBData%AngleType=APERIODIC
    else
      CBData%AngleType=TIME_INDEPENDENT
      write(*,*) 'Unknown angle type defaulting to TIME INDEPENDENT'
    end if
    if (.not.FlightPath) then
      if (trim(TranslationType) == "TIMEINDEPENDENT") then
        CBData%TranslationType=TIME_INDEPENDENT
      else if (trim(TranslationType) == "KNOWNFUNCTION") then
        CBData%TranslationType=KNOWN_FUNCTION
      else if (trim(TranslationType) == "PERIODIC") then
        CBData%TranslationType=PERIODIC
      else if (trim(TranslationType) == "NONPERIODIC" .or. trim(TranslationType) == "APERIODIC") then
        CBData%TranslationType=APERIODIC
      else
        CBData%TranslationType=TIME_INDEPENDENT
        write(*,*) 'Unknown translation type defaulting to TIME INDEPENDENT'
      end if
    else
      if (trim(TranslationType) == "KNOWNFUNCTION") then
        CBData%TranslationType=KNOWN_FUNCTION
      else        
        CBData%TranslationType=KNOWN_FUNCTION
        write(*,*) 'Unknown translation type defaulting to KNOWN FUNCTION'
      end if
      
      if (trim(flightPathAccel) == "CONSTANTVELOCITY") then
        CBData%flightPathAccel=CONSTANT_VELOCITY
      else if (trim(flightPathAccel) == "CONSTANTACCELERATION") then
        CBData%flightPathAccel=CONSTANT_ACCELERATION
      else        
        CBData%flightPathAccel=CONSTANT_VELOCITY
        write(*,*) 'Unknown flight path type defaulting to CONSTANT VELOCITY'
      end if
      
      if (trim(FlightPathProfile) == "TAKEOFF") then
        CBData%flightPathProfile = TAKEOFF
      else if (trim(FlightPathProfile) == "OVERFLIGHT") then
        CBData%flightPathProfile = OVERFLIGHT
      else if (trim(FlightPathProfile) == "APPROACH") then
        CBData%flightPathProfile = APPROACH
      else
        CBData%flightPathProfile = CUSTOM
        write(*,*) 'Unknown flight path profile defaulting to CUSTOM'
      end if    
    end if
     
    CBData%funcData%AH      = vectorSetValue(AH)
    CBData%funcData%VH      = vectorSetValue(VH)
    CBData%funcData%Y0      = vectorSetValue(Y0)
    CBData%funcData%OMEGA   = OMEGA
    CBData%funcData%PSI0    = PSI0
    CBData%funcData%A0      = A0
    CBData%funcData%A       = A
    CBData%funcData%B       = B
    if (A1 /= 0)  CBData%funcData%A(1) = A1
    if (A2 /= 0)  CBData%funcData%A(2) = A2
    if (B1 /= 0)  CBData%funcData%B(1) = B1
    if (B2 /= 0)  CBData%funcData%B(2) = B2
    CBData%filename         = filename
    CBData%axisvalue        = vectorSetValue(axisvalue)
    CBData%translationvalue = vectorSetValue(translationvalue)
    CBData%anglevalue       = anglevalue
    nullify(CBData%next)

  end subroutine readInCBData

  !*******************************************************
  !SUBROUTINE destroyCBList(CBList)
  !  This routine will remove a CB from memory as well as all of it's
  ! children, it does this by calling itself recursively.
  !ARGUEMENTS:
  ! - CBList: a pointer to the CB we want to destroy 
  !SUBROUTINES USED:
  ! destroyCBList(CBList) which is this routine
  !******************************************************
  recursive subroutine destroyCBList(CBList)
    type(CBStructure),pointer::CBList

    if(associated(CBList)) then
      if (associated(CBList%next)) then
        call destroyCBList(CBList%next)
      end if
      deallocate(CBList)
      nullify(CBList)
    end if

  end subroutine destroyCBList

  !******************************************************
  !SUBROUTINE appendCB(CBList, nextCB)
  !  Because the order of a CBList is important, we need to add
  ! a CB to the end of the list. We do this by transversing the 
  ! list until we get to the last structure then creating a CB
  ! and append it onto the last one.
  !ARGUEMENTS:
  ! - CBList: The list that we want to append the CB onto.
  ! - nextCB: Thet new CB that will be inserted onto the end of the CBList.
  !SUBROUTINES USED:
  ! copyCB(CBStructure, CBStructure) in this file, appendCB(CBList, CBStructure)
  ! which is this file.
  !*********************************************************
  recursive subroutine appendCB(CBList,nextCB)
    type(CBStructure), pointer::CBList
    type(CBStructure), intent(in)::nextCB
    if(.not.associated(CBList)) then
      allocate(CBList)
      nullify(CBList%next)
      call copyCB(CBList,nextCB)
    else 
      call appendCB(CBList%next,nextCB)
    end if

  end subroutine appendCB

  !***************************************************
  !SUBROUTINE copyCB(new, CB)
  !  Because we can not simply change pointer locations
  ! to add a new CB to a CBList we need a routine which will
  ! physically copy data members and the data members of it's children
  ! from one CB to another, this is that routine.
  !ARGUEMENTS:
  ! - new: The CB that we want to create.
  ! - CB: The old CB that we want to copy from.
  !SUBROUTINES USED:
  ! copyCB(CBStructure, CBStructure) which is this file.
  !*****************************************************
  subroutine copyOneCB(new,CB)
    type(CBStructure)::new
    type(CBStructure),intent(in)::CB
  
    new%Title            = CB%Title
    new%iB               = CB%iB
    call CopyNonPeriodicCB(CB%nonPeriodicCB, new%nonPeriodicCB)
    !new%nonPeriodicCB    = CB%nonPeriodicCB
    new%rotation         = CB%rotation
    new%windframe        = CB%windFrame
    new%axistype         = CB%axistype
    new%translationtype  = CB%translationtype
    new%angletype        = CB%angletype  
    new%funcData         = CB%funcData
    new%filename         = CB%filename
    new%axisvalue        = CB%axisvalue
    new%translationvalue = CB%translationvalue
    new%angleValue       = CB%anglevalue
    nullify(new%next)

  end subroutine copyOneCB

  !****************************************************
  !SUBROUTINE copyCB(new, CB)
  !  Because we can not simply change pointer locations
  ! to add a new CB to a CBList we need a routine which will
  ! physically copy data members and the data members of it's children
  ! from one CB to another, this is that routine.
  !ARGUEMENTS:
  ! - new: The CB that we want to create.
  ! - CB: The old CB that we want to copy from.
  !SUBROUTINES USED:
  ! copyCB(CBStructure, CBStructure) which is this file.
  !*****************************************************
  recursive subroutine copyCB(new,CB)
    type(CBStructure),pointer::new
    type(CBStructure),intent(in)::CB
    if (.not. associated(new)) then
      return
    end if

    new%Title            = CB%Title
    new%iB               = CB%iB
    call CopyNonPeriodicCB(CB%nonPeriodicCB, new%nonPeriodicCB)
!    new%nonPeriodicCB    = CB%nonPeriodicCB
    new%rotation         = CB%rotation
    new%windframe        = CB%windFrame
    new%flightPath       = CB%flightPath
    new%axistype         = CB%axistype
    new%translationtype  = CB%translationtype
    new%angletype        = CB%angletype  
    new%funcData         = CB%funcData
    new%storedData       = CB%storedData
    new%filename         = CB%filename
    new%axisvalue        = CB%axisvalue
    new%translationvalue = CB%translationvalue
    new%angleValue       = CB%anglevalue
    nullify(new%next)
    if(associated(CB%next)) then
       allocate(new%next)
       call copyCB(new%next,CB%next)
    end if
  end subroutine copyCB

  subroutine CopyNonPeriodicCB(old, new)
    implicit none
    type(NonPeriodicData), intent(in)::old
    type(NonPeriodicData), intent(inout)::new
    new%iB = old%iB
    new%nTau = old%nTau
    nullify(new%tau, new%angle, new%dangle, new%ddangle)
    nullify(new%axis, new%trans, new%dtrans, new%ddtrans)
    if(.not.associated(old%tau)) return
    allocate(new%tau(size(old%tau))); new%tau = old%tau
    allocate(new%angle(size(old%angle))); new%angle = old%angle
    allocate(new%dangle(size(old%dangle))); new%dangle = old%dangle
    allocate(new%ddangle(size(old%ddangle))); new%ddangle = old%ddangle
    allocate(new%axis(size(old%axis))); new%axis = old%axis
    allocate(new%trans(size(old%trans))); new%trans = old%trans
    allocate(new%dtrans(size(old%dtrans))); new%dtrans = old%dtrans
    allocate(new%ddtrans(size(old%ddtrans))); new%ddtrans = old%ddtrans
  end subroutine CopyNonPeriodicCB

  !********************************************
  !SUBROUTINE replaceCB(CBNum, CBList, newCB)
  !  This routine will take a number, a CBList, and a newCB
  ! and insert the newCB into the CBList at the CBNum location
  ! erasing whatever CB was there before. It does this by 
  ! transversing the list until CBNum is zero then copy the data
  ! members into the CB that is already there.
  !ARGUEMENTS:
  ! - CBNum: The number location along the list that the newCB will
  !           be inserted.
  ! - CBList: The CBList that we are inserting into.
  ! - newCB: The new CB we are inserting.
  !SUBROUTINES USED:
  ! replaceCB(int, CBList, CBStructure) which is this routine.
  !**************************************************
  recursive subroutine replaceCB(CBNum,CBlist,newCB)
    integer, intent(in)::CBNum
    type(CBStructure), pointer::CBList
    type(CBStructure), intent(inout)::newCB
    
    if (.not. associated(CBList)) then
      return
    end if

    if(CBNum/=1) then
       call replaceCB(CBNum-1,CBList%next,newCB)
    else
       CBList%Title            = newCB%Title
       CBList%iB               = newCB%iB  
       call CopyNonPeriodicCB(CBList%nonPeriodicCB, newCB%nonPeriodicCB)
!       CBList%nonPeriodicCB    = newCB%nonPeriodicCB
       CBList%rotation         = newCB%rotation
       CBList%windframe        = newCB%windFrame
       CBList%axistype         = newCB%axistype
       CBList%translationtype  = newCB%translationtype
       CBList%angletype        = newCB%angletype  
       CBList%funcData         = newCB%funcData
       CBList%filename         = newCB%filename
       CBList%axisvalue        = newCB%axisvalue
       CBList%translationvalue = newCB%translationvalue
       CBList%angleValue       = newCB%anglevalue
    end if

  end subroutine replaceCB

  !*****************************************************
  !FUNCTION getCB(CBNum, CBList) result (res)
  !  This function will transverse the list until CBNum is 
  ! zero then return resulting CB at that location.
  !ARGUEMENTS:
  ! - CBNum: This is the number location along the list that 
  !           cooresponds to the CB we want.
  ! - CBList: The CBList that we are transversing.
  ! - res: The resulting CBStructure at CBNum location along CBList.
  !SUBROUTINES USED:
  ! getCB(integer, CBStructure) which is this routine.
  !******************************************************
  recursive subroutine getCB(CBNum,CBList,res)
    integer::CBNum
    type(CBStructure),pointer::CBList

    type(CBStructure)::res
    if (.not. associated(CBList)) then
      return
    end if

    if(CBNum==1) then
       res=CBList
    else
       call getCB(CBNum-1,CBList%next,res)
    end if

  end subroutine getCB

  !**************************************************
  !SUBROUTINE getSize(CBList) result(res)
  !  This routine will transverse the list, adding one to 
  ! the result until the end is reached, then return the result.
  !ARGUEMENTS:
  ! - CBList: The CBList that we want the size of.
  ! - res: The resulting size of the CBList.
  !*****************************************************
  function getSize(CBList) result(res)
    type(CBStructure), pointer::CBList

    type(CBStructure), pointer::temp
    integer::res

    if (.not. associated(CBList)) then
      res = 0
      return
    end if
    
    res=0
    if(associated(CBList)) then
      temp=>CBList
      do 
        if (.not.associated(temp%next)) then
          exit
        end if
        res=res+1
        temp=>temp%next
      end do
      res=res+1
    end if

  end function getSize

  !****************************************************
  !SUBROUTINE PrintCB(CBData)
  !  This routine will Print the data in CBData.
  !ARGUEMENTS:
  ! - CBData: The CB that we are Printing the information for.
  !**************************************************
  recursive subroutine PrintCB(CBData)
    type(CBStructure), intent(in)::CBData
    write(*,*)'/'
    write(*,*)'CB INFORMATION'
    write(*,*)'Title:  ',       trim(CBData%Title)
    write(*,*)'iB:  ',               CBData%iB
    write(*,*)'rotation:  ',         CBData%rotation
    write(*,*)'windframe:  ',        CBData%windFrame
    write(*,*)'axistype:  ',         CBData%axistype
    write(*,*)'translationtype:  ',  CBData%translationtype
    write(*,*)'angletype:  ',        CBData%angletype  
    write(*,*)'filename:  ',    trim(CBData%filename)
    write(*,*)'axisvalue:  ',        CBData%axisvalue
    write(*,*)'translationvalue:  ', CBData%translationvalue
    write(*,*)'angleValue:  ',       CBData%anglevalue
    write(*,*)'next pointer:  ',     associated(CBData%next)
    write(*,*)'/'
    if(associated(CBData%next))call PrintCB(CBData%next)

  end subroutine PrintCB

  !**************************************************************
  !SUBROUTINE createNonPeriodicCB(CBData,unitNumber)
  !   This routine reads from the file contianed in CBData, and initializes all
  !  of the non periodic values
  !ARGUEMENTS:
  !  - nonPerData: The non periodic data that we are reading in.
  !  - fileName: The file name that has this data.
  !  - unitNumber: The unitNumber that is associated with namelist, but only a number here.
  !SUBROUTINES USED:
  !  DifferentiateArrayReal(array,arraytau,darray) and DifferentiateArrayVector(array,arraytau,darray)
  !  in integrationModule in constants.f90
  !************************************************************************************************
  subroutine createNonPeriodicCB(nonPerData, fileName, unitNumber)
    implicit none
    type(nonPeriodicData), intent(inout)::nonPerData
    character(len=4096), intent(in)::fileName
    integer, intent(in)::unitNumber

    integer::j,k,angleFlag, axisFlag, translationFlag
    ! A flag is an integer which takes the value 0 when the component is not nonperiodic
    ! and takes the value 1 when it's nonperiodic
    integer :: nTau, unitNum
    unitNum = 66666
    open(unitNum, file=trim(globalFolderName)//trim(fileName), status='old')
    ! first we read the CB to which corresponds the data
    read(unitNum,*) k
    !CBData%nonPeriodicCB%iB=CBData%
    read(unitNum,*) nonPerData%nTau
    nTau=nonPerData%nTau
    read(unitNum,*) angleFlag, translationFlag, axisFlag
    allocate(nonPerData%tau(nTau))
    allocate(nonPerData%angle(nTau))
    allocate(nonPerData%dangle(nTau))
    allocate(nonPerData%ddangle(nTau))  
    allocate(nonPerData%axis(nTau))
    allocate(nonPerData%trans(nTau))
    allocate(nonPerData%dtrans(nTau))
    allocate(nonPerData%ddtrans(nTau))
    read(unitNum,*) (nonPerData%tau(j), j=1, nTau)
    if (angleFlag == 1) then
      read(unitNum,*) (nonPerData%angle(j), j=1, nTau)
    end if
    if (translationFlag == 1) then
      read(unitNum,*) (nonPerData%trans(j)%A(1), j=1, nTau)
      read(unitNum,*) (nonPerData%trans(j)%A(2), j=1, nTau) 
      read(unitNum,*) (nonPerData%trans(j)%A(3), j=1, nTau)
    end if
    if (axisFlag == 1) then
      read(unitNum,*) (nonPerData%axis(j)%A(1), j=1, nTau)
      read(unitNum,*) (nonPerData%axis(j)%A(2), j=1, nTau) 
      read(unitNum,*) (nonPerData%axis(j)%A(3), j=1, nTau)
    end if
    close(unitNum)     
    
    !Now that we have built the different arrays corresponding to the angle, the axis and the translation,
    !we can differentiate those quantities
    if (angleFlag==1) then
       call differentiateArrayReal(nonPerData%angle,nonPerData%tau,nonPerData%dangle)
       call differentiateArrayReal(nonPerData%dangle,nonPerData%tau,nonPerData%ddangle)  
    end if
    if (translationFlag==1) then
       call differentiateArrayVector(nonPerData%trans,nonPerData%tau,nonPerData%dtrans)
       call differentiateArrayVector(nonPerData%dtrans,nonPerData%tau,nonPerData%ddtrans)
    end if
  end subroutine createNonPeriodicCB
  
  
  subroutine BuildFlightPathProfile(nonPerData, unitNumber, flightPathProfile, iB, nbWayPoints)
    implicit none
    type(nonPeriodicData):: nonPerData    
    integer:: unitNumber, flightPathProfile, iB, nbWayPoints

    integer:: i
    real:: tauMin
    type(vector):: yStart, VH, AH    
    
    namelist / WP /  tauMin, yStart, VH    
    
    nonPerData%iB = iB   
    select case (flightPathProfile)
      case(CUSTOM)
        nonPerData%nTau = nbWaypoints
      case(TAKEOFF)
        nonPerData%nTau = 4
      case(OVERFLIGHT)
        nonPerData%nTau = 3
      case(APPROACH)
        nonPerData%nTau = 3
    end select
    
    if (nonPerData%nTau .ne. nbWayPoints) then
      call Error('The number of input waypoints is not the same as the number',&
                 'of waypoints in the requested flight path profile. Stopping')
    end if
    
    allocate(nonPerData%tau(nonPerData%nTau))     ;nonPerData%tau     = -huge(tauMin)
    allocate(nonPerData%angle(nonPerData%nTau))   ;nonPerData%angle   = 0.0
    allocate(nonPerData%dangle(nonPerData%nTau))  ;nonPerData%dangle  = 0.0
    allocate(nonPerData%ddangle(nonPerData%nTau)) ;nonPerData%angle   = 0.0
    allocate(nonPerData%axis(nonPerData%nTau))    ;nonPerData%angle   = 0.0
    allocate(nonPerData%trans(nonPerData%nTau))   ;nonPerData%trans   = vectorSetCoordinates(0.0, 0.0, 0.0)
    allocate(nonPerData%dtrans(nonPerData%nTau))  ;nonPerData%dtrans  = vectorSetCoordinates(0.0, 0.0, 0.0)
    allocate(nonPerData%ddtrans(nonPerData%nTau)) ;nonPerData%ddtrans = vectorSetCoordinates(0.0, 0.0, 0.0)
      
    select case (flightPathProfile)
      case(TAKEOFF)
        nonPerData%tau(1)     = 0.0      
        nonPerData%trans(1)   = vectorSetCoordinates(2000.0  , 0.0, 20.0)
        nonPerData%dtrans(1)  = vectorSetCoordinates(30.867  , 0.0, 0.0)    !60kts
        nonPerData%ddtrans(1) = vectorSetCoordinates(0.0     , 0.0, 0.0)
        
        nonPerData%trans(2)   = vectorSetCoordinates(500.0   , 0.0, 20.0)   !Level flight
        nonPerData%dtrans(2)  = vectorSetCoordinates(30.867  , 0.0, 0.0)    !60kts
        nonPerData%ddtrans(2) = vectorSetCoordinates(0.0     , 0.0, 0.0)
        
        nonPerData%trans(3)   = vectorSetCoordinates(0.0     , 0.0, 153.9746) !15deg angle of climb
        nonPerData%dtrans(3)  = vectorSetCoordinates(30.867  , 0.0, 0.0)      !60kts
        nonPerData%ddtrans(3) = vectorSetCoordinates(0.0     , 0.0, 0.0)
        
        nonPerData%trans(4)   = vectorSetCoordinates(-1000.0 , 0.0, 421.9238) !15deg angle of climb
        nonPerData%dtrans(4)  = vectorSetCoordinates(30.867  , 0.0, 0.0)      !60kts
        nonPerData%ddtrans(4) = vectorSetCoordinates(0.0     , 0.0, 0.0)
      case(OVERFLIGHT)
        nonPerData%tau(1)     = 0.0
        nonPerData%trans(1)   = vectorSetCoordinates(2000.0  , 0.0, 150.0)
        nonPerData%dtrans(1)  = vectorSetCoordinates(30.867  , 0.0, 0.0)      !60kts
        nonPerData%ddtrans(1) = vectorSetCoordinates(0.0     , 0.0, 0.0)
      
        nonPerData%trans(2)   = vectorSetCoordinates(0.0     , 0.0, 150.0)    !Level flight
        nonPerData%dtrans(2)  = vectorSetCoordinates(30.867  , 0.0, 0.0)      !60kts
        nonPerData%ddtrans(2) = vectorSetCoordinates(0.0     , 0.0, 0.0)
      
        nonPerData%trans(3)   = vectorSetCoordinates(-1000.0 , 0.0, 150.0)    !Level fligh
        nonPerData%dtrans(3)  = vectorSetCoordinates(30.867  , 0.0, 0.0)      !60kts
        nonPerData%ddtrans(3) = vectorSetCoordinates(0.0     , 0.0, 0.0)
      case(APPROACH)
        nonPerData%tau(1)     = 0.0
        nonPerData%trans(1)   = vectorSetCoordinates(2000.0  , 0.0, 330.208)
        nonPerData%dtrans(1)  = vectorSetCoordinates(30.867  , 0.0, 0.0)      !60kts
        nonPerData%ddtrans(1) = vectorSetCoordinates(0.0     , 0.0, 0.0)
      
        nonPerData%trans(2)   = vectorSetCoordinates(0.0     , 0.0, 120.0)    !6deg angle of descent
        nonPerData%dtrans(2)  = vectorSetCoordinates(30.867  , 0.0, 0.0)      !60kts
        nonPerData%ddtrans(2) = vectorSetCoordinates(0.0     , 0.0, 0.0)
      
        nonPerData%trans(3)   = vectorSetCoordinates(-1141.72, 0.0, 0.0)      !6deg angle of descent
        nonPerData%dtrans(3)  = vectorSetCoordinates(30.867  , 0.0, 0.0)      !60kts
        nonPerData%ddtrans(3) = vectorSetCoordinates(0.0     , 0.0, 0.0)
    end select      
    
    do i=1,nbWayPoints
      ! default values
      tauMin   = -huge(tauMin)
      yStart%A = -huge(tauMin)
      VH%A     = -huge(tauMin)
      read(unitNumber, nml=WP)
      if (tauMin.ne.-huge(tauMin)) then
        nonPerData%tau(i) = tauMin
      else if (i.eq.1) then
        nonperData%tau(i) = 0.0
      end if
      if (yStart%A(1).ne.-huge(tauMin)) then
        nonPerData%trans(i) = yStart
      end if
      if (VH%A(1).ne.-huge(tauMin)) then
        nonPerData%dtrans(i)%A  = VH%A
      end if
      nonPerData%ddtrans(i) = vectorSetCoordinates(0.0, 0.0, 0.0)
    end do
    do i=1,nbWaypoints-1
      if (vectorAbsolute(nonPerData%dtrans(i)).eq.0.0) then
        if(vectorAbsolute(nonPerData%trans(i+1)-nonPerData%trans(i)).ne.0.0) then
          call Error('The vehicle has zero velocity and thus cannot move between waypoints #' &
              //trim(integertostring(i))//' and #'//trim(integertostring(i+1))//'.')
        else if (nonPerData%tau(i+1).eq.-huge(tauMin)) then
          call Error('The vehicle is stationary between waypoints #'//trim(integertostring(i))// &
              ' and #'//trim(integertostring(i+1))//' but no time information was given.')
        end if
      else if(vectorAbsolute(nonPerData%trans(i+1)-nonPerData%trans(i)).eq.0.0) then
        call Error('The vehicle has non-zero velocity but the positions of waypoints #'// &
            trim(integertostring(i))//' and #'//trim(integertostring(i+1))//' are the same.')
      end if
    end do
    
  end subroutine BuildFlightPathProfile
  
  
  
  !**************************************************************
  !SUBROUTINE BuildFlightPath(Title, iB, wayPoints)
  !  This subroutine creates a non-periodic CB with the same 
  !  result as createNonPeriodicCB except that this CB is built
  !  within PSU-WOPWOP from waypoint data given by the user.
  !**************************************************************
  
  subroutine buildFlightPath(nonPerData, flightPathMotion)
    implicit none
    type(nonPeriodicData)::nonPerData
    type(functionData):: funcData
    integer, intent(in)::flightPathMotion
    
    integer:: i, j
    real:: dt
    real, dimension(:), pointer:: speed, userTau
    type(vector):: yHat, vHat, distance, dx, vel
      
    allocate(speed(nonPerData%nTau), userTau(nonPerData%nTau))
    userTau = nonPerData%tau
    
    distance = nonPerData%trans(nonPerData%nTau) - nonPerData%trans(1)
    do i=1,nonPerData%nTau-1
      dt = userTau(i+1) - UserTau(i)
      dx = nonPerData%trans(i+1)-nonPerData%trans(i)
      speed(i) = vectorAbsolute(nonPerData%dtrans(i))
      if (vectorAbsolute(dx).eq.0.0 .and. speed(i).eq.0.0 .and. dt.ne.0.0 .and. nonPerData%tau(i+1).ne.-huge(dt)) then
          if (userTau(i).gt.1.001*nonPerData%tau(i) .or. userTau(i).lt.0.999*nonPerData%tau(i)) then
            ! The time between waypoints is not correct. Use userTau(i+1)-userTau(i)
            call Notice('The vehicle is stationary.  The time between waypoints #'//trim(integertostring(i))// &
                ' and #'//trim(integertostring(i+1))//'',&
                'will not be modified but the actual start times will be set to ',&
                't('//trim(integertostring(i))//') = '//trim(realtostring(nonPerData%tau(i)))//'',&
                'and t('//trim(integertostring(i+1))//') = '// &
                trim(realtostring(nonPerData%tau(i)+userTau(i+1)-userTau(i)))//'')
          else
            call Notice('The vehicle is stationary between waypoints #'//trim(integertostring(i))//' and #'// &
                trim(integertostring(i+1))//'.')
          end if
      else if (vectorAbsolute(dx).eq.0.0 .and. speed(i).eq.0.0 .and. dt.eq.0.0) then
        call Error('The vehicle is stationary between waypoints #'//trim(integertostring(i))//' and #'// &
            trim(integertostring(i+1))//' and improper time information was given.')
      else
        vel = nonPerData%dtrans(i)
        ! To be sure that the velocity is correct we're going to use the 
        ! magnitude of the given velocity, multiplied by the unit vector between each waypoint.
        vHat = vectorRealDivide(dx, vectorAbsolute(dx))
        nonPerData%dtrans(i) = realVectorMultiply(vectorAbsolute(vel), vHat)
        ! We need to let the user know if we've changed the velocity vector.
        do j=1,3
          if (abs(nonPerData%dtrans(i)%A(j)-vel%A(j)).gt. 0.0001) then
            yHat = vectorRealDivide(vel, speed(i))
            call Notice('Travel from waypoint '//trim(integertostring(i))//' to waypoint '// &
                trim(integertostring(i+1))//' is not', &
                'possible with the given velocity.  The speed' ,&
                'will be kept the same but will be applied along the' ,&
                'unit vector between the two waypoints.')
            exit        
          end if
        end do      
        do j=1,3
          if (flightPathMotion.eq.CONSTANT_ACCELERATION .and. dx%A(j).ne.0.0) then
            nonPerData%ddtrans(i)%A(j) = (speed(i+1)**2.0 - speed(i)**2.0)*(vHat%A(j)**2.0)/(4.0*dx%A(j))
          end if
        end do
        dt = solveQuadratic(vectorAbsolute(nonPerData%ddtrans(i)), speed(i), -vectorAbsolute(dx))
      end if
      nonPerData%tau(i+1) = nonPerData%tau(i)+dt
    end do    
    deallocate(speed, userTau)
  
  end subroutine buildFlightPath
  
  
  function SolveQuadratic(A, B, C) result(dt)
    implicit none
    real:: A, B, C, temp, dt
  
    dt = -huge(dt)
    if (A.ne.0.0) then
      temp = B**2.0-4.0*A*C    
      if (temp.ge.0.0) then
        dt = (-B-sqrt(temp))/(2.0*A)
        if (dt.le.0.0) then
          dt  = (-B+sqrt(temp))/(2.0*A)
        end if
      else
        call Error('The given velocity and acceleration do not produce',&
                   'a real solution of time from the quadratic equation.')
      end if
    else if (A.eq.0.0 .and. B.ne.0.0) then
      dt = abs(C/B)
    end if
  
  end function SolveQuadratic
  

  !************************************************************************************
  !SUBROUTINE buildWindFrame(CBData,baseNum, unitNumber)
  !   The subroutine WindFrame is the input file that gives the rotation axis
  !  with respect to the horizon axis and the angle of rotation.
  !  We create this frame such as xw is tangent to the velocity of the helicopter.
  !ARGUEMENTS:
  !  - CBData: The CB that contains the non periodic data that cooresponds to 
  !             the wind base change
  !  - baseNum: The number of the base that cooresponds to the wind base change.
  !  - unitNumber: The number for the namelist, only used as a number.
  !************************************************************************************
  subroutine buildWindFrame(CBData, previousCB) 
    type(CBStructure), intent(inout)::CBData, previousCB

    type(vector)::xw,xh,yh
    integer::i,nTau
    type(vector),dimension(:),allocatable::zw 
    real,dimension(:),allocatable::Bangle  !KSB -angle used for array & function.
    real::abs

    nTau=CBData%nonPeriodicCB%nTau
    allocate(zw(nTau))
    allocate(Bangle(nTau))
    !open(unitNumber+2, file=CBData%filename)
    !write(unitNumber+2,*) baseNum
    !write(unitNumber+2,*) nTau
    !write(unitNumber+2,*) 1, 0, 1
    xh%A(1)=1.0
    xh%A(2)=0.0
    xh%A(3)=0.0
    yh%A(1)=0.0
    yh%A(2)=1.0
    yh%A(3)=0.0
    ! We have decided that the ground base has an axis
    ! z directed toward the center of the earth, such as
    ! the horizon base. Thus xw and xh are expresssed in the
    ! same base.
    ! If we want to define it otherwise we would have to change
    ! the definition of xw so that we put it in the horizon base.
    do i=1, nTau
       abs=(previousCB%nonPeriodicCB%dtrans(i)*previousCB%nonPeriodicCB%dtrans(i))**0.5
       abs=1.0/abs
       xw=abs*previousCB%nonPeriodicCB%dtrans(i)
       if (xw*xh/=1.0) then
          zw(i)=xh.cross.xw
          !     zw(i)=1.0/((zw(i)*zw(i))**0.5)*zw(i) 
          !     if(xw*yh>=0) then
          Bangle(i)=acos(xh*xw)
          !     else
          !       Bangle(i)=-acos(xh*xw)
          !     end if
          !     if (i<10) then
          !       call PrintVector(xw)
          !     end if  
       end if
       if (xw*xh==1.0) then    
          zw(i)%A(1)=0.0
          zw(i)%A(2)=0.0
          zw(i)%A(3)=1.0
          Bangle(i)=0.0   
       end if
    end do
    !write(unitNumber+2,*) (previousCB%nonPeriodicCB%tau(j), j=1, nTau)
    !write(unitNumber+2,*) (Bangle(j), j=1,nTau)
    !write(unitNumber+2,*) (zw(j)%A(1), j=1, nTau)
    !ERROR write(unitNumber+2,*) (zw(i)%A(1), j=1, nTau)
    !write(unitNumber+2,*) (zw(j)%A(2), j=1, nTau)
    !write(unitNumber+2,*) (zw(j)%A(3), j=1, nTau)         
    !close(unitNumber+2) 
    CBData%nonPeriodicCB%tau  =previousCB%nonPeriodicCB%tau
    CBData%nonPeriodicCB%angle=Bangle
    CBData%nonPeriodicCB%axis =zw
    

  end subroutine BuildWindFrame


  !*******************************************************
  ! Here follows the functions we use in matrot to give the values
  ! of a change of base at the time tau.
  ! These functions depend on the input. We have used an
  ! object approach such that the input data carries with itself
  ! the information necessary to realize the different functions.
  ! Thus we won't have to change those functions when we change the 
  ! type of CB in the input.
  !********************************************************


  !********************************************************
  ! The function angle aims at reading the angle of the
  ! rotation of the frame i into the frame i-1 at the 
  ! time tau.
  !********************************************************
  recursive function angle(CBData, tau, index) result(res)
    implicit none
    integer :: index, pos
    real :: tau, res, psi
    type(CBStructure)::CBData

    select case(CBData%angletype)
    case(APERIODIC)
       res=interpolateReal(tau, CBData%nonPeriodicCB%tau, &
                           CBData%nonPeriodicCB%angle, index)  
    case(TIME_INDEPENDENT)
       res=CBData%anglevalue
    case(KNOWN_FUNCTION)
       res=CBData%funcData%OMEGA*tau+CBData%funcData%PSI0
    case(PERIODIC)
       psi=angle(rotationCB, tau, index)+CBData%funcData%PSI0
       res=(CBData%funcData%A0)
       do pos=1,size(CBData%funcData%A)
         res = res - &
            CBData%funcData%A(pos)*cos(pos*psi)  - CBData%funcData%B(pos)*sin(pos*psi)
       end do
    case default
       write (*,*) "ERROR: ", CBData%angletype, " is not a valid angle type."
       stop
    end select

  end function angle

  !*************************************************************
  ! The function axis aims at reading what is the axis of the
  ! rotation of the frame i into the frame i-1 at the 
  ! time tau.     
  !**************************************************************
  function axis(CBData, tau, index) result(res)
    implicit none
    integer :: index
    real :: tau
    type(vector) :: res
    type(CBStructure)::CBData

    select case(CBData%axistype)
      case(APERIODIC)
         res=interpolateVector(tau, CBData%nonPeriodicCB%tau, &
                               CBData%nonPeriodicCB%axis, index)
      case(TIME_INDEPENDENT)
       res=CBData%axisvalue
      case(KNOWN_FUNCTION)
       write(*,*) 'ERROR: KnownFunction COB for axis not written yet.'  
       stop
      case(PERIODIC)
       write(*,*) 'ERROR: Periodic COB for axis not written yet.'  
       stop
      case default
         write (*,*) "ERROR: ", CBData%axistype, " is not a valid axis type."
         stop
      end select     

    end function axis

    !***********************************************************
    ! The function trans aims at reading what is the  
    ! translation of the frame i into the frame i-1 at the 
    ! time tau.  
    !***********************************************************
    function readTranslation(CBData, tau, tau0, index) result(res)
      implicit none
      integer :: index
      real :: tau, tau0
      type(vector) :: res
      type(CBStructure) :: CBData

      select case(CBData%translationtype)
      case(APERIODIC)
        res=interpolateVector(tau, CBData%nonPeriodicCB%tau, &
                             CBData%nonPeriodicCB%trans, index)
      case(TIME_INDEPENDENT)
        res=CBData%translationvalue
      case(KNOWN_FUNCTION)
        res=(tau-tau0)**2*CBData%funcData%AH+(tau-tau0)*CBData%funcData%VH+CBData%funcData%Y0
      case(PERIODIC) 
        write(*,*) 'ERROR: Periodic COB for translation not written yet.'  
        stop
      case default
        write (*,*) "ERROR: ", CBData%translationtype, " is not a valid translation type."
        stop
      end select

    end function readTranslation

    !*****************************************************
    ! The function dangle aims at reading what is the  
    ! angular speed of the frame i into the frame i-1 at the 
    ! time tau.  
    !******************************************************
    recursive function dangle(CBData, tau, index) result(res)
      implicit none
      type(CBStructure)::CBData
      integer :: index,pos
      real :: tau, res 
      real :: psi
      real :: dpsi

      select case(CBData%angletype)
      case(APERIODIC)
         res=interpolateReal(tau, CBData%nonPeriodicCB%tau, &
                             CBData%nonPeriodicCB%dangle, index)  
      case(TIME_INDEPENDENT)
         res=0
      case(KNOWN_FUNCTION)
         res=CBData%funcData%OMEGA
      case(PERIODIC)
         psi=angle(rotationCB, tau, index)+CBData%funcData%PSI0
         dpsi=dangle(rotationCB, tau, index)
         res = 0.0
         do pos=1,size(CBData%funcData%A)
           res=res + (dpsi*pos*(CBData%funcData%A(pos)*sin(pos*psi)- &
                                CBData%funcData%B(pos)*cos(pos*psi)))
         end do
      case default
         write (*,*) "ERROR: ", CBData%angletype, " is not a valid angle type."
         stop 
      end select

    end function dangle

    !************************************************************
    ! The function ReadVelocity aims at reading what is the  
    ! velocity of the origin of the frame i into the frame i-1 at the 
    ! time tau.
    !**************************************************************
    function readVelocity(CBData,tau, tau0,index) result(res)
      implicit none
      type(CBStructure)::CBData
      integer :: index
      real :: tau,tau0
      type(vector) :: res

      select case(CBData%translationtype)
      case(APERIODIC)
         res=interpolateVector(tau, CBData%nonPeriodicCB%tau, &
                               CBData%nonPeriodicCB%dtrans, index)   
      case(TIME_INDEPENDENT)
         res%A(1)=0.
         res%A(2)=0.
         res%A(3)=0.
      case(KNOWN_FUNCTION)
         res=(tau-tau0)*CBData%funcData%AH+CBData%funcData%VH
      case(PERIODIC)
         write(*,*) 'ERROR: Periodic COB for translation not written yet.'  
         stop
      case default
         write (*,*) "ERROR: ", CBData%translationtype, " is not a valid translation type."
         stop    
      end select

    end function readVelocity

    !***********************************************************
    ! The function ddangle aims at reading what is the  
    ! angular acceleration of the frame i into the frame i-1 at the 
    ! time tau.
    !*************************************************************
    function ddangle(CBData,tau, index) result(res)
      implicit none
      integer :: index, pos
      real :: tau, res, psi, dpsi 
      type(CBStructure)::CBData 

      select case(CBData%angletype)
      case(APERIODIC)
         res=interpolateReal(tau,CBData%nonPeriodicCB%tau, &
                             CBData%nonPeriodicCB%ddangle, index)
      case(TIME_INDEPENDENT)
         res=0
      case(KNOWN_FUNCTION)
         res=0
      case(PERIODIC)
         psi=angle(rotationCB, tau, index)+CBData%funcData%PSI0
         dpsi=dangle(rotationCB, tau, index)
         res = 0.0
         do pos=1,size(CBData%funcData%A)
           res=res + (dpsi**2*pos**2*&
               (CBData%funcData%A(pos)*cos(pos*psi)+&
                CBData%funcData%B(pos)*sin(pos*psi)))
         end do
      case default
         write (*,*) "ERROR: ", CBData%angletype, " is not a valid angle type."
         stop
      end select

    end function ddangle

    !**************************************************************
    ! The function accelerationOi aims at reading what is the  
    ! acceleration of the origin of the frame i into the frame i-1 at the 
    ! time tau.
    !***************************************************************
    function readacceleration(CBData,tau, index) result(res)
      implicit none
      integer :: index
      real :: tau
      type (vector) :: res  
      type(CBStructure)::CBData 

      select case(CBData%translationtype)
      case(APERIODIC)
         res=interpolateVector(tau, CBData%nonPeriodicCB%tau, &
                               CBData%nonPeriodicCB%ddtrans, index)
      case(TIME_INDEPENDENT)
         res%A(1)=0
         res%A(2)=0
         res%A(3)=0
      case(KNOWN_FUNCTION)
         res=CBData%funcData%AH
      case(PERIODIC)
         write(*,*) 'ERROR: Periodic COB for translation not written yet.'  
         stop
      case default
         write (*,*) "ERROR: ", CBData%translationtype, " is not a valid translation type."
         stop
      end select

    end function readacceleration


    !! A function to determine if the only motion is uniform translation
    !function OnlyUniformlyTranslating (CBData) result(Result)
    !  implicit none
    !  type(CBStructure), intent(in) :: CBData
    !  logical :: Result ! OnlyUniformlyTranslating
    !
    !  Result = .false.
    !  ! Organize these from most likely to least likely
    !  if (associated(CBData%next)) then
    !    ! If there is more than one CB, we call the result false. Someday maybe
    !    ! allow multiple translations to be catenated together. For now, call it
    !    ! a maneuvering observer.
    !    Result = .false.
    !  else if (CBdata%angleValue /= 0.0) then
    !    Result = .false.
    !  else if (associated(CBData%nonPeriodicCB%tau)) then
    !    Result = .false.
    !  else if (any(CBData%funcData%AH%A /= 0)) then
    !    Result = .false.
    !  else if (any(CBData%funcData%VH%A /= 0)) then
    !    Result = .true.
    !  end if
    !end function OnlyUniformlyTranslating
    
    ! A function to determine if the only motion is uniform translation
    ! Modified 10/2/2017 KSB.  The intent is that even if multiple COB are used for the observer
    ! but only translation (no acceleration or rotation) is used, then there can be more than one
    ! CB.  This is true if there are TimeIndependent rotations.  This seems to work, but we should
    ! keep our eye out for potential problems in the future, just to be sure.
    !
    recursive function OnlyUniformlyTranslating (CBData) result(Result)
      implicit none
      type(CBStructure), intent(in) :: CBData
      logical :: Result ! OnlyUniformlyTranslating

      Result = .false.

      if (associated(CBData%next)) then
        ! If there is more than one CB, will search through the CB list and 
        ! make sure all the CB's are only uniformly translating (which likely
        ! means they might be TimeIndependent and not translating at all, but
        ! this is zero translation).
        Result = OnlyUniformlyTranslating(CBData%next)
      else if ((CBdata%axisType == NO_VALUE .or. CBdata%axisType == TIME_INDEPENDENT) .and. &
               (CBdata%angleType== NO_VALUE .or. CBdata%angleType== TIME_INDEPENDENT) ) then
        if ( (CBdata%translationType == TIME_INDEPENDENT) .or. &
             (CBdata%translationType == KNOWN_FUNCTION) ) then
          if (all(CBData%funcData%AH%A == 0)) then
            Result = .true.
          end if
        end if
      end if
    end function OnlyUniformlyTranslating


    subroutine GetUniformTranslation (CBData, pos, vel)
      implicit none
      type(CBStructure), intent(in) :: CBData
      type(vector), intent(out) :: pos, vel

      pos = CBData%funcData%Y0
      vel = CBData%funcData%VH
      
    end subroutine GetUniformTranslation

    recursive function TrueRotation(CBData) result(Rotation)
      implicit none
      type(CBStructure), intent(in):: CBData
      logical:: Rotation
      
      Rotation=.false.
      if (CBData%rotation) then
         Rotation=.true.
         return
      else if (associated(CBData%next)) then
        Rotation = TrueRotation(CBData%next)
      end if
    end function TrueRotation

  end module COBObject