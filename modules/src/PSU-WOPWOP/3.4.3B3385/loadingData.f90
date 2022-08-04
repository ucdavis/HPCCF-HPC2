! PSU-WOPWOP
! $Id: loadingData.f90 3366 2016-10-31 17:16:11Z brentner $
! $LastChangedDate: 2016-10-31 13:16:11 -0400 (Mon, 31 Oct 2016) $
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


!****************************************************************************
!
! MODULE 
!   loadingObject
! 
! DESCRIPTION
!   This module defines a class called "loading" with its constructors,
!   destructors and access functions. This class reads in the loading from a
!   file and calculates its derivatives. It functions as a black box to the
!   patch object, hiding the type of loading from it (i.e. constant, periodic
!   and aperiodic). The box isn't perfect - sometimes patch needs to know
!   anyway, but for the most part it is opaque.
!
! TYPES
!   loading
!   SingleTimeLoad
!
! PUBLIC SUBROUTINES
!  CreateLoading (newLoading, loadFileUnit, dataFormat, timeOffset)
!  CopyLoading (newLoading, oldLoading, timeOffset)
!  DestroyLoading (oldLoading)
!  DestroySingleLoad (nKey, nKeyD)
!  CreateLoadingData(newLoading, loadFileUnit, loadInfo(LOAD_TIMETYPE))
!  ReadLoadingDimensions(newLoading, loadFileUnit, loadInfo(LOAD_TIMETYPE))
!  ReadConstantLoadingAtTime(newLoading,streamNumber,eLoadFlag)
!  ReadPeriodicLoadingAtTime(newLoading,streamNumber,keyOffset,eLoadFlag, iTau)
!  ReadAperiodicLoadingAtTime(newLoading,streamNumber,eLoadFlag, iTau)
! 
! PUBLIC FUNCTIONS
!  HasImplicitTau (ld) result (implied)
!  GetLoading (loadObject, time) result (returnedLoad)
!  GetLoadingDerivative (loadObject, time) result (returnedDerivative)
!  GetSourceTime (loadObject, iTau) result (sourceTime)
!  GetNTau (loadObject) result(nTau)
!  GetLoadType (loads) result (loadType)
!  GetLoadingSizeInMemory(loads) result(res)
!  GetLoadSpecies(loadObject) result(loadSpecies)
!  GetLoadnkey(loadObject)
!
! OTHER PUBLIC DATA
!  Integer parameters (constants) describing the type of loading:
!    NO_LOADING
!    CONSTANT_LOADING
!    PERIODIC_LOADING
!    APERIODIC_LOADING
! 
! HISTORY
!   Created 7/10/03 by Chris Hennes
!
!****************************************************************************

module loadingObject
  use constantsModule, only: LOAD_GRIDTYPE, RHO, C, GAMMA, VERSION_1, &
    LOAD_TIMETYPE, LOAD_TIMETYPE_CONSTANT, LOAD_TIMETYPE_PERIODIC, LOAD_TIMETYPE_APERIODIC, &
    LOAD_DATATYPE, LOAD_DATATYPE_PRESSURE, LOAD_DATATYPE_LOADING, LOAD_DATATYPE_FLOW, &
    APERIODICARRAYSIZE, LOAD_TIMETYPE_QPERIODIC
  use mathModule, only: vector, VectorSetValue, operator(*), operator(-), operator(+), vectorSetCoordinates, vectorSubtract
  use IOModule, only: isGood, ReadBinaryReal, ReadBinaryInteger, ReadBinaryReal1DArray, ReadBinaryReal2DArray
  use debugModule
  use MPIModule
  use strings, only: RealtoString, IntegertoString
  implicit none

  private
  public:: loading, SingleTimeLoad, ReadLoadingAtTime
  public:: ReadConstantLoadingAtTime,ReadPeriodicLoadingAtTime,ReadAperiodicLoadingAtTime
  public:: CreateLoading, CopyLoading, GetLoading, DestroyLoading, DestroySingleLoad
  public:: GetLoadingDerivative, GetLoadType, GetLoadingTauRange
  public:: GetLoadingDt, LoadingHasImplicitTau, GetLoadingNTau
  public:: GetLoadingSourceTime, GetLoadingOffset, GetLoadPeriod
  public:: CreateLoadingTauArray, WriteLoadDebugInfo
  public:: GetLoadSpecies, CreateLoadingData, ReadLoadingDimensions
  public:: GetLoadnKey, RotateLoadingArray, ResetLoadingKeyCount, SetLoadingNKey
  public:: NO_LOADING, CONSTANT_LOADING, PERIODIC_LOADING, APERIODIC_LOADING
  
  ! A type to store the complete loading vector at a single time
  type SingleTimeLoad
    ! The time for the load
    real::time
    real, dimension(:), pointer:: pressVals =>NULL()
    real, dimension(:), pointer:: densVals =>NULL()
    type(vector), dimension(:), pointer::loadVectors =>NULL()
    type(vector), dimension(:), pointer::momVectors =>NULL()
    ! These are calculated from the data when it is read
    type(vector), dimension(:), pointer::derivatives =>NULL()    
  end type SingleTimeLoad

  ! A type to store all of the loading at all times
  type loading
    private
    ! The load type is one of CONSTANT_LOADING, PERIODIC_LOADING or 
    ! APERIODIC_LOADING, as defined by the constants below. It could also be 
    ! NO_LOADING, but that is an error status.
    
    integer::loadType, loadSpecies
    
    integer::keycount,imax,jmax,nbnodes,gridType
    
    ! The period of the data, in seconds.
    real::period
    ! The first and last key values (for periodic data)
    real::keyMin, keyMax
    ! The differential time segment for the loading
    real::dt
    ! The time offset
    integer::timeStepOffset
    ! This is what is actually read in from the file:
    type(SingleTimeLoad), dimension(:), pointer::loads
    ! This is a storage space for a zero-load vector to return when an error
    ! occurs or the data ends.
    type(SingleTimeLoad),pointer :: zeroLoad =>NULL()
    !real, dimension(:,:),pointer::temp2DArray 
    ! The three variables are used for MPI communication
    integer::nkey
    logical::isCopy
  end type loading
  
  ! Some constants:
  integer, parameter:: NO_LOADING          = 0, &
                       CONSTANT_LOADING    = 1, &
                       PERIODIC_LOADING    = 2, &
                       APERIODIC_LOADING   = 3, &
                       CYLINDER_OBJECT     = 4, &
                       TRAILINGEDGE_OBJECT = 5

  integer, parameter:: STRUCTURED_GRID     = 1, &
                       UNSTRUCTURED_GRID   = 2
contains

  subroutine WriteSingleTimeLoadDebugInfo(loads,unitnum)
    type(SingleTimeLoad), intent(in)::loads
    integer::unitnum
    write(unitnum,*) '*** SingleTimeLoadDebugInfo ***'
    write(unitnum,*) 'time= ', trim(realtostring(loads%time))
    write(unitnum,*) 'associated pressure values? ', associated(loads%pressVals)
    if(associated(loads%pressVals)) then
      write(unitnum,*) 'Size of vectors', trim(integertostring(size(loads%pressVals)))
      write(unitnum,*) 'First value=',loads%pressVals(1)
      write(unitnum,*) 'Last value=',loads%pressVals(size(loads%pressVals))
    end if  
    write(unitnum,*) 'associated density values? ', associated(loads%densVals)
    if(associated(loads%densvals)) then
      write(unitnum,*) 'Size of vectors', trim(integertostring(size(loads%densVals)))
      write(unitnum,*) 'First value=',loads%densVals(1)
      write(unitnum,*) 'Last value=',loads%densVals(size(loads%densVals))
    end if      
    write(unitnum,*) 'associated load vectors? ', associated(loads%loadVectors)
    if(associated(loads%loadVectors)) then
      write(unitnum,*) 'Size of vectors', trim(integertostring(size(loads%loadVectors)))
      write(unitnum,*) 'First vector=',loads%loadVectors(1)
      write(unitnum,*) 'Last vector=',loads%loadVectors(size(loads%loadVectors))
    end if
    write(unitnum,*) 'associated momentum vectors? ', associated(loads%momVectors)
    if(associated(loads%momVectors)) then
      write(unitnum,*) 'Size of vectors', trim(integertostring(size(loads%momVectors)))
      write(unitnum,*) 'First vector=',loads%momVectors(1)
      write(unitnum,*) 'Last vector=',loads%momVectors(size(loads%momVectors))
    end if    
    write(unitnum,*) ' associated derivatives? ', associated(loads%derivatives)
    if(associated(loads%derivatives)) then
      write(unitnum,*) 'size of derivatives ', trIM(integertostring(size(loads%derivatives)))
      write(unitnum,*) 'First derivatives follows'
      write(unitnum,*) loads%derivatives(1)
      write(unitnum,*) 'Last derivatives follows'
      write(unitnum,*) loads%derivatives(size(loads%derivatives))
    end if    
    write(unitnum,*) '--- End SingleTimeLoadDebugInfo ---'

  end subroutine WriteSingleTimeLoadDebugInfo

  subroutine WriteLoadDebugInfo(load,unitnum)
    type(loading), intent(in)::load
    integer::unitnum
    write(unitnum,*) '*** LoadingDebugInfo ***'
    write(unitnum,*) 'loadType= ', trim(integertostring(load%loadType))
    write(unitnum,*) 'keyCount= ', trim(integertostring(load%keyCount))
    write(unitnum,*) 'iMax= ', trim(integertostring(load%iMax))
    write(unitnum,*) 'jMax= ', trim(integertostring(load%jMax))
    write(unitnum,*) 'nbNodes= ', trim(integertostring(load%nbNodes))
    write(unitnum,*) 'gridType= ', trim(integertostring(load%gridType))
    write(unitnum,*) 'period= ', trim(realtostring(load%period))
    write(unitnum,*) 'keyMax= ', trim(realtostring(load%keyMax))
    write(unitnum,*) 'dt= ', trim(realtostring(load%dt))
    write(unitnum,*) 'timeStepOffset= ', trim(integertostring(load%timeStepOffset))
    write(unitnum,*) 'nKey=', trim(integertostring(load%nkey))
    write(unitnum,*) 'isCopy=', load%isCopy
    write(unitnum,*) 'associated zeroLoad? ', associated(load%zeroLoad)
    if(associated(load%zeroLoad)) then
      write(unitnum,*) 'zero load'
      call WriteSingleTimeLoadDebugInfo(load%zeroLoad,unitnum)
    end if
    write(unitnum,*) 'loads associated? ', associated(load%loads)
    if(associated(load%loads)) then
      write(unitnum,*) 'size of loads ', trim(integertostring(size(load%loads)))
      write(unitnum,*) 'First loads follows'
      call WriteSingleTimeLoadDebugInfo(load%loads(1),unitnum)
      write(unitnum,*) 'Last loads follows'
      call WriteSingleTimeLoadDebugInfo(load%loads(size(load%loads)),unitnum)
    end if
    write(unitnum,*) '--- End LoadDebugInfo ---'
    
  end subroutine WriteLoadDebugInfo
  
  

  ! A constructor. Creates the loading by reading it in from the file.
  ! Parameters:
  !   newLoading - The new loading object to create. Note that the memory must
  !                already exist for this object.
  !   loadFileUnit - The unit number of the file which contains the load
  !                  information. The file must already be open and the header
  !                  read in prior to reaching this function.
  !   loadInfo  - information from the loading file fixed header
  subroutine CreateLoading (newLoading, loadFileUnit,loadInfo)
    type(loading), intent(inout):: newLoading
    integer, intent(in) :: loadFileUnit
    integer,dimension(:),optional::loadInfo
    integer:: keyIndex
    
    newLoading%isCopy = .false.
    newLoading%keycount = 0    
    newLoading%gridtype = loadInfo(LOAD_GRIDTYPE)
    ! ReadLoadingDimensions read the period, nKey, and the size of the loading 
    ! data dimensions (iMax, jMax, nbNodes).  It saves these values in the 
    ! loading data structure.
    call ReadLoadingDimensions(newLoading, loadFileUnit, loadInfo)
    ! Allocate the memory and read the data
    if (loadInfo(LOAD_TIMETYPE) .ne. APERIODIC_LOADING) then
      allocate (newLoading%loads(newLoading%nkey))
    else
    ! The loading information is stored for aperiodicArraySize (nominally 4) time 
    ! segments: [L, C, R, RR] when working with aperiodic data.  We store one point 
    ! back in time for use with gradient calculations, while also storing two points
    ! forward in time for use in central difference equations and loading pressure 
    ! gradient calculations.
      allocate (newLoading%loads(aperiodicArraySize), stat=ierr)
      if (ierr >0) call error("Error in allocating single time loads arrays")
    end if
    
    newLoading%loadSpecies = loadInfo(LOAD_DATATYPE)
    allocate (newLoading%zeroLoad)
    select case (loadInfo(LOAD_DATATYPE))
      case (LOAD_DATATYPE_LOADING)
        ! Allocate the memory and read the data
        do keyIndex=1, size(newLoading%loads)
          allocate (newLoading%loads(keyIndex)%loadVectors(newLoading%nbNodes))
          allocate (newLoading%loads(keyIndex)%derivatives(newLoading%nbNodes)) 
        end do
        ! Create the zero-load vector 
        allocate (newLoading%zeroLoad%loadVectors(newLoading%nbNodes))
        newLoading%zeroLoad%loadVectors(:) = vectorSetCoordinates(0.,0.,0.)
      case (LOAD_DATATYPE_FLOW)
        ! Allocate the memory and read the data
        do keyIndex=1,size(newLoading%loads)
          allocate (newLoading%loads(keyIndex)%densVals(newLoading%nbNodes))
          allocate (newLoading%loads(keyIndex)%momVectors(newLoading%nbNodes))
          allocate (newLoading%loads(keyIndex)%pressVals(newLoading%nbNodes))
        end do
        ! Create the zero-load vector
        allocate (newLoading%zeroLoad%densVals(newLoading%nbNodes))
        allocate (newLoading%zeroLoad%momVectors(newLoading%nbNodes))
        allocate (newLoading%zeroLoad%pressVals(newLoading%nbNodes))
        newLoading%zeroLoad%densVals(:) = 1.
        newLoading%zeroLoad%momVectors(:) = vectorSetCoordinates(0.,0.,0.)
        newLoading%zeroLoad%pressVals(:) = 0.
      case (LOAD_DATATYPE_PRESSURE)
        ! Allocate the memory and read the data
        do keyIndex=1,size(newLoading%loads)
          allocate (newLoading%loads(keyIndex)%pressVals(newLoading%nbNodes))
        end do
        ! Create the zero-load vector
        allocate (newLoading%zeroLoad%pressVals(newLoading%nbNodes))
        newLoading%zeroLoad%pressVals(:) = 0.
    end select

    if (debuglevel > 14) then
      ! write out loading information for general code debugging.
      write(*,*) '*** Loading Data Information ***'
      write(*,*) 'Loadtype: ',newLoading%loadType, '(1 = constant 2 = periodic 3 = aperiodic)'
      write(*,*) 'Imax: ',newLoading%iMax
      write(*,*) 'Jmax: ',newLoading%jMax
      write(*,*) 'nbNodes: ',newLoading%nbNodes
      write(*,*) 'Grid type: ',newLoading%gridType , '(1 = structured 2 = unstructured)'
      
      select case (newLoading%loadType)
      case (PERIODIC_LOADING)
        write(*,*) 'Periodic Loading'
        write(*,*) 'period: ',newLoading%period
        write(*,*) 'nkey: ',newLoading%nKey
        write(*,*) 'dTau, ',newLoading%dt
        write(*,*) 'timeStepOffset: ',newLoading%timeStepOffset
      case (APERIODIC_LOADING)
        write(*,*) 'Aperiodic Loading'
        write(*,*) 'nkey: ',newLoading%nKey
        write(*,*) 'dTau not know yet for Aperiodic loading'
      case default  !Constant Loading
        write(*,*) 'Constant Loading'
      end select 
      write(*,*) '--- End Loading Data Information ---'
    end if
  
 end subroutine CreateLoading


 subroutine ReadLoadingDimensions (newLoading, streamNumber, loadInfo)
   type(loading), intent(inout):: newLoading
   integer,dimension(:),optional::loadInfo
   integer, intent(in):: streamNumber
   
   real(kind=4):: temp
   integer:: status
   
   select case (loadInfo(LOAD_TIMETYPE))
   case (CONSTANT_LOADING)
     newLoading%loadType = CONSTANT_LOADING       
     newLoading%period = 0.
     newLoading%nkey=1
     newLoading%dt = 0.0              
   case (PERIODIC_LOADING, LOAD_TIMETYPE_QPERIODIC)
     newLoading%loadType = PERIODIC_LOADING
     ! read in the period
     call ReadBinaryReal(streamNumber, newLoading%period, status)
     if (newLoading%period <= 0 .or. status /= 0) then
       call Error ("The period of the loading file is invalid.", &
                   "Period = "//trim(RealToString(newLoading%period)))
     end if

     if (newLoading%period > 1.e15) then
       call Error ("The period of the loading file is too large.", &
                   "Period = "//trim(RealToString(newLoading%period)))
     end if

     ! Read in the number of Key locations
     call ReadBinaryInteger(streamNumber,newLoading%nKey, status)
     if (newLoading%nKey < 0 .or. status /= 0) then
       call Error("ERROR: Invalid number of key locations in periodic loading file.")
     end if
     newLoading%dt = newLoading%Period/(newLoading%nKey-1)
   case (APERIODIC_LOADING)
     newLoading%loadType = APERIODIC_LOADING
     newLoading%period = 0.
     ! Read in the number of time steps (nKey)
     call ReadBinaryInteger(streamNumber, newLoading%nKey, status)
     
     if (newLoading%nKey < 0 .or. status /= 0) then
       call Error("Invalid number of key locations in aperiodic loading file.")
     end if
     ! We don't know newLoading%dt yet.  It should be computed in after the
     ! second aperiodic key (time) has been read in.
   case default
      write (*,*) "ERROR: An unknown internal error has occurred."
      write (*,*) "       Loading noise will be incorrect."
      newLoading%loadType = NO_LOADING
      newLoading%period = 0.
      newLoading%nKey = 0
      newLoading%dt = 0.
   end select
 
   call  ReadDimensions(newLoading, streamNumber) 
 end subroutine ReadLoadingDimensions     
    

  ! First step in getting the loading information.  Once within 
  ! this subroutine the load type (const/ periodic/ aperiodic)
  ! is determined at which time another subroutine is called
  ! which is built for that specific load type.
  subroutine ReadLoadingAtTime(newLoading, loadFileUnit, keyOffset, timeType)
    type(loading):: newLoading
    integer,intent(in):: loadFileUnit, timeType
    real::keyOffset

    newLoading%keyCount = newLoading%keyCount + 1
    if( newLoading%keyCount <= newLoading%nKey ) then
      select case (timeType)
        case (LOAD_TIMETYPE_CONSTANT)
          call ReadConstantLoadingAtTime(newLoading,loadFileUnit)
        case (LOAD_TIMETYPE_PERIODIC,LOAD_TIMETYPE_QPERIODIC)
          call ReadPeriodicLoadingAtTime(newLoading,loadFileUnit, keyOffset)
        case (LOAD_TIMETYPE_APERIODIC)
          call ReadAperiodicLoadingAtTime(newLoading,loadFileUnit)
        case default
          call Error ("Unknown time type in loading ("//&
                      trim(IntegerToString(timeType))//")")
      end select
    end if
  end subroutine ReadLoadingAtTime
  

  !Grabs the constant loading data
  subroutine ReadConstantLoadingAtTime(newLoading,streamNumber)
    type(loading), intent(inout):: newLoading
    integer, intent(in):: streamNumber
    integer:: i, status, loadSpecies
    real :: patm
    real, pointer::r
    real(kind=4), dimension(:), allocatable:: tempVals    
    real,         dimension(:,:), allocatable::temp2DArray
    real(kind=4), dimension(:,:), allocatable::temp2DArrayKind4
    ! In the case of constant loading this is only done once
    loadSpecies = newLoading%loadSpecies
    select case (loadSpecies)
      case (LOAD_DATATYPE_PRESSURE)
        patm = (rho*c*c)/gamma
        allocate(tempVals(size(newLoading%loads(1)%pressVals)))
        ! Read in the data:
        call ReadBinaryReal1DArray(streamNumber, tempVals, status)
        newLoading%loads(1)%pressVals = tempVals
        if (status.ne.0) then
          call Error("Unable to read constant pressVals data file.", &
                     "       Error in pressVals specification.")
        end if
        if (status.ne.0) then
          call Error("Unable to read pressVals data file.")
        end if
        deallocate(tempVals)
      case (LOAD_DATATYPE_LOADING)
        allocate (temp2DArray(newLoading%nbNodes,3))
        allocate (temp2DArrayKind4(newLoading%nbNodes,3))
        ! Read in the data:
        call ReadBinaryReal2DArray(streamNumber, temp2DArrayKind4, status)
        temp2DArray = temp2DArrayKind4
        if (status.ne.0) then
          if (status == 2) then
            write (*,*) "       The file ended prematurely."
          end if
          call Error("Unable to read constant loading data file. Error in loading specification.")
        end if
        do i=1, newLoading%nbNodes
          if (isGood(temp2DArray(i,1)) .and. &
              isGood(temp2DArray(i,2)) .and. &
              isGood(temp2DArray(i,3))) then
            ! loading vector set to -1 since internally we need force on the fluid
            ! NOT force on the surface (like in aerodynamics).
            newLoading%loads(1)%loadVectors(i)=-1.0*VectorSetValue(temp2DArray(i,:))
          else
            call Error ("Loading value at i="//trim(IntegerToString(i))//" is not a number ("//&
                        trim(RealToString(temp2DArray(i,1)))//", "//&
                        trim(RealToString(temp2DArray(i,2)))//",G4 "//&
                        trim(RealToString(temp2DArray(i,3)))//").")
          end if
        end do
        ! Create the time derivative (which is obviously zero in the constant loading
        ! case).
        newLoading%loads(1)%derivatives(:) = vectorSetCoordinates(0.,0.,0.)
        newLoading%dt = 0.0
        deallocate(temp2DArray, temp2DArrayKind4)
      case (LOAD_DATATYPE_FLOW)
        allocate (temp2DArray(newLoading%nbNodes,3))
        allocate (temp2DArrayKind4(newLoading%nbNodes,3))
        allocate (tempVals(size(newLoading%loads(1)%densVals)))          
        patm = (rho*c*c)/gamma
        ! Read in the data:
        call ReadBinaryReal1DArray(streamNumber, tempVals, status)
        newLoading%loads(1)%densVals = tempVals
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read constant load data file."
          write (*,*) "       Error in densVals specification."
          stop
        end if
        call ReadBinaryReal2DArray(streamNumber, temp2DArrayKind4, status)
        temp2DArray = temp2DArrayKind4
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read constant load data file."
          write (*,*) "       Error in momVectors specification."
          stop
        end if
        call ReadBinaryReal1DArray(streamNumber, tempVals, status)
        newLoading%loads(1)%pressVals = tempVals
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read constant load data file."
          write (*,*) "       Error in pressVals specification."
          stop
        end if
        do i=1,newLoading%nbNodes
          newLoading%loads(1)%momVectors(i)=VectorSetValue(temp2DArray(i,:))
          ! Assign some pointers to make the equation more understandable
          r  => newLoading%loads(1)%densVals(i)
          if (r <= 0) then
            write (*,*) "ERROR2: densVals must be greater than zero. (",i,") -> ", r
            stop
          end if
        end do
        ! I'm not sure this is helpful at this point, but I'll leave it (KSB)  
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read load data file."
          stop
        end if
        deallocate (temp2DArray, temp2DArrayKind4, tempVals)
        newLoading%timeStepOffset = 0.0
      case default
        write (*,*)  "ERROR:  An unknown internal error has occured."
     end select
  end subroutine ReadConstantLoadingAtTime
   
  subroutine ReadPeriodicLoadingAtTime(newLoading,streamNumber,keyOffset)
    type(loading), intent(inout):: newLoading
    integer, intent(in):: streamNumber
    real, intent(in) :: keyOffset
    integer:: i, keyIndex, status
    real(kind=4):: temp
    real :: keyMin, keyMax, dt, timeOffset, patm
    real, pointer::r
    real(kind=4), dimension(:), allocatable:: tempVals
    real(kind=4), dimension(:,:), allocatable:: temp2DArrayKind4
    real,         dimension(:,:), allocatable:: temp2DArray
    
    keyIndex = newLoading%keycount
    select case (newLoading%loadSpecies)
      case (LOAD_DATATYPE_PRESSURE)
        allocate(tempVals(size(newLoading%loads(1)%pressVals)))
        patm = (rho*c*c)/gamma
        call ReadBinaryReal(streamNumber, temp, status)
        newLoading%loads(keyIndex)%time = temp
        if (keyIndex > 1) then
          if (newLoading%loads(keyIndex)%time <= newLoading%loads(keyIndex-1)%time) then
            call Error("ERROR: Keys in periodic pressVals data file must be increasing.")
          end if
        end if
        call ReadBinaryReal1DArray(streamNumber, tempVals, status)
        newLoading%loads(keyIndex)%pressVals = tempVals
        if (status.ne.0) then
          call Error("ERROR: Unable to read pressVals data file.", &
                     "       Failed after reading "//trim(integertostring(keyIndex-1))//" values successfully,", &
                     "       while attempting to read the next set of pressVals data.")
        end if
        if (status.ne.0) then
          call Error("ERROR: Unable to read pressVals data file.")
        end if
        deallocate(tempVals)      
      case (LOAD_DATATYPE_LOADING)
        allocate (temp2DArray(newLoading%nbNodes,3), temp2DArrayKind4(newLoading%nbNodes,3))
        call ReadBinaryReal(streamNumber, temp, status)
        newLoading%loads(keyIndex)%time = temp
        if (keyIndex > 1) then
          if (newLoading%loads(keyIndex)%time <= newLoading%loads(keyIndex-1)%time) then
            call Error("Keys in periodic loading file must be increasing.")
          end if
        end if
        call ReadBinaryReal2DArray(streamNumber, temp2DArrayKind4, status)
        temp2DArray = temp2DArrayKind4
        if (status.ne.0) then
          if (status == 2) then
            write (*,*) "       The file ended prematurely."
          end if
          call Error("ERROR: Unable to read the loading file.", &
                     "       Failed after reading "//trim(integertostring(keyIndex-1))//" values successfully,", &
                     "       while attempting to read the next set of densVals data.")
        end if
        do i=1, newLoading%nbNodes
          if (isGood(temp2DArray(i,1)) .and. &
              isGood(temp2DArray(i,2)) .and. &
              isGood(temp2DArray(i,3))) then
            newLoading%loads(keyIndex)%loadVectors(i)= &
                     -1.0*VectorSetValue(temp2DArray(i,:))
          else
            call Error ("Loading value at i="//trim(IntegerToString(i))//" is not a number ("//&
                         trim(RealToString(temp2DArray(i,1)))//", "//&
                         trim(RealToString(temp2DArray(i,2)))//", "//&
                         trim(RealToString(temp2DArray(i,3)))//").")
          end if
        end do
        deallocate (temp2DArray, temp2DArrayKind4)
      case (LOAD_DATATYPE_FLOW)
        allocate (temp2DArray(newLoading%nbNodes,3), temp2DArrayKind4(newLoading%nbNodes,3))
        allocate (tempVals(size(newLoading%loads(1)%densVals)))
        patm = (rho*c*c)/gamma
        call ReadBinaryReal(streamNumber, temp, status)
        newLoading%loads(keyIndex)%time = temp
        if (keyIndex > 1) then
          if (newLoading%loads(keyIndex)%time <= newLoading%loads(keyIndex-1)%time) then
            write (*,*) "ERROR: Keys in periodic load data file must be increasing."
            stop
          end if
        end if
        call ReadBinaryReal1DArray(streamNumber, tempVals, status)
        newLoading%loads(keyIndex)%densVals = tempVals
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read load data file."
          write (*,*) "       Failed after reading ", keyIndex-1, " values successfully,"
          write (*,*) "       while attempting to read the next set of densVals data."
          stop
        end if
        call ReadBinaryReal2DArray(streamNumber, temp2DArrayKind4, status)
        temp2DArray = temp2DArrayKind4
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read load data file."
          write (*,*) "       Failed after reading ", keyIndex-1, " values successfully,"
          write (*,*) "       while attempting to read the next set of momVectors data."
          stop
        end if
        call ReadBinaryReal1DArray(streamNumber, tempVals, status)
        newLoading%loads(keyIndex)%pressVals = tempVals
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read load data file."
          write (*,*) "       Failed after reading ", keyIndex-1, " values successfully,"
          write (*,*) "       while attempting to read the next set of pressVals data."
          stop
        end if
        do i=1, newLoading%nbNodes
          newLoading%loads(keyIndex)%momVectors(i)=VectorSetValue(temp2DArray(i,:))
          r  => newLoading%loads(keyIndex)%densVals(i)
          if (r <= 0) then
            write (*,*) "ERROR3: densVals must be greater than zero. (",i,") -> ", r
            stop
          end if
        end do
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read load data file."
          stop
        end if
        deallocate (temp2DArray, temp2DArrayKind4, tempVals)
      end select
    if (newLoading%keycount == newLoading%nkey) then
      ! First, convert to seconds (it might be in some other units right now)
      keyMax = newLoading%loads(newLoading%nKey)%time
      keyMin = newLoading%loads(1)%time
      newLoading%keyMax = keyMax
      newLoading%keyMin = keyMin
      do keyIndex=1, newLoading%nKey
        newLoading%loads(keyIndex)%time = (newLoading%loads(keyIndex)%time - keyMin) / &
                                     (keyMax - keyMin)*newLoading%period
      end do          
      ! Use the new keyOffset
      if (keyOffset == 0) then
        newLoading%timeStepOffset = 0
      else
        timeOffset = (keyOffset - keyMin)/(keyMax - keyMin)*newLoading%period
        if( debuglevel > 6 ) then
          write (*,*) "Key offset is ",keyOffset
          write (*,*) "Key max is ",newLoading%keyMax
          write (*,*) "Period is ",newLoading%period
          write (*,*) "Time offset is ", timeOffset
        end if
        ! Next, find dt
        dt = newLoading%loads(2)%time - newLoading%loads(1)%time
        ! Finally, calculate the offsetting index (rounded to the nearest integer):
        newLoading%timeStepOffset = floor (timeOffset/dt + 0.5)
        newLoading%dt = dt
      end if 
      if( newLoading%loadSpecies==LOAD_DATATYPE_LOADING ) then
        ! Calculate the derivatives
        ! This is done with a central difference method: start with the middle
        ! f'(t)=(f(t-dt)-f(t+dt))/2dt
          
        ! NOTE: In the middle case, the variable "dt" actually stores 2dt
        ! allocate (newLoading%loads(1)%derivatives(newLoading%nKey-1))
        do keyIndex=2,newLoading%nKey-1
          dt = newLoading%loads(keyIndex+1)%time - newLoading%loads(keyIndex-1)%time
          do i=1,3
            newLoading%loads(keyIndex)%derivatives(:)%A(i) = &
              (newLoading%loads(keyIndex+1)%loadVectors(:)%A(i) - &
              newLoading%loads(keyIndex-1)%loadVectors(:)%A(i)) / dt
          end do
        end do

        ! Now we have to do the first one (we can ignore the last one since we are
        ! assuming that the first and last points are the same, and the code later
        ! ignores the last one)
        dt = newLoading%loads(2)%time - newLoading%loads(newLoading%nKey-1)%time + newLoading%period
        do i=1,3
          newLoading%loads(1)%derivatives(:)%A(i) = &
              (newLoading%loads(2)%loadVectors(:)%A(i) - &
               newLoading%loads(newLoading%nKey-1)%loadVectors(:)%A(i)) / dt
        end do
      end if
    end if
          
  end subroutine ReadPeriodicLoadingAtTime
   

  subroutine ReadAperiodicLoadingAtTime(newLoading,streamNumber)
    type(loading), intent(inout):: newLoading
    integer, intent(in):: streamNumber
    integer:: i, status, loadSpecies, keyIndex
    real(kind=4):: temp
    real :: dt, patm
    real, pointer :: r
    real(kind=4), dimension(:),   allocatable:: tempVals
    real(kind=4), dimension(:,:), allocatable:: temp2DArrayKind4
    real,         dimension(:,:), allocatable:: temp2DArray
    
    loadSpecies = GetLoadSpecies(newLoading)
    ! For aperiodic data we want this counter to stop at 4
    ! so the arrays do not overflow      
    if (newLoading%keycount < aperiodicArraySize) then
      ! ksb debug:  This is a problem when I have multi-file aperiodic files.  
      ! I want to go on just like it was one big file.  Maybe I could set a 
      ! new %nKey is the sum of the nKey for each file? (2/15/2015- Not sure if this 
      ! comment is still true since I made major changes to the code.)
      keyIndex = newLoading%keycount
    else
      keyIndex = aperiodicArraySize
    end if  
      
    call ReadBinaryReal(streamNumber, temp, status)
    if (status.ne.0) then
      call Error("Unable to read loading file.", &
       "       Failed after reading ",trim(integertostring(newLoading%keycount)), &
       " values successfully.")
    end if
    newLoading%loads(keyIndex)%time = temp
      if (keyIndex == 2) newLoading%dt = (newLoading%loads(keyIndex)%time-newLoading%loads(keyIndex-1)%time)        
      
      select case (loadSpecies)
        case (LOAD_DATATYPE_LOADING)
          allocate (temp2DArray(newLoading%nbNodes,3), temp2DArrayKind4(newLoading%nbNodes,3))
          if (debugLevel >= 14) then
            write (*,*) "  Reading loading at t=",newLoading%loads(keyIndex)%time
          end if
          if (keyIndex > 1) then
            if (newLoading%loads(keyIndex)%time <= newLoading%loads(keyIndex-1)%time) then
              call Error("Times in aperiodic loading file must be increasing.")
          end if
        end if
        call ReadBinaryReal2DArray(streamNumber,temp2DArrayKind4, status)
        temp2DArray = temp2DArrayKind4
        if (status.ne.0) then
          call Error("Unable to read loading file.", &
                     "       Failed after reading ",trim(integertostring(newLoading%keycount)), " values successfully,", &
                     "       while attempting to read the next set of loading data.")
        end if
        do i=1, newLoading%nbNodes
          if (isGood(temp2DArray(i,1)) .and. isGood(temp2DArray(i,2)) .and. &
              isGood(temp2DArray(i,3))) then
            ! The loading vector in put is force acting on the surface.  
            ! Formulation 1A needs force acting on the fluid
            newLoading%loads(keyIndex)%loadVectors(i)= -1.0*VectorSetValue(temp2DArray(i,:))
          else
            call Error ("Loading value at t="//trim(IntegerToString(2))//&
                        ", i="//trim(IntegerToString(i))//" is not a number ("//&
                        trim(RealToString(temp2DArray(i,1)))//", "//&
                        trim(RealToString(temp2DArray(i,2)))//", "//&
                        trim(RealToString(temp2DArray(i,3)))//").")
          end if
        end do          
    !ksb debug: 2/6/2015
        if (keyIndex >= 3) then
          if (keyIndex == 3) then
            dt = newLoading%loads(keyIndex)%time - newLoading%loads(keyIndex-1)%time
            ! Now do the first point using a forward difference
            do i=1,3
               newLoading%loads(1)%derivatives(:)%A(i) = &
                    (newLoading%loads(2)%loadVectors(:)%A(i) - &
                     newLoading%loads(1)%loadVectors(:)%A(i)) / dt
            end do  
          end if
          ! NOTE: In the middle case, the variable "dt" actually stores 2dt
          dt = newLoading%loads(keyIndex)%time - newLoading%loads(keyIndex-2)%time
          do i=1,3
            newLoading%loads(keyIndex-1)%derivatives(:)%A(i) = &
              (newLoading%loads(keyIndex)%loadVectors(:)%A(i) - &
               newLoading%loads(keyIndex-2)%loadVectors(:)%A(i)) / dt
          !!ksb debug: backward difference
          !  newLoading%loads(keyIndex)%derivatives(:)%A(i) = &
          !    (3*newLoading%loads(keyIndex  )%loadVectors(:)%A(i) - &
          !     4*newLoading%loads(keyIndex-1)%loadVectors(:)%A(i) + &
          !       newLoading%loads(keyIndex-2)%loadVectors(:)%A(i)) / dt
          end do
          if (newLoading%keycount .eq. newLoading%nkey) then
            ! Do the last point with a backwards difference 
            dt = newLoading%loads(keyIndex)%time - newLoading%loads(keyIndex-1)%time          
            do i=1,3              
              newLoading%loads(keyIndex)%derivatives(:)%A(i) = &
                (newLoading%loads(keyIndex)%loadVectors(:)%A(i) - &
                 newLoading%loads(keyIndex-1)%loadVectors(:)%A(i)) / dt
            end do
          end if
        end if     
        deallocate (temp2DArray, temp2DArrayKind4)

      case (LOAD_DATATYPE_FLOW)
        allocate (temp2DArray(newLoading%nbNodes,3), temp2DArrayKind4(newLoading%nbNodes,3))
        allocate(tempVals(newLoading%nbNodes))
        patm = (rho*c*c)/gamma
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read load data file."
          write (*,*) "       Failed after reading ", newLoading%keycount, " values successfully."
          stop
        end if
        if (keyIndex > 1) then
          if (newLoading%loads(keyIndex)%time <= newLoading%loads(keyIndex-1)%time) then
            write (*,*) "ERROR: Times in aperiodic load data file must be increasing."
            write (*,*) "       Timestep "//trim(IntegerToString(1))//": "//&
                        trim(RealToString(newLoading%loads(keyIndex-1)%time))//&
                        ", Timestep "//trim(IntegerToString(keyIndex))//": "//&
                        trim(RealToString(newLoading%loads(keyIndex)%time))
            stop
          end if
        end if
        call ReadBinaryReal1DArray(streamNumber, tempVals, status)
        newLoading%loads(keyIndex)%densVals = tempVals
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read load data file."
          write (*,*) "       Failed after reading ", newLoading%keycount, " values successfully,"
          write (*,*) "       while attempting to read the next set of densVals data."
          stop
        end if
        call ReadBinaryReal2DArray(streamNumber,  temp2DArrayKind4, status)
        temp2DArray = temp2DArrayKind4
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read load data file."
          write (*,*) "       Failed after reading ", newLoading%keycount, " values successfully,"
          write (*,*) "       while attempting to read the next set of momVectors data."
          stop
        end if
        call ReadBinaryReal1DArray(streamNumber, tempVals, status)
        newLoading%loads(keyIndex)%pressVals = tempVals
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read load data file."
          write (*,*) "       Failed after reading ", newLoading%keycount, " values successfully,"
          write (*,*) "       while attempting to read the next set of pressVals data."
          stop
        end if
        do i=1, newLoading%nbNodes
          newLoading%loads(keyIndex)%momVectors(i)=VectorSetValue(temp2DArray(i,:))
          r  => newLoading%loads(keyIndex)%densVals(i)
          if (r <= 0) then
            write (*,'(A)') "ERROR1: densVals must be greater than zero."
            write (*,'(A,E14.7,A,I3,A,E14.7)') "       t=",&
                    newLoading%loads(2)%time,&
                    ", location=(",i,") densVals=", r
            stop
          end if
        end do
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read load data file."
          stop
        end if
        deallocate(temp2DArray, temp2DArrayKind4, tempVals)
  
      case (LOAD_DATATYPE_PRESSURE)
        allocate (tempVals(size(newLoading%loads(1)%pressVals)))
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read pressure data file."
          write (*,*) "       Failed after reading ", newLoading%keycount, " values successfully."
          stop
        end if
        if (keyIndex > 1) then
          if (newLoading%loads(keyIndex)%time <= newLoading%loads(keyIndex-1)%time) then
            write (*,*) "ERROR: Times in aperiodic pressure data file must be increasing."
            stop
          end if
        end if
        call ReadBinaryReal1DArray(streamNumber, tempVals, status)
        newLoading%loads(keyIndex)%pressVals = tempVals
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read pressure data file."
          write (*,*) "       Failed after reading ", newLoading%keycount, " values successfully,"
          write (*,*) "       while attempting to read the next set of pressure data."
          stop
        end if
        if (status.ne.0) then
          write (*,*) "ERROR: Unable to read pressure data file."
          stop
        end if
        deallocate(tempVals)
      
      case default
        return
    end select

  end subroutine ReadAperiodicLoadingAtTime


!!This subroutine is used by all Create( )Loading subroutines to allocate data
 subroutine CreateLoadingData(newLoading, loadInfo)
    type(loading), intent(inout):: newLoading
    integer:: keyIndex
    integer,dimension(:)::loadInfo
    
    ! Allocate the memory and read the data
    if (loadInfo(LOAD_TIMETYPE) .ne. APERIODIC_LOADING) then
      allocate (newLoading%loads(newLoading%nkey))
    else
    ! The loading information is stored for 4 time segments: [L, C, R, RR]
    ! when working with aperiodic data.  We store one point back in time for 
    ! use with gradient calculations, while also storing two points
    ! forward in time for use in central difference equations and loading 
    ! pressure gradient calculations.
      allocate (newLoading%loads(4), stat=ierr)
      if (ierr >0) call error("Error in allocating single time loads arrays")
    end if
    
    newLoading%loadSpecies = loadInfo(LOAD_DATATYPE)
    allocate (newLoading%zeroLoad)
    select case (loadInfo(LOAD_DATATYPE))
      case (LOAD_DATATYPE_LOADING)
        ! Allocate the memory and read the data
        do keyIndex=1,newLoading%NKey ! ksb debug: size(newLoading%loads)
          allocate (newLoading%loads(keyIndex)%loadVectors(newLoading%nbNodes))
          allocate (newLoading%loads(keyIndex)%derivatives(newLoading%nbNodes)) 
        end do

        ! Create the zero-load vector 
        allocate (newLoading%zeroLoad%loadVectors(newLoading%nbNodes))
        newLoading%zeroLoad%loadVectors(:) = vectorSetCoordinates(0.,0.,0.)
      case (LOAD_DATATYPE_FLOW)
        ! Allocate the memory and read the data
!        allocate (temp2DArray(newLoading%nbNodes,3))
        do keyIndex=1,size(newLoading%loads)
          allocate (newLoading%loads(keyIndex)%densVals(newLoading%nbNodes))
          allocate (newLoading%loads(keyIndex)%momVectors(newLoading%nbNodes))
          allocate (newLoading%loads(keyIndex)%pressVals(newLoading%nbNodes))
        end do

        ! Create the zero-load vector
        allocate (newLoading%zeroLoad%densVals(newLoading%nbNodes))
        allocate (newLoading%zeroLoad%momVectors(newLoading%nbNodes))
        allocate (newLoading%zeroLoad%pressVals(newLoading%nbNodes))
        newLoading%zeroLoad%densVals(:) = 1
        newLoading%zeroLoad%momVectors(:) = vectorSetCoordinates(0.,0.,0.)
        newLoading%zeroLoad%pressVals(:) = 0
      case (LOAD_DATATYPE_PRESSURE)
        ! Allocate the memory and read the data
        do keyIndex=1,size(newLoading%loads)
          allocate (newLoading%loads(keyIndex)%pressVals(newLoading%nbNodes))
        end do

        ! Create the zero-load vector
        allocate (newLoading%zeroLoad%pressVals(newLoading%nbNodes))
        newLoading%zeroLoad%pressVals(:) = 0
    end select

  end subroutine CreateLoadingData


  ! A subroutine to read in the dimensions of a structured or unstructured
  ! loading file.
  subroutine ReadDimensions (newLoading, streamNumber)
    type(loading), intent(inout):: newLoading
    integer, intent(in):: streamNumber
    integer :: status
    
    ! Check to see it the file is structured or unstructured:
    if (newLoading%gridType == STRUCTURED_GRID) then
      call ReadBinaryInteger(streamNumber, newLoading%iMax, status)
      call ReadBinaryInteger(streamNumber, newLoading%jMax, status)
      if (newLoading%iMax <= 0 .or. newLoading%jMax <= 0 .or. status /= 0) then
        call Error ("Invalid structured dimensions in loading data.", &
                    "Reading "//IntegerToString(newLoading%iMax)//" by " &
                    //IntegerToString(newLoading%jMax)//&
                    " with status "//IntegerToString(status))
      end if
      newLoading%nbNodes = newLoading%iMax * newLoading%jMax 
      if (debugLevel >= 4) then
        call Message ("Structured loading dimensions: ["//trim(IntegerToString(newLoading%iMax)) &
          //", "//trim(IntegerToString(newLoading%jMax))//"].")
      end if
    else  !Unstructured Grid
      call ReadBinaryInteger(streamNumber, newLoading%nbNodes, status)
      if (newLoading%nbNodes <= 0) then
        call Error ("Invalid number of nodes: must be greater than zero.")
      else if (status /= 0) then
        call Error ("Error reading patch file.")
      end if
      if (debugLevel >= 4) then
        call Message ("Unstructured loading dimensions: [nodes] ["// &
          trim(IntegerToString(newLoading%nbNodes))//"].")
      end if
    end if
  end subroutine ReadDimensions
  
  ! This subroutine creates the loading for a cylinder object
  ! The loading is modeled by cylinder.f90
 
  ! A quick check of the load type.  If it is aperiodic
  ! then the loading has an implicit tau.
  function LoadingHasImplicitTau (ld) result (implied)
    type(loading), intent(in) :: ld
    logical :: implied

    if (ld%loadType == APERIODIC_LOADING) then
      implied = .true.
    else
      implied = .false.
    end if
    
  end function LoadingHasImplicitTau

  ! A quick subroutine to grab the timestep
  function GetLoadingDt(load) result(dt)
    type(loading), intent(in) :: load
    real::dt

    dt = load%dt
    
  end function GetLoadingDt
  
 ! A type to store all of the loading at all times
    ! A copy constructor (only valuable for periodic loading when key
  ! is not equal to zero).
  subroutine CopyLoading (newLoading, oldLoading, keyOffset)
    type(loading), intent(inout):: newLoading
    type(loading), intent(in):: oldLoading
    real, intent (in):: keyOffset
    real :: timeOffset, keyMin, keyMax
    
    newLoading%isCopy = .true.
    ! Copy the data that needs to be copied
    newLoading%loadType = oldLoading%loadType
    newLoading%iMax = oldLoading%iMax
    newLoading%jMax = oldLoading%jMax
    newLoading%nkey = oldLoading%nkey
    newLoading%gridType = oldLoading%gridType
    newLoading%nbNodes = oldLoading%nbNodes
    newLoading%period = oldLoading%period
    newLoading%keyMax = oldLoading%keyMax
    newLoading%keyMin = oldLoading%keyMin
    newLoading%dt = oldLoading%dt
    keyMax = newLoading%keyMax
    keyMin = newLoading%keyMin   
    newLoading%loadSpecies = oldLoading%loadSpecies    
    ! Use the same loading data (just store pointers to it)
    newLoading%loads => oldLoading%loads
    newLoading%zeroLoad => oldLoading%zeroLoad

    ! Use the new keyOffset
    if (keyOffset == 0.0 .or. newLoading%loadType /= PERIODIC_LOADING) then
      newLoading%timeStepOffset = 0
    else
      ! First, convert to seconds (it might be in some other units right now)
      !timeOffset = (keyOffset / newLoading%keyMax) * newLoading%period

      timeOffset = (keyOffset - keyMin)/(keyMax - keyMin)*newLoading%period
      !print*,"newLoading%loadType=",newLoading%loadtype
      !write (*,*) "Key offset is ",keyOffset
      !write (*,*) "Key max is ",newLoading%keyMax
      !write (*,*) "Period is ",newLoading%period
      !write (*,*) "Time offset is ", timeOffset
      ! Finally, calculate the offsetting index (rounded to the nearest integer):
      newLoading%timeStepOffset = floor (timeOffset/newLoading%dt + 0.5)
    end if
       
  end subroutine CopyLoading


  ! A destructor
  subroutine DestroyLoading (oldLoading)
    type(loading), intent(inout):: oldLoading
    integer:: nKey
!ksb debug:  This does not seem to have as much "zeroing" as DestroySurface.  I might need
! to check that sometime to ensure this does not result in a memory leak.
    if(.not. oldLoading%isCopy) then
      nKey = size(oldLoading%loads)
!      if (oldLoading%loadType == LOAD_TIMETYPE_APERIODIC) nKey = nKey - 1

      if (associated (oldLoading%loads)) then
        !    Deallocate the inner memory:
        call DestroySingleLoad(oldLoading, nKey)
          
        ! Deallocate the outer memory:
        deallocate (oldLoading%loads)
        nullify(oldLoading%loads)
      end if    
      if (associated (oldLoading%zeroLoad)) then
        deallocate (oldLoading%zeroLoad)
      end if
    end if
  end subroutine DestroyLoading
  
  !This subroutine deals only with the SingleTimeLoad deallocation.
  subroutine DestroySingleLoad(oldLoading, nKey)
    type(loading), intent(inout):: oldLoading
    integer:: k, nKey !, loadSpecies

    do k=1,nKey
      if (associated(oldLoading%loads(k)%loadVectors)) then
        deallocate (oldLoading%loads(k)%loadVectors)
      end if
      if (associated(oldLoading%loads(k)%densVals)) then
        deallocate (oldLoading%loads(k)%densVals)
      end if
      if (associated(oldLoading%loads(k)%momVectors)) then
        deallocate (oldLoading%loads(k)%momVectors)
      end if        
      if (associated(oldLoading%loads(k)%pressVals)) then    
        deallocate (oldLoading%loads(k)%pressVals)
      end if        
      if (associated(oldLoading%loads(k)%derivatives)) then    
        deallocate (oldLoading%loads(k)%derivatives)
      end if        
    end do

    if (associated (oldLoading%zeroLoad%loadVectors)) then
      deallocate (oldLoading%zeroLoad%loadVectors)
    end if
  end subroutine DestroySingleLoad

  
  ! An access function for the load data.
  function GetLoading(loadObject, time) result (returnedLoad)
    integer, intent(in):: time
    type(loading), intent(in)::loadObject
    type(SingleTimeLoad), pointer:: returnedLoad
    integer:: realTime
  
    select case (loadObject%loadType)
      case (CONSTANT_LOADING)
        ! For constant loading, nothing needs to be calculated - just return the
        ! loads
        returnedLoad => loadObject%loads(1)
        return  
      case (PERIODIC_LOADING)
        ! We are accessing the data here, so we use the offset.
        realTime = GetLoadingOffset(loadObject, time)
        returnedLoad => loadObject%loads(realTime)
        return      
      case (APERIODIC_LOADING)
        if (time <= size(loadObject%loads) .and. &
            time >= 1) then
          returnedLoad => loadObject%loads(time)
        else
          call Error ("Timestep is outside range of data in GetLoading.", &
                      "Aperiodic loading exists from 1 to "&
                      //trim(IntegerToString(size(loadObject%loads))),&
                      "Request was for t="//trim(IntegerToString(time)))           
        end if      
        return
      case default
        write (*,*) "ERROR: An unknown internal error has occurred."
        return
    end select
    
  end function GetLoading

  ! An access function for the loading derivative data.
  function GetLoadingDerivative (loadObject, inode, itime) result (returnDeriv)
    integer, intent(in):: itime, inode
    type(loading), intent(in)::loadObject
    type(vector), pointer, dimension(:):: returnedDerivative
    type(vector), allocatable, dimension(:), target:: derivative
    type(vector):: returnDeriv
    integer:: ip1, im1, i
    integer:: realTime

    ! Determine which type of loading we have:
    select case (loadObject%loadType) 
      case (CONSTANT_LOADING)
        ! For constant loading, nothing needs to be calculated 
        returnDeriv = vectorSetCoordinates(0.,0.,0.)
        return
        
      case (PERIODIC_LOADING)     
        ! We are accessing the data here, so we use the offset.
        realTime = GetLoadingOffset(loadObject, itime)
        returnDeriv = loadObject%loads(realTime)%derivatives(inode)
        return
        
      case (APERIODIC_LOADING)
        if( itime >0 .and. itime <= aperiodicArraySize ) then
          returnDeriv%A = loadObject%loads(itime)%derivatives(inode)%A
        else
          call Error ("Timestep is outside range of data in GetLoadingDerivative.", &
                      "Aperiodic loading exists from 1 to "&
                      //trim(IntegerToString(size(loadObject%loads))),&
                      "Request was for itime="//trim(IntegerToString(itime)))
        end if
        return
        
      case default
        call Error("ERROR: An unknown internal error has occurred in GetLoadingDerivatve.")
        return
    end select
    return
    
  end function GetLoadingDerivative
  
  ! Function for determining the loading source time.  Needed if
  ! the source time is not explicitly given and necessary to 
  ! determine the observer time of a given signal.
  function GetLoadingSourceTime (loadObject, iTau) result (sourceTime)
    type(loading), intent(in) :: loadObject
    integer, intent(in) :: iTau
    real :: sourceTime, addTime, p
    integer :: realTime
    
    
    select case (loadObject%loadType)
      case (CONSTANT_LOADING)
        ! For constant loading, nothing needs to be calculated - just return the
        ! loads
        sourceTime = 0
        return
      case (PERIODIC_LOADING)
        ! Note that the timestep offset IS NOT used in this function. It is only
        ! applied when actually accessing the load data. See note in
        ! function GetIndexAtTime ().

        ! If iTau is less than 1, we shift the whole range over until it is
        ! positive.
        if (iTau <= 0) then
          ! Calculate the number of periods to shift
          p = (iTau)/(loadObject%nKey-1)
          realTime = iTau + (abs(p)+1)*(loadObject%nKey-1)
          addTime = (p-1) * loadObject%period
        else
          realTime = modulo (iTau, loadObject%nKey-1)
          addTime = int((iTau-1) / (loadObject%nKey-1)) * loadObject%period
        end if 
        if (realTime == 0) realTime = loadObject%nKey-1
        if (realTime == loadObject%nKey) realTime = 1
        sourceTime = addTime + &
                     loadObject%loads(realTime)%time
        return
 
      case (APERIODIC_LOADING)

        ! Make sure the time is in the range:
        if (iTau > loadObject%nKey) then
          call Error ("The index iTau is greater than the size of loadObject%load")
        else
          sourceTime = loadObject%loads(iTau)%time
        end if
        return

      case default
        write (*,*) "ERROR: An unknown internal error has occurred."
        write (*,*) "       Noise will be incorrect."
        
        ! Return a zero time
        sourceTime = 0
        return
    end select
    
  end function GetLoadingSourceTime

  function GetLoadingOffset(loadObject, time) result(realTime)
    implicit none
    type(loading), intent(in):: loadObject
    integer, intent(in):: time
    
    real:: p
    integer:: nbLoadSteps, offsetTime, realTime
    
    nbLoadSteps = size(loadObject%loads)
    offsetTime = time + loadObject%timeStepOffset
    
    if (offsetTime <= 0) then
      p = (offsetTime-1)/(nbLoadSteps-1)
      realTime = offsetTime + (abs(p)+1)*(nbLoadSteps-1)
    else
      realTime = modulo (offsetTime, nbLoadSteps-1)
    end if
    if (realTime == 0) realTime = nbLoadSteps-1
    if (realTime == nbLoadSteps) realTime = 1
  
  end function GetLoadingOffset
  
  ! Access function for the number of source timesteps.
  !  Really only useful/necessary for aperiodic cases.
  function GetLoadingNTau (loadObject) result(nTau)
    type(loading), intent(in)::loadObject
    integer :: nTau
    
    ! Determine which type of loading we have:
    select case (loadObject%loadType)
      case (CONSTANT_LOADING)
        ! This one is sort of meaningless - a more appropriate answer might be
        ! infinity. I don't forsee this case being useful.
        nTau = 0
        return
        
      case (PERIODIC_LOADING)
        ! This one is sort of meaningless - a more appropriate answer might be
        ! infinity. I don't forsee this case being useful.
        nTau = loadObject%nKey-1
        return
        
      case (APERIODIC_LOADING)
        nTau = loadObject%nKey
        return 
        
      case default
        nTau = 0
        return
    end select
    
  end function GetLoadingNTau

  ! Gets the source time range, the bounds of the time range
  ! needed to calculate the signal for a given observer time range.
  subroutine GetLoadingTauRange (loadObject, tMin, tMax, tMinIndex, tMaxIndex)
    type(loading), intent(in)::loadObject
    real, intent(inout) :: tMin, tMax
    integer, intent(out) :: tMinIndex, tMaxIndex

    integer :: t

    if (tMin >= tMax) then
      call Error("Invalid time range in getTauRange.")
    end if

    ! If this object has constant loading, this function is basically
    ! meaningless, and should probably not have been called.
    if (loadObject%loadType == CONSTANT_LOADING) then
      if (debugLevel >= 1) then
        write (*,*) "WARNING: GetTauRange called for a constant loading object."
      end if
      tMinIndex = 1
      tMaxIndex = 1
      return
    end if

    ! If the loading is periodic we need to figure out what direction the
    ! searches need to go. Then we perform a linear search for the index of the
    ! timestep just below the minimum time. This is the slowest possible
    ! algorithm, but it is very easy to understand. If time here becomes an
    ! issue it can easily be replaced.
    t = 1
    if (loadObject%loadType == PERIODIC_LOADING .and. &
        GetLoadingSourceTime (loadObject, t) > tMin) then
      do  
        t = t - 1
        if (GetLoadingSourceTime (loadObject, t) < tMin) then
          exit
        end if
      end do
    else
      do
        t = t + 1
        if (GetLoadingSourceTime (loadObject, t) > tMin) then
          exit
        end if
      end do
      t = t - 1
    end if
    tMinIndex = t
    tMin = GetLoadingSourceTime (loadObject, tMinIndex)

    ! Do basically the same thing for the maximum, but start at tMinIndex, so we
    ! always want to increment t.
    t = tMinIndex + 1
    do
      if (GetLoadingSourceTime(loadObject, t) > tMax) then
        exit
      end if
      t = t + 1
    end do
    if (t > size(loadObject%loads) .and. &
        loadObject%loadType == APERIODIC_LOADING) then
      t = size(loadObject%loads)
    end if
    tMaxIndex = t
    tMax = GetLoadingSourceTime (loadObject, tMaxIndex)

  end subroutine GetLoadingTauRange


  ! This is an internal utility function, not called by external code.
  function GetIndexAtTime (S, time) result (i)
    type(loading), intent(in)::S
    real, intent(inout) :: time
    integer :: i, loadSpecies

    real :: realtime
    if (S%loadType == CONSTANT_LOADING) then
      i = 1
    else if (S%loadType == PERIODIC_LOADING) then
      ! Don't use the time offset here: that is ONLY used when accessing the
      ! loading itself. It really offsets the DATA, not the tau array, and we
      ! are looking for tau here. Even when the data is offset, the source time
      ! returned at iTau = 1 should be the first time in the array. This ensures
      ! that all objects return the same source time at a given integer time
      ! index (even though they will return different loading data at that time
      ! index).
      if (time >= 0) then
        realtime = time - (floor(time / S%period)*S%period)
      else
        realtime = time + real(1+abs(ceiling(time / S%period)))*S%period
      end if
      i = floor((realtime-S%loads(1)%time) / S%dt)+1
      i = i + ((floor(time / S%period)*(size(S%loads)-1)))
    else !if (S%loadType == APERIODIC_LOADING)
      if( debugLevel.gt.13 ) then
        call Warning('Requesting an index into the source time array ', &
        'at a time BEFORE the first source time in the aperiodic surface', &
        'Time requested: '//RealToString(time), &
        'Time set to surface start time: '//RealToString(S%loads(1)%time))
      end if
      if( debugLevel.gt.13 ) then
        call Warning('Requesting an index into the source time array ', &
        'at a time AFTER the first source time in the aperiodic surface', &
        'Time requested: '//RealToString(time), &
        'Time set to surface last time: '//RealToString(S%loads(S%nkey)%time))
      end if    
      loadSpecies = S%loadSpecies
      select case (loadSpecies)
        case (LOAD_DATATYPE_LOADING)
          i = floor((time-S%loads(1)%time) / S%dt)+1
          if (i>S%nkey) then
            i = S%nkey
          end if
          return
        case (LOAD_DATATYPE_FLOW)
          if( time < S%loads(1)%time ) then
            time = S%loads(1)%time
            i = 1
          else if( time > S%loads(S%nkey)%time ) then
            time = S%loads(S%nkey)%time
            i = S%nkey 
          else
            i = floor((time-S%loads(1)%time) / S%dt)+1
          end if
          return
        case (LOAD_DATATYPE_PRESSURE)
          if( time < S%loads(1)%time ) then
            time = S%loads(1)%time
            i = 1
          else if( time > S%loads(S%nkey)%time ) then
            time = s%loads(S%nkey)%time
            i = S%nkey 
          else
            i = floor((time-S%loads(1)%time) / S%dt)+1
          end if
          return
        case default
          write (*,*) "ERROR: An unknown internal error has occurred."
          return
      end select
    end if

  end function GetIndexAtTime

  ! Generates the source time array given nTau, tauMin and tauMax.
  subroutine CreateLoadingTauArray (loads, tauMin, tauMax, tau, iTauMin)
    type(loading), intent(in)::loads
    real, intent(inout):: tauMin, tauMax
    real, dimension (:), pointer :: tau
    integer, intent(out) :: iTauMin
    
    integer :: i, iTauMax
    if (loads%loadType == CONSTANT_LOADING) then
      call Error ("Cannot create tau array from constant loading.")
    else   
      if(debugLevel.gt.13) then
        write(*,*) 'Creating surface tau array.' !ksb - I think this should say loading tau array
        write(*,*) 'Time range exists from ', tauMin, ' to ', tauMax
      end if
      itauMin = GetIndexAtTime (loads, tauMin) - 1
      itauMax = GetIndexAtTime (loads, tauMax) + 1
      !Check if we are at the very end of the data and correct itauMax if 
      !necessary
      if (itauMax == GetLoadingNTau(loads)+1) then
        itauMax = GetLoadingNTau(loads)
      end if  
      if (loads%loadType == APERIODIC_LOADING) then
        if (itauMin < 1) then
          itauMin = 1
        else if (iTauMax > GetLoadingNTau(loads)) then
          itauMax = GetLoadingNTau(loads)
        end if
        if (iTauMin > GetLoadingNTau(loads)) then
          call Error("The minimum source time to calculate acoustic pressVals is larger than", &
                     "the maximum aperiodic loading data.", "Reset your observer time.", &
                     "Minumum source time = "//RealToString(tauMin), &
                     "Maximum aperiodic loading time = "// &
                     RealToString(GetLoadingSourceTime(loads, GetLoadingNTau(loads))))
        end if
        if (iTauMax < 1) then
          call Error("The maximum source time to calculate acoustic pressVals is smaller than", &
                     "the minimum aperiodic loading data.", "Reset your observer time.", &
                     "Maximum source time = "//RealToString(tauMax), &
                     "Minimum aperiodic loading time = "// &
                     RealToString(GetLoadingSourceTime(loads, 1)))
        end if
      end if
      if(debugLevel.gt.13) then
        write(*,*) 'Indices of surface array are ', itauMin, ' and ', itauMax
      end if
      allocate (tau(itauMax-iTaumin+1))
      do i=1,size(tau)
        tau(i) = GetLoadingSourceTime (loads, iTauMin+i-1)
      end do
    end if

  end subroutine CreateLoadingTauArray
  
  subroutine RotateLoadingArray(loadObject)
    type(loading):: loadObject
    type(singleTimeLoad):: tempLoad
    integer:: i
    ! Because PSU-WOPWOP employs pointers for its arrays we can't simply rotate 
    ! the 'loads' arrays. We must reassign the arrays within loads at the lowest
    ! level of pointers.
    tempLoad              = loadObject%loads(1)
    do i=1,aperiodicArraySize-1 
      loadObject%loads(i) = loadObject%loads(i+1)
    end do
    loadObject%loads(aperiodicArraySize) = tempLoad

    return
  end subroutine RotateLoadingArray
  
  subroutine ResetLoadingKeyCount(loadObject, val)
    type(loading):: loadObject
    integer, intent(in), optional:: val
    if( present(val) ) then
      loadObject%keycount = val
    else
      loadObject%keycount = 0
    end if
  end subroutine ResetLoadingKeyCount

  subroutine SetLoadingNKey(loadObject, val)
    type(loading):: loadObject
    integer, intent(in), optional:: val
    if( present(val) ) then
      loadObject%nKey = val
    end if
  end subroutine SetLoadingNKey
   
  ! function to determine load type (const./periodic/aperiodic)
  function GetLoadType (loads) result (loadType)
    type(loading), intent(in)::loads
    integer :: loadType

    loadType = loads%loadType
    
    return
    
  end function GetLoadType
  
  ! function to determine load type (const./periodic/aperiodic)
  function GetLoadSpecies (loadObject) result (loadSpecies)
    type(loading), intent(in)::loadObject
    integer :: loadSpecies

    loadSpecies = loadObject%loadSpecies
    return
    
  end function GetLoadSpecies
  
  function GetLoadnKey (loadObject) result(nKey)
    type(loading), intent(in):: loadObject
    integer:: nKey
    
    nKey = loadObject%nKey
  end function GetLoadnkey
  
  function GetLoadPeriod(loadObject) result(period)
    type(loading), intent(in):: loadObject
    real:: period
    
    period = loadObject%period
  end function GetLoadPeriod
  
end module loadingObject
