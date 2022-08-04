! PSU-WOPWOP
! $Id: surface.f90 3366 2016-10-31 17:16:11Z brentner $
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
! 3/1/2003 - 2/28/2006 NASA Langlreadconnectey Grant NAG-03025 
! (NASA Langley Research Center, D. P. Lockard, Technical Officer)
!
!
! Written by Guillaume Bres, Guillaume Perez, Leonard Lopes, Hsuan-Nien Chen, 
! Christopher Hennes, Rui Cheng and Benjamin Goldman. 
! Faculty advisor Dr. Kenneth S. Brentner.


module surfaceObject
!****************************************************************************
!
! MODULE 
!   surfaceObject
! 
! DESCRIPTION
!   This module defines a class called "surface" with its constructors,
!   destructors and access functions. This class reads in the surface from a
!   file. It functions as a black box to the patch object, hiding the type of 
!   surface from it (i.e. constant, periodic and aperiodic). The box isn't 
!   perfect - sometimes patch needs to know anyway, but for the most part it 
!   is opaque.
!
! TYPES
!   surface
!   SingleSurface
!
! PUBLIC SUBROUTINES
!  CreateSurface (newSurface, surfaceFileUnit, dataFormat, timeOffset)
!  CopySurface (newSurface, oldSurface, timeOffset)
!  DestroySurface (oldSurface)
!  ReadSurfaceDimensions(newSurface, surfaceFileUnit)
!
! PRIVATE SUBROUTINES
!  CreateSurfaceData (newSurface)
!  ReadSurfaceDimensions(newSurface, streamNumber)
!  ReadConstantSurfaceAtTime(newSurface,streamNumber,endFlag)
!  ReadPeriodicSurfaceAtTime(newSurface,streamNumber,keyOffset,endFlag)
!  ReadAperiodicSurfaceAtTime(newSurface,streamNumber,endFlag)
! 
! PUBLIC FUNCTIONS
!  HasImplicitTau (ld) result (implied)
!  GetSurface (surfaceObject, time) result (returnedSurface)
!  GetSurfaceDerivative (surfaceObject, time) result (returnedDerivative)
!  GetSurfaceSourceTime (surfaceObject, iTau) result (sourceTime)
!  GetNTau (surfaceObject) result(nTau)
!  GetLoadType (grid) result (surfaceType)
!  GetSurfaceSizeInMemory(grid) result(res)
!  GetSurfacenKey(surfaceObject)
!
!  Integer parameters (constants) describing the type of surface:
!    CONSTANT_SURFACE
!    PERIODIC_SURFACE
!    APERIODIC_SURFACE
! 
! HISTORY
!   Created 12/17/04 by Chris Hennes
!
!****************************************************************************
  use constantsModule, only: GEO_GRIDTYPE, GEO_GRIDTYPE_UNSTRUCTURED, GEO_TIMETYPE, &
    GEO_TIMETYPE_APERIODIC, GEO_TIMETYPE_PERIODIC, GEO_TIMETYPE_CONSTANT, GEO_IBLANK, &
    GEO_NORMTYPE, APERIODICARRAYSIZE, GEO_NORMTYPE_NODES, GEO_NORMTYPE_CELLS, &
    GEO_TIMETYPE_QPERIODIC
  use mathModule, only: vector, VectorSetValue, operator(*), operator(-), operator(+), &
    vectorSetCoordinates, vectorSubtract, vectorAbsolute, operator(.cross.), operator (/)
  use IOModule, only: isGood, ReadBinaryReal, ReadBinaryInteger, ReadBinaryReal1DArray, &
    ReadBinaryInteger2DArray, ReadBinaryReal2DArray
  use debugModule
  use MPIModule

  implicit none

  private
  public:: surface, SingleSurface, ReadSurfaceAtTime
  public:: CreateSurface, CopySurface, GetSurface, DestroySurface
  public:: GetSurfaceDimensions,ReadSurfaceDimensions
  public:: readDimensions, readConnectivity
  public:: GetSurfaceVelAccel, GetSurfaceType, GetSurfaceGridType
  public:: GetSurfaceTauRange, GetSurfaceUnstructuredFaces
  public:: GetBoundingBoxAtTime, WriteSurfaceDebugInfo
  public:: GetSurfaceDt, SurfaceHasImplicitTau, GetSurfaceNTau, GetSurfaceSourceTime
  public:: CreateSurfaceTauArray, SetSurfaceNKey
  public:: IsSurfaceCompact, GetSurfaceNbNodes, GetPositionAtTime
  public:: CONSTANT_SURFACE, PERIODIC_SURFACE, APERIODIC_SURFACE
  public:: STRUCTURED_GRID, UNSTRUCTURED_GRID
  public:: IsNodeBased, GetSurfaceOffset
  public:: GetSurfaceMaxTime, GetSurfaceMinTime, GetSurfacenKey, GetSurfaceGridSize, GetSurfaceIblanking
  public:: RotateSurfaceArray, ResetSurfaceKeyCount, GetSurfaceTimeStepOffset
  public:: GetSurfaceFirstNode, GetSurfaceLastNode, SetSurfaceNodeRange

  ! A type to store the complete surface vector at a single time
  type SingleSurface
    ! The time the surface is for
    real::time
    ! The surface data
    type(vector) :: minCorner, maxCorner ! the bounding box
    real, pointer, dimension(:)::areas => NULL()
    type(vector), pointer, dimension(:)::vectors => NULL()
    type(vector), pointer, dimension(:)::normals => NULL()
    integer, pointer, dimension(:)::iblanks => NULL()
    ! If we are recording sigma data and are using cell-centered data, we need
    ! to store off the node positions seperately, since WOPWOP doesn't use them
    ! for anything.
    type(vector), pointer, dimension(:)::nodeCoord => NULL()
  end type SingleSurface

  
  ! A type to store all of the surface at all times
  type Surface
    private
    character (len=32) :: title
    ! The surface type is one of CONSTANT_SURFACE, PERIODIC_SURFACE or 
    ! APERIODIC_SURFACE, as defined by the constants below.
    integer::surfaceType
    ! gridType is one of STRUCTURED_GRID, UNSTRUCTURED_GRID. If it is
    ! structured, iMax and jMax must be specified. If it is unstructured,
    ! nbNodes and nbNodesPerFace. nbNodes must always be set by the read-in
    ! routines.
    integer::gridType
    integer::iMax,jMax,nbNodes,nbFaces,nbData
    integer::normType
    integer, dimension(:,:),pointer:: connectivity =>NULL()
    !
    ! logical to keep track if this surface uses iblanking or not
    logical:: iblanking
    integer:: FirstNode, LastNode
    ! The period of the data, in seconds.
    real::period
    ! The starting and ending key values (for periodic data)
    real::keyMin, keyMax
    ! The differential time segment for the surface
    real::dt
    ! The timestep offset
    integer::timeStepOffset
    ! This is what is actually read in from the file:
    type(SingleSurface), dimension(:), pointer::grid =>NULL()
    ! These are calculated from the data when it is read
    type(SingleSurface), dimension(:), pointer::velocity =>NULL()
    type(SingleSurface), dimension(:), pointer::acceleration =>NULL()
    ! This is a storage space for a zero vector to return when an error
    ! occurrs or the data ends.
    type(SingleSurface),pointer :: zero
    ! The variables used for MPI communication
    integer::nkey
    !Array for storing coord data
    integer::keycount
    !real,dimension(:,:),pointer::tempPositionArray
    !real,dimension(:,:),pointer::tempDataArray
    !Variable to determine whether to delete any data:
    logical::isCopy
    type(vector),dimension(:),pointer::nodeCoord
  end type Surface
  
  
  ! Some constants:
  integer, parameter:: CONSTANT_SURFACE  = 001, &
                       PERIODIC_SURFACE  = 002, &
                       APERIODIC_SURFACE = 003, &
                       CYLINDER_MODEL    = 004

  integer, parameter:: STRUCTURED_GRID   = 101, &
                       UNSTRUCTURED_GRID = 102
  
  integer, parameter:: IBLANKED_GRID     = 1
  integer, parameter:: MAX_NODES_PER_FACE= 10

contains

  subroutine WriteSurfaceDebugInfo(surf,unitnum)
    type(surface), intent(in)::surf
    integer:: unitnum
    write(unitnum,*) '*** SurfaceDebugInfo ***'
    write(unitnum,*) 'title= ', trim(surf%title)
    write(unitnum,*) 'surfaceType= ', trim(integertostring(surf%surfaceType))
    write(unitnum,*) 'gridType= ', trim(integertostring(surf%gridType))
    write(unitnum,*) 'iMax= ', trim(integertostring(surf%iMax))
    write(unitnum,*) 'jMax= ', trim(integertostring(surf%jMax))
    write(unitnum,*) 'nbNodes= ', trim(integertostring(surf%nbNodes))
    write(unitnum,*) 'nbFaces= ', trim(integertostring(surf%nbFaces))
    write(unitnum,*) 'nbData= ', trim(integertostring(surf%nbData))
    write(unitnum,*) 'normType= ',trim(integertostring( surf%normType))
    write(unitnum,*) 'period= ', trim(realtostring(surf%period))
    write(unitnum,*) 'keyMax= ', trim(realtostring(surf%keyMax))
    write(unitnum,*) 'dt= ', trim(realtostring(surf%dt))
    write(unitnum,*) 'nKey= ', trim(integertostring(surf%nKey))
    write(unitnum,*) 'keyCount= ', trim(integertostring(surf%keyCount))
    write(unitnum,*) 'isCopy= ', surf%isCopy
    write(unitnum,*) 'timeStepOffset= ', surf%timeStepOffset
    write(unitnum,*) 'connectivity associated? ', associated(surf%connectivity)
    if(associated(surf%connectivity)) then
      write(unitnum,*) 'Size of sonnectivity ', trim(integertostring(size(surf%connectivity)))
    end if
    write(unitnum,*) 'zero associated? ', associated(surf%zero)
    if(associated(surf%zero)) then
      call WriteSingleSurfaceDebugInfo(surf%zero,unitnum)
    end if
    write(unitnum,*)'grid associated? ', associated(surf%grid)
    if(associated(surf%grid)) then
      write(unitnum,*) 'Size of grid', trim(integertostring(size(surf%grid)))
      write(unitnum,*) 'First grid follows'
      call WriteSingleSurfaceDebugInfo(surf%grid(1),unitnum)
      write(unitnum,*) 'Last grid follows'
      call WriteSingleSurfaceDebugInfo(surf%grid(size(surf%grid)),unitnum)
    end if
    write(unitnum,*)'Velocity associated? ', associated(surf%velocity)
    if(associated(surf%velocity)) then
      write(unitnum,*) 'Size of velocity', trim(integertostring(size(surf%velocity)))
      write(unitnum,*) 'First velocity follows'
      call WriteSingleSurfaceDebugInfo(surf%velocity(1),unitnum)
      write(unitnum,*) 'Last velocity follows'
      call WriteSingleSurfaceDebugInfo(surf%velocity(size(surf%velocity)),unitnum)
    end if
    write(unitnum,*) 'Acceleration associated? ', associated(surf%acceleration)
    if(associated(surf%acceleration)) then
      write(unitnum,*) 'Size of acceleration ', trim(integertostring(size(surf%acceleration)))
      write(unitnum,*) 'First acceleration follows'
      call WriteSingleSurfaceDebugInfo(surf%acceleration(1),unitnum)
      write(unitnum,*) 'Last acceleration follows'
      call WriteSingleSurfaceDebugInfo(surf%acceleration(size(surf%acceleration)),unitnum)
    end if
    write(unitnum,*) '--- End SurfaceDebugInfo ---'

  end subroutine WriteSurfaceDebugInfo
  
  subroutine WriteSingleSurfaceDebugInfo(surf,unitnum)
    type(SingleSurface), intent(in)::surf
    integer::unitnum
    write(unitnum,*) '*** SingleSurfaceDebugInfo ***'
    write(unitnum,*) 'time=', trim(realtostring(surf%time))
    write(unitnum,*) 'minCorner=', surf%minCorner
    write(unitnum,*) 'maxCorner=', surf%maxCorner
    write(unitnum,*) 'associated vectors? ', associated(surf%vectors)
    if(associated(surf%vectors)) then
      write(unitnum,*) 'Size of vectors ', trim(integertostring(size(surf%vectors)))
    end if
    write(unitnum,*) 'associated normals? ', associated(surf%normals)
    if(associated(surf%normals)) then
      write(unitnum,*) 'Size of normals ', trim(integertostring(size(surf%normals)))
    end if
    write(unitnum,*) 'associated areas? ', associated(surf%areas)
    if(associated(surf%areas)) then
      write(unitnum,*) 'Size of areas ', trim(integertostring(size(surf%areas)))
    end if
    write(unitnum,*)'--- End SingleSurfaceDebugInfo ---'

  end subroutine WriteSingleSurfaceDebugInfo

  ! A constructor. Creates the surface by reading it in from the file.
  ! Parameters:
  !   newSurface - The new surface object to create. Note that the memory must
  !                already exist for this object.
  !   surfaceFileUnit - The unit number of the file which contains the load
  !                     information. The file must already be open and the header
  !                     read in prior to reaching this function.
  !   keyOffset - The amount to offset the key value. For example, if the data
  !               is periodic and this blade starts at an azimuth of 90 degrees
  !               relative to the first data point in the surface file then the
  !               key offset is 90 degrees.
  !   geoInfo - information from the patch file fixed header
  !   title   - surface title
  subroutine CreateSurface (newSurface, surfaceFileUnit,&
                            keyOffset,geoInfo,title)
    type(surface), intent(inout):: newSurface
    integer, intent(in) :: surfaceFileUnit
    real, intent(in) :: keyOffset
    character(len=32),optional :: title
    integer,dimension(:)::geoInfo
    character(len=30):: temp
    !nullify all pointers before we begin
    nullify(newSurface%connectivity,newSurface%grid,newSurface%velocity,&
            newSurface%acceleration,newSurface%zero,newSurface%nodeCoord)
            
    if (present(title)) then
      newSurface%title = title
    end if
    newSurface%isCopy = .false.

    if (geoInfo(GEO_GRIDTYPE) == GEO_GRIDTYPE_UNSTRUCTURED) then
      newSurface%gridType = UNSTRUCTURED_GRID
    else
      ! Defaults to structured if not specified
      newSurface%gridType = STRUCTURED_GRID
    end if

    select case(geoInfo(GEO_TIMETYPE))
      case (GEO_TIMETYPE_CONSTANT)
        newSurface%surfaceType = CONSTANT_SURFACE
        newSurface%dt = 0.0
      case (GEO_TIMETYPE_PERIODIC, GEO_TIMETYPE_QPERIODIC)
        newSurface%surfaceType = PERIODIC_SURFACE        
      case default !(GEO_TIMETYPE_APERIODIC or GEO_TIMETYPE_MTFAPERIODIC)
        newSurface%surfaceType = APERIODIC_SURFACE
    end select
    if( geoInfo(GEO_IBLANK) == IBLANKED_GRID ) then
      newSurface%iblanking = .true.
    else
      newSurface%iblanking = .false.
    end if
    newSurface%keycount = 0
    newSurface%normType = geoInfo(GEO_NORMTYPE)
    newSurface%nKey     = 1
    ! read the dimensions from the dynamic header (just for this patch)
    call ReadSurfaceDimensions(newSurface, surfaceFileUnit)
    !ksb debug: set defaults for first and last node in the patch (for patch subsetting - for David Lockard; 3/20/2015)
    newSurface%FirstNode = 1
    newSurface%LastNode = newSurface%nbNodes
    call CreateSurfaceData (newSurface)
    
    if (debuglevel > 14) then
      ! write out surface information for general code debugging.
      write(*,*) '#### Surface Information for, ',trim(newSurface%title),' ###'
      write(*,*) 'Surfacetype is ',newSurface%surfaceType,' (1= constant, 2=periodic, 3=aperiodic)'
      write(*,*) 'Gridypte is ',newSurface%gridtype,' (101 = structured, 102= unstructured)'
      write(*,*) 'IsCopy: ',newSurface%IsCopy
      
      if (newSurface%gridtype == STRUCTURED_GRID) then
        write(*,*) 'Imax: ', newSurface%imax
        write(*,*) 'Jmax: ', newSurface%jmax
        write(*,*) 'Associated surf%nodeCoord: ',associated(newSurface%nodeCoord)
        write(*,*) 'Associated grid(1)%nodeCoord: ',associated(newSurface%grid(1)%nodeCoord)
        write(*,*) 'Associated grid(1)%vectors: ',associated(newSurface%grid(1)%vectors)
      else
        write(*,*) 'nbNodes: ',newSurface%nbnodes
        write(*,*) 'nbFaces: ',newSurface%nbFaces
      end if
      
      write(*,*) 'nbData: ', newSurface%nbData
      write(*,*) 'Normal type: ', newSurface%normtype,' (1= node centered, 2 = face centered)'
      if (newSurface%surfacetype == 2 ) then
        write(*,*) 'Periodic Surface'
        write(*,*) 'Period: ',newSurface%period
        write(*,*) 'nkey: ', newSurface%nkey
        write(*,*) 'timeStepOffset: ',newSurface%timeStepOffset
      elseif (newSurface%surfacetype ==3) then
        write(*,*) 'Aperiodic Surface'
        write(*,*) 'nkey:',newSurface%nkey
      else
        write(*,*) 'Constant Surface'
      end if
      
      write(*,*) 'dTau: ',newSurface%dt

      write(*,*) '### End Surface Information ###'
    end if
 end subroutine CreateSurface


 subroutine ReadSurfaceDimensions(newSurface, streamNumber)
   type(surface), intent(inout):: newSurface
   integer, intent(in):: streamNumber
   
   integer:: status
   real(kind=4):: temp
   
    select case(newSurface%surfaceType)
      case (PERIODIC_SURFACE)
        call ReadBinaryReal(streamNumber, newSurface%period, status)
        if (debugLevel >= 4) then
          call Message ("Period: "//trim(RealToString(newSurface%period))//" s.")
        end if
        if (newSurface%period < 0 .or. status /= 0) then
          call Error("ERROR: Invalid period in periodic surface file.")
        end if

        if (newSurface%period > 1.e15) then
          call Error ("The period of the surface file is too large.", &
                      "Period = "//trim(RealToString(newSurface%period)))
        end if

        ! Read in the number of Key locations
        call ReadBinaryInteger(streamNumber, newSurface%nKey, status)
        if (debugLevel >= 4) then
          call Message ("Number of keys: "//trim(IntegerToString(newSurface%nKey)))
        end if


        if (newSurface%nKey < 0 .or. status /= 0) then
          call Error("ERROR: Invalid number of key locations in periodic surface file.")
        end if
        newSurface%dt = newSurface%Period/(newSurface%nKey-1)
      
      case (APERIODIC_SURFACE)
      ! Read in the number of times
        call ReadBinaryInteger(streamNumber, newSurface%nKey, status)
        if (newSurface%nKey < 0 .or. status /= 0) then
          call Error ("Invalid number of key locations in aperiodic flow datafile.", &
                      "Surface '"//trim(newSurface%title)//"'")
        end if
        if (debugLevel >= 4) then
          call Message ("Number of keys: "//trim(IntegerToString(newSurface%nKey)))
        end if
        if (newSurface%nKey==1 )then
          !
          ! Even though this surface was specified as Aperiodic, it onlhy has one time step,
          ! hence it is really constant.
          !
          newSurface%surfaceType = CONSTANT_SURFACE
        end if
    end select
          
    ! Read in the dimensions (they get stored in the object:
    call ReadDimensions (newSurface, streamNumber)
    if (debugLevel >= 4) then
      if (newSurface%gridType == STRUCTURED_GRID) then
        write (*,*) "Structured grid dimensions: [",&
                        trim(IntegerToString(newSurface%iMax)),", ", &
                        trim(IntegerToString(newSurface%jMax)),"]"
      else
        write (*,*) "Unstructured grid dimensions: [nodes,faces] [", &
                        trim(IntegerToString(newSurface%nbNodes)), ", ", &
                        trim(IntegerToString(newSurface%nbFaces)), "]"
      end if
    end if
    ! Read in the connectivity (also stored in the object):
    call ReadConnectivity (newSurface, streamNumber)    
  end subroutine ReadSurfaceDimensions


  !  Allocates all the necessary arrays for each type of loading.
  subroutine CreateSurfaceData(newSurface)
    type(surface), intent(inout):: newSurface
    integer:: nKey, keyIndex
    
    ! Allocate the memory and read the data
    if (newSurface%surfaceType .ne. APERIODIC_SURFACE) then
      ! for constant or periodic surfaces, we allocate the entire time history (if 
      ! there is one).  Setting nKey here specifies how long it is.
      nKey = newSurface%nKey
    else 
    ! The surface information is stored for aperiodicArraySize (4) time segments: [L, C, R, RR]
    ! when working with aperiodic data.  We store one point back in time for 
    ! use with gradient calculations, while also storing two points
    ! forward in time for use in central difference equations and loading 
    ! pressure gradient calculations.  (Nominally aperiodicArraySize is a constant = 4)
      nKey = aperiodicArraySize
    end if 
    
    allocate (newSurface%grid(nKey))
    
    do keyIndex=1,nKey
      allocate (newSurface%grid(keyIndex)%vectors(newSurface%nbData))
      allocate (newSurface%grid(keyIndex)%normals(newSurface%nbData))
      allocate (newSurface%grid(keyIndex)%iblanks(newSurface%nbData))
      allocate (newSurface%grid(keyIndex)%areas(newSurface%nbData))
      allocate (newSurface%grid(keyIndex)%nodeCoord(newSurface%nbNodes))
    end do
 
    ! Create the zero-load vector
    allocate (newSurface%zero)
    allocate (newSurface%zero%vectors(newSurface%nbData))
    allocate (newSurface%zero%normals(newSurface%nbData))
    allocate (newSurface%zero%areas(newSurface%nbData)) 
    
    if (newSurface%normType==2) then  !ksb debug - this section needs to be fixed for face centered normals
      allocate (newSurface%nodeCoord(newSurface%nbNodes))
!      allocate (newSurface%grid(1)%nodecoord(newSurface%nbNodes))
    end if
    newSurface%zero%vectors(:) = vectorSetCoordinates(0.,0.,0.)
    newSurface%zero%normals(:) = vectorSetCoordinates(0.,0.,0.)
    newSurface%zero%areas(:) = 0
    
  end subroutine CreateSurfaceData

  subroutine ReadSurfaceAtTime(newSurface, surfaceFileUnit, keyOffset, timeType) 
    type(surface), intent(inout):: newSurface
    integer, intent(in):: surfaceFileUnit,timeType
    real::keyOffset
    
    newSurface%keyCount = newSurface%keyCount + 1
    if( newSurface%keyCount <= newSurface%nKey ) then
      select case (timeType)
        case (GEO_TIMETYPE_CONSTANT)
          call ReadConstantSurfaceAtTime(newSurface,SurfaceFileUnit)
        case (GEO_TIMETYPE_PERIODIC, GEO_TIMETYPE_QPERIODIC)
          call ReadPeriodicSurfaceAtTime(newSurface,SurfaceFileUnit,keyOffset)
        case (GEO_TIMETYPE_APERIODIC)
          call ReadAperiodicSurfaceAtTime(newSurface,SurfaceFileUnit)
        case default
          call Error ("Unknown time type in geometry ("//&
                      trim(IntegerToString(timeType))//")")
        end select
    end if
  end subroutine ReadSurfaceAtTime
  
  subroutine ReadConstantSurfaceAtTime(newSurface,streamNumber)
    type(Surface),intent(inout)::newSurface
    integer::streamNumber
    ! for constant surface this is only done once
    call ReadSurfaceData (newsurFace,streamNumber,1)
    call CalculateSurfaceAreas (newSurface, 1)
    call createVelAccelArrays(newSurface) 
    call CalculateVelocity (newSurface)
    call CalculateAcceleration (newSurface)
    !In the case of constant surface the dt is unknown so it is set to zero
    newSurface%dt = 0.0
  end subroutine ReadConstantSurfaceAtTime

  subroutine ReadPeriodicSurfaceAtTime(newSurface,streamNumber,keyOffset)
    type(Surface),intent(inout)::newSurface
    integer::streamNumber
    real,intent(in)::keyOffset
    integer::stat, keyIndex, TimeSeg, index
    real(kind=4)::temp
    real::keyMin,keyMax,dt,timeOffset    
    
    keyIndex = newSurface%keycount
    ! read in the keyvalue (time step)
    call ReadBinaryReal(streamNumber, temp, stat)
    newSurface%grid(keyIndex)%time = temp
      
    if (debugLevel >= 14) then
      write (*,*) "Reading key = ",newSurface%grid(keyIndex)%time
    end if
   
    call ReadSurfaceData (newSurface,streamNumber,keyIndex)
    
    if (keyIndex == newSurface%nKey) then 
      ! We just read the last time step in the periodic data.
      ! Convert all key values to be time-based:
      keyMax = newSurface%grid(newsurface%nKey)%time
      keyMin = newSurface%grid(1)%time
      newSurface%keyMax = keyMax
      newSurface%keyMin = keyMin
        
      do keyIndex=1,newsurface%nKey
        newSurface%grid(keyIndex)%time = (newSurface%grid(keyIndex)%time-keyMin)/ &
                                         (keyMax - keyMin)*newSurface%period
      end do

      ! Use the new keyOffset
      if (keyOffset == 0) then
        newSurface%timeStepOffset = 0
      else
        ! First, convert to seconds (it might be in some other units right now)
        timeOffset = (keyOffset - keyMin)/(keyMax - keyMin)*newSurface%period
        ! Next, find dt
        dt = newSurface%grid(2)%time - newSurface%grid(1)%time
        ! Finally, calculate the offsetting index (rounded to the nearest integer):
        newSurface%timeStepOffset = floor (timeOffset/dt + 0.5)
        newSurface%dt = dt
      end if
    
      ! For periodic surface the differential time step is the difference between
      ! two time steps
 !ksb debug:  this is where dt gets changed. 5/14/2016 ksb     
      newSurface%dt = newSurface%grid(2)%time - newSurface%grid(1)%time
      dt=newSurface%grid(3)%time - newSurface%grid(2)%time
      call CalculateSurfaceAreas (newSurface, 1)
      !Create velocity and acceleration vectors
      call createVelAccelArrays(newSurface)
      call CalculateVelocity (newSurface)
      call CalculateAcceleration (newSurface)                         
    end if

  end subroutine ReadPeriodicSurfaceAtTime


  subroutine ReadAperiodicSurfaceAtTime(newSurface,streamNumber)
    type(surface), intent(inout):: newSurface
    integer, intent(in):: streamNumber
    integer:: t,stat, keyIndex
    real(kind=4):: temp
    real :: dt

    ! For aperiodic data we want this counter to stop at aperiodicArraySize (nominally 4)
    ! so the arrays do not overflow      
    if (newSurface%keycount < aperiodicArraySize) then
      keyIndex = newSurface%keycount
    else
      keyIndex = aperiodicArraySize
    end if  
            
    call ReadBinaryReal(streamNumber, temp, stat)
    newSurface%grid(keyIndex)%time = temp
    if (debugLevel >= 14) then
      write (*,*) "  Reading surface data at t=",newSurface%grid(keyIndex)%time
    end if
    if (keyIndex == 2) then
      newSurface%dt = newSurface%grid(keyIndex)%time-newSurface%grid(keyIndex-1)%time
    end if
      
    call ReadSurfaceData (newSurface,streamNumber,keyIndex)
    call CalculateSurfaceAreas (newSurface, keyIndex)
    if (keyIndex == 1) then
      !Create velocity and acceleration vectors
      call createVelAccelArrays(newSurface)
    else !(keyIndex > 1)
      if (newSurface%grid(keyIndex)%time <= newSurface%grid(keyIndex-1)%time) then
        call Error ("Keys in aperiodic surface file must be increasing.",&
                    "Surface '"//trim(newSurface%title)//"'",&
                    "Key #"//trim(IntegerToString(keyIndex-1))//" = "&
                           //trim(RealToString(newSurface%grid(keyIndex-1)%time))&
                    //", Key #"//trim(IntegerToString(keyIndex))//" = "&
                               //trim(RealToString(newSurface%grid(keyIndex)%time)))
      end if      
!      if (dt > 1.001 * newSurface%dt .or. dt < .999 * newSurface%dt) then
!        call Error ("dt changes between surface timesteps: this is not allowed",&
!                    "First dt: "//trim(RealToString(newSurface%dt)),&
!                    "Found dt="//trim(RealToString(dt))//&
!                    " at step "//trim(IntegerToString(t)) )
!        end if
    end if
    if (keyIndex >= 3) then
      call CalculateVelocity (newSurface, keyIndex)
      call CalculateAcceleration (newSurface, keyIndex)
    end if
      
  end subroutine ReadAperiodicSurfaceAtTime

  subroutine CreateVelAccelArrays(newSurface)
    type(surface), intent(inout):: newSurface
    integer:: t, i
    
    allocate (newSurface%velocity(size(newSurface%grid)))        
    allocate (newSurface%acceleration(size(newSurface%grid)))                
    do t = 1,size(newSurface%grid)
      allocate(newSurface%velocity(t)%vectors(newSurface%nbData))
      allocate(newSurface%velocity(t)%normals(newSurface%nbData))
      allocate(newSurface%velocity(t)%areas(newSurface%nbData))   
      
      allocate(newSurface%acceleration(t)%vectors(newSurface%nbData))
      allocate(newSurface%acceleration(t)%normals(newSurface%nbData))
      allocate(newSurface%acceleration(t)%areas(newSurface%nbData))    
      do i=1,newSurface%nbData
        newSurface%velocity(t)%vectors(i)%A = 0.0
        newSurface%velocity(t)%normals(i)%A = 0.0
        newSurface%velocity(t)%areas(i) = 0.0
        
        newSurface%acceleration(t)%vectors(i)%A = 0.0
        newSurface%acceleration(t)%normals(i)%A = 0.0
        newSurface%acceleration(t)%areas(i) = 0.0            
      end do
    end do
  end subroutine CreateVelAccelArrays

  subroutine ReadSurfaceData (newSurface, streamNumber, pos)
    implicit none
    type(Surface),intent(inout)::newSurface
    integer,intent(in)::streamNumber, pos
    integer::stat, i, n
    integer, dimension(:,:), allocatable:: tempIntArray
    real(kind=4), dimension(:,:), allocatable:: tempPositionArrayKind4, tempDataArrayKind4
    real,         dimension(:,:), allocatable:: tempPositionArray, tempDataArray
    
    allocate (tempPositionArray(newSurface%nbNodes,3), tempPositionArrayKind4(newSurface%nbNodes,3))
    allocate (tempDataArray(newSurface%nbData,3), tempDataArrayKind4(newSurface%nbData,3))
     
    call ReadBinaryReal2DArray(streamNumber,tempPositionArrayKind4, stat)
    tempPositionArray = tempPositionArrayKind4
    if (stat .ne. 0) then
      call Error('Unable to read surface position information from data file at step '//IntegerToString(pos),&
                 'Error in surface specification', &
                 'Stat came back as '//IntegerToString(stat))
    end if
      
    call ReadBinaryReal2DArray(streamNumber,tempDataArrayKind4, stat)
    tempDataArray = tempDataArrayKind4
  
    if (stat .ne. 0) then
      call Error('Unable to read surface normal vector information from data file at step '//IntegerToString(pos),&
                 'Error in surface specification', &
                 'Stat came back as '//IntegerToString(stat), &
                 'Expected array size was '//IntegerToString(newSurface%nbData))
    end if
    
    if(newSurface%iblanking) then
      allocate (tempIntArray(newSurface%nbData,1)) 
      call ReadBinaryInteger2DArray(streamNumber,tempIntArray,stat)
      if (stat /= 0) then
        call Error('Unable to read iblank information ',&
                   'Error in surface specification', &
                   'Stat came back as '//IntegerToString(stat))
      end if
    end if

    ! Grab the bounding box:
    newSurface%grid(pos)%minCorner%A(1) = &
      minval(tempPositionArray(:,1))
    newSurface%grid(pos)%minCorner%A(2) = &
      minval(tempPositionArray(:,2))
    newSurface%grid(pos)%minCorner%A(3) = &
      minval(tempPositionArray(:,3))
    newSurface%grid(pos)%maxCorner%A(1) = &
      maxval(tempPositionArray(:,1))
    newSurface%grid(pos)%maxCorner%A(2) = &
      maxval(tempPositionArray(:,2))
    newSurface%grid(pos)%maxCorner%A(3) = &
      maxval(tempPositionArray(:,3))

    !if (debugLevel >= 14) then
    !  write (*,*) newSurface%grid(pos)%minCorner%A, " --> ", &
    !              newSurface%grid(pos)%maxCorner%A
    !end if
    select case (newSurface%normType)
      case (GEO_NORMTYPE_NODES)
        ! Node centered data, no conversion necessary
        do i=1, newSurface%nbNodes
          newSurface%grid(pos)%vectors(i) = VectorSetValue(tempPositionArray(i,1:3))
          newSurface%grid(pos)%normals(i) = VectorSetValue(tempDataArray(i,1:3))
          if (vectorAbsolute(newSurface%grid(pos)%normals(i)) > 1.0e8) then
            call Error ("Normal vectors seem excessively large - check your file format.")
          end if 
          ! iblank information
          if( newSurface%iblanking ) then
            newSurface%grid(pos)%iblanks(i) = tempIntArray(i,1)
          else
            newSurface%grid(pos)%iblanks(i) = 1
          end if
          ! Use just a subset of the nodes - for Dave Lockard (3/20/2015) - may change this
          ! when we can actually do source parallel.
          if( i < newSurface%FirstNode .or. i > newSurface%LastNode ) then
              newSurface%grid(pos)%iblanks(i) = 0
          end if
        end do
      case (GEO_NORMTYPE_CELLS)
        ! We really are only interested in the cell centers, so we need to
        ! calculate them: we do this by adding up all the nodes at each face and
        ! averaging their values.
        do i=1,newSurface%nbNodes
          newSurface%grid(pos)%nodecoord(i) = vectorSetValue(temppositionarray(i,1:3))
        end do
        do i=1,newSurface%nbFaces
          newSurface%grid(pos)%vectors(i)%A = 0.0
          ! n=1 is how many connected cells, from 2 on is the connectivity
          do n=2,(newSurface%connectivity(i,1)+1)
            newSurface%grid(pos)%vectors(i)%A = newSurface%grid(pos)%vectors(i)%A + &
               tempPositionArray(newSurface%connectivity(i,n),1:3)
          end do
          newSurface%grid(pos)%vectors(i)%A = newSurface%grid(pos)%vectors(i)%A / &
            (newSurface%connectivity(i,1))
          ! The normals are specified at the cell centers already, so no
          ! conversion is neccessary. 
          newSurface%grid(pos)%normals(i) =  VectorSetValue(tempDataArray(i,1:3))
          if (vectorAbsolute(newSurface%grid(pos)%normals(i)) > 1.0e8) then
            call Error ("Normal vectors seem excessively large - check your file format.")
          end if 
          ! iblank information
          if( newSurface%iblanking ) then
            newSurface%grid(pos)%iblanks(i) = tempIntArray(i,1)
          else
            newSurface%grid(pos)%iblanks(i) = 1
          end if 
          ! Use just a subset of the nodes - for Dave Lockard (3/20/2015) - may change this
          ! when we can actually do source parallel.
          if( i < newSurface%FirstNode .or. i > newSurface%LastNode ) then
              newSurface%grid(pos)%iblanks(i) = 0
          end if
        end do     
    end select
    
    deallocate (tempPositionArray, tempPositionArrayKind4)
    deallocate (tempDataArray, tempDataArrayKind4)

  end subroutine ReadSurfaceData


  ! A subroutine to read in the dimensions of a structured or unstructured
  ! surface file.
  subroutine ReadDimensions (newSurface, streamNumber)
    type(surface), intent(inout):: newSurface
    integer, intent(in):: streamNumber

    integer :: status
    
    ! Check to see it the file is structured or unstructured:
    if (newSurface%gridType == STRUCTURED_GRID) then
      call ReadBinaryInteger(streamNumber, newSurface%iMax, status)
      call ReadBinaryInteger(streamNumber, newSurface%jMax, status)
      
      if (newSurface%iMax <= 0 .or. newSurface%jMax <= 0 .or. status /= 0) then
        call Error ("Invalid structured dimensions in surface file.", &
                    "Reading "//IntegerToString(newSurface%iMax)//" by " &
                    //IntegerToString(newSurface%jMax)//&
                    " with status "//IntegerToString(status))
      end if
      newSurface%nbNodes = newSurface%iMax * newSurface%jMax
      newSurface%nbFaces = (newSurface%iMax-1) * (newSurface%jMax-1)
      if (newSurface%iMax > 10000 .or. newSurface%jMax > 10000) then
        call Warning ("The surface dimensions are very large. This may indicate",&
                      "a file format problem.", &
                      "imax = "//trim(IntegerToString(newSurface%iMax))//&
                      ", jmax = "//trim(IntegerToString(newSurface%jMax)))
      end if
    else  !Unstructured grid
      call ReadBinaryInteger(streamNumber, newSurface%nbNodes, status)
      call ReadBinaryInteger(streamNumber, newSurface%nbFaces, status)
      !ksb debug: (9/3/15) I don't think this check is really appropriate.  The number of faces can be 
      ! larger than the number of nodes.  Think of a cube with each side divided into 2 triangles.
      ! This will have 8 nodes and 12 faces.  This comment and this if block can be deleted in the
      ! future 
      !if(newSurface%nbNodes < newSurface%nbFaces)then
      !  call Warning ("Number of nodes is less than the number of faces.")
      !end if
      if (newSurface%nbNodes <= 0) then
        call Error ("Invalid number of nodes: must be greater than zero.")
      else if (status /= 0) then
        call Error ("Error reading patch file.")
      end if
    end if
    
    if (newSurface%normType == GEO_NORMTYPE_NODES) then
      newSurface%nbData = newSurface%nbNodes 
    else
      newSurface%nbData = newSurface%nbFaces
    end if
    
  end subroutine ReadDimensions

  
  ! A subroutine to read in the connectivity list for an unstructured surface
  ! file. When called on a structured file the function does nothing and returns
  ! cleanlu (i.e. - not an error).
  subroutine ReadConnectivity (newSurface, streamNumber)
    type(surface), intent(inout):: newSurface
    integer, intent(in):: streamNumber

    integer :: stat,nodesPerFace,node,i,j
    stat = 0
    select case (newSurface%gridType)
      case (UNSTRUCTURED_GRID)
        ! Allocate the necessary memory
        allocate (newSurface%connectivity(newSurface%nbFaces, MAX_NODES_PER_FACE+1))
        ! Read the connectivity list:
        do i=1,newSurface%nbFaces
          call ReadBinaryInteger(streamNumber,nodesPerFace,stat)
          if (debugLevel >= 5) then
            write (*,*) "Unstructured nodes per face was ",nodesPerFace
          end if
          if (nodesPerFace >= 10) then
            call Error("Maximum number of supported nodes per face ("     &
                //trim(IntegerToString(MAX_NODES_PER_FACE))//") exceeded.",&
                       "Check your file format if this is not correct.")      
          end if
          newSurface%connectivity(i,1)=nodesPerFace
          do j=1,nodesPerFace
            call ReadBinaryInteger(streamNumber,node,stat)
            newSurface%connectivity(i,j+1)=node
          end do
        end do
        if (stat/= 0) then
          call Error ("Failed reading the unstructured connectivity list.")
        end if
      case default !(STRUCTURED_GRID)
        ! Do nothing if the grid is structured: there is no connectivity list
        return        
    end select
  end subroutine ReadConnectivity


  ! A subroutine to calculate the surface area for a single surface (i.e. at a
  ! single time). It assumes that all necessary memory has been allocated and
  ! that the appropriate data has been read in.
  subroutine CalculateSurfaceAreas (newSurface, keyIndex)
    type(Surface), intent(inout) :: newSurface
    real :: aux
    integer :: i,timestep, minTimestep, maxTimestep, surf=1
    logical :: unitNormals
    integer:: keyIndex

    ! There are two possibilities here: first, if the normal vectors are of
    ! unit length then the area will be calculated using a cross-product
    ! methodology. If the normals are not of unit length, it is assumed that
    ! their length represents the area. The area is stored off and the vectors
    ! are normalized.

    ! First, find out if the vectors are of unit length: we will only check a
    ! small sample, assuming that if these are unit length, so is everything
    ! else.
    unitNormals = .true.
    do i=1,newSurface%nbData, 5 !size(newSurface%grid(keyIndex)%normals),5
      !ksb debug - this looks wrong - it only checks the first normal.  That might
      ! be OK, but then there is no need for the do loop.  Really, I just wonder if
      ! I should always calculate structured surface areas and use normals for 
      ! unstructured cases.
      aux = vectorAbsolute(newSurface%grid(keyIndex)%normals(1))
      if (aux < 1E-15) then
        ! We're going to treat this as zero length, which is a reasonable value
        ! for a degenerate geometry: don't fail on this. Note that 1E-15 is
        ! entirely heuristic.
        cycle
      else if (aux > 1.0e6) then
        call Error ("Surface normals seem to be abnormally large - check your file format")
      end if
      ! Do a fuzzy match since it will rarely be exactly 1.0000      
      if (aux > 1.001 .or. aux < .999) then
        unitNormals = .false.
        exit ! Bail out of the loop early
      end if
    end do
    if (.not. unitNormals .or. newSurface%gridType == UNSTRUCTURED_GRID ) then
      ! We just pull them from the normals:
      ! The length of unit vector represents area, so area can be determined by normal vector easily.
      if (newSurface%surfaceType==APERIODIC_SURFACE) then 
        minTimestep = keyIndex
        maxTimestep = keyIndex
      else
        minTimestep = 1
        maxTimestep = newSurface%nKey !size(newSurface%grid)
      end if
      
      do timestep = minTimestep,maxTimestep
        do i=1,size(newSurface%grid(timestep)%normals)
          aux = VectorAbsolute(newSurface%grid(timestep)%normals(i))
          if (aux > 0.0) then
            newSurface%grid(timestep)%normals(i) = &
              newSurface%grid(timestep)%normals(i) / aux
          else
            newSurface%grid(timestep)%normals(i)%A = 0
          end if
          newSurface%grid(timestep)%areas(i) = aux
        end do
      end do
    else ! if (newSurface%gridType == STRUCTURED_GRID) then
      ! We calculate the areas
      call CalculateStructuredSurfaceAreas (newSurface, keyIndex)
    !else
    !  call Error (&
    !  "This unstructured grid has unit normal vectors. For unstructured", &
    !  "grids, the normal vectors must be scaled by the area they represent.")
    end if
      
  end subroutine CalculateSurfaceAreas


  ! A subroutine to calculate, based on the positions of the vertices, the area
  ! associated with each node in a structured mesh. If the mesh is a compact
  ! mesh then this routine calculates the length instead.
  subroutine CalculateStructuredSurfaceAreas (newSurface, keyIndex)
    type(Surface), intent(inout) :: newSurface
    type(vector), dimension(:,:), allocatable :: coordinates
    real, dimension(:,:), allocatable :: areas
    integer :: i, j, iMax, jMax, timestep,l, pos
    real :: aux, aux2
    integer:: keyIndex

    iMax = newSurface%iMax
    jMax = newSurface%jMax
    if (iMax <=0 .or. jMax <=0) then
      call Error (&
      "CalculateStructuredSurfaceAreas was called for a mesh that does not", &
      "appear to be structured.")
    end if

    ! Set up the temp memory
    allocate (coordinates(iMax,jMax), areas(iMax,jMax))

    do timestep = 1, size(newSurface%grid)
      ! The existance of keyIndex in this routine means that it has been called
      ! from readAperiodicSurfaceAtTime.  Thus we need some special logic to 
      ! determine which grid we need calculate the coordinates of.
      if (.false.) then !ksb debug: 3/8/15: (keyIndex .ne. -1) then
        pos = keyIndex
      else
        pos = timestep
      end if
      ! Recalling that everything is stored in flat arrays, which we want to
      ! index using 2 coordinates (i,j), we reshape the original array so that
      ! we can access it this way. This makes the following code much easier to
      ! understand.
      l=0
      do j=1,jMax
        do i=1,iMax
          l=l+1
          coordinates(i,j)=newSurface%grid(pos)%vectors(l)
        enddo
      enddo
     
      if(iMax==1 .and. jMax==1) then
        !This is a point compact patch 
        areas(1,1)=0.1 
      else if (iMax == 1) then
        ! This is a compact patch and we don't really want area, we want length.
        ! ( increment in j direction for this compact patch )
        aux=vectorAbsolute(coordinates(1,2)-coordinates(1,1))
        areas(1,1)=0.5*aux
        do j=2,jMax-1
           aux2=vectorAbsolute(coordinates(1,j+1)-coordinates(1,j))
           areas(1,j)=0.5*(aux+aux2)
           aux=aux2
        end do
        areas(1,jMax)=0.5*aux

      else if (jMax == 1) then
        ! This is a compact patch and we don't really want area, we want length.
        ! ( increment in i direction for this compact patch )
        aux=vectorAbsolute(coordinates(2,1)-coordinates(1,1))
        areas(1,1)=0.5*aux
        do i=2,iMax-1
           aux2=vectorAbsolute(coordinates(i+1,1)-coordinates(i,1))
           areas(i,1)=0.5*(aux+aux2)
           aux=aux2
        end do
        areas(iMax,1)=0.5*aux
      else
        ! This is a real surface patch: calculate the cell areas:
        do i=2,iMax-1
           areas(i,1)=0.25*(CellArea(coordinates,i,1)+&
                            CellArea(coordinates,i-1,1))
           do j=2,jMax-1
              areas(i,j)=0.25*(CellArea(coordinates,i,j)+&
                               CellArea(coordinates,i-1,j)+&
                               CellArea(coordinates,i,j-1)+&
                               CellArea(coordinates,i-1,j-1))                
           end do
           areas(i,jMax)=0.25*(CellArea(coordinates,i,jMax-1)+&
                               CellArea(coordinates,i-1,jMax-1))
        end do
        do j=2,jMax-1
           areas(1,j)=0.25*(CellArea(coordinates,1,j)+&
                            CellArea(coordinates,1,j-1))
           areas(iMax,j)=0.25*(CellArea(coordinates,iMax-1,j)+&
                               CellArea(coordinates,iMax-1,j-1))
        end do
        areas(1,1)=0.25*CellArea(coordinates,1,1)
        areas(1,jMax)=0.25*CellArea(coordinates,1,jMax-1)
        areas(iMax,1)=0.25*CellArea(coordinates,iMax-1,1)
        areas(iMax,jMax)=0.25*CellArea(coordinates,iMax-1,jMax-1)
      end if
      ! Pack the 2D area array back onto the 1D storage array:
      newSurface%grid(pos)%areas = PACK (areas, .true.)
    end do ! End of loop over timesteps

    ! Deallocate the temp memory
    deallocate (coordinates, areas)
    
  end subroutine CalculateStructuredSurfaceAreas

  
  ! This subroutine calculates the area of a single cell.
  ! First, the area of the 4 triangles defined by the 4 points are
  ! computed, using the cross product of 2 of the 3 vectors of a triangle 
  ! (Recall: 0.5*abs(OM cross product ON) is the area of the triangle OMN). Then
  ! the average of the 4 areas is computed to get the surface of the all cell. 
  function CellArea(coordinates,i,j) result(dScell)
    type (vector), dimension(:,:), intent(in) :: coordinates
    integer, intent(in)::i,j
    real::dScell

    dScell=.25*&
      (vectorAbsolute((coordinates(i,  j+1) - coordinates(i,  j))   .cross. &
                      (coordinates(i+1,j)   - coordinates(i,  j)))   +      &
       vectorAbsolute((coordinates(i+1,j+1) - coordinates(i,  j+1)) .cross. &
                      (coordinates(i,  j)   - coordinates(i,  j+1))) +      &
       vectorAbsolute((coordinates(i,  j+1) - coordinates(i+1,j+1)) .cross. &
                      (coordinates(i+1,j)   - coordinates(i+1,j+1))) +      &
       vectorAbsolute((coordinates(i+1,j+1) - coordinates(i+1,j))   .cross. &
                      (coordinates(i,  j)   - coordinates(i+1,j))))
    
  end function CellArea


  function SurfaceHasImplicitTau (ld) result (implied)
    type(surface), intent(in) :: ld
    logical :: implied

    if (ld%surfaceType == APERIODIC_SURFACE) then
      implied = .true.
    else
      implied = .false.
    end if
    
  end function SurfaceHasImplicitTau

  ! A copy constructor (only valuable for periodic surface when key
  ! is not equal to zero).
  subroutine CopySurface (newSurface, oldSurface, keyOffset)
    type(surface), intent(inout):: newSurface
    type(surface), intent(in):: oldSurface
    real, intent(in):: keyOffset

    real :: timeOffset, keyMin, keyMax
    newSurface%isCopy = .true.
    
    ! Copy the data that needs to be copied
    newSurface%surfaceType = oldSurface%surfaceType
    newSurface%gridType = oldSurface%gridType
    if(newSurface%gridType == STRUCTURED_GRID) then
      newSurface%iMax = oldSurface%iMax
      newSurface%jMax = oldSurface%jMax
      nullify(newsurface%connectivity)
    else
      newSurface%nbNodes = oldSurface%nbNodes
      newSurface%nbFaces = oldSurface%nbFaces
      newSurface%connectivity => oldSurface%connectivity
    end if
    newSurface%nbData = oldSurface%iMax*oldSurface%jMax
    newSurface%period = oldSurface%period
    newSurface%keyMax = oldSurface%keyMax
    newSurface%keyMin = oldSurface%keyMin
    newSurface%nkey = oldSurface%nkey
    ! Use the same surface data (just store pointers to it)
    newSurface%grid => oldSurface%grid
    newSurface%velocity => oldSurface%velocity
    newSurface%acceleration => oldSurface%acceleration
    newSurface%zero => oldSurface%zero
    ! Use the new keyOffset
    if (keyOffset == 0 .or. newSurface%surfaceType /= PERIODIC_SURFACE) then
      newSurface%timeStepOffset = 0
    else
      ! First, convert to seconds (it might be in some other units right now)
      !timeOffset = (keyOffset / newSurface%keyMax) * newSurface%period
      keyMax = newSurface%keyMax
      keyMin = newSurface%keyMin 
      timeOffset = (keyOffset - keyMin)/(keyMax - keyMin)*newSurface%period
      ! Finally, calculate the offsetting index (rounded to the nearest integer):
      newSurface%timeStepOffset = floor (timeOffset/newSurface%dt + 0.5)
    end if
       
  end subroutine CopySurface
  
  ! A destructor for Surface type
  subroutine DestroySurface (oldSurface)
    type(surface), intent(inout):: oldSurface
    integer:: k

    if (.not. oldSurface%isCopy) then
      oldSurface%title =''
      oldSurface%surfaceType = 0
      oldSurface%gridType = 0
      oldSurface%iMax=0; oldSurface%jMax=0; oldSurface%nbNodes=0; oldSurface%nbFaces=0
      oldSurface%nbData = 0
      oldSurface%normType = 0
      oldSurface%period = 0.0
      oldSurface%keyMin=0.0; oldSurface%keyMax=0.0
      oldSurface%dt = 0.0
      oldSurface%timeStepOffset = 0
      oldSurface%nkey = 0
      oldSurface%keycount = 0
      
      if (associated(oldSurface%connectivity)) then
        deallocate(oldSurface%connectivity)
      end if
      if (associated(oldSurface%grid)) then
        do k=lbound(oldSurface%grid,1),ubound(oldSurface%grid,1)
          call DestroySingleSurface(oldSurface%grid(k))
        end do
        deallocate(oldSurface%grid)
      end if
      if (associated(oldSurface%velocity)) then
        do k=lbound(oldSurface%velocity,1),ubound(oldSurface%velocity,1)
          call DestroySingleSurface(oldSurface%velocity(k))
        end do
        deallocate(oldSurface%velocity)
      end if
      if (associated(oldSurface%acceleration)) then
        do k=lbound(oldSurface%acceleration,1),ubound(oldSurface%acceleration,1)
          call DestroySingleSurface(oldSurface%acceleration(k))
        end do
        deallocate(oldSurface%acceleration)
      end if
      if (associated(oldSurface%zero)) then
        call DestroySingleSurface(oldSurface%zero)
        deallocate(oldSurface%zero)
      end if
      if (associated(oldSurface%nodeCoord)) then
        deallocate(oldSurface%nodeCoord)
      end if
    end if
  end subroutine DestroySurface

  ! A destructor for Single Surface type
  subroutine DestroySingleSurface(oldSingleSurface)
    type(SingleSurface), intent(inout):: oldSingleSurface
    
    oldSingleSurface%time=0.0
    if (associated(oldSingleSurface%areas)) then   
      deallocate(oldSingleSurface%areas)
    end if
    if (associated(oldSingleSurface%vectors)) then   
      deallocate(oldSingleSurface%vectors)
    end if
    if (associated(oldSingleSurface%normals)) then      
      deallocate(oldSingleSurface%normals)
    end if
    if (associated(oldSingleSurface%iblanks)) then                   
      deallocate(oldSingleSurface%iblanks)
    end if
    if (associated(oldSingleSurface%nodeCoord)) then                   
      deallocate(oldSingleSurface%nodeCoord)
    end if 
  end subroutine DestroySingleSurface
  
  ! An access function
  function GetSurface (surfaceObject, time) result (returnedSurface)
    integer, intent(in):: time
    type(surface), intent(in)::surfaceObject
    type(SingleSurface), pointer:: returnedSurface
    integer:: realTime, offsetTime
    real :: p

    ! Determine which type of surface we have:
    select case (surfaceObject%surfaceType)
    
      case (CONSTANT_SURFACE)
        ! For constant surface, nothing needs to be calculated - just return the
        ! grid
        returnedSurface => surfaceObject%grid(lbound(surfaceObject%grid,1))
        return
 
      case (PERIODIC_SURFACE)
        realTime = GetSurfaceOffset(surfaceObject, time)
        returnedSurface => surfaceObject%grid(realTime)
        return
 
      case (APERIODIC_SURFACE)
        if( time >= 1 .and. time <= size(surfaceObject%grid) ) then
          returnedSurface => surfaceObject%grid(time)
        else
          !returnedSurface => surfaceObject%zero
          call Error ("Timestep is outside range of data in GetSurface.", &
                      "Aperiodic surface exists from 1 to "&
                      //trim(IntegerToString(size(surfaceObject%grid))),&
                      "Request was for t="//trim(IntegerToString(time)))           
        end if
        
        return
 
      case default
        ! There is no graceful way of handling this, it's an internal error.
        ! Die.
        call Error ("In internal error occurred in GetSurface (surface.f90)")
        return

    end select

    return
    
  end function GetSurface


  ! An access function
  subroutine GetSurfaceVelAccel (surfaceObject, time, vel, accel, keyIndexSurf)
    integer, intent(in):: time
    type(surface), intent(in)::surfaceObject
    type(SingleSurface), pointer:: vel, accel
    integer:: realTime, offsetTime, keyIndexSurf
    real :: p   
    
    ! Do a little error-checking
    if (time < 1 .and. surfaceObject%surfaceType == APERIODIC_SURFACE) then
      call Error("ERROR: Timestep number must be greater than zero.")
    end if
    ! Determine which type of surface we have:
    select case (surfaceObject%surfaceType)
    
      case (CONSTANT_SURFACE)
        ! For constant surface, nothing needs to be calculated - just return
        ! zero
        vel => surfaceObject%zero
        accel => surfaceObject%zero
        return        
      case (PERIODIC_SURFACE)
        realTime = GetSurfaceOffset(surfaceObject, time)
        vel => surfaceObject%velocity(realTime)
        accel => surfaceObject%acceleration(realTime)
        return
        
      case (APERIODIC_SURFACE)
        if (time <= size(surfaceObject%grid) .and. &
            time >= 1) then
          vel   => surfaceObject%velocity(keyIndexSurf)
          accel => surfaceObject%acceleration(keyIndexSurf)
        else
          vel   => surfaceObject%zero
          accel => surfaceObject%zero
        end if
        return

      case (CYLINDER_MODEL)
        vel => surfaceObject%zero
        accel => surfaceObject%zero
        return
        
      case default
        write (*,*) "ERROR: An unknown internal error has occurred."
        stop

    end select

    return
    
  end subroutine GetSurfaceVelAccel
  
  function GetSurfaceOffset(surfaceObject, time) result(realTime)
    implicit none
    type(surface), intent(in):: surfaceObject
    integer, intent(in):: time
    
    real:: p
    integer:: nbSurfSteps, offsetTime, realTime
    
    nbSurfSteps = size(surfaceObject%grid)
    offsetTime = time + surfaceObject%timeStepOffset
    if (offsetTime <= 0) then
      p = (offsetTime-1)/(nbSurfSteps-1)
      realTime = offsetTime + (abs(p)+1)*(nbSurfSteps-1)
    else
      realTime = modulo (offsetTime, nbSurfSteps-1)
    end if
    if (realTime == 0) realTime = nbSurfSteps-1
    if (realTime == nbSurfSteps) realTime = 1  
  
  end function GetSurfaceOffset


  function GetIndexAtTime (S, time) result (i)
    type(surface), intent(in)::S
    real, intent(inout) :: time
    integer :: i, maxKeyIndex

    real :: realtime

    select case (S%surfaceType)
      case (CONSTANT_SURFACE)
        i = 1
      case (PERIODIC_SURFACE)
        ! Normalize the incoming time so it falls withing a period:
        if (time > 0) then
          realtime = time - (floor(time / S%period)*S%period)
        else
          realtime = time + real(1+abs(ceiling(time / S%period)))*S%period
        end if
        i = floor((realtime-S%grid(1)%time) / S%dt)+1
        i = i + ((floor(time / S%period)*(size(S%grid)-1)))
      case (APERIODIC_SURFACE)
        !
        ! if time is outside of range, set it to the nearest endpoint of time range
        !
        maxKeyIndex = size(S%grid)-3
        if( time < S%grid(1)%time ) then
          if( debugLevel.gt.13 ) then
            call Warning('Requesting an index into the source time array ', &
              'at a time BEFORE the first source time in the aperiodic surface', &
              'Time requested: '//RealToString(time), &
              'Time set to surface start time: '//RealToString(S%grid(1)%time))
          end if
          time = s%grid(1)%time
          i = 1
        else if( time > S%grid(maxKeyIndex)%time ) then
          if( debugLevel.gt.13 ) then
            call Warning('Requesting an index into the source time array ', &
              'at a time AFTER the first source time in the aperiodic surface', &
              'Time requested: '//RealToString(time), &
              'Time set to surface last time: '//RealToString(S%grid(maxKeyIndex)%time))
          end if
          time = s%grid(maxKeyIndex)%time
          i = maxKeyIndex 
        else
          i = floor((time-S%grid(1)%time) / S%dt)+1
        end if
    end select

  end function GetIndexAtTime


  function GetPositionAtTime (S, time, i) result (pos)
    type(surface), intent(in)::S
    real, intent(inout) :: time
    integer, intent(in) :: i
    type(vector) :: pos
    integer :: timestep, index, maxKeyIndex
    real :: interp
    type(SingleSurface), pointer:: surfAtT, surfAtTplus1

    select case (S%surfaceType)
      case (CONSTANT_SURFACE)
        ! The constant surface case is easy -- try it first'
        pos = S%grid(lbound(S%grid,1))%vectors(i)
      case (PERIODIC_SURFACE)
        ! If we don't have a constant surface then we will need to do an
        ! interpolation.
        timestep = floor(time/S%dt)
        interp = (time - (timestep * S%dt)) / S%dt
        surfAtT => GetSurface (S, timestep)
        surfAtTplus1 => GetSurface (S, timestep+1)
        pos = (1.0-interp) * surfAtT%vectors(i) + &
                   interp  * surfAtTplus1%vectors(i)
      case (APERIODIC_SURFACE)
        if(debugLevel.gt.13) then
          write(*,*) 'Getting position on surface at given time ', time
          write(*,*) 'Index of location on surface is ', i, &
                     ' of ', GetSurfaceNbNodes(S)
        end if 
        maxKeyIndex = size(S%grid)-3
        if(time.lt.S%grid(1)%time) then
          if(debugLevel.gt.13) then
            write(*,*) 'Adjusting source time to first time of surface'
          end if
          time = S%grid(1)%time
          surfAtT => GetSurface(S, 1)
          pos = surfAtT%vectors(i)
        else if(time.gt.S%grid(maxKeyIndex)%time) then
          if(debugLevel.gt.13) then
            write(*,*) 'Adjusting source time to last time of surface'
          end if
          time = S%grid(maxKeyIndex)%time
          surfAtT => GetSurface(S, maxKeyIndex)
          pos = surfAtT%vectors(i)
        else
          if(debugLevel.gt.13) then
            write(*,*) 'Source time is within the surface time range.'
            write(*,*) 'Time wanted is ', time, ', surface delta T is ', S%dt
          end if
          index = GetIndexAtTime(S,time) 
          interp = (time - S%grid(index)%time)/S%dt
          if(debugLevel.gt.13) then
            write(*,*) 'Determined that index is ', index
          end if
          if( index < maxKeyIndex .and. index > 0 ) then
            surfAtT => GetSurface (S, index)
            surfAtTplus1 => GetSurface (S, index+1)
            pos = (1.0-interp) * surfAtT%vectors(i) + &
                   interp  * surfAtTplus1%vectors(i)
          else if( index == maxKeyIndex ) then
            surfAtT => GetSurface (S, index)
            pos = surfAtT%vectors(i) 
          else      
            print*,'Not sure how this can happen, index = ', index
          end if 
        end if
        if(debugLevel.gt.13) then
          write(*,*) 'Determined position as ', pos 
        end if
    end select
  end function GetPositionAtTime
  

  subroutine GetBoundingBoxAtTime (S, time, minCorner, maxCorner)
    type(surface), intent(in)::S
    real, intent(in) :: time
    type(vector), intent(out) :: minCorner, maxCorner    
    integer :: start
    real :: interp, realtime

    select case (S%surfaceType)
      case (CONSTANT_SURFACE)
        minCorner = S%grid(1)%minCorner
        maxCorner = S%grid(1)%maxCorner
      case (PERIODIC_SURFACE)
        ! Normalize the incoming time so it falls withing a period:
        if (time > 0) then
          realtime = time - (floor(time / S%period)*S%period)
        else
          realtime = time + real(1+abs(ceiling(time / S%period)))*S%period
        end if
        
        ! Get the lower bound for the interpolation:
        start = floor((realtime-S%grid(1)%time) / S%dt) + 1
        ! Get the interpolation position
        interp = (realtime-S%grid(start)%time) / S%dt
  
        if (start .ge. size(S%grid)) then
          minCorner = S%grid(size(S%grid))%minCorner
          maxCorner = S%grid(size(S%grid))%maxCorner
        else
          maxCorner = (1.0-interp)*S%grid(start  )%maxCorner + &
                           interp *S%grid(start+1)%maxCorner
          maxCorner = (1.0-interp)*S%grid(start  )%maxCorner + &
                           interp *S%grid(start+1)%maxCorner
        end if
      
      case (APERIODIC_SURFACE)
    
        ! Get the lower bound for the interpolation:
        start = floor((time-S%grid(1)%time) / S%dt) + 1
        if (start < 1) then
          minCorner = S%grid(1)%minCorner
          maxCorner = S%grid(1)%maxCorner
          !call Warning ("Outside time bounds for data. Not extrapolating.")
        else if (start >= size(S%grid)) then
          minCorner = S%grid(size(S%grid))%minCorner
          maxCorner = S%grid(size(S%grid))%maxCorner
          !call Warning ("Outside time bounds for data. Not extrapolating.")
        else
          ! Get the interpolation position
          interp = (time-S%grid(start)%time) / S%dt

          maxCorner = (1.0-interp)*S%grid(start  )%maxCorner + &
                           interp *S%grid(start+1)%maxCorner
          maxCorner = (1.0-interp)*S%grid(start  )%maxCorner + &
                           interp *S%grid(start+1)%maxCorner
        end if
    end select

  end subroutine GetBoundingBoxAtTime
  

  function GetSurfaceSourceTime (surfaceObject, iTau) result (sourceTime)
    type(surface), intent(in) :: surfaceObject
    integer, intent(in) :: iTau
    real :: sourceTime, addTime, p
    integer :: realTime
    
! HENNES NOTE: We can't die here, the source time routine sometimes goes out of
! range when it's bracketing and we need to handle that situation gracefully.
! What we'd really like to do is throw an exception here that the bracketing
! routine could catch. Since obviously we can't do that, this should probably 
! be modified to return a logical that indicates success or failure, so that 
! external code can decide if it cares or not. Notice that the code below
! already handles these cases by returning the "zero" surface.
!    ! Do a little error-checking
!    if (iTau < 1 .and. surfaceObject%surfaceType == APERIODIC_SURFACE) then
!      call Error ("Timestep number must be greater than zero "//&
!                  "in GetSurfaceSourceTime.")
!    else if (iTau > size(surfaceObject%grid) .and.&
!             surfaceObject%surfaceType == APERIODIC_SURFACE) then
!      call Error ("Timestep number must be less than the array size "//&
!                  "in GetSurfaceSourceTime.")
!    end if

    ! Determine which type of surface we have:
    select case (surfaceObject%surfaceType)

      case (CONSTANT_SURFACE)
        ! For constant surface, nothing needs to be calculated - just return the
        ! grid
        sourceTime = 0
        return
 
      case (PERIODIC_SURFACE) 
        ! Note that the timestep offset IS NOT used in this function. It is only
        ! applied when actually accessing the load data. This makes it appear
        ! from the outside like the surface is different for different copied
        ! patches, but that the time range is the same (as it should be).
        
        ! If iTau is less than 1, we shift the whole range over until it is
        ! positive.
        if (iTau <= 0) then
          ! Calculate the number of periods to shift
          p = (iTau)/(size(surfaceObject%grid)-1)
          realTime = iTau + (abs(p)+1)*(size(surfaceObject%grid)-1)
          addTime = (p-1) * surfaceObject%period
        else
          realTime = modulo (iTau, size(surfaceObject%grid)-1)
          addTime = ((iTau-1) / (size(surfaceObject%grid)-1)) * surfaceObject%period
        end if 
        if (realTime == 0) realTime = size(surfaceObject%grid)-1
        if (realTime == size(surfaceObject%grid)) realTime = 1
        sourceTime = addTime + surfaceObject%grid(realTime)%time
        return
 
      case (APERIODIC_SURFACE)

        ! Make sure the time is in the range:
        if (iTau > size(surfaceObject%grid)) then
          sourceTime = HUGE(sourceTime)
        else
          sourceTime = surfaceObject%grid(iTau)%time
        end if

        return

      case default
        call Error("ERROR: An unknown internal error has occurred.")

    end select
    
  end function GetSurfaceSourceTime


  function GetSurfaceNTau (surfaceObject) result(nTau)
    type(surface), intent(in)::surfaceObject
    integer :: nTau

    ! Determine which type of surface we have:
    select case (surfaceObject%surfaceType)
    
      case (CONSTANT_SURFACE)
        ! This one is sort of meaningless - a more appropriate answer might be
        ! infinity. I don't forsee this case being useful.
        nTau = 0
        return
        
      case (PERIODIC_SURFACE)
        ! This one is sort of meaningless - a more appropriate answer might be
        ! infinity. I don't forsee this case being useful.
        nTau = size(surfaceObject%grid)-1
        return
        
      case (APERIODIC_SURFACE)
        nTau = size(surfaceObject%grid)
        return
        
      case default
        nTau = 0
        return

    end select
    
  end function GetSurfaceNTau


  subroutine GetSurfaceTauRange (surfaceObject, tMin, tMax, tMinIndex, tMaxIndex)
    type(surface), intent(in)::surfaceObject
    real, intent(inout) :: tMin, tMax
    integer, intent(out) :: tMinIndex, tMaxIndex

    integer :: t

    if (tMin >= tMax) then
      call Error("ERROR: Invalid time range in getTauRange.")
    end if

    ! If this object has constant surface, this function is basically
    ! meaningless, and should probably not have been called.
    if (surfaceObject%surfaceType == CONSTANT_SURFACE) then
      if (debugLevel >= 1) then
        write (*,*) "WARNING: GetTauRange called for a constant surface object."
      end if
      tMinIndex = 1
      tMaxIndex = 1
      return
    end if

    ! If the surface is periodic we need to figure out what direction the
    ! searches need to go. Then we perform a linear search for the index of the
    ! timestep just below the minimum time. This is the slowest possible
    ! algorithm, but it is very easy to understand. If time here becomes an
    ! issue it can easily be replaced.
    t = 1
    if (surfaceObject%surfaceType == PERIODIC_SURFACE .and. &
        GetSurfaceSourceTime (surfaceObject, t) > tMin) then
      do  
        t = t - 1
        if (GetSurfaceSourceTime (surfaceObject, t) < tMin) then
          exit
        end if
      end do
    else
      do
        t = t + 1
        if (GetSurfaceSourceTime (surfaceObject, t) > tMin) then
          exit
        end if
      end do
      t = t - 1
    end if
    tMinIndex = t
    tMin = GetSurfaceSourceTime (surfaceObject, tMinIndex)

    ! Do basically the same thing for the maximum, but start at tMinIndex, so we
    ! always want to increment t.
    t = tMinIndex + 1
    do
      if (GetSurfaceSourceTime(surfaceObject, t) > tMax) then
        exit
      end if
      t = t + 1
    end do
    if (t > size(surfaceObject%grid) .and. &
        surfaceObject%surfaceType == APERIODIC_SURFACE) then
      t = size(surfaceObject%grid)
    end if
    tMaxIndex = t
    tMax = GetSurfaceSourceTime (surfaceObject, tMaxIndex)

  end subroutine GetSurfaceTauRange


  subroutine CreateSurfaceTauArray (newSurface, tauMin, tauMax, tau, iTauMin)
    type(surface), intent(in)::newSurface
    real, intent(inout):: tauMin, tauMax
    real, dimension (:), pointer :: tau
    integer, intent(out) :: iTauMin
    
    integer :: i, iTauMax

    if (newSurface%surfaceType == CONSTANT_SURFACE) then
      call Error ("Cannot create tau array from constant surface.")
    else
      if(debugLevel.gt.13) then
        write(*,*) 'Creating surface tau array.'
        write(*,*) 'Time range exists from ', tauMin, ' to ', tauMax
      end if
      itauMin = GetIndexAtTime (newSurface, tauMin) - 1  !ksb - need to check why -1 here
      itauMax = GetIndexAtTime (newSurface, tauMax) + 1  !ksb - need to check why +1 here
      if (newSurface%surfaceType == APERIODIC_SURFACE) then
        if (itauMin < 1) then
          itauMin = 1
        end if
        if (iTauMax > GetSurfaceNTau(newSurface)) then
          iTauMax = GetSurfaceNTau(newSurface)
        end if
      end if
      if(debugLevel.gt.13) then
        write(*,*) 'Indices of surface array are ', itauMin, ' and ', itauMax
      end if
      allocate (tau(itauMax-iTaumin+1))
      do i=1, size(tau)
        tau(i) = GetSurfaceSourceTime (newSurface, iTauMin+i-1)
      end do
    end if
  end subroutine CreateSurfaceTauArray
  
  !! This function gets the minimum time that this surface exists 
  function GetSurfaceMinTime(S) result (time)
    implicit none
    !! The surface
    type(surface), intent(in)::S
    !! The resultant max time that it exists
    real :: time
    select case (S%surfaceType)
      case (CONSTANT_SURFACE)
        time = -HUGE(time)
      case (PERIODIC_SURFACE)
        time = 0.0
      case (APERIODIC_SURFACE)
        time = S%grid(1)%time 
      case default
        call Error('Unknown surface type in GetSurfaceMaximumTime')
    end select
  end function GetSurfaceMinTime
  
  !! This function gets the maximum time that this surface exists 
  function GetSurfaceMaxTime(S) result (time)
    implicit none
    !! The surface
    type(surface), intent(in)::S
    !! The resultant max time that it exists
    real :: time
    select case (S%surfaceType)
      case (CONSTANT_SURFACE)
        time = HUGE(time)
      case (PERIODIC_SURFACE)
        time = S%period
      case (APERIODIC_SURFACE)
        time = S%grid(size(S%grid))%time 
      case default
       call Error('Unknown surface type in GetSurfaceMaximumTime')
    end select
  end function GetSurfaceMaxTime
  
  
  function GetSurfaceType (newSurface) result (surfaceType)
    type(surface), intent(in)::newSurface
    integer :: surfaceType

    surfaceType = newSurface%surfaceType
    return
    
  end function GetSurfaceType
  
  
  function GetSurfaceGridType (newSurface) result (gridType)
    type(surface), intent(in)::newSurface
    integer :: gridType

    gridType = newSurface%gridType
    return
    
  end function GetSurfaceGridType


  function IsSurfaceCompact (newSurface) result (isCompact)
    type(surface), intent(in)::newSurface
    logical :: isCompact

    if ((newSurface%gridType == STRUCTURED_GRID) .and. &
            (newSurface%iMax == 1 .or. newSurface%jMax == 1)) then
      isCompact = .true.
    else if (newSurface%gridType == UNSTRUCTURED_GRID .and. &
             newSurface%nbFaces == 2) then
      isCompact = .true.
    else
      isCompact = .false.
    end if
    
  end function IsSurfaceCompact


  function GetSurfaceNbNodes (newSurface) result (nbNodes)
    type(surface), intent(in)::newSurface
    integer :: nbNodes
    if (newSurface%gridtype == STRUCTURED_GRID) then
      nbNodes = newSurface%iMax*newSurface%jMax
    else
      !ksb debug:  If the surface is an unstructured grid, then use
      ! nbData instead of nbNodes.  nbData should be the higher number
      ! of nbFaces and nbNodes.  Not sure if this will break anything
      ! because I don't have cell centered cases.  But the one case I 
      ! have from Rick Gaeta's student at OK State it seems to work.
      ! KSB 4/15/2015
      ! nbNodes = newSurface%nbData !ksb debug: newSurface%nbnodes
      !ksb debug: removed previous change - it failed checksuite case 3_P_1
      nbNodes = newSurface%nbnodes 
    end if

  end function GetSurfaceNbNodes
  
  function GetSurfaceFirstNode (newSurface) result (firstNode)
    type(surface), intent(in)::newSurface
    integer:: firstNode
    firstNode = newSurface%FirstNode
  end function GetSurfaceFirstNode
  
  function GetSurfaceLastNode (newSurface) result (lastNode)
    type(surface), intent(in)::newSurface
    integer:: lastNode
    lastNode = newSurface%LastNode
  end function GetSurfaceLastNode
  
  subroutine SetSurfaceNodeRange (newSurface,firstNode,lastNode)
    type(surface), intent(inout)::newSurface
    integer, intent(in):: firstNode, lastNode
    newSurface%FirstNode = firstNode
    newSurface%LastNode  = lastNode
  end subroutine SetSurfaceNodeRange

  subroutine GetSurfaceDimensions (newSurface, dims)
    type(surface), intent(in)::newSurface
    integer, dimension(:) :: dims

    if (newSurface%gridType == STRUCTURED_GRID) then
      dims(1) = newSurface%iMax
      dims(2) = newSurface%jMax
    else
      dims(1) = newSurface%nbNodes
      dims(2) = newSurface%nbFaces
    end if
  end subroutine GetSurfaceDimensions


  subroutine CalculateVelocity (newSurface, keyIndex)
    type(surface), intent(inout)::newSurface
    integer:: i, t, imax
    integer, optional:: keyIndex
    ! Note that for aperiodic data we store 4 points in time:  [L C R RR].  Because of this limited
    ! storage we must add a few exceptional cases to our set of calculations.  All of the equations 
    ! employed are second order accurate.    
     
    imax=size(newSurface%grid)
    select case (newSurface%surfaceType)
      case (CONSTANT_SURFACE)
        do i=1,newSurface%nbData
          newSurface%velocity(1)%vectors(i)%A = 0
          newSurface%velocity(1)%normals(i)%A = 0
          newSurface%velocity(1)%areas(i) = 0
        end do
      case (PERIODIC_SURFACE)
        ! Do the first point:
        ! We can still central difference, but we have to wrap around. Recall that
        ! the last point and first point are the same point, so we only use one of
        ! them:
        do i=1,newSurface%nbData
          newSurface%velocity(1)%vectors(i) = &
               (newSurface%grid(2)%vectors(i)-&
                newSurface%grid(size(newSurface%grid)-1)%vectors(i))/(2.0*newSurface%dT)
          newSurface%velocity(1)%normals(i) = &
               (newSurface%grid(2)%normals(i)-&
                newSurface%grid(size(newSurface%grid)-1)%normals(i))/(2.0*newSurface%dT)
          newSurface%velocity(1)%areas(i) = &
               (newSurface%grid(2)%areas(i)-&
                newSurface%grid(size(newSurface%grid)-1)%areas(i))/(2.0*newSurface%dT)
        end do
        
        do t=2,size(newSurface%velocity)-1
          do i=1,newSurface%nbData
            newSurface%velocity(t)%vectors(i) = &
               (newSurface%grid(t+1)%vectors(i) - newSurface%grid(t-1)%vectors(i))/(2.0*newSurface%dT)
            newSurface%velocity(t)%normals(i) = &
               (newSurface%grid(t+1)%normals(i) - newSurface%grid(t-1)%normals(i))/(2.0*newSurface%dT)
            newSurface%velocity(t)%areas(i) = &
               (newSurface%grid(t+1)%areas(i) - newSurface%grid(t-1)%areas(i))/(2.0*newSurface%dT)
          end do
        end do        
        
        ! The last point will be the same as the first point, which we already
        ! calculated above
        newSurface%velocity(imax)%vectors = newSurface%velocity(1)%vectors
        newSurface%velocity(imax)%normals = newSurface%velocity(1)%normals
        newSurface%velocity(imax)%areas = newSurface%velocity(1)%areas        
      case (APERIODIC_SURFACE)
        if (keyIndex .ge. 3) then
          if (keyIndex .eq. 3) then
            ! Do the first point:
            ! We cannot central difference, so we use a 2nd order forward difference
            do i=1,newSurface%nbData
              newSurface%velocity(1)%vectors(i)%A = &
                   -3*newSurface%grid(1)%vectors(i)%A+ &
                    4*newSurface%grid(2)%vectors(i)%A- &
                    newSurface%grid(3)%vectors(i)%A / (2.0*newSurface%dt)
              newSurface%velocity(1)%normals(i)%A = &
                   -3*newSurface%grid(1)%normals(i)%A+ &
                    4*newSurface%grid(2)%normals(i)%A- &
                    newSurface%grid(3)%normals(i)%A / (2.0*newSurface%dt)
              newSurface%velocity(1)%areas(i) = &
                   -3*newSurface%grid(1)%areas(i)+ &
                    4*newSurface%grid(2)%areas(i)- &
                    newSurface%grid(3)%areas(i) / (2.0*newSurface%dt)
            end do
          end if
           
          ! In an aperiodic case we only store 4 grid points in time so we can only perform the central difference
          ! for iTau-1, as it is the only point in time with the necessary information.
          do i=1,newSurface%nbData
            newSurface%velocity(keyIndex-1)%vectors(i) = &
               (newSurface%grid(keyIndex)%vectors(i) - newSurface%grid(keyIndex-2)%vectors(i))/(2.0*newSurface%dT)
            newSurface%velocity(keyIndex-1)%normals(i) = &
               (newSurface%grid(keyIndex)%normals(i) - newSurface%grid(keyIndex-2)%normals(i))/(2.0*newSurface%dT)
            newSurface%velocity(keyIndex-1)%areas(i) = &
              (newSurface%grid(keyIndex)%areas(i) - newSurface%grid(keyIndex-2)%areas(i))/(2.0*newSurface%dT)
          end do
          
          if (newSurface%keycount .eq. newSurface%nkey) then
            ! We cannot central difference, so we use a 2nd order backward difference
            do i=1,newSurface%nbData
              newSurface%velocity(imax)%vectors(i)%A = &
                    3*newSurface%grid(imax)%vectors(i)%A-&
                    4*newSurface%grid(imax-1)%vectors(i)%A+&
                    newSurface%grid(imax-2)%vectors(i)%A / (2.0*newSurface%dt)
              newSurface%velocity(imax)%normals(i)%A = &
                    3*newSurface%grid(imax)%normals(i)%A-&
                    4*newSurface%grid(imax-1)%normals(i)%A+&
                    newSurface%grid(imax-2)%normals(i)%A / (2.0*newSurface%dt)
              newSurface%velocity(imax)%areas(i) = &
                    3*newSurface%grid(imax)%areas(i)-&
                    4*newSurface%grid(imax-1)%areas(i)+&
                    newSurface%grid(imax-2)%areas(i) / (2.0*newSurface%dt)
            end do             
          end if
        end if 
    end select

  end subroutine CalculateVelocity


  subroutine CalculateAcceleration (newSurface, keyIndex)
    type(surface), intent(inout)::newSurface
    integer:: i, t, imax, pos
    integer, optional:: keyIndex
    ! Note that for aperiodic data we store 4 points in time:  [L C R RR].  Because of this limited
    ! storage we must add a few exceptional cases to our set of calculations.  All of the equations 
    ! employed are second order accurate.

    imax=size(newSurface%acceleration)
    select case (newSurface%surfaceType)
      case (CONSTANT_SURFACE)
        do i=1,newSurface%nbData
          newSurface%acceleration(1)%vectors(i)%A = 0.0
          newSurface%acceleration(1)%normals(i)%A = 0.0
        end do
      case (PERIODIC_SURFACE)
        ! Do the first point:
        ! We can still central difference, but we have to wrap around. Recall that
        ! the last point and first point are the same point, so we only use one of
        ! them:
        do i=1,newSurface%nbData
          newSurface%acceleration(1)%vectors(i)%A = &
             (newSurface%grid(2)%vectors(i)%A - &
            2*newSurface%grid(1)%vectors(i)%A +&
              newSurface%grid(size(newSurface%grid)-1)%vectors(i)%A) / (newSurface%dT**2)
          newSurface%acceleration(1)%normals(i)%A = &
             (newSurface%grid(2)%normals(i)%A - &
            2*newSurface%grid(1)%normals(i)%A +&
              newSurface%grid(size(newSurface%grid)-1)%normals(i)%A) / (newSurface%dT**2)
        end do

        ! This function uses second order schemes at all points.        
        do t=2,size(newSurface%acceleration)-1
          do i=1,newSurface%nbData
            newSurface%acceleration(t)%vectors(i)%A = &
               (newSurface%grid(t+1)%vectors(i)%A - &
              2*newSurface%grid(t  )%vectors(i)%A +&
                newSurface%grid(t-1)%vectors(i)%A) / (newSurface%dT**2)
            newSurface%acceleration(t)%normals(i)%A = &
               (newSurface%grid(t+1)%normals(i)%A - &
              2*newSurface%grid(t  )%normals(i)%A +&
                newSurface%grid(t-1)%normals(i)%A) / (newSurface%dT**2)
          end do
        end do   
        
        ! Do the last point:
        ! The last point will be the same as the first point, whcih we already
        ! calculated above
        newSurface%acceleration(imax)%vectors = &
           newSurface%acceleration(1)%vectors
        newSurface%acceleration(imax)%normals = &
           newSurface%acceleration(1)%normals             
      case (APERIODIC_SURFACE)
        if (newSurface%keycount .ge. 4) then
          if (newSurface%keycount .eq. 4) then
            ! Do the first point:
            ! We cannot central difference, so we use a 2nd order forward difference
            do i=1,newSurface%nbData
              newSurface%acceleration(1)%vectors(i)%A = &
                     -newSurface%grid(4)%vectors(i)%A+&
                    4*newSurface%grid(3)%vectors(i)%A-&
                    5*newSurface%grid(2)%vectors(i)%A+&
                    2*newSurface%grid(1)%vectors(i)%A/ (newSurface%dt**2)
              newSurface%acceleration(1)%normals(i)%A = &
                     -newSurface%grid(4)%normals(i)%A+&
                    4*newSurface%grid(3)%normals(i)%A-&
                    5*newSurface%grid(2)%normals(i)%A+&
                    2*newSurface%grid(1)%normals(i)%A/ (newSurface%dt**2)
            end do
            
            ! Do the second point (the next eqn will handle the third and so on)
            do i=1,newSurface%nbData
              newSurface%acceleration(2)%vectors(i)%A = &
                 (newSurface%grid(3)%vectors(i)%A - &
                2*newSurface%grid(2)%vectors(i)%A +&
                  newSurface%grid(1)%vectors(i)%A) / (newSurface%dT**2)
              newSurface%acceleration(2)%normals(i)%A = &
                 (newSurface%grid(3)%normals(i)%A - &
                2*newSurface%grid(2)%normals(i)%A +&
                  newSurface%grid(1)%normals(i)%A) / (newSurface%dT**2)
            end do            
          end if

          ! In an aperiodic case we store 4 grids in time but we only perform the central difference
          ! for iTau, as we've already calulated iTau-1 and we can't yet calculate for iTau+1    
          do i=1,newSurface%nbData
            newSurface%acceleration(keyIndex-1)%vectors(i)%A = &
               (newSurface%grid(keyIndex  )%vectors(i)%A - &
              2*newSurface%grid(keyIndex-1)%vectors(i)%A +&
                newSurface%grid(keyIndex-2)%vectors(i)%A) / (newSurface%dT**2)
            newSurface%acceleration(keyIndex-1)%normals(i)%A = &
               (newSurface%grid(keyIndex  )%normals(i)%A - &
              2*newSurface%grid(keyIndex-1)%normals(i)%A +&
                newSurface%grid(keyIndex-2)%normals(i)%A) / (newSurface%dT**2)
          end do     
          if (newSurface%keyCount .eq. newSurface%nkey) then
            ! We cannot central difference, so we use a 2nd order backward difference
            do i=1,newSurface%nbData
              newSurface%acceleration(imax)%vectors(i)%A = &
                     -newSurface%grid(imax-3)%vectors(i)%A+&
                    4*newSurface%grid(imax-2)%vectors(i)%A-&
                    5*newSurface%grid(imax-1)%vectors(i)%A+&
                    2*newSurface%grid(imax )%vectors(i)%A/ newSurface%dt**2
              newSurface%acceleration(imax)%normals(i)%A = &
                     -newSurface%grid(imax-3)%normals(i)%A+&
                    4*newSurface%grid(imax-2)%normals(i)%A-&
                    5*newSurface%grid(imax-1)%normals(i)%A+&
                    2*newSurface%grid(imax )%normals(i)%A/ newSurface%dt**2
             end do
          end if        
        end if
    end select

  end subroutine CalculateAcceleration

!!! Currently Unused
  subroutine CalculateAcceleration4thOrder (newSurface)
    type(surface), intent(inout)::newSurface
    integer :: t, i, imax
    ! This function uses fourth order schemes at all points except the ends,
    ! which use second order.
    ! First, allocate the necessary memory for the number of timesteps that we
    ! have.
    allocate (newSurface%acceleration(size(newSurface%grid)))
    do t=3,size(newSurface%acceleration)-2
      allocate(newSurface%acceleration(t)%vectors(newSurface%nbData))
      allocate(newSurface%acceleration(t)%normals(newSurface%nbData))
      do i=1,newSurface%nbData
        newSurface%acceleration(t)%vectors(i)%A =  &
           (  -newSurface%grid(t+2)%vectors(i)%A + &
            16*newSurface%grid(t+1)%vectors(i)%A - &
            30*newSurface%grid(t  )%vectors(i)%A + &
            16*newSurface%grid(t-1)%vectors(i)%A - &
               newSurface%grid(t-2)%vectors(i)%A) / (12*newSurface%dT**2)
        newSurface%acceleration(t)%normals(i)%A =  &
           (  -newSurface%grid(t+2)%normals(i)%A + &
            16*newSurface%grid(t+1)%normals(i)%A - &
            30*newSurface%grid(t  )%normals(i)%A + &
            16*newSurface%grid(t-1)%normals(i)%A - &
               newSurface%grid(t-2)%normals(i)%A) / (12*newSurface%dT**2)
      end do
    end do

    ! Do the second and second-to-last points as second order
    if (size(newSurface%grid) > 3) then
      allocate(newSurface%acceleration(2)%vectors(newSurface%nbData))
      allocate(newSurface%acceleration(size(newSurface%grid)-1)%vectors(newSurface%nbData))
      allocate(newSurface%acceleration(2)%normals(newSurface%nbData))
      allocate(newSurface%acceleration(size(newSurface%grid)-1)%normals(newSurface%nbData))
    
      t = 2
      do i=1,newSurface%nbData
        newSurface%acceleration(t)%vectors(i)%A = &
           (newSurface%grid(t+1)%vectors(i)%A - &
          2*newSurface%grid(t  )%vectors(i)%A +&
            newSurface%grid(t-1)%vectors(i)%A) / (newSurface%dT**2)
        newSurface%acceleration(t)%normals(i)%A = &
           (newSurface%grid(t+1)%normals(i)%A - &
          2*newSurface%grid(t  )%normals(i)%A +&
            newSurface%grid(t-1)%normals(i)%A) / (newSurface%dT**2)
      end do
      
      t = size(newSurface%grid)-1
      do i=1,newSurface%nbData
        newSurface%acceleration(t)%vectors(i)%A = &
           (newSurface%grid(t+1)%vectors(i)%A - &
          2*newSurface%grid(t  )%vectors(i)%A +&
            newSurface%grid(t-1)%vectors(i)%A) / (newSurface%dT**2)
        newSurface%acceleration(t)%normals(i)%A = &
           (newSurface%grid(t+1)%normals(i)%A - &
          2*newSurface%grid(t  )%normals(i)%A +&
            newSurface%grid(t-1)%normals(i)%A) / (newSurface%dT**2)
      end do
      
    end if
    
    ! Do the first point:
    allocate(newSurface%acceleration(1)%vectors(newSurface%nbData))
    allocate(newSurface%acceleration(1)%normals(newSurface%nbData))
    select case (newSurface%surfaceType)
      case (CONSTANT_SURFACE)
        do i=1,newSurface%nbData
          newSurface%acceleration(1)%vectors(i)%A = 0
          newSurface%acceleration(1)%normals(i)%A = 0
        end do
      case (PERIODIC_SURFACE)
        ! We can still central difference, but we have to wrap around. Recall that
        ! the last point and first point are the same point, so we only use one of
        ! them:
        do i=1,newSurface%nbData
          newSurface%acceleration(1)%vectors(i)%A = &
             (newSurface%grid(2)%vectors(i)%A - &
            2*newSurface%grid(1  )%vectors(i)%A +&
              newSurface%grid(size(newSurface%grid)-1)%vectors(i)%A) / (newSurface%dT**2)
          newSurface%acceleration(1)%normals(i)%A = &
             (newSurface%grid(2)%normals(i)%A - &
            2*newSurface%grid(1  )%normals(i)%A +&
              newSurface%grid(size(newSurface%grid)-1)%normals(i)%A) / (newSurface%dT**2)
        end do
      case (APERIODIC_SURFACE)
        ! We cannot central difference, so we use a 2nd order forward difference
        do i=1,newSurface%nbData
          newSurface%acceleration(1)%vectors(i)%A = &
                 -newSurface%grid(4)%vectors(i)%A+&
                4*newSurface%grid(3)%vectors(i)%A-&
                5*newSurface%grid(2)%vectors(i)%A+&
                2*newSurface%grid(1)%vectors(i)%A/ (newSurface%dt**2)
          newSurface%acceleration(1)%normals(i)%A = &
                 -newSurface%grid(4)%normals(i)%A+&
                4*newSurface%grid(3)%normals(i)%A-&
                5*newSurface%grid(2)%normals(i)%A+&
                2*newSurface%grid(1)%normals(i)%A/ (newSurface%dt**2)
        end do
    end select

    ! Do the last point:
    imax=size(newSurface%acceleration)
    select case (newSurface%surfaceType)
      case (PERIODIC_SURFACE)
        ! The last point will be the same as the first point, which we already
        ! calculated above
        allocate(newSurface%acceleration(imax)%vectors(newSurface%nbData))
        allocate(newSurface%acceleration(imax)%normals(newSurface%nbData))
        newSurface%acceleration(imax)%vectors = &
           newSurface%acceleration(1)%vectors
        newSurface%acceleration(imax)%normals = &
           newSurface%acceleration(1)%normals
      case (APERIODIC_SURFACE)
        allocate(newSurface%acceleration(imax)%vectors(newSurface%nbData))
        allocate(newSurface%acceleration(imax)%normals(newSurface%nbData))
        ! We cannot central difference, so we use a 2nd order backward difference
        do i=1,newSurface%nbData
          newSurface%acceleration(1)%vectors(i)%A = &
                 -newSurface%grid(imax-3)%vectors(i)%A+&
                4*newSurface%grid(imax-2)%vectors(i)%A-&
                5*newSurface%grid(imax-1)%vectors(i)%A+&
                2*newSurface%grid(imax-0)%vectors(i)%A/ newSurface%dt**2
          newSurface%acceleration(1)%normals(i)%A = &
                 -newSurface%grid(imax-3)%normals(i)%A+&
                4*newSurface%grid(imax-2)%normals(i)%A-&
                5*newSurface%grid(imax-1)%normals(i)%A+&
                2*newSurface%grid(imax-0)%normals(i)%A/ newSurface%dt**2
        end do
    end select

  end subroutine CalculateAcceleration4thOrder
 
 
  subroutine RotateSurfaceArray(S)
    type(Surface), pointer:: S
    type(singleSurface):: tempGrid, tempVel, tempAccel
    integer:: i
    ! Because PSU-WOPWOP employs pointers for its arrays we can't simply rotate 
    ! the 'surface' arrays. We must reassign the arrays within Surface at the lowest
    ! level of pointers. 
   
    tempGrid  = S%grid(1)
    do i=1,size(S%grid)-1
      S%grid(i) = S%grid(i+1)
    end do
    S%grid(size(S%grid)) = tempGrid
      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tempVel       = S%velocity(1)
    do i=1,size(S%velocity)-1
      S%velocity(i) = S%velocity(i+1)
    end do
    S%velocity(size(S%velocity)) = tempVel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tempAccel         = S%acceleration(1)
    do i=1,size(S%acceleration)-1
      S%acceleration(i) = S%acceleration(i+1)
    end do
    S%acceleration(size(S%acceleration)) = tempAccel
      
    return
  end subroutine RotateSurfaceArray
  
  subroutine ResetSurfaceKeyCount(S,val)
    type(surface):: S
    integer, intent(in), optional:: val
    if( present(val) ) then
      s%keycount = val
    else
      S%keycount = 0
    end if
  end subroutine ResetSurfaceKeyCount
  
  subroutine SetSurfaceNKey(S,val)
    type(surface):: S
    integer, intent(in), optional:: val
    if( present(val) ) then
      S%NKey = val
    end if
  end subroutine SetSurfaceNKey
  
  subroutine GetSurfaceUnstructuredFaces (S, faces)
    type(surface), intent(inout)::S
    integer, dimension(:,:), pointer :: faces

    faces => S%connectivity
  end subroutine GetSurfaceUnstructuredFaces
  
  logical function IsNodeBased(S) result(val)
    type(surface):: S
  
    if( S%normType == GEO_NORMTYPE_NODES ) then
      val = .true.
    else
      val = .false.
    end if
    return
  end function IsNodeBased
  
  function GetSurfacenKey(surfaceObject) result(nKey)
    type(surface):: surfaceObject
    integer:: nKey
    
    nKey = surfaceObject%nKey
    return
  end function GetSurfacenkey   

  function GetSurfaceDt(surfaceObject) result(dt)
    type(surface), intent(in) :: surfaceObject
    real::dt

    dt = surfaceObject%dt
    return
  end function GetSurfaceDt
  
  function GetSurfaceTimeStepOffset(surfaceObject) result(timeStepOffset)
    type(surface), intent(in):: surfaceObject
    integer:: timeStepOffset
    
    timeStepOffset = surfaceObject%timeStepOffset
  end function GetSurfaceTimeStepOffset
  
  function GetSurfaceGridSize(surfaceObject) result(gridSize)
    type(surface), intent(in):: surfaceObject
    integer:: gridSize
    
    gridSize = size(surfaceObject%grid)
  end function GetSurfaceGridSize
  
  function GetSurfaceIblanking(surfaceObject) result(iblanking)
    type(surface), intent(in):: surfaceObject
    logical:: iblanking
    
    iblanking = surfaceObject%iblanking
  end function GetSurfaceIblanking
  
end module surfaceObject