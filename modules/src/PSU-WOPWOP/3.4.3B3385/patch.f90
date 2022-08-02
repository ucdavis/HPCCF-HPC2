module patchObject
! PSU-WOPWOP
! $Id: patch.f90 3375 2017-10-10 14:54:20Z brentner $
! $LastChangedDate: 2017-10-10 10:54:20 -0400 (Tue, 10 Oct 2017) $
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
! DARPA Helicopter Quieting Program
! Awarding Agency: Army Research Office
! Prime Award No: W911NF-04-1-0419 (to Georgia Tech)
! Subcontract Award No: E-16-X20-G1
!
!
! Written by Guillaume Bres, Guillaume Perez, Leonard Lopes, Hsuan-Nien Chen, 
! Christopher Hennes, Rui Cheng and Benjamin Goldman. 
! Faculty advisor Dr. Kenneth S. Brentner.


   
!**********************************************************************
!MODULE patch
!  This module contains all the patch definitions and subroutines.  A patch 
! is the lowest level of object we contain and all the integration happens
! here.  
!**********************************************************************

!****** Patch
! NAME
!   Patch -- A module containing the Patch object and associated functions.
! DESCRIPTION
!   The Patch modules defines the patch object and its associated functions.
! METHODS
!   DestroyPatch
!   CreatePatch
!   CopyPatch
!   InitializeSigma
!   IntegratePatch
!   ObserverTime
!   TimeError
!   SourceTime
!   TauError
!   CalculatepressureParameters
!   NodeIntegration
!   GenerateThicknessNoise
!   GenerateLoadingNoise
!   TimeInterpol
!   AddSignalByFrequency
!   ReadFlowData
!   ReadLoading
!   ReadSourceTime
!   CreateSourceTimeArray
!   DerivativeL
!   WritePatch
!   RecordSigma
!   GetPatchDimensions
!   GetPatchImplicitTimeRange
!   RotatePressureParameters
!   DestroyDataParameters
! MODULE GLOBAL VARIABLES
!   mg_obsCBList
!   mg_srcCBList
!   mg_srcInitialVector
!   mg_srcVector
!   mg_obsInitialVector
!   mg_radiationVector
!   mg_obsVector
!   mg_tauLocal 
!   mg_tLocal
!   mg_flowParametersL
!   mg_flowParametersC
!   mg_flowParametersR
!   mg_initialized
!***
  use COBObject
  use constantsModule
  use frequencyAnalysisModule , only: ConvertToSpectrum, SumSignalTogether, AddSignal !, SumSignalTogetherVector 
  use surfaceObject
  use loadingObject
  use wallObject
  use IOModule
  use MPIModule
  use debugModule
  use mathModule
  use broadbandModule
  !ksb debug: 10/10/2017
  use observerObject, only: firstObsSetup

  implicit none

  private
  public:: ReadPatchAtTime, ReadPatchGeomAtTime, ReadPatchLoadAtTime
  public:: patch, CreatePatch, CopyPatch, DestroyPatch, DeterminePatchDtau
  public:: InitializeSigma, IsStructured
  public:: GetPatchImplicitTimeRange, WritePatchPlot3DStructuredSigma
  public:: GetPatchDimensions, WritePatchFVUNSHeader, WritePatchFVUNSData
  public:: WritePatchDebugInfo, IntegratePatch, GetPatchNbNodes
  public:: GetAssociatedLoad, GetPatchSurfacenKey, GetPatchLoadnKey
  public:: GetPatchNTau, GetPatchDTau, RotatePatchTimeData, GetPatchTimeArrays
  public:: ResetGlobals, InitializeStoredData, SetPatchNkey, SetPatchNodeRange
  public:: GetGeoTitle, GetLoadTitle, GetDimensions, GetLoading, GetPatchITauMin
  public:: ResetPatchKeyCount, GetPatchTauArray, CheckPatchNSect
  public:: GetNKey, GetPatchMinObsTimeInData, GetPatchMaxObsTimeInData, GetPatchGeoNKey
  public:: CreatePatchObsTimeArrays, CreatePatchTimeRangeArrays, GetPatchMaxSourceTime
  public:: ComputePatchThrust, ComputePatchArea, ComputePatchSpan, ComputePatchSectLength
  public:: PredictPeggBB, CompareDataTimeInfo, PatchITauInRange, PatchLastITau
  
  !****************************************************************************
  !TYPE PATCH
  ! - iMax: Integer, 1st dimension of the patch (FROM INPUT DATA).
  ! - jMax: Integer, 2nd dimension of the patch (FROM INPUT DATA).
  ! - nTau: Integer, number of source time points (FROM INPUT DATA).
  ! - coordinate: Array of vectors, dimension (iMax, jMax), coordinates 
  !          of each node on the patch (FROM INPUT DATA).
  ! - normal: Array of vectors, dimension (iMax,jMax), normal vectors 
  !      at each node of the patch (COMPUTED).
  ! - tangent: Array of vectors, dimension (iMax,jMax), tangent vectors 
  !         at each node of the patch (COMPUTED).
  ! - loading: Array of vectors, dimension (iMax, jMax, nTau,2), normal and
  !         tangent components of the loading at each node on the patch, 
  !         for each source time (FROM INPUT DATA or COMPUTED).
  ! - tau: array of real, dimension(nTau), source time history (FROM 
  !        INPUT DATA or COMPUTED).
  ! - dS: array of real, dimension(iMax,jMax), corresponding area or 
  !       linear contribution for each node (COMPUTED).
  ! - computeThicknessFlag, computeLoadingFlag: logical, flag to know if the thickness or
  !                       loading noise, or both, need to be computed
  ! - numberVariables: Number of sigma surface default variables, usually 5(x, y, z, t, tau) 
  !         where x, y, z are the spacial coordinates, t is the observer time, 
  !         and tau is source time.
  ! - numberFunctions: Number of extra functions that we want to save, a few are tangent, normal,
  !        total noise.
  ! - sigmaSurface: A large array dimension(numberVariables+numberFunctions, nbNodes, ntau) that will
  !                 store all information we want to see later in fieldview or another
  !                 program, see corresponding papers for information. 
  !**********************************************************************
  type patch
     private  
     character(len=32)::geoTitle,loadTitle, prefix
     integer:: loadingOutput ! this variable is added to put the loading noise somewhere
                             ! other than loading; e.g., in thickness noise.
                             ! when this is 1 put loadin into loading; if 2 put it in thickness
     integer:: nTau, numberVariables, numberFunctions, iTau
     integer:: iTauMin,iMax,jMax
     integer :: prefixLength, refFrame
     integer:: nbNodes, firstNode, lastNode, loadSpecies, nSourceBase, srcIndex
     logical::isCopy, flowData
     logical::computeThicknessFlag, computeLoadingFlag, trailingEdgeFlag, computePressGradFlag
     logical::surfaceAperiodic, loadAperiodic
     real:: dtau, tauMin, tauMax, tauLocal
     real, dimension(:), pointer::tau =>NULL()
     real, dimension(:,:), pointer:: ObsTimeHistory =>NULL(), qTdS =>NULL(), qLdS =>NULL()
     real, dimension(:,:,:), pointer::sigmaSurface =>NULL()
     type (vector):: srcVector     
     type (surface), pointer :: surfaceObject =>NULL()
     type (loading), pointer :: loadObject =>NULL()
     type (StoredData), dimension(:), pointer:: StoredDataArray =>NULL()
     type (CBStructure), pointer :: srcCBList => NULL(), obsCBList => NULL()
     type(vector), dimension(:), pointer :: accelerationsaveLL =>NULL(), accelerationsaveL =>NULL(), accelerationsaveC =>NULL()
     type(vector), dimension(:), pointer :: nhatdotsaveLL =>NULL(),      nhatdotsaveL =>NULL(),      nhatdotsaveC =>NULL()
     type(obsTimeRange), dimension(:), pointer:: ObsTimeInRange =>NULL()
  end type patch
  
  !*************************
  !TYPE ObsTimeRange
  ! tMin stores the minimum observer time from each patch at each source time
  ! tMax stores the maximum observer time from each patch at each source time
  type ObsTimeRange
    real(kind=8), dimension(:), pointer:: tMin =>NULL(), tMax =>NULL()
  end type ObsTimeRange

  !****************************************************************************
  !TYPE PressureParameters
  ! - U: Ui is defined as [1-(rho/rho0)]vi + (rho*ui/rho0)
  ! - L: Li is defined as Pij * nhatj + rho * ui * (un - vn)
  !****************************************************************************
  type DataParameters
    real :: time
    type(vector), dimension(:), pointer :: U =>NULL()
    type(vector), dimension(:), pointer :: L =>NULL()
    type(vector), dimension(:), pointer :: nhat =>NULL()
  end type DataParameters

  !**************************************************************************
  !TYPE: StoredData
  !  This data type stores the information required for calculating
  !  derivatives and noise terms for a single patch using only a few
  !  time steps.
  !***************************************************************************     
  type StoredData
    private
    real:: tLocal, initialr
    real(kind=8), dimension(:), pointer :: qTdSL =>NULL(),             qTdSC =>NULL()
    real(kind=8), dimension(:), pointer :: qLdSL =>NULL(),             qLdSC =>NULL()
    real(kind=8), dimension(:), pointer :: obsTimeHistoryL =>NULL(),   obsTimeHistoryC =>NULL()
    type(vector) :: radiationVector, obsVector 
    type(vectorKind8), dimension(:), pointer :: dPdTL =>NULL(),   dPdTC =>NULL()
    type(vectorKind8), dimension(:), pointer :: dPdL1L =>NULL(), dPdL1C =>NULL()
    type(DataParameters), pointer :: parameterL =>NULL(), parameterC =>NULL(), parameterR =>NULL()
    integer:: radIndex
  end type StoredData
 
  !**********************************************************************************
  !  We are using some outside routines for estimating observer time for a moving observer
  ! and it became necessary to store some outside values here that are used in obstime
  ! routine.
  !**********************************************************************************    
  integer                       :: mg_srcIndex, mg_patchNum, mg_radIndex
  real                          :: mg_tauLocal, mg_tLocal 
  type(vector)                  :: mg_srcVector, mg_obsInitialVector, mg_radiationVector, mg_obsVector 
  type(Surface)       , pointer :: mg_srcSurface => NULL()
  type(CBStructure)   , pointer :: mg_obsCBList => NULL() , mg_srcCBList => NULL()
  type(DataParameters), pointer :: mg_flowParametersL =>NULL() ,mg_flowParametersC =>NULL(), &
                                   mg_flowParametersR =>NULL()

  integer, parameter :: NO_SURFACE_DATA = 0,       &
                        PRESSURE_SURFACE_DATA = 1, &
                        LOADING_SURFACE_DATA = 2,  &
                        FLOW_SURFACE_DATA = 3
contains 

  subroutine WritePatchDebugInfo(P,unitnum)
    type(patch), intent(in)::P
    integer::unitnum
    write(unitnum,*) '*** PatchDebugInfo ***'
    write(unitnum,*) 'Ntau= ', trim(integertostring(P%ntau))
    write(unitnum,*) 'NumberVariables= ', trim(integertostring(P%numberVariables))
    write(unitnum,*) 'NumberFunctions= ', trim(integertostring(P%numberFunctions))
    write(unitnum,*) 'iTauMin= ', trim(integertostring(P%iTauMin))
    write(unitnum,*) 'iMax= ', trim(integertostring(P%iMax))
    write(unitnum,*) 'jMax= ', trim(integertostring(P%jMax))
    write(unitnum,*) 'dtau= ', trim(realtostring(P%dtau))
    write(unitnum,*) 'isCopy= ', P%isCopy
    write(unitnum,*) 'computeThicknessFlag= ', P%computeThicknessFlag
    write(unitnum,*) 'computeLoadingFlag= ', P%computeLoadingFlag
    write(unitnum,*) 'trailingEdgeFlag= ', P%trailingEdgeFlag
    write(unitnum,*) 'computePressGradFlag= ', P%computePressGradFlag
    write(unitnum,*) 'geoTitle= ', trim(P%geoTitle)
    write(unitnum,*) 'loadTitle= ', trim(P%loadTitle)
    write(unitnum,*) 'prefix= ', trim(P%prefix)
    write(unitnum,*) 'prefixLength= ', trim(integertostring(P%prefixLength))
    write(unitnum,*) 'refFrame= ', P%refFrame
    write(unitnum,*) 'surfaceAperioc= ',P%surfaceAperiodic
    write(unitnum,*) 'loadAperiodic= ',P%loadAperiodic
    write(unitnum,*) 'associated surface? ', associated(P%surfaceObject)
    if(associated(P%surfaceObject)) then
      call WriteSurfaceDebugInfo(P%surfaceObject,unitnum)
    end if  
    write(unitnum,*) 'associated loading? ', associated(P%loadObject)
    if(associated(P%loadObject)) then
      call WriteLoadDebugInfo(P%loadObject,unitnum)
    end if  
    write(unitnum,*) 'associated tau? ', associated(P%tau)
    if (associated(P%tau)) then
      write(unitnum,*) 'size tau=', trim(integertostring(size(P%tau)))
    end if
    write(unitnum,*) 'associated sigmaSurface? ', associated(P%sigmaSurface)
    if (associated(P%sigmaSurface)) then
      write(unitnum,*) 'size sigmaSurface= ', trim(integertostring(size(P%sigmaSurface)))
    end if
    write(unitnum,*) '--- End Patch Debug Info ---'

  end subroutine WritePatchDebugInfo
  
  !**********************************************************************************
  !SUBROUTINE destroyPatch(P)
  !  This routime will just go through the patch P and destroy all of it's allocated 
  ! array, neccessary for clean up.
  !ARGUEMENTS:
  ! - P: I patch we want to destroy.
  !**********************************************************************************
  subroutine destroyPatch(P)
    type(patch), intent(inout)::P
    integer:: i
    
    call DestroySurface(P%surfaceObject)
    deallocate (P%surfaceObject)
    
    if(associated(P%loadObject)) then
      call DestroyLoading (P%loadObject)
      deallocate (P%loadObject)
    end if
    

    if(associated(P%storedDataArray)) then
      do i=0, mg_NbWall
        call DestroyStoredData(P, P%storedDataArray(i))
      end do
      deallocate(P%storedDataArray)
    end if
    !if(associated(P%/)) then
    !  call DestroyStoredData(P, P%storedDataArray) 
    !  deallocate(P%storedDataArray)
    !end if    
    !if (mg_NbWall.gt.0) then
    !  do i=1,mg_NbWall
    !    call DestroyStoredData(P, P%storedWallArray(i))
    !  end do
    !  deallocate(P%StoredWallArray)
    !end if
  
    if(associated(P%sigmaSurface)) deallocate(P%sigmaSurface)
    if(associated(P%tau))          deallocate (P%tau)
    !if(associated(P%ObsTimeInRange)) then !ksb debug: result of merge - not sure what to do
    !  if(associated(P%ObsTimeInRange(1)%tMin)) then
    !    do i=1,size(P%obsTimeInRange)
    !      deallocate(P%obsTimeInRange(i)%tMin)
    !      deallocate(P%obsTimeInRange(i)%tMax)
    !    end do
    !  end if
    !  deallocate(P%ObsTimeInRange)
    !end if
    call ResetGlobals()
  end subroutine destroyPatch
  
  !ksb16 debug:  These comments are out of data 2/11/2014
  !**********************************************************************************
  !SUBROUTINE createPatch(P,file1,file2)  !ksb debug: the comments are out of date
  !  The subroutine reads the grid coordinates in the file "name1", and the source
  ! time history in the file "name2" (if present), and computes the others components 
  ! of the patch.
  !ARGUMENTS:
  ! - P: Patch we are creating.
  ! - file1: String of character, name of the grid geometry file.
  ! - file2: String of character, OPTIONAL, name of the source time history file.
  !REMARKS:
  ! - If file2 is present, the source time history is read from the file.
  !   "name2" by the subroutine "readSourceTime". Otherwise, the source time
  !   history is computed by the subroutine "setsourcetime".
  ! - once the patch coordinates are read, the normal vectors and
  !   integration area are computed by the subroutines "normalVector" and "area".
  !SUBROUTINES USED:
  ! - "readCoordinates","createNormalVectors","area","unitNormal","readSourceTime","setsourcetime",
  !    createSigma(P) in patchmod.f90"
  !**********************************************************************************
  subroutine CreatePatch(newPatch, surfaceFileUnit, surfMajorVer, surfMinorVer, geoInfo, &
                         loadFileUnit, loadMajorVer, loadMinorVer,loadInfo, &
                         hasLoading, compThickness, prefix) 
    type(patch), intent(inout):: newPatch
    character(len=*), intent(in) :: prefix
    integer:: surfaceFileUnit, loadFileUnit, surfType, stat
    integer,intent(in):: surfMajorVer, surfMinorVer, loadMajorVer, loadMinorVer
    logical:: hasLoading, compThickness
    integer,dimension(:):: geoInfo,loadInfo
    integer, save :: patchIndex = 0  ! ksb debug:  I think this should be defined in the module
                                     ! and treated as a module global variable.  But it looks like
                                     ! it is only used for the default patch title. 
    ! This is going to create a new patch
    ! Make sure that all the pointers are nullified before we start

    nullify (newPatch%surfaceObject, newPatch%loadObject, &
             newPatch%tau, newPatch%sigmaSurface)

    ! Reset the member global variables
    newPatch%ntau=0
    newPatch%dTau=0.0
    newPatch%prefix=prefix
    newPatch%prefixLength=len(prefix)
    newPatch%isCopy = .false.
    newPatch%flowData = .false.
    if (loadMajorVer==VERSION_0) then
      newPatch%refFrame = LOAD_REFFRAME_MIXED
    else
      newPatch%refFrame = loadInfo(LOAD_REFFRAME)
    end if

    ! Create the default geometry name:
    patchIndex = patchIndex + 1
    newPatch%geoTitle  = "Patch "//trim(IntegerToString(patchIndex))
    newPatch%loadTitle = "Patch "//trim(IntegerToString(patchIndex))
    
    ! Set default loadingOutput
    newPatch%loadingOutput = 2  !normal loading output
   
    call ResetGlobals()
    ! Read in the patch name from the geometry (and function data) file dynamic header 
    ! (i.e., in the "Structured" or "Unstructured" header) for this patch (zone).
    if (surfMajorVer==VERSION_1) then
      call ReadBinaryString(surfaceFileUnit,newPatch%geoTitle, stat)
      call Message('  Geometry: '//prefix//trim(newPatch%geoTitle))
      if (hasLoading) then
        call ReadBinaryString(loadFileUnit,newPatch%loadTitle, stat)
        call Message('  Loading:  '//prefix//trim(newPatch%loadTitle))
        if(loadMinorVer==VERSION_1) then
          call ReadBinaryInteger(loadFileUnit,newPatch%loadingOutput, stat)
        end if
      end if 
    end if

    allocate (newPatch%surfaceObject)
    call CreateSurface (newPatch%surfaceObject,surfaceFileUnit,0.0, &
                       geoInfo, newPatch%geoTitle)
    ! provide the patch with key information (that is determined by the surface)
    newPatch%nTau      = getSurfaceNkey(newPatch%surfaceObject) 
    newPatch%nbNodes   = GetSurfaceNbNodes (newPatch%surfaceObject)
    newPatch%firstNode = GetSurfaceFirstNode (newPatch%surfaceObject)
    newPatch%lastNode  = GetSurfaceLastNode (newPatch%surfaceObject)
   !ksb debug: this is redundant: newPatch%nTau      = getSurfaceNkey(newPatch%surfaceObject) 
    
    select case(GetSurfaceType(newPatch%surfaceObject))
    case (APERIODIC_SURFACE)
      newPatch%surfaceAperiodic=.true.
    case default
      newPatch%surfaceAperiodic=.false.
    end select
    if (surfMajorVer==VERSION_0) then
        surfType = DetermineAcousticSurfaceType (loadFileUnit)
        if (surfType .ne. NO_SURFACE_DATA) then
          hasLoading = .TRUE.
        end if
    else
      if (hasLoading) then
        select case (loadInfo(LOAD_DATATYPE))
        case (LOAD_DATATYPE_PRESSURE)
          surfType = PRESSURE_SURFACE_DATA
        case (LOAD_DATATYPE_LOADING) 
          surfType = LOADING_SURFACE_DATA
        case (LOAD_DATATYPE_FLOW) 
          surfType = FLOW_SURFACE_DATA
        case default
          call Error ("Invalid load data type ("//&
                      trim(IntegerToString(loadInfo(LOAD_DATATYPE)))//")")
        end select
      else
        surfType = NO_SURFACE_DATA
      end if
    end if
    ! Next, figure out what to do with the loading

    if ((surfType.ne.NO_SURFACE_DATA).and.(loadingNoiseFlag .or. totalNoiseFlag)) then
      newPatch%computeLoadingFlag = .true.
      allocate (newPatch%loadObject)
      call CreateLoading (newPatch%loadObject,loadFileUnit,loadinfo)
      newPatch%nTau = max(getLoadNkey(newPatch%loadObject),newPatch%nTau)       
      select case(GetLoadType(newPatch%loadObject))
      case (APERIODIC_SURFACE)
        newPatch%loadAperiodic=.true.
      case default
        newPatch%loadAperiodic=.false.
      end select
    else
      newPatch%computeLoadingFlag = .false.
      newPatch%loadAperiodic = .false.
    end if

    ! If this is not a compact patch and the user has said to calculate either
    ! thickness or total noise, set that flag to true.   
    if (.not. IsSurfaceCompact (newPatch%surfaceObject) .and.&
          (thicknessNoiseFlag .or. totalNoiseFlag) .and. compThickness) then
        newPatch%computeThicknessFlag = .true.
    else
        newPatch%computeThicknessFlag = .false.
    end if

    !Determine pressure gradient logic
    if (PressureGradientFlag) then
      newPatch%computePressGradFlag = .true.
    else
      newPatch%computePressGradFlag = .false.
    end if

    newPatch%trailingEdgeFlag=.false.
    newPatch%itau=0

  end subroutine CreatePatch


  subroutine ReadPatchAtTime(newPatch, surfaceFileUnit,loadFileUnit,& 
                             geoInfo, loadInfo, initialAzimuth,&
                             hasLoading, vFlag)
                           
    type(Patch), intent(inout)::newPatch
    integer::surfaceFileUnit,loadFileUnit,i
    integer,intent(in)::vFlag
    integer,dimension(:)::loadInfo,geoInfo
    logical:: hasLoading
    real::initialAzimuth
    
    integer::SurfnKey, LoadnKey
    
    ! V0_0 files are organized such that all of the data in time for each patch is 
    ! listed immediately after the respective patch dimension information, thus it 
    ! must be gathered all at once (where SurfnKey  and LoadnKey are the number of 
    ! timesteps).  V1_0 files are organized such that the data in time for
    ! each patch is listed together and after all of the patch dimension information, 
    ! thus the Read( )AtTime subroutines need only be called once.  (The outer loop 
    ! within container is used in place of SurfnKey and LoadnKey for V1_0).
    
    if (vFlag == VERSION_0) then
      SurfnKey = GetSurfacenKey(newPatch%surfaceObject)
      if (hasloading) then
        Loadnkey = GetLoadnKey(newPatch%loadObject)
      end if
    else
      SurfnKey = 1
      LoadnKey = 1
    end if
    
    do i = 1,SurfnKey
      ! ReadSurfaceAtTime reads the surface (geometry) data AND
      ! if this is the last time step, it compute other things, like the times, 
      ! velocity and acceleration vectors, etc.  The variable eflag helps control that.
      call ReadSurfaceAtTime(newPatch%surfaceObject,SurfaceFileUnit,&
                             initialAzimuth, geoInfo(GEO_TIMETYPE))
    end do
    
    if (hasloading .and. (loadingNoiseFlag .or. totalNoiseFlag) ) then    
      do i = 1,LoadnKey
        ! ReadLoadingAtTime reads the loading (function data) AND
        ! if this is the last time step, it compute other things, like dt (that may be all)
        ! surface and loading have different dt (they should be the same value though)
        call ReadLoadingAtTime(newPatch%loadObject,loadFileUnit,&
                        initialAzimuth, loadInfo(LOAD_TIMETYPE))
      end do
    end if

  end subroutine ReadPatchAtTime

  subroutine ReadPatchGeomAtTime(newPatch, surfaceFileUnit,geoInfo, &
                                 initialAzimuth, vFlag)
    type(Patch), intent(inout)::newPatch
    integer, intent(in)::surfaceFileUnit, vflag
    integer, intent(in), dimension(:)::geoInfo  
    real, intent(in)::initialAzimuth    
    integer:: i, SurfnKey
    
    ! V0_0 files are organized such that all of the data in time for each patch is 
    ! listed immediately after the respective patch dimension information, thus it 
    ! must be gathered all at once (where SurfnKey  and LoadnKey are the number of 
    ! timesteps).  V1_0 files are organized such that the data in time for
    ! each patch is listed together and after all of the patch dimension information, 
    ! thus the Read( )AtTime subroutines need only be called once.  (The outer loop 
    ! within container is used in place of SurfnKey and LoadnKey for V1_0).
    
    if (vFlag == VERSION_0) then
      SurfnKey = GetSurfacenKey(newPatch%surfaceObject)
    else
      SurfnKey = 1
    end if
    do i = 1,SurfnKey
      call ReadSurfaceAtTime(newPatch%surfaceObject,SurfaceFileUnit,&
                             initialAzimuth, geoInfo(GEO_TIMETYPE))
    end do
  end subroutine ReadPatchGeomAtTime

  subroutine ReadPatchLoadAtTime(newPatch, loadFileUnit, loadInfo, &
                                 initialAzimuth, vFlag)
    type(Patch), intent(inout)::newPatch
    integer, intent(in):: loadFileUnit, vFlag
    integer, intent(in), dimension(:)::loadInfo
    real, intent(in):: initialAzimuth    
    integer:: i, LoadnKey
    
    ! V0_0 files are organized such that all of the data in time for each patch is 
    ! listed immediately after the respective patch dimension information, thus it 
    ! must be gathered all at once (where SurfnKey  and LoadnKey are the number of 
    ! timesteps).  V1_0 files are organized such that the data in time for
    ! each patch is listed together and after all of the patch dimension information, 
    ! thus the Read( )AtTime subroutines need only be called once.  (The outer loop 
    ! within container is used in place of SurfnKey and LoadnKey for V1_0).
    
    if (vFlag == VERSION_0) then
      Loadnkey = GetLoadnKey(newPatch%loadObject)
    else
      LoadnKey = 1
    end if
    if (loadingNoiseFlag .or. totalNoiseFlag) then    
      do i = 1,LoadnKey
        call ReadLoadingAtTime(newPatch%loadObject,loadFileUnit,&
                        initialAzimuth, loadInfo(LOAD_TIMETYPE))
      end do
    end if
  end subroutine ReadPatchLoadAtTime                             
  
  function GetNKey(newPatch) result(nKey)
    type(patch):: newPatch
    integer:: nKey
    
    if (associated(newPatch%loadObject)) then
      nKey = GetLoadnKey(newPatch%loadObject)
    else
      nKey = GetSurfacenKey(newPatch%surfaceObject)
    end if
    
  end function GetNKey
  
   subroutine ResetPatchKeyCount(P,geokeyval,loadkeyval,forceGeo,forceLoad)
    type(Patch):: P
    integer, intent(in), optional:: geokeyval, loadkeyval
    logical, intent(in), optional:: forceGeo, forceLoad
    logical:: forceGeoReset, forceLoadReset
    
    forceGeoReset = .false.
    if( present(forceGeo) ) then
      if( forceGeo ) forceGeoReset=.true.
    end if
  
    if (getSurfaceType(P%surfaceObject) == APERIODIC_SURFACE .or. forceGeoReset) then
      if( present(geokeyval) ) then
        call ResetSurfaceKeyCount(P%surfaceObject,geokeyval)
      else
        call ResetSurfaceKeyCount(P%surfaceObject)
      end if
    end if

    forceLoadReset = .false.
    if( present(forceLoad) ) then
      if( forceLoad ) forceLoadReset=.true.
    end if
  
    if (associated(P%loadObject) ) then
      if (getLoadType(P%loadObject) == APERIODIC_LOADING .or. forceLoadReset) then
        if( present(loadkeyval) ) then
          call ResetLoadingKeyCount(P%loadObject,loadkeyval)
        else
          call ResetLoadingKeyCount(P%loadObject)
        end if
      end if
    end if
  
  end subroutine ResetPatchKeyCount
   
  subroutine SetPatchNKey(P,val, geoNKeyval,loadNKeyval)
    type(Patch):: P
    integer, intent(in):: val
    integer, intent(in), optional:: geoNKeyval, loadNKeyval
  
    P%nTau = val
    if (getSurfaceType(P%surfaceObject) == APERIODIC_SURFACE) then
      if( present(geoNKeyval) ) then
        call SetSurfaceNKey(P%surfaceObject,geoNKeyval)
      else
        call SetSurfaceNKey(P%surfaceObject,val)   
      end if
    end if
    if (associated(P%loadObject) ) then
      if (getLoadType(P%loadObject) == APERIODIC_LOADING) then
        if( present(loadNKeyval) ) then
          call SetLoadingNKey(P%loadObject,loadNKeyval)
        else
          call SetLoadingNKey(P%loadObject,val)
        end if
      end if
    end if
  
  end subroutine SetPatchNkey
  
  subroutine SetPatchNodeRange(P,firstNode,lastNode)
    type(Patch), intent(inout):: P
    integer, intent(in):: firstNode, lastNode
    P%FirstNode = firstNode
    P%LastNode  = lastNode
    call SetSurfaceNodeRange(P%surfaceObject,firstNode,lastNode)
  end subroutine SetPatchNodeRange
  
    
  subroutine DeterminePatchDtau(newpatch,tauMin,tauMax,ntau,dtau)
    type(patch),intent(inout)::newpatch
    integer::i,ntau
    real:: epsilon
    real::tempDTau,taumin,taumax,dtau,MinimumTau

    ! We need to throw this check in here because we bring in dtau from C%dtau
    ! but for a copy this information only exists in the patch.
    if (newpatch%isCopy) then
      dtau = newpatch%dtau
      return  !ksb debug
    end if
    ! need to make sure that the input parameters given by the container make
    ! sense    
    if (dtau.eq.0.0.and.(ntau.gt.1)) then
      if(tauMin.eq.tauMax) then
        call Error("tauMin and tauMax given in namelist are the same", &
                   "ntau was specified so assuming this is input error")
      else
        dtau = (tauMax-tauMin)/real(nTau)
      end if
    end if
    ! Get dTau. Defer to that in the surface and loading files if it exists,
    ! otherwise take the one that is input. Note that the surface and loading
    ! must have matching dTaus if they are non-constant.
    newPatch%dTau = GetSurfaceDt (newPatch%surfaceObject)
    ! Loading noise time calculations
    if (newPatch%computeLoadingFlag) then
      if (associated(newPatch%loadObject)) then
        tempDtau = GetLoadingDt(newPatch%loadObject)
      end if
      if (tempDtau /= 0.0 .and. newPatch%dTau /= 0.0 .and. &
          abs(tempDtau-newPatch%dTau) > tempDTau*epsilon(tempDtau) ) then  !epsilon doesn't work Intel Fortran 9.1 (Windows)
        call Error ("The timestep in the surface file does not match that in the",&
                    "loading file")
      else if (tempDtau == 0.0 .and. dtau /= 0.0) then ! take on the parent containers dtau
        newPatch%dTau = dtau
      else
        newPatch%dTau = tempDtau  !this seems redundant - first if should catch this - ksb
      end if
      if(GetNbSlaves().gt.0) then
        newPatch%iTauMin   = 1       
        if      (associated(newPatch%loadObject)) then
          minimumTau    = GetLoadingSourceTime(newPatch%loadObject,1)
          newPatch%nTau = GetLoadingNTau      (newPatch%loadObject)
        end if
        return
      end if
      !End of loading dTau calculation
    end if
    ! Thickness noise time calculations
    !If patch loading file has not previously set a dTau, we need to calculate
    !it here.
    if (newPatch%dTau == 0.0) then
      if (dTau==0.0.and.(tauMin==tauMax.or.ntau==0)) then
        call Error('A (dTau) or (tauMin, tauMax, ntau) must be set for this patch.')
      else
        if (dTau==0.0) then
          newPatch%nTau = nTau
          allocate(newPatch%tau(nTau))
          newPatch%dtau = (tauMax-tauMin)/ntau
          do i=1, newPatch%nTau
            newPatch%tau(i) = ((i-1) * newPatch%dTau) + tauMin
          end do
        else
          newPatch%dTau=dTau
        end if
      end if
    end if
!    ! This is a parallel temporary check and still needs to be addressed
!    if ((.not.associated(newPatch%tau)).and.(GetNbSlaves().gt.0)) then
!      call Error('Temporary Bug','Please specify a tauMin, tauMax and nTau')
!      stop
!    end if
    newPatch%iTauMin = 1
    newPatch%trailingEdgeFlag=.false. 
  end subroutine DeterminePatchDtau
 

  function DetermineAcousticSurfaceType(loadFileUnit) result (surfType)
    integer, intent(in) :: loadFileUnit
    integer :: surfType, status  
    character(len=32) :: charType

    if (loadFileUnit <= 0) then
      surfType = NO_SURFACE_DATA
      return
    end if
    
          
    call ReadBinaryString (loadFileUnit, charType, status)
    if (status /= 0) then
      if (status == 1) then    
        surfType = NO_SURFACE_DATA
        return
      else if (status == 2) then
        call Error ("Loading file ended prematurely.")
      else 
        call Error ("Failed while reading loading file header.")
      end if
    end if

    ! Parse the string:
    call strToUpper (charType, 32)
    if (index(charType, "NO") /= 0) then
      surfType = NO_SURFACE_DATA
    else if (index(charType, "FLOW") /= 0) then
      surfType = FLOW_SURFACE_DATA
    else if (index(charType, "LOAD") /= 0) then
      surfType = LOADING_SURFACE_DATA      
    else
      call Error ("Unrecognized data type: "//trim(charType))
    end if
          
    
  end function DetermineAcousticSurfaceType


  function IsStructured (P)
    type(patch), intent(in) :: P
    logical :: IsStructured
    IsStructured = (GetSurfaceGridType (P%surfaceObject) == STRUCTURED_GRID)
    
  end function IsStructured


  subroutine CopyPatch (newPatch, oldPatch, initialAzimuth)
    type(patch), intent(inout) :: newPatch
    type(patch), intent(in) :: oldPatch
    real, intent(in) :: initialAzimuth

    ! Store a pointer to the same data as the original patch.
    newPatch%geoTitle="Copy of "//trim(oldPatch%geoTitle)
    newPatch%loadTitle="Copy of "//trim(oldPatch%loadTitle)
    newPatch%prefix               =  oldPatch%prefix 
    newPatch%loadingOutput        =  oldPatch%loadingOutput
    newPatch%ntau                 =  oldPatch%ntau 
    newPatch%numberVariables      =  oldPatch%numberVariables
    newPatch%numberFunctions      =  oldPatch%numberFunctions
    newPatch%iTau                 =  oldPatch%iTau  
    newPatch%iTauMin              =  oldPatch%iTauMin
    newPatch%iMax                 =  oldPatch%iMax
    newPatch%jMax                 =  oldPatch%jMax      
    newPatch%prefixLength         =  oldPatch%prefixLength
    newPatch%refFrame             =  oldPatch%refFrame
    newPatch%nbNodes              =  oldPatch%nbNodes  
    newPatch%loadSpecies          =  oldPatch%loadSpecies
    newPatch%nSourceBase          =  oldPatch%nSourceBase
    newPatch%srcIndex             =  oldPatch%srcIndex          
    newPatch%nbNodes = oldPatch%nbNodes
    newPatch%isCopy = .true.
    newPatch%flowData             =  oldPatch%flowData
    newPatch%computeThicknessFlag =  oldPatch%computeThicknessFlag
    newPatch%computeLoadingFlag   =  oldPatch%computeLoadingFlag    
    newPatch%trailingEdgeFlag     =  oldPatch%trailingEdgeFlag
    newPatch%computePressGradFlag =  oldPatch%computePressGradFlag
    newPatch%surfaceAperiodic     =  oldPatch%surfaceAperiodic
    newPatch%loadAperiodic        =  oldPatch%loadAperiodic    
    newPatch%dtau                 =  oldPatch%dtau  
    newPatch%tauMin               =  oldPatch%tauMin
    newPatch%tauMax               =  oldPatch%tauMax
    newPatch%tauLocal             =  oldPatch%tauLocal
    if(associated(oldPatch%tau)) then
      allocate (newPatch%tau(size(oldPatch%tau)))
      newPatch%tau=oldPatch%tau
    else
      nullify(newPatch%tau)
    end if
    if(associated(oldPatch%ObsTimeHistory)) then
      allocate (newPatch%ObsTimeHistory(size(oldPatch%ObsTimeHistory,1),size(oldPatch%ObsTimeHistory,2)))
      newPatch%ObsTimeHistory=oldPatch%ObsTimeHistory
    else
      nullify(newPatch%ObsTimeHistory)
    end if
    if(associated(oldPatch%qTdS)) then
      allocate (newPatch%qTdS(size(oldPatch%qTdS,1),size(oldPatch%qTdS,2)))
      newPatch%qTdS=oldPatch%qTdS
    else
      nullify(newPatch%qTdS)
    end if
    if(associated(oldPatch%qLdS)) then
      allocate (newPatch%qLdS(size(oldPatch%qLdS,1),size(oldPatch%qLdS,2)))
      newPatch%qLdS=oldPatch%qLdS
    else
      nullify(newPatch%qLdS)
    end if   
    if(associated(oldPatch%sigmaSurface)) then
      allocate (newPatch%sigmaSurface(size(oldPatch%sigmaSurface,1), &
        size(oldPatch%sigmaSurface,2),size(oldPatch%sigmaSurface,2)))
      newPatch%sigmaSurface=oldPatch%sigmaSurface
    else
      nullify(newPatch%sigmaSurface)
    end if    
    newPatch%srcVector   =  oldPatch%srcVector
    ! Copy the surface. This doesn't copy all of the memory, just stores the
    ! azimuthal offset and moves some pointers around.
    allocate (newPatch%surfaceObject)
    call CopySurface (newPatch%surfaceObject, oldPatch%surfaceObject, initialAzimuth)
    ! Not sure this is really used but copy anyway for now (ksb 7/31/2014)
    !!allocate (newPatch%srcSurface)
    !!call CopySurface (newPatch%srcSurface, oldPatch%srcSurface, initialAzimuth)    
    ! Copy the loading (which doesn't actually copy the data, just stores the
    ! offset, so there's no memory hit)
    nullify (newPatch%loadObject)
    if (associated (oldPatch%loadObject)) then
      allocate (newPatch%loadObject)
      call CopyLoading (newPatch%loadObject, oldPatch%loadObject, initialAzimuth)
    end if    
    if(associated(oldPatch%StoredDataArray)) then
      allocate (newPatch%StoredDataArray(size(oldPatch%StoredDataArray)))
      newPatch%StoredDataArray=oldPatch%StoredDataArray
    else
      nullify(newPatch%StoredDataArray)
    end if
    if(associated(oldPatch%srcCBList)) then
      allocate (newPatch%srcCBList)
      newPatch%srcCBList=oldPatch%srcCBList
    else
      nullify(newPatch%srcCBList)
    end if
    if(associated(oldPatch%obsCBList)) then
      allocate (newPatch%obsCBList)
      newPatch%obsCBList=oldPatch%obsCBList
    else
      nullify(newPatch%obsCBList)
    end if
    if(associated(oldPatch%accelerationsaveLL)) then
      allocate (newPatch%accelerationsaveLL(size(oldPatch%accelerationsaveLL)))
      newPatch%accelerationsaveLL=oldPatch%accelerationsaveLL
    else
      nullify(newPatch%accelerationsaveLL)
    end if
    if(associated(oldPatch%accelerationsaveL)) then
      allocate (newPatch%accelerationsaveL(size(oldPatch%accelerationsaveL)))
      newPatch%accelerationsaveL=oldPatch%accelerationsaveL
    else
      nullify(newPatch%accelerationsaveL)
    end if
    if(associated(oldPatch%accelerationsaveC)) then
      allocate (newPatch%accelerationsaveC(size(oldPatch%accelerationsaveC)))
      newPatch%accelerationsaveC=oldPatch%accelerationsaveC
    else
      nullify(newPatch%accelerationsaveC)
    end if
    if(associated(oldPatch%nhatdotsaveLL)) then
      allocate (newPatch%nhatdotsaveLL(size(oldPatch%nhatdotsaveLL)))
      newPatch%nhatdotsaveLL=oldPatch%nhatdotsaveLL
    else
      nullify(newPatch%nhatdotsaveLL)
    end if
    if(associated(oldPatch%nhatdotsaveL)) then
      allocate (newPatch%nhatdotsaveL(size(oldPatch%nhatdotsaveL)))
      newPatch%nhatdotsaveL=oldPatch%nhatdotsaveL
    else
      nullify(newPatch%nhatdotsaveL)
    end if
    if(associated(oldPatch%nhatdotsaveC)) then
      allocate (newPatch%nhatdotsaveC(size(oldPatch%nhatdotsaveC)))
      newPatch%nhatdotsaveC=oldPatch%nhatdotsaveC
    else
      nullify(newPatch%nhatdotsaveC)
    end if
    if(associated(oldPatch%ObsTimeInRange)) then
      allocate (newPatch%ObsTimeInRange(size(oldPatch%ObsTimeInRange)))
      newPatch%ObsTimeInRange=oldPatch%ObsTimeInRange
    else
      nullify(newPatch%ObsTimeInRange)
    end if

  end subroutine CopyPatch


  !*************************************************************************
  !SUBROUTINE InitializeSigma(P)
  !  This routine will create a sigma array for the patch depending on which
  ! sigma flags are true. The sigma flags are defined in the global module 
  ! constants.
  !ARGUEMENTS:
  ! - P: The patch that we are creating the sigmaSurface for.
  !*************************************************************************
  subroutine InitializeSigma(P)
    implicit none
    type(patch), intent(inout)::P

    ! A patch has 5 variables (X, Y, Z, t, tau) and can have a number of
    ! functions associated with it.  Vectors are recorded as having 3 functions
    ! numbers and are defined as vectors in the sigma surface namelist.
    if( GetSurfaceIblanking(P%surfaceObject) ) then
      P%numberVariables = 6
    else
      P%numberVariables = 5
    end if
    P%numberFunctions = 0
    if(loadingNoiseSigmaFlag)   P%numberFunctions = P%numberFunctions + 1
    if(thicknessNoiseSigmaFlag) P%numberFunctions = P%numberFunctions + 1
    if(totalNoiseSigmaFlag)     P%numberFunctions = P%numberFunctions + 1
    if(normalSigmaFlag)         P%numberFunctions = P%numberFunctions + 3
    if(machSigmaFlag)           P%numberFunctions = P%numberFunctions + 1
    if(velocitySigmaFlag)       P%numberFunctions = P%numberFunctions + 3
    if(accelerationSigmaFlag)   P%numberFunctions = P%numberFunctions + 3
    if(loadingSigmaFlag)        P%numberFunctions = P%numberFunctions + 3
    if(densitySigmaFlag)        P%numberFunctions = P%numberFunctions + 1
    if(momentumSigmaFlag)       P%numberFunctions = P%numberFunctions + 3
    if(pressureSigmaFlag)       P%numberFunctions = P%numberFunctions + 1
    if(areaSigmaFlag)           P%numberFunctions = P%numberFunctions + 1
    if(MdotrSigmaFlag)          P%numberFunctions = P%numberFunctions + 1
    if(iblankSigmaFlag)         P%numberFunctions = P%numberFunctions + 1

  end subroutine InitializeSigma


  !**************************************************************************
  !SUBROUTINE IntegratePatch(P, obsVector, obsCBList, CBlist, &
  !                      thicknessNoise, loadingNoise, iBlankArray, &
  !                      time, PTGradient, PL1Gradient,PL2Gradient, &
  !                      keyIndex, nTau, tauMin, tauMax, initialAzimuth)
  !  The subroutine computes the thickness and loading noise for a patch for
  ! an observer position obsInitVector, and returns the result in thicknessNoise and loadingNoise. 
  !ARGUMENTS:
  ! - P: The patch we are integrating over.
  ! - obsInitVector: The observer initial vector, if the observer moves then
  !                  this value will change with time, tau.
  ! - obsCBList: The observer CB list, the actual location of the observer will
  !              we computing within the time loop.
  ! - CBList: This is the CBList that is defines the position of the observer 
  !           with respect to the initial frame.
  ! - thicknessNoise, loadingNoise: array of real, dimension nt, arrays to store the 
  !                          overall thickness and loading noise  
  !REMARKS:
  !  - This is the main subroutine of the code.
  !**************************************************************************
  subroutine IntegratePatch(P, BPMData, obsInitVector, obsCBList, CBList, &
                            iBlankArray, thicknessNoise, loadingNoise, PTGradient, &
                            PL1Gradient, observerTimeArray, &   
                            dt, iTau, tauMax, keyIndexSurf, keyIndexLoad, &
                            patchNum, nbZones, obsNT, obsCount, lastGeoZone, lastBlade, radDistance, segmentSize)
    type(patch), intent(inout):: P
    type(BPM), pointer:: BPMData    
    type(Vector), intent(in):: obsInitVector
    type(CBStructure), pointer::obsCBList, CBList
    integer, dimension(:)::iBlankArray
    real(kind=8), dimension(:), intent(inout)::thicknessNoise, loadingNoise    
    type(vectorKind8), dimension(:), intent(inout)::PTGradient,PL1Gradient
    real(kind=8), dimension(:), intent(in)::observerTimeArray
    real(kind=8), intent(inout):: dt
    real:: dtKind4
    real:: segmentSize
    real, dimension(:), pointer:: radDistance
    integer, intent(in)::iTau, patchNum, nbZones, keyIndexSurf, keyIndexLoad, obsNT, obsCount
    integer, intent(inout):: tauMax
    integer::i, j, iTau1, iTau2
    logical:: Initializing, inRange, dataLimit, lastGeoZone, lastBlade, storePeggData, storeBPMData

    if (iTau.eq.1) P%iTau = 0
    ! We only need to worry about the data being in range of the observer time
    ! if it is aperiodic.  Otherwise there should be no problem.
    inRange = .false.
    if (.not. (P%surfaceAperiodic .or. P%loadAperiodic) ) then
      inRange = .true.
    else
      if (P%ObsTimeInRange(obsCount)%tMin(iTau).ne.huge(dt)) then
        if (iTau.eq.1) then
          iTau1 = iTau
          itau2 = iTau+1
        else
          iTau1 = iTau - 1
          iTau2 = iTau
        end if
        inRange = sourceTimeInRange(P%ObsTimeInRange(obsCount)%tMin(iTau-1), &
                                    P%ObsTimeInRange(obsCount)%tMax(iTau),   &
                                    observerTimeArray(1), obsNT, dt)
      else
        inRange = .true.
      end if
    end if
  
    if (inRange) then
      P%iTau = P%iTau + 1
      mg_radIndex = 1

      ! Usually this is the same as iTau.eq.1, but it won't necessarily be with aperiodic data
      ! if we're not starting the computations until the data is in range of the observer time array.
      if (P%iTau.eq.1) then
        ! Reset the static member variables
        interpolatePosition = 1
        if (patchNum.eq. 1) then
          call ResetGlobals
        end if
        mg_obsInitialVector = obsInitVector
        call initializeInterpolateIndex ()
        
        !If we have a moving observer we need to initalize some things
        if(associated(obsCBList)) then 
          !We set the radiation vector just to be the first position of the observer
          mg_radiationVector=mg_obsInitialVector
          P%storedDataArray%initialr=vectorAbsolute(mg_obsInitialVector)
          !Because we need the observer CBlist to be a global parameter in patch we need to copy it over
          mg_obsCBList=>obsCBList
        end if     
      
        ! We need to add this check in case of an aperiodic case which will not have
        ! the observerTimeArray built at this point.  In the case of an aperiodic case
        ! we will already have dt.
        if (observerTimeArray(1) /= observerTimeArray(size(observerTimeArray))) then
          dt = observerTimeArray(2)-observerTimeArray(1)
        end if
        ! Allocate the memory for the noise:
        call CreateStoredData(P,CBList)

        !If we have a moving observer we need to initalize some things
        if(associated(obsCBList)) then 
          !We set the radiation vector just to be the first position of the observer
          mg_radiationVector=mg_obsInitialVector
          P%storedDataArray(0)%initialr=vectorAbsolute(mg_obsInitialVector)
          !Because we need the observer CBlist to be a global parameter in patch we need to copy it over
          mg_obsCBList=>obsCBList
        end if
      
        if(P%computeThicknessFlag.or.P%computeLoadingFlag.or.&
           isomsThicknessNoiseFlag.or.sigmaFlag) then
          
          if(debugLevel.gt.3 .and. IsMaster()) then
            write (*,'(A,$)') P%prefix(1:P%prefixLength)//&
              ' Integrating '//trim(P%geoTitle)//&
              ' with nTau = '     
            write (*,'(I6)') P%nTau
            write (*,*) P%prefix(1:P%prefixLength)//&
              "  Source time range = ", P%tau(1), " --> ", P%tau(P%ntau)
          end if
        end if
      end if    
    
    ! This replaces the current (and likely incorrect) module globals with the correct information for
    ! current patch.
      call patch2MG(P, P%storedDataArray(mg_CurrentWallNb), iTau)
      mg_CurrentWallNb = 0
   
      initializing = .FALSE.
      mg_patchNum = patchNum
      if(P%computeThicknessFlag.or.P%computeLoadingFlag.or.&
        &  isomsThicknessNoiseFlag.or.sigmaFlag) then
        !call CalculateDifferentialNoise(P, CBList, obsCBList, &
        !                                thicknessNoise, loadingNoise, iBlankArray, &
        !                                PTGradient,PL1Gradient, &
        !                                observerTimeArray, dt, iTau, tauMax, &
        !                                keyIndexSurf, keyIndexLoad, &
        !                                nbZones, .FALSE., .FALSE.)
        call CalculateDifferentialNoise (P, BPMData, CBList, obsCBList, iBlankArray, thicknessNoise, &
        &   loadingNoise, PTGradient,PL1Gradient, observerTimeArray, dt, iTau, tauMax, keyIndexSurf, &
        &   keyIndexLoad, obsCount, lastGeoZone, lastBlade, nbZones, .FALSE., .FALSE., radDistance, segmentSize) 
        if (mg_NbWall.gt.0) then
          call IntegrateOverWall(P, BPMData, CBList, obsCBList, dt, thicknessNoise, loadingNoise, &
            PTGradient, PL1Gradient, iBlankArray, observerTimeArray, iTau, tauMax, keyIndexSurf,  &
            keyIndexLoad, nbZones, obsCount, lastGeoZone, lastBlade, radDistance,segmentSize)
        end if                                    
      end if
      dtKind4 = dt
      if (iTau /= tauMax) then
        do i=0,mg_NbWall
          call mg2Patch(P, P%storedDataArray(i), Initializing)
        end do
      else
        deallocate(P%tau)      
      end if
    end if
  end subroutine IntegratePatch

  subroutine CreateStoredData(P, CBList)
    type(patch):: P
    type(CBStructure), pointer:: CBList
    type(storedData), pointer:: Z
    integer::nbNodes
    integer, dimension(2):: surfDims
    real(kind=8), parameter:: zero=0.

    
    P%nSourceBase = GetSize(CBList)
    nbNodes = GetSurfaceNbNodes(P%surfaceObject)
    P%nbNodes = nbNodes
    call GetSurfaceDimensions (P%surfaceObject, surfDims)
    if( GetSurfaceGridType(P%surfaceObject) == STRUCTURED_GRID) then
      P%iMax = surfDims(1)
      P%jMax = surfDims(2)
    end if
    
    if (associated(P%loadObject)) then
      P%loadSpecies = GetLoadSpecies(P%loadObject)
      if (P%loadSpecies == LOAD_DATATYPE_FLOW) then
        P%flowData = .TRUE.
      end if
    endif
 

    if(PressureGradient1AFlag) then
      ! allocate and initialize the acceleration and nhatdot arrays for the patch
      allocate(P%accelerationsaveLL(nbNodes), P%accelerationsaveL(nbNodes), &
      &        P%accelerationsaveC(nbNodes))
      P%accelerationsaveLL = vectorSetCoordinates(0.,0.,0.)
      P%accelerationsaveL  = vectorSetCoordinates(0.,0.,0.)
      P%accelerationsaveC  = vectorSetCoordinates(0.,0.,0.)
      allocate(P%nhatdotsaveLL(nbNodes), P%nhatdotsaveL(nbNodes), &
      &        P%nhatdotsaveC(nbNodes))  
      P%nhatdotsaveLL = vectorSetCoordinates(0.,0.,0.)
      P%nhatdotsaveL  = vectorSetCoordinates(0.,0.,0.)
      P%nhatdotsaveC  = vectorSetCoordinates(0.,0.,0.)
    end if
    
    do mg_CurrentWallNb=0,mg_NbWall
      Z=>P%storedDataArray(mg_CurrentWallNb)

      ! Allocate the time history data:
      allocate (Z%obsTimeHistoryL(P%nbNodes), Z%obsTimeHistoryC(P%nbNodes))
      Z%obsTimeHistoryL = zero
      Z%obsTimeHistoryC = zero

      ! If needed, allocate the flow parameter arrays.
      if (P%loadSpecies == LOAD_DATATYPE_FLOW) then
        allocate (Z%parameterC, Z%parameterR)      
        call CreatePermeableParameters(P, Z%parameterC, nbNodes)
        call CreatePermeableParameters(P, Z%parameterR, nbNodes)
      else if (P%loadSpecies == LOAD_DATATYPE_PRESSURE) then
        allocate (Z%parameterC, Z%parameterR)      
        call CreatePressureParameters(P, Z%parameterC, nbNodes)
        call CreatePressureParameters(P, Z%parameterR, nbNodes)     
      end if

    
      if(P%computeThicknessFlag) then
        ! allocate memory for thickness noise
        allocate(Z%qTdSL(nbNodes), Z%qTdSC(nbNodes)) 
        Z%qTdSL = zero
        Z%qTdSC = zero
        if(PressureGradient1AFlag) then
          ! Allocate the memory for pressure gradient
          allocate(Z%dPdTL(nbNodes), Z%dPdTC(nbNodes))
          Z%dPdTL  = vectorSetKind8Coordinates(zero,zero,zero)
          Z%dPdTC  = vectorSetKind8Coordinates(zero,zero,zero)
        end if

      end if
      if(P%computeLoadingFlag.or.IsomsThicknessNoiseFlag) then
        ! allocate memory for loading noise (or Isom's thickness noise)
        allocate(Z%qLdSL(nbNodes), Z%qLdSC(nbNodes))
        Z%qLdSL = zero
        Z%qLdSC = zero
        if(PressureGradient1AFlag) then
          ! Allocate the memory for pressure gradient
          allocate(Z%dPdL1L(nbNodes), Z%dPdL1C(nbNodes))
          Z%dPdL1L  = vectorSetKind8Coordinates(zero,zero,zero)
          Z%dPdL1C  = vectorSetKind8Coordinates(zero,zero,zero)
       end if
      end if
      if(sigmaFlag) then
        call InitializeSigma (P)
        allocate(P%sigmaSurface(P%numberFunctions+P%numberVariables,&
                                nbNodes,P%nTau))
      end if 
      
    end do
    mg_CurrentWallNb=0   !set this back to zero.
  end subroutine CreateStoredData
  

  subroutine CalculateDifferentialNoise (P, BPMData, CBList, obsCBList, iBlankArray, thicknessNoise, &
  &   loadingNoise, PTGradient,PL1Gradient, observerTimeArray, dt, iTau, tauMax, keyIndexSurf,&
  &   keyIndexLoad, obsCount, lastGeoZone, lastBlade, nbZones, isWall, lastWall, radDistance, segmentSize)
    type(patch)::P
    type(StoredData), pointer:: Z => NULL()
    type(CBStructure), pointer::CBList, obsCBList
    integer, dimension(:)::iBlankArray
    type(vectorKind8), dimension(:), intent(inout)::PTGradient,PL1Gradient
    type(BPM), pointer:: BPMData
    real(kind=8), dimension(:), intent(in)::observerTimeArray
    real(kind=8), dimension(:), intent(inout)::thicknessNoise, loadingNoise    
    integer, intent(in) :: iTau, keyIndexSurf, keyIndexLoad, nbZones, obsCount
    integer, intent(inout) :: tauMax
    real(kind=8)::dt
    logical, intent(in):: isWall, lastWall, lastGeoZone, lastBlade
    
    real::sourceTime, sourceTimeL, obsTimeL
    integer::i, j, BBkey, index
    type(vector)::velocity, acceleration, nhat
    type(vector)::nhatdot, velocityL, accelerationL
    type(vector)::nhatL, nhatdotL, velocityLL, accelerationLL
    type(vector)::nhatLL, nhatdotLL
    type(vector):: srcVectorL, srcVectorLL
    real :: surfAreaL, surfAreaLL, dtKind4, observerTimeKind4
    real:: delt1, delt2, alpha, segmentSize
    type(SingleSurface), pointer:: currentSurface =>NULL(), surfV =>NULL(), surfA =>NULL()
    type(vector) :: accelerationdot,nhat2dot    
    type(vector), dimension(:), pointer :: VectorTemp =>NULL()
    real        , dimension(:), pointer :: RealTemp =>NULL(), radDistance
    logical:: unknownSourceTimeFlag = .FALSE., doOpen = .TRUE.   
    real(kind=8),parameter:: zero=0.
    
    Z => P%StoredDataArray(mg_CurrentWallnb)  !ksb debug:  check to see if I really need Z here (or just use P)
    sourceTime = P%tau(iTau)
    call CalculateResultantChangeOfBase(CBList, P%nSourceBase, sourceTime)
    dtKind4 = dt  
    if (associated(P%loadObject)) then
      if (iTau == 1) then
        ! Calculate the value at the first time and store it twice to take
        ! care of the initialization. Not that efficient, but the speed hit is
        ! inconsequential.
        if (P%loadSpecies == LOAD_DATATYPE_FLOW) then 
          call CalculatePermeableParameters (P, CBList, P%nSourceBase, mg_flowParametersR, 1, keyIndexSurf, keyIndexLoad)
          call CalculatePermeableParameters (P, CBList, P%nSourceBase, mg_flowParametersC, 1, keyIndexSurf, keyIndexLoad)
        else if (iTau == 1 .and. P%loadSpecies == LOAD_DATATYPE_PRESSURE) then
          call CalculatePressureParameters (P, mg_flowParametersR, 1, keyIndexSurf, keyIndexLoad)
          call CalculatePressureParameters (P, mg_flowParametersC, 1, keyIndexSurf, keyIndexLoad)
        end if
      end if
     
      ! Rotate the data parameters (i.e., circular shift of module global variables: LCR -> CRL ), then
      ! overwrite the right value (R) with the new data
      if (P%loadSpecies == LOAD_DATATYPE_FLOW .or. P%loadSpecies == LOAD_DATATYPE_PRESSURE) then
        if (iTau < P%nTau) then  ! don't shift the iTau===P%nTau - just use it again.
          call Rotate_mg_flowParameters
          if (P%flowData) then
            if (P%iTau == 1) then
              allocate (mg_flowParametersR)           
              call CreatePermeableParameters(P, mg_flowParametersR, P%nbNodes)
            end if
            call CalculatePermeableParameters (P, CBList, P%nSourceBase, & 
                                                 mg_flowParametersR, P%iTau+1, keyIndexSurf+1, keyIndexLoad+1)
          else  ! pressure data
            if (P%iTau.eq.1) then
              allocate (mg_flowParametersR)
              call CreatePressureParameters(P, mg_flowParametersR, P%nbNodes)
            end if
            call CalculatePressureParameters (P, mg_flowParametersR, P%iTau+1, keyIndexSurf+1, keyIndexLoad+1)
          end if
        end if
      end if          
    end if

    if (globalBPMNoiseFlag .and. isSurfaceCompact(P%surfaceObject)) then    
      BBkey = 1
      if (BPMData%timeType.eq.LOAD_TIMETYPE_PERIODIC) BBkey = getBBoffset(size(BPMData%tau), &
          BPMData%timeOffset, BPMData%period, keyIndexLoad+P%iTauMin-1)
    end if

    ! Do the integration loop
    currentSurface => GetSurface (P%surfaceObject, keyIndexSurf+P%iTauMin-1)
    call GetSurfaceVelAccel (P%surfaceObject, keyIndexSurf+P%iTauMin-1, surfV, surfA, keyIndexSurf)

    do i=1, P%nbNodes
      if( currentSurface%iblanks(i) /=1 .and. .not. sigmaFlag) cycle !don't do all this if iblank = 0 
      ! Calculate the velocity and acceleration of this point at emission time
      call CalculateVelAccel (currentSurface%vectors(i), velocity, acceleration)
      
      ! Tack on the velocity and acceleration due to the flexing surface:
      ! First term: due to the COB specified in namelist, 
      ! Second term: calculated by motion of point read from input data in velocity and acceleration
      velocity     = velocity     + VectorBaseChange(surfV%vectors(i))
      acceleration = acceleration + VectorBaseChange(surfA%vectors(i))
      nhat         = VectorBaseChange(currentSurface%normals(i))
      nhatdot      = VectorDerivative(currentSurface%normals(i))+&
                     VectorBaseChange(surfV%normals(i))
     
      ! Calculate the derivative of acceleration and second derivative of nhat for pressure gradient 1A
      if(PressureGradient1AFlag) then 
        P%accelerationsaveC(i)=acceleration  !ksb debug: Should I save accelerationdot like this and nhat2dot?
        P%nhatdotsaveC(i)=nhatdot
            
        if (P%iTau == 1) then
          accelerationdot%A=0.0
          nhat2dot%A=0.0
        else if(P%iTau==2) then
        else if(iTau==2) then
          accelerationdot = (P%accelerationsaveC(i)-P%accelerationsaveL(i))/P%dtau
          nhat2dot        = (P%nhatdotsaveC(i)     -P%nhatdotsaveL(i))     /P%dtau
        else ! (iTau>=3)
          accelerationdot = (P%accelerationsaveLL(i) - 4.0*P%accelerationsaveL(i) + &
                             3.0*P%accelerationsaveC(i))/(2.0*P%dtau)
          nhat2dot        = (P%nhatdotsaveLL(i) - 4.0*P%nhatdotsaveL(i) + &
                             3.0*P%nhatdotsaveC(i))/(2.0*P%dtau)  
        end if
      end if
      
      ! Calculate the observer time
      call ObserverTime (currentSurface%vectors(i), sourceTime, observerTimeKind4, &
                        dtKind4, P%dTau, Z%initialr, obsCBList)  
      Z%obsTimeHistoryC(i) = observerTimeKind4

      if (AtmAtten) then
        if (.not.isWall .and. mg_radIndex.le.size(radDistance)) then
          ! Grab the radiation distance at the middle of the observer time segment.
          if (observerTimeKind4.gt.(observerTimeArray(1)+(mg_radIndex-1)*segmentSize+segmentSize/2.0) .and. &
                    radDistance(mg_radIndex).eq.0.0) then
            radDistance(mg_radIndex) = vectorAbsolute(mg_radiationVector)
            mg_radIndex = mg_radIndex + 1
          end if
        end if
      end if

      if (debugLevel >= 14) then
        write(*,*) "Time data: ", i, iTau, dtKind4, P%dTau, Z%initialr, &
                   associated(obsCBList), sourceTime,  Z%obsTimeHistoryC(i)
      end if
      ! Do the actual integration calculation
      call NodeIntegration (P, i, P%iTau, Z%qTdSC, Z%qLdSC, velocity, acceleration,&
                            nhat, nhatdot, accelerationdot, nhat2dot, surfV%areas(i), &
                            Z%dPdTC, Z%dPdL1C, CBList, P%nSourceBase, &
                            sourceTime, keyIndexSurf, keyIndexLoad)

      ! BPM's method includes contributions from individual segments, thus terms need to be calculated
      ! within the node loop
      !if (associated(BPMData) .and. isSurfaceCompact(P%surfaceObject)) then
      if( globalBPMNoiseFlag .and. isSurfaceCompact(P%surfaceObject) ) then
        if (BPMData%iTau.le.size(BPMData%obsTimeArray)) then
          call ComputeBPMBroadbandTerms(BPMdata, currentSurface%normals(i), mg_radiationVector, velocity, obsCBList, i, &
                getBPMtermKey(BPMData, BPM_SECTAOA, BBkey), getBPMtermKey(BPMData, BPM_U, BBkey), sourceTime)
        end if
        ! Put everything back.
        if (associated(ObsCBList)) call CalculateResultantChangeOfBase(CBList, P%nSourceBase, sourceTime)        
      end if

      ! If we're supposed to store the sigma info, do that now
      if(sigmaFlag) then
        if( currentSurface%iblanks(i) /=1 ) then
          if( P%computeThicknessFlag ) Z%qTdSC(i) = 0.0
          if( P%computeLoadingFlag )   Z%qLdSC(i) = 0.0
          if( PressureGradient1AFlag ) then
              if( P%computeThicknessFlag ) Z%dPdTC(i) = vectorSetKind8Coordinates(zero,zero,zero)
              if( P%computeLoadingFlag )   Z%dPdL1C(i)= vectorSetKind8Coordinates(zero,zero,zero)
          end if
        end if
        call RecordSigma (P, i, iTau, Z%qLdSC, Z%qTdSC, velocity, &
                          acceleration, Z%obsTimeHistoryC(i), mg_flowParametersC, &
                          keyIndexSurf, keyIndexLoad)
      end if      
    end do
    
    !ksb debug:  this is where Ben does average for BB noise (over a period - it actually it does the whole window size now. )
    !ksb debug: if (associated(BPMData) .and. isSurfaceCompact(P%surfaceObject) .and. .not.isWall) then
    if (globalBPMNoiseFlag .and. isSurfaceCompact(P%surfaceObject) .and. .not.isWall) then
      if (BPMData%iTau.le.size(BPMData%obsTimeArray)) then
        if (minval(Z%obsTimeHistoryC).ge. BPMData%obsTimeArray(BPMData%iTau) .and. &
            maxval(Z%obsTimeHistoryC).lt.(BPMData%obsTimeArray(BPMData%iTau)+segmentSize)) then
          do i=1,P%nbNodes
            call ComputeBPMintensity(BPMData, BPMData%Intensity(:,:,BPMData%iTau), BBkey, i)
          end do
          BPMdata%nbContributions(BPMData%iTau) = BPMdata%nbContributions(BPMData%iTau) + 1
        elseif (maxval(Z%obsTimeHistoryC).ge.(BPMData%obsTimeArray(BPMData%iTau)+segmentSize)) then
          BPMData%iTau = BPMData%iTau + 1
        end if
      end if
    end if  
    
    ! Then, check the observerTime array to make sure it is monotonically increasing
    if (P%iTau > 1) then
      do i=1,P%nbNodes
        if( Z%obsTimeHistoryC(i) < Z%obsTimeHistoryL(i) ) then
          call Warning('The computed observerTimeHistory is not monotonically increasing.')
          print*,'i=',i, ' in CalculateDifferentialNoise, patch.f90'
          exit
        end if
      end do

      if ( P%surfaceAperiodic .or. P%loadAperiodic ) then
        if (iTau.eq.2) then
          P%ObsTimeInRange(obsCount)%tMin(1) = min(P%ObsTimeInRange(obsCount)%tMin(1),minval(Z%obsTimeHistoryL))
          P%ObsTimeInRange(obsCount)%tMax(1) = max(P%ObsTimeInRange(obsCount)%tMax(1),maxval(Z%obsTimeHistoryL))
        end if
        P%ObsTimeInRange(obsCount)%tMin(iTau) = min(P%ObsTimeInRange(obsCount)%tMin(iTau),minval(Z%obsTimeHistoryC))
        P%ObsTimeInRange(obsCount)%tMax(iTau) = max(P%ObsTimeInRange(obsCount)%tMax(iTau),maxval(Z%obsTimeHistoryC))
      end if
 
      if ( SourceTimeInRange(minval(Z%obsTimeHistoryL), maxval(Z%obsTimeHistoryC), &
           observerTimeArray(1), size(observerTimeArray), dt)) then

        if (integrationType==TIME_DOMAIN) then
          if(P%computeThicknessFlag) then
          
            call precisetinterpolate(Z%obsTimeHistoryL, Z%obsTimeHistoryC, &
                                     Z%qTdSL, Z%qTdSC, observerTimeArray, thicknessNoise, dt, iBlankArray)
          end if
          if(P%computeLoadingFlag.or.isomsThicknessNoiseFlag) then
            select case (P%LoadingOutput) !ksb debug
            case (THICKNESS_OUT)
              call precisetinterpolate(Z%obsTimeHistoryL, Z%obsTimeHistoryC, &
                                       Z%qLdSL, Z%qLdSC, observerTimeArray, thicknessNoise, dt, iBlankArray)
            case default ! P%LoadingOutout==LOADING_OUT
              call precisetinterpolate(Z%obsTimeHistoryL, Z%obsTimeHistoryC, &
                                       Z%qLdSL, Z%qLdSC, observerTimeArray, loadingNoise, dt, iBlankArray)
            end select
          end if
        end if
        ! Do the interpolation into observer time for pressure gradient
        if(PressureGradient1AFlag) then
          if(mg_CurrentWallNb == mg_NbWall) then
            if(P%computeThicknessFlag) then
              do i=1,3
                if (P%iTau.eq.2) then
                  Z%dPdTL(:)%A(i) = Z%dPdTC(:)%A(i)
                end if
                call precisetinterpolate(Z%obsTimeHistoryL, Z%obsTimeHistoryC, Z%dPdTL(:)%A(i), & 
                  Z%dPdTC(:)%A(i), observerTimeArray, PTGradient(:)%A(i), dt, iBlankArray) 
              end do              
            end if
            if(P%computeLoadingFlag) then
              do i=1,3   
                if (P%iTau.eq.2) then
                  Z%dPdL1L(:)%A(i) = Z%dPdL1C(:)%A(i)
                end if
                select case (P%LoadingOutput) !ksb debug
                case (THICKNESS_OUT)
                  call precisetinterpolate(Z%obsTimeHistoryL, Z%obsTimeHistoryC, Z%dPdL1L(:)%A(i), &
                    Z%dPdL1C(:)%A(i), observerTimeArray, PTGradient(:)%A(i), dt, iBlankArray)
                case default ! P%LoadingOutout==LOADING_OUT
                  call precisetinterpolate(Z%obsTimeHistoryL, Z%obsTimeHistoryC, Z%dPdL1L(:)%A(i), &
                    Z%dPdL1C(:)%A(i), observerTimeArray, PL1Gradient(:)%A(i), dt, iBlankArray)
                end select
              end do    
            end if
          end if          
        end if
      else if (minval(Z%obsTimeHistoryL).gt.observerTimeArray(size(observerTimeArray))) then
        ! This check saves time for aperiodic data by exiting this routine once the data
        ! has passed the observer time range.
        ! However, we need to be careful not to stop the calculations until we've also
        ! finished computing the contribution from any walls being use
        if( (mg_CurrentWallNb == mg_nbWall) .and. (P%surfaceAperiodic .or. P%loadAperiodic) ) then
          tauMax = iTau
        end if
      end if
    end if

!!!!!!!!!!!!!!!!   
    !ksb debug: if (iTau /= tauMax) then
    ! I don't think there is any harm in rotating the data even if this is the last time step
    ! in the file.  It might help make it all work with multiple time file (MTF) aperiodic.  
    ! I think all this does is dispense with data rotating if it is the last time step. KSB 2/14/2015
    if( .true. ) then
      !Rotate the necessary data from previous time steps
      call RotateStoredReal2DKind8(Z%obsTimeHistoryL, Z%obsTimeHistoryC)
      
      if(P%computeThicknessFlag) then
        call RotateStoredReal2DKind8(Z%qTdSL,  Z%qTdSC)
      end if
      if(P%computeLoadingFlag .or. isomsThicknessNoiseFlag) then
        call RotateStoredReal2DKind8(Z%qLdSL,  Z%qLdSC)
      end if

      if(PressureGradient1AFlag) then
        if(P%computeThicknessFlag) call RotateStoredVectorKind8(Z%dPdTL,  Z%dPdTC)
        if(P%computeLoadingFlag)   call RotateStoredVectorKind8(Z%dPdL1L, Z%dPdL1C)
        if(mg_CurrentWallNb == mg_NbWall) then
          call RotateStoredVector(P%accelerationsaveLL, P%accelerationsaveL, P%accelerationsaveC)
          call RotateStoredVector(P%nhatdotsaveLL,      P%nhatdotsaveL,      P%nhatdotsaveC)
        end if
      end if
      
    end if

  end subroutine CalculateDifferentialNoise

  subroutine RotateStoredReal2DKind8 (arrayL, arrayC)
    real(kind=8), dimension(:), pointer:: tempArray, arrayL, arrayC
    
    tempArray => arrayL
    arrayL    => arrayC
    arrayC    => tempArray
    
  end subroutine RotateStoredReal2DKind8

  subroutine RotateStoredReal3D (arrayLL, arrayL, arrayC)
    real, dimension(:), pointer:: tempArray, arrayLL, arrayL, arrayC
    
    tempArray => arrayLL
    arrayLL   => arrayL
    arrayL    => arrayC
    arrayC    => tempArray
    
  end subroutine RotateStoredReal3D


  subroutine RotateStoredReal3DKind8 (tempArray, arrayLL, arrayL, arrayC)
    real(kind=8), dimension(:), pointer:: tempArray, arrayLL, arrayL, arrayC
    
    tempArray => arrayLL
    arrayLL   => arrayL
    arrayL    => arrayC
    arrayC    => tempArray
    
  end subroutine RotateStoredReal3DKind8

  subroutine RotateStoredVector(arrayLL, arrayL, arrayC)
    type(vector), dimension(:), pointer:: tempArray, arrayLL, arrayL, arrayC
    
    tempArray => arrayLL
    arrayLL   => arrayL
    arrayL    => arrayC
    arrayC    => tempArray
    
  end subroutine RotateStoredVector

  subroutine RotateStoredVectorKind8(arrayL, arrayC)
    type(vectorKind8), dimension(:), pointer:: tempArray, arrayL, arrayC
    
    tempArray => arrayL
    arrayL    => arrayC
    arrayC    => tempArray
    
  end subroutine RotateStoredVectorKind8
  
  subroutine PredictPeggBB(PeggData, initialR, obsCBList, sourceTime, dt, dTau, time, segmentSize)
    implicit none
    type(CBStructure), pointer:: obsCBList
    type(PEGG), pointer:: PeggData
    real:: initialR, sourceTime, dTau, segmentSize
    real(kind=8):: dt
    integer, intent(in):: time
    
    real(kind=8), dimension(:), pointer:: currentIntensity
    type(vector):: tempVect
    real:: dtKind4, ObsTime
    integer:: BBkey
    
    tempVect = vectorSetCoordinates(0.0, 0.0, 0.0)
    dtKind4  = dt
    call ObserverTime (tempVect, sourceTime, ObsTime, &
                        dtKind4, dTau, initialR, obsCBList)

    ! Pegg's BB predictions should be taken as an 
    ! an average over each 0.5s window.
    if (ObsTime.ge.PeggData%obsTimeArray(PeggData%iTau) .and. &
      ObsTime.lt.(PeggData%obsTimeArray(PeggData%iTau)+segmentSize)) then
!	if (PeggData%maxSrctime.eq.huge(PeggData%maxSrcTime)) PeggData%maxSrcTime = sourceTime+0.5  
      BBkey = 1
      if (PeggData%timetype.eq.LOAD_TIMETYPE_PERIODIC) BBkey = getBBoffset(size(PeggData%tau), &
          PeggData%timeOffset, Peggdata%period, time)
      ! Pegg's method considers terms related to the whole blade (eg. blade area), 
      ! thus this subroutine should be called after all the data for the nodes has been computed.
      call ComputePeggBroadbandTerms(PeggData, mg_radiationVector, BBkey)
        currentIntensity => PeggData%Intensity(:,PeggData%iTau)
       call ComputePeggIntensity(PeggData, currentIntensity, BBkey)
      PeggData%nbContributions(PeggData%iTau) = PeggData%nbContributions(PeggData%iTau) + 1 
!    else if (PeggData%nbContributions(PeggData%iTau).gt.0) then
    elseif (ObsTime.ge.(PeggData%obsTimeArray(PeggData%iTau)+segmentSize)) then
      PeggData%iTau = PeggData%iTau + 1
    end if
    
  end subroutine PredictPeggBB

  !function GetObsTimeArray (tMin, dtLocal, tLocal) result(tGlobal)
  !  real::tMin, dtLocal, tempTime
  !  real, dimension(:):: tLocal
  !  real(kind=8), dimension(:), allocatable:: tGlobal
  !  integer:: i, it, ntLocal, ntGlobal
  !  tempTime = 0.0
  !  ntLocal = size(tLocal)
  !  do i=1, ntLocal-1  
  !    if (tLocal(i+1) > tMin) then
  !      it = ceiling((tLocal(i)-tMin)/dtLocal)
  !      tempTime =  tMin + dtLocal*it
  !      exit
  !    end if
  !  end do
  !  ntGlobal = floor((tLocal(ntLocal)-tempTime)/dtLocal)
  !  allocate(tGlobal(ntGlobal))
  !  do i = 1, size(tGlobal)
  !    tGlobal(i) = tempTime + dtLocal*(i-1)
  !  end do
  !  
  !end function GetObsTimeArray


  recursive subroutine IntegrateOverWall(P, BPMData, CBList, obsCBList, dt, thicknessNoise, &
    loadingNoise, PTGradient, PL1Gradient, iBlankArray, observerTimeArray, iTau, tauMax,  &
    keyIndexSurf, keyIndexLoad, nbZones, obsCount, lastGeoZone, lastBlade, radDistance,segmentSize) 

       
    type(patch)::P
    real(kind=8), intent(in)::dt
    type(CBStructure), pointer::CBList, obsCBList
    real(kind=8), dimension(:), intent(inout)::thicknessNoise, loadingNoise
    type(vectorKind8), dimension(:), intent(inout)::PTGradient,PL1Gradient
    type(BPM), pointer:: BPMData
    integer, dimension(:)::iBlankArray
    real(kind=8), dimension(:):: observerTimeArray
    integer, intent(in):: itau, keyIndexSurf, keyIndexLoad, nbZones
    integer, intent(inout):: tauMax
    integer:: obsCount
    real, dimension(:), pointer:: radDistance
    real:: segmentSize
         
    type(vector)::tempObsPosition
    integer::i, n
    logical::  lastWall, lastGeoZone, lastBlade
  
    tempObsPosition = mg_obsInitialVector
    do mg_CurrentWallNb=1, mg_nbWall
      if (mg_CurrentWallNb.eq.mg_NbWall) then
        lastWall = .true.
      else
        lastWall = .false.
      end if
      !call CalculateDifferentialNoise(P, CBList, obsCBList, thicknessNoise, loadingNoise, &
      !                                iBlankArray, PTGradient,PL1Gradient, observerTimeArray, dt, &
      !                                iTau, tauMax, keyIndexSurf, keyIndexLoad, nbZones, .TRUE., lastWall)
      
      call CalculateDifferentialNoise (P, BPMData, CBList, obsCBList, iBlankArray, thicknessNoise, &
      &   loadingNoise, PTGradient,PL1Gradient, observerTimeArray, dt, iTau, tauMax, keyIndexSurf, &
      &   keyIndexLoad, obsCount, lastGeoZone, lastBlade, nbZones, .TRUE., lastWall, radDistance, segmentSize)

       mg_obsInitialVector = tempObsPosition
    end do
    mg_CurrentWallNb = 0
  end subroutine IntegrateOverWall
  
  subroutine Rotate_mg_flowParameters
    type(DataParameters), pointer:: pressureParametersTemp
    
    ! Rotate the storage so that the center is at the current timestep and
    ! the left memory gets reused to calculate the new right data:
    ! We are starting with (from the previous time in this loop):
    !      L          C       R
    ! [itau - 2] [itau - 1] [itau]
    ! So we shift the C and R to the left and take the memory from L and
    ! re-use it to calculate the new R (this avoids having to do any
    ! reallocations)
    !      L        C         R
    ! [itau - 1]  [itau]  [itau +1] <-- Typically the memory from the
    !                                   previous L, overwritten with new
    !                                   data 
    ! We can now take a central difference across L and R to get the time
    ! derivatives at C.
    ! For pressure gradient problem,
    ! additional shifting is needed as following
    !    LLL       LL        L
    ! [itau-3]  [itau-2]  [itau-1]

    pressureParametersTemp => mg_flowParametersL
    
    mg_flowParametersL     => mg_flowParametersC
    mg_flowParametersC     => mg_flowParametersR
    mg_flowParametersR     => pressureParametersTemp

  end subroutine Rotate_mg_flowParameters

  !this subroutine is not currently used. ksb 3/26/2013
  subroutine CounterRotateData(P, CBlist, nSourceBase, iTau, &
                                    pressureParametersTemp, flowData, keyIndexLoad, keyIndexSurf, patchNum)
    type(patch), intent(inout)::p
    integer, intent(in)::iTau, nSourceBase, keyIndexSurf, keyIndexLoad, patchNum
    type(DataParameters), pointer:: pressureParametersTemp
    type(CBStructure), pointer::CBList
    logical:: flowData
    ! Rotate the storage so that the center is at the current timestep and
    ! the left memory gets reused to calculate the new right data:
    ! We are starting with (from the previous time in this loop):
    !      L          C       R
    !   [itau]   [itau + 1] [itau + 2]
    ! So we shift the C and R to the left and take the memory from L and
    ! re-use it to calculate the new R (this avoids having to do any
    ! reallocations)
    !      L        C         R
    ! [itau - 1]  [itau]  [itau +1] <-- Really the memory from the
    !                                   previous L, overwritten with new
    !                                   data 
    ! We can now take a central difference across L and R to get the time
    ! derivatives at C.
    pressureParametersTemp => mg_flowParametersR
    mg_flowParametersR => mg_flowParametersC
    mg_flowParametersC => mg_flowParametersL
    mg_flowParametersL => pressureParametersTemp
    if (iTau > 1) then
      if (flowData) then
        call CalculatePermeableParameters (P, CBList, nSourceBase, & 
                                           mg_flowParametersL, iTau-1, keyIndexSurf, keyIndexLoad)
      else
        call CalculatePressureParameters (P, mg_flowParametersL, iTau-1, keyIndexSurf, keyIndexLoad)
      end if
    else
      ! On the last time through we need to delete the memory associated
      ! with the right memory location, and point it at the center
      ! instead..[dp][al]t
      call DestroyDataParameters (mg_flowParametersL)
      mg_flowParametersL = mg_flowParametersC
    end if
  end subroutine CounterRotateData

    
  !**************************************************************************
  !SUBROUTINE obstime (P,i,iTau,tau,obsTime)
  !  The subroutine finds what time the noise generated by the blade will 
  ! reach the observer.  Many definitions used here are defined in this file
  ! globally.  The problem is trying to find out what time the observer will
  ! hear the noise if the observer is moving. That is where the helper files
  ! zbrac and zbrent come in.
  !**************************************************************************
  subroutine ObserverTime (eta, tau, obsTime, dt, dTau, initialr, observerCBList)

    real, intent(in)::tau, initialr, dTau 
    real,intent(in)::dt
    real,intent(out)::obsTime
    type(CBStructure), pointer::observerCBList
    type(vector), intent(in)::eta
    
    type(vector) :: pos, V, ray, tempObsVec
    real::errorUpper, errorLower, temp
    logical::success
    real::tLow, tHigh, tTolerance
    double precision, save::a1, b1, c1, r2, rM, M2
    real, external:: zbrent
    
    real:: obsTime2
    type(vector):: tempVec
    type(matrix4)::T
    
    mg_srcVector = position(eta)
    mg_obsCBList => observerCBList
    if(.not.associated(observerCBList)) then  ! stationary observer case
      !
      ! If the observer is not moving then the observer time is just defined by the 
      ! speed of sound
      !
      mg_obsVector = mg_obsInitialVector
      if( mg_currentWallNb /= 0 ) then
        mg_obsVector = GetImagePosition(mg_wallArray(mg_currentWallNb),mg_obsVector)
      end if
      mg_radiationVector=mg_obsVector-mg_srcVector
      obsTime=tau+vectorAbsolute(mg_radiationVector)/c
      !ksb debug: !If the observer is not moving then the observer time is just defined by the speed of sound
      !mg_radiationVector=mg_obsInitialVector-mg_srcVector
      !obsTime=tau+vectorAbsolute(mg_radiationVector)/c
    else if (observerCBList%translationType==KNOWN_FUNCTION .and. OnlyUniformlyTranslating(observerCBList)) then  ! the COB with the constant velocity needs to be first right now !ksb 10/8/2017
      ! This is a wind-tunnel type case: don't need to do root-finding here
    !ksb debug: 10/10/2017
    if( firstObsSetup ) then
      call GetUniformTranslation (observerCBList, pos, V)  ! the COB with the constant velocity needs to be first right now !ksb 10/8/2017
      !ksb debug: 10/8/2017
      T=observerTransformationMatrix(mg_obsCBList, 0.)      
      tempObsVec = T*mg_obsInitialVector
      M2 = vectorDotProductDouble(V,V)/(c**2.0)
      a1 = 1.0-M2
      firstObsSetup=.false.
    end if
      !tempObsVec = mg_obsInitialVector + pos
      if( mg_currentWallNb /= 0 ) then
        tempObsVec = GetImagePosition(mg_wallArray(mg_currentWallNb),tempObsVec)
      end if
      ! Calculate the observer time, t:
      ray = tempObsvec - mg_srcVector
      r2 = vectorDotProductDouble(ray,ray)
      rM = vectorDotProductDouble(ray,V)/c
    !ksb debug: 10/10/2017  M2 = vectorDotProductDouble(V,V)/(c**2.0)
      b1 = rM/c+tau  !removed factor of -2.0, LVL
    !ksb debug: 10/10/2017  a1 = 1.0-M2
      c1 = tau**2.0-r2/(c**2.0)
      obsTime =(b1 + sqrt(b1**2.0-a1*c1))/a1  !removed constant coefficients that cancelled out, LVL
      !ksb debug: 10/8/2017
      T=observerTransformationMatrix(mg_obsCBList, obsTime)
      mg_obsVector = T*mg_obsInitialVector
      if( mg_currentWallNb /= 0 ) then
        mg_obsVector = T*GetImagePosition(mg_wallArray(mg_currentWallNb),mg_obsVector)
      end if
     ! mg_obsVector = tempObsVec + obsTime*V

      mg_radiationVector = mg_obsVector-mg_srcVector
    else
      !If the observer is moving arbitrarily it gets much more complicated
      !First off we need to calculate the patch global tau so that the timeError routine can use it
      mg_tauLocal=tau
      !r is the radiation distance guess based on the previous value
      !r=vectorAbsolute(mg_radiationVector)
      !We need to guess an initial range of t based on the guessed radiation distance
      tLow =tau+initialr/c-dtau
      tHigh=tau+initialr/c+dtau
      !This tolerance value is set to dt because anything smaller would be lost at the end
      ! when we interpolate back into the observer time array
      tTolerance=dt ! check this - does it need to be smaller KSB
      !zBrac brakets the results between tLow and tHigh
      !errorUpper and errorLower are the errors in the routine timeError for tLow and tHigh
      !We have them here because we dont want to do the same calculation twice
      call zbrac(timeError, tLow, tHigh, success, errorLower, errorUpper)
      !------------------------------------------------------------------------
      ! Workaround for bug in Windows version of Intel compiler
      ! Added 7/10/2006 by Chris Hennes 
      ! temp = timeError (tLow, tMin)
      !------------------------------------------------------------------------
      !If zbrac was successful then we can go ahead and find the actual observertime
      !If zbrac was not successfull then zbrent will fail and the code will stop
      !zbrent will call timeError routine until the tolerance is met
      if (success) then
        obsTime=zbrent(timeError, tLow, tHigh, tTolerance, errorLower, errorUpper)
      else
        call Error('zbrent could not bracket the observer time. ',                 &
                   'The observer is moving with source time: '//RealToString(tau), &
                   'The initial range was '//RealToString(tau+initialr/c)//' +- '//&
                   RealToString(dTau))
        stop
      end if
    end if
  end subroutine ObserverTime

  !***********************************************************
  !SUBROUTINE timeError(observerTime)
  !  This routine is used by zbrac and zbrent.  Because those files
  ! can only take in a function name, there needs to be global variables
  ! that are used by this function.
  !ARGUMENTS:
  ! - observerTime: The observer time we are guessing.
  ! - tError: The error associated with the guess
  !**************************************************************
  function timeError(observerTime) result(tError)
    real,intent(in)::observerTime

    real::tError
    type(matrix4)::T

    !We first find the transformation matrix for the observer at this observerTime
    T=observerTransformationMatrix(mg_obsCBList, observerTime)
    !We next calculate the vector position of the observer"
    !ksb debug:  I think we need this fix for walls - but it has not been checked yet
    mg_obsVector = T*mg_obsInitialVector
    if( mg_currentWallNb /= 0 ) then
      mg_obsVector = T*GetImagePosition(mg_wallArray(mg_currentWallNb),mg_obsVector)
    end if
    !ksb debug: mg_obsVector=T*mg_obsInitialVector
    !Calculate the radiation vector based on the position of the source
    mg_radiationVector=mg_obsVector-mg_srcVector
    !The error will be the tau time plus the radiation vector divided by the speed of 
    ! sound minus our guess
    tError=(mg_tauLocal+vectorAbsolute(mg_radiationVector)/c)-observerTime 
    !zBrac and zBrent will then take the error value and change their guess depending on that
  end function timeError


  subroutine SourceTime (x, t, dt, tau, sourceCBList)
    type(vector), intent(in) :: x
    real, intent(in)::t, dt
    real, intent(out)::tau
    type(CBStructure), pointer::sourceCBList

    real::errorUpper,errorLower
    logical::success
    real::tauLow, tauHigh, tauTolerance, initialr
    real, external:: zbrent
    ! Store a local copy of the neccessary variables
    nullify(mg_srcCBList)  !ksb debug: check for memory leak here.
    if (associated(sourceCBList)) then
      mg_srcCBList=>sourceCBList
    end if
    mg_tLocal=t
    mg_obsVector = x
    ! Calculate the approximate t range; tauLow and tauHigh don't have to
    ! be inside the tau range of the surface (but they both can't be low
    ! or high) because GetPositionAtTime adjusts the source time if a request
    ! is made out of the range.  I hope that works everywhere! ksb
    initialr = vectorAbsolute(mg_radiationVector)
    tauLow = t-(initialr/c)-dt  
    tauHigh= t-(initialr/c)+dt

    tauTolerance= 0.001*dt  !ksb: 0.1% of dt seems small enough !.000001*dt
    call zbrac(tauError,tauLow,tauHigh,success,errorLower,errorUpper)
    if( .not.success ) then
      !
      !ksb:  It seem that if the initial guess is at the end of the actual tau range
      ! (and GetPositionAtTime will not let it be outside the tau range now)
      ! then zbrac might not be able to bracket the tau range, but it gets very
      ! close and doesn't know it.  We can check that here, by checking the error
      ! returned.  If the error is less than the tauTolerance, then change the sign
      ! of the error term (Lower or Upper - it depends on which one is very close
      ! to zero) so that zbrent will think zbrac was successful.
      !
      if( abs(errorLower)< tauTolerance ) then
        errorLower = -errorLower
        success = .true.
      else if( abs(errorUpper)<tauTolerance ) then
        errorUpper = -errorUpper
        success = .true.
      else
        !
        !ksb: Another reason zbrac can fail is if both initial guesses are outside the
        ! tau range on the same end (both high or both low).  Because GetPositionAtTime
        ! will ensure we can't go outside the tau range of the aperiod data surface,
        ! try expanding the range of the initial guess greately to see if we can bracket
        ! it.
        tauLow = -huge(dt)
        tauHigh=  huge(dt)
        call zbrac(tauError,tauLow,tauHigh,success,errorLower,errorUpper)
      end if
    end if

    !--------------------------------------------------------------------------
    ! Workaround for bug in Windows version of Intel compiler
    ! Added 7/10/2006 by Chris Hennes 
  !   temp = tauError (tauLow)  !ksb: I'm not sure why this was needed.
    !--------------------------------------------------------------------------
    if(success) then
      tau = zbrent(tauError,tauLow,tauHigh,tauTolerance,errorLower,errorUpper)
    else
      call Error('Could not bracket source time.',&
                 'The source is moving with observer time: '//RealToString(t), &
                 'The initial range was '//&
                 trim(RealToString(t-(initialr/c)-dt))//' --> '//&
                 trim(RealToString(t-(initialr/c)+dt)))
      stop
    end if

  end subroutine SourceTime

  
  function tauError (sourceTime) result (tauErrorValue)
    real,intent(inout)::sourceTime
    
    real::tauErrorValue
    type(matrix4)::T
    type(vector):: pos
    
    !ksb:  Need to call GetPostionAtTime first because if the sourceTime
    ! is not in the surface source time range, GetPositionAtTime will adjust
    ! the source time to the nearest endpoint.  sourceTime is used in the
    ! transformation maxtrix - so we would like it to be correct.
    !
    pos = GetPositionAtTime(mg_srcSurface,sourceTime,mg_srcIndex)
    !We first find the transformation matrix for the source at this 
    !observerTime (pay no attention to the function name -- legacy code)
    T=observerTransformationMatrix(mg_srcCBList, sourceTime)
    !We next calculate the vector postion of the observer
    mg_srcVector=T*pos !GetPositionAtTime(mg_srcSurface,sourceTime,mg_srcIndex)
    !Calculate the radiation vector based on the position of the source
    mg_radiationVector=mg_obsVector-mg_srcVector
    !The error will be the tau time plus the radiation vector divided by the 
    !speed of sound minus our guess
    tauErrorValue=(sourceTime+vectorAbsolute(mg_radiationVector)/c)-mg_tLocal
    !zBrac and zBrent will then take the error value and change there guess 
    !depending on that
    
  end function TauError 
  
  
  subroutine CreatePermeableParameters(P, flowParameters, nbNodes)
    type(patch):: P
    type(DataParameters):: flowParameters
    integer:: nbNodes
    
    nullify(flowParameters%L,flowParameters%U,flowParameters%nhat)
    if (P%computeLoadingFlag) then
      ! If the memory has not been allocated yet, do so now
      if (.not. associated(flowParameters%L)) then
        allocate (flowParameters%L(nbNodes))
        flowParameters%L = vectorSetCoordinates(0.,0.,0.)
      end if
    end if
    if (P%computeThicknessFlag) then
      ! If the memory has not been allocated yet, do so now
      if (.not. associated(flowParameters%U)) then
        allocate (flowParameters%U(nbNodes))
        flowParameters%U = vectorSetCoordinates(0.,0.,0.)
      end if
    end if
    if (.not. associated(flowParameters%nhat))then
      allocate (flowParameters%nhat(nbNodes))
      flowParameters%nhat = vectorSetCoordinates(0.,0.,0.)
    end if
  end subroutine CreatePermeableParameters
  
  
  subroutine CreatePressureParameters(P, flowParameters, nbNodes)
    type(patch):: P
    type(DataParameters):: flowParameters
    integer:: nbNodes
    real:: sourceTime
    
    nullify(flowParameters%L,flowParameters%nhat)
    if (P%computeLoadingFlag) then
      ! If the memory has not been allocated yet, do so now
      if (.not. associated(flowParameters%L)) then
        allocate (flowParameters%L(nbNodes))
        flowParameters%L = vectorSetCoordinates(0.,0.,0.)
      end if
    end if    
    if (.not. associated(flowParameters%nhat))then
      allocate (flowParameters%nhat(nbNodes))
      flowParameters%nhat = vectorSetCoordinates(0.,0.,0.)
    end if 
  end subroutine CreatePressureParameters
 

  !**************************************************************************
  !SUBROUTINE CalulatepressureParameters (pressureParameters, iTau)
  !  This subroutine calculates the U and L parameters for a permeable surface
  !  at the given timestep and stores them into the given pressureParameters
  !  object.
  !**************************************************************************
  subroutine CalculatePermeableParameters (P, cbList, nBase, pressureParameters, iTau, keyIndexSurf, keyIndexLoad)
    type(patch), intent(inout)::P
    type(DataParameters), pointer:: pressureParameters
    integer, intent(in) :: iTau
    integer :: i, nBase, nbNodes, keyIndexSurf, keyIndexLoad
    type(SingleTimeLoad), pointer:: flow
    type(SingleSurface), pointer:: surf, surfV, surfA
    type(vector) :: velocity, acceleration, momentum
    type(vector), dimension(:), pointer :: U, L
    real :: vn, un, sourceTime, vinf
    type(CBStructure), pointer::CBList

    nbNodes = GetSurfaceNbNodes (P%surfaceObject)
    
    U => pressureParameters%U
    L => pressureParameters%L
     
    sourceTime = P%tau(iTau)
    
    pressureParameters%time = sourceTime 
    call CalculateResultantChangeOfBase(CBList, nBase,sourceTime)
    flow => GetLoading(P%loadObject, keyIndexLoad+P%iTauMin-1)

    surf => GetSurface (P%surfaceObject, keyIndexSurf+P%iTauMin-1)
    call GetSurfaceVelAccel (P%surfaceObject, keyIndexSurf+P%iTauMin-1, surfV, surfA, keyIndexSurf)
    
    do i=1,size(surf%vectors)
      ! Calculate the velocity and acceleration of this point at emission time
      call CalculateVelAccel (surf%vectors(i),  velocity, acceleration)
      ! Tack on the velocity and acceleration due to the flexing surface:
      velocity = velocity + VectorBaseChange(surfV%vectors(i))
      acceleration = acceleration + VectorBaseChange(surfA%vectors(i))
      ! Get the normal vector at this time:
      pressureParameters%nhat(i)=vectorBaseChange(surf%normals(i))
      ! Get the density and momentum at this point and time
      if (P%refFrame == LOAD_REFFRAME_GROUND) then
        momentum = flow%momVectors(i)
      else if (P%refFrame == LOAD_REFFRAME_MIXED) then
        momentum = vectorBaseChange(flow%momVectors(i))
      else if (P%refFrame == LOAD_REFFRAME_BLADE) then
        momentum = vectorBaseChange(flow%momVectors(i)) + flow%densVals(i)*velocity
      else if (P%refFrame == LOAD_REFFRAME_BLADE_SCALEV) then
        vinf = sqrt(velocity*velocity)
        momentum = vinf*vectorBaseChange(flow%momVectors(i)) + flow%densVals(i)*velocity
      end if
      if (associated(U)) then
        U(i) = (1.0 - (flow%densVals(i)/rho))*velocity + (momentum/rho)
      end if
      if (associated(L)) then
        ! Note that the vector multiplication here is really a dot product:
        ! calculate the normal components of the surface and fluid velocities
        un = (momentum/(flow%densVals(i))) * pressureParameters%nhat(i)
        vn = velocity * pressureParameters%nhat(i)
        ! Calculate L
        L(i) = flow%pressvals(i)*pressureParameters%nhat(i) + &
               (un - vn) * momentum
      end if
    end do
  end subroutine CalculatePermeableParameters
  
  subroutine CalculatePressureParameters(P, PressureParameters, iTau, keyIndexSurf, keyIndexLoad)
    type(patch), intent(inout)::P
    type(DataParameters), pointer:: pressureParameters
    integer, intent(in) :: iTau, keyIndexSurf, keyIndexLoad
    integer :: i, nbNodes
    type(SingleTimeLoad), pointer:: loads
    type(SingleSurface), pointer:: surf
    type(vector), dimension(:), pointer :: L
    real :: sourceTime
      
    nbNodes = GetSurfaceNbNodes (P%surfaceObject)

    L => pressureParameters%L

    sourceTime = P%tau(iTau)
    pressureParameters%time = sourceTime

    loads => GetLoading(P%loadObject, keyIndexLoad+P%iTauMin-1)
    surf => GetSurface (P%surfaceObject, keyIndexSurf+P%iTauMin-1)
    do i=1,size(surf%vectors)
      ! Get the normal vector at this time:
      pressureParameters%nhat(i)=vectorBaseChange(surf%normals(i))
      if (associated(L)) then
        ! Calculate L
        L(i)%A(:) = loads%pressVals(i)*pressureParameters%nhat(i)%A(:)
      end if
    end do
    
  end subroutine CalculatePressureParameters
  

  subroutine DestroyDataParameters (p)
    type(DataParameters), pointer:: p
    
    if (associated(p)) then
      if (associated(p%L)) then
        deallocate(p%L)
      end if
      if (associated(p%U)) then
        deallocate(p%U)
      end if
      if (associated(p%nhat)) then
        deallocate(p%nhat)
      end if
    end if
    
  end subroutine DestroyDataParameters


  subroutine DestroyStoredData(P, Z)
    type(patch):: P
    type(storedData):: Z  !ksb cleanup:  this and all cleanup routines need to be reviewed.
                          !ksb cleanup:  they should check if elements exist.

    if(associated(Z%obsTimeHistoryC)) then
      deallocate(Z%obsTimeHistoryL, Z%obsTimeHistoryC)
      nullify   (Z%obsTimeHistoryL, Z%obsTimeHistoryC)
    end if

    if(PressureGradient1AFlag .and.associated(P%accelerationsaveLL)) then !ksb debug
      deallocate(P%accelerationsaveLL, P%accelerationsaveL, P%accelerationsaveC)
      nullify   (P%accelerationsaveLL, P%accelerationsaveL, P%accelerationsaveC)
      deallocate(P%nhatdotsaveLL,      P%nhatdotsaveL,      P%nhatdotsaveC)
      nullify   (P%nhatdotsaveLL,      P%nhatdotsaveL,      P%nhatdotsaveC)
    end if
    if(P%computeThicknessFlag) then
      deallocate(Z%qTdSL,             Z%qTdSC)
      nullify   (Z%qTdSL,             Z%qTdSC)
    end if
    if(P%computeLoadingFlag .or. isomsThicknessNoiseFlag) then
      deallocate(Z%qLdSL,             Z%qLdSC)
      nullify   (Z%qLdSL,             Z%qLdSC)    
    end if
    if(PressureGradient1AFlag.and.P%computeThicknessFlag) then
      deallocate(Z%dPdTL,             Z%dPdTC)
      nullify   (Z%dPdTL,             Z%dPdTC)
    end if
    if(PressureGradient1AFlag.and.P%computeLoadingFlag) then
      deallocate(Z%dPdL1L,            Z%dPdL1C)
      nullify   (Z%dPdL1L,            Z%dPdL1C)
    end if
    !ksb debug: merge result
    !if(PressureGradient1AFlag) then 
    !  deallocate(Z%accelerationsaveLL, Z%accelerationsaveL, Z%accelerationsaveC)
    !  nullify   (Z%accelerationsaveLL, Z%accelerationsaveL, Z%accelerationsaveC)
    !  deallocate(Z%nhatdotsaveLL,      Z%nhatdotsaveL,      Z%nhatdotsaveC)
    !  nullify   (Z%nhatdotsaveLL,      Z%nhatdotsaveL,      Z%nhatdotsaveC)
    !end if    
    
    call DestroyDataParameters(Z%parameterL)
    call DestroyDataParameters(Z%parameterC)
    call DestroyDataParameters(Z%parameterR)
  end subroutine DestroyStoredData   

  
  !**************************************************************************
  !SUBROUTINE NodeIntegration (P,i,time,qTdS,qLdS, velocity, acceleration)
  !  The subroutine computes the both integral contribution at point (i) on 
  ! the patch P, at the source time 'time', for a radiation vector 'radiationVector'.
  !ARGUMENTS:
  ! - P: patch.
  ! - i: integer, position on the patch.
  ! - time: real, source time tau.
  ! - radiationVector: vector, radiation vector.
  ! - qTdS,qLdS: array of real, dimension(P%imaobsVect,P%jMax,P%nTau), Thickness
  !              and loading integrant.
  !REMARKS:
  ! - The radiation vector is computed by the subroutine "retardedtime": 
  !   therefore, the observer position is not needed.
  ! - All the vectors must be expressed in the same frame of reference
  !   (observer frame).
  ! - The area used S corresponds to the exact surface contribution of the 
  !   point (i) of the patch.
  !**************************************************************************
  subroutine NodeIntegration(P, i, time, qTdS, qLdS, velocity, acceleration,&
                             nhat, nhatdot,accelerationdot, nhat2dot, dAdTau, &
                             dPdT, dPdL1, CBList, nSourceBase, &
                             sourcetime, keyIndexSurf, keyIndexLoad)
    use constantsModule, only: THICKNESS_OUT, LOADING_OUT
    type(patch), intent(inout)::P
    integer, intent(in)::i,time
    real(kind=8),dimension(:), intent(out)::qTdS,qLdS
    type(vectorKind8),dimension(:), intent(out)::dPdT,dPdL1
    type(vector)::velocity,acceleration, nhat, nhatdot,dPdTtemp, dPdL1temp, dPdL2temp, Ldot, L2dot
    type(vector), pointer :: L
    type(SingleSurface), pointer :: surf
    real,intent(in) :: dAdTau, sourcetime
    integer :: loadSpecies
    type(vector)::M,rhat,Mdot,accelerationdot,nhat2dot,Udot,U2dot
    real::r,M2,Mr,Mdotr,S,qT,qL, Un, Lr,Udotn,Undot,U2dotn,M2dotr
    type(vector),pointer :: U
    type(CBStructure), pointer :: CBList
    integer :: nSourceBase, keyIndexSurf, keyIndexLoad
 
    surf => GetSurface (P%surfaceObject, keyIndexSurf+P%iTauMin-1)
    !radiation vector:unit vector
    r=vectorAbsolute(mg_radiationVector)
    rhat=(1/r)*mg_radiationVector
    !Mach number:
    M=(1/c)*velocity
    M2=M*M
    Mr=M*rhat
    Mdotr=(1/c)*(acceleration*rhat)
    !area
    S=surf%areas(i) 
    
    if (associated(P%loadObject)) then
      loadSpecies = GetLoadSpecies(P%loadObject)
    endif
    if(P%computeThicknessFlag) then  !ksb debug: I think this should be done when PressureGradient1AFlag is true too.
      call GenerateThicknessNoise(P, i, time, nhat, nhatdot, r, Mr, Mdotr, M2, &
                                  S, velocity, acceleration,accelerationdot,   &
                                  dAdTau,qT,Un,U,Udot,U2dot)
      qTdS(i)=qT/(4*pi)
    end if
    if (isomsThicknessNoiseFlag) then
      call GenerateIsomsThicknessNoise(rhat, nhat, nhatdot, M, r,&
                                       Mr, Mdotr, M2, S, qL)
      qLdS(i)=qL/(4*pi)
    else if((P%ComputeLoadingFlag.and.(.not.P%trailingEdgeFlag)) .and.      &
      (associated(P%loadObject))) then
      call GenerateLoadingNoise(P, i, time, rhat, M, r, Mr, Mdotr,     &
                                M2, S, dAdtau, qL, L, Ldot, L2dot, Lr, &
                                CBList, nSourceBase, sourcetime, keyIndexLoad)
      qLdS(i)=qL/(4*pi)
    else if((P%ComputeLoadingFlag.and.(P%trailingEdgeFlag)) .and.      &
       (associated(P%loadObject) .and. loadSpecies .ne. LOAD_DATATYPE_PRESSURE)) then
      qL = GenerateTrailingEdgeNoise(P, i, Time, r, S, velocity) 
      qLdS(i)=qL/(4*pi)
    end if
    if(PressureGradient1AFlag) then
      M2dotr=(1/c)*(accelerationdot*rhat)
      Mdot=(1/c)*acceleration
      if(P%computeThicknessFlag) then
        Udotn=Udot*nhat
        Undot=U*nhatdot
        U2dotn=U2dot*nhat
        call GenerateThicknessPGrad1A(Un,Udotn,Undot,U2dotn,Udot,nhatdot,U,nhat2dot,r,M,Mdotr,Mr,M2,M2dotr,Mdot,rhat,dPdTTemp,S)
        dPdT(i)%A(:)=dPdTtemp%A(:)/(4*pi)
      end if
      if(P%ComputeLoadingFlag .and. &
        (associated(P%loadObject))) then
        call GenerateLoadingPGrad1A(L,Lr,Ldot,L2dot,r,M,Mdotr,Mr,M2,M2dotr,Mdot,rhat,dPdL1Temp,S)
        dPdL1(i)%A(:) = dPdL1temp%A(:)/(4*pi)
      end if
    end if
  end subroutine NodeIntegration

  !***************************************************************************
  !SUBROUTINE GenerateThicknessNoise(P,i,time,nhat,r,Mr,Mdotr,M2,S,qT,Un,Index)
  !  The subroutine computes the integral contribution qT corresponding to 
  ! Thickness noise at point (i) on the patch P, at the source time 'time'. 
  !ARGUMENTS:
  ! - P: patch
  ! - i: integer, position on the patch
  ! - time: real, source time tau
  ! - nhat: vector, unit normal vector
  ! - r: real, magnitude of the radiation vector 
  ! - Mr: real, dot product of Mach number vector and unit radiation vector
  ! - Mdotr: real: dot product of the derivative of the Mach number vector 
  !                and unit radiation vector
  ! - M2: real, dot product of Mach number vector with itself
  ! - S: real, area
  ! - qT: real: integral value
  ! - Index : integer value to determine the time step for pressure gradient formulation
  !REMARKS:
  ! - The flow velocity U can be either read from an input file or 
  !   computed in some particular situations
  ! - all the vectors must be expressed in the same frame of reference
  !   (observer frame). The function vectbc and vectderiv performed the change of
  !   base (along others things)
  !***************************************************************************
  subroutine GenerateThicknessNoise(P, i, time, nhat, nhatdot, r, Mr, Mdotr, M2, &
                                    S, velocity, acceleration, accelerationdot,  &
                                    dAdtau, qT, Un, U,Udot,U2dot)
    type(patch), intent(inout)::P
    integer, intent(in)::i, time
    type (vector), intent(in)::nhat, nhatdot
    real, intent(in)::r,Mr,Mdotr,M2,S,dAdtau
    real, intent(out)::qT
    type(vector), intent(in) :: acceleration,accelerationdot
    type(vector), intent(in),target ::velocity
    type(vector) :: Udot, U2dot
    type(vector), pointer :: U

    real::Udotn,Undot,Un, qT1, qT2,timedeno
    integer :: loadSpecies

    ! Check to see if we have a permeable surface
    if (associated(P%loadObject)) then
      loadSpecies = GetLoadSpecies(P%loadObject)
    else
      loadSpecies = 0
    end if
      
    if (loadSpecies == LOAD_DATATYPE_FLOW) then
      ! Permeable surface calculations
      U       => mg_flowParametersC%U(i)
      Udot    = (mg_flowParametersR%U(i) - mg_flowParametersL%U(i)) / &
              (mg_flowParametersR%time - mg_flowParametersL%time)
      
      if(PressureGradient1AFlag) then
        if(time==1) then
          timedeno= mg_flowParametersR%time - mg_flowParametersC%time
        else 
          timedeno= mg_flowParametersC%time - mg_flowParametersL%time
        end if
        U2dot=(mg_flowParametersR%U(i) -2.0*mg_flowParametersC%U(i) + &
                 mg_flowParametersL%U(i)) / timedeno**2
      end if
     
      ! Get the vectors that we need
      Udotn = Udot * nhat
      Undot = U * nhatdot
      Un    = U * nhat
    else

      ! Impermeable surface calculations
      Undot  = velocity*nhatdot
      Un     = velocity*nhat
      Udotn  = acceleration*nhat

      U      => velocity
      Udot   = acceleration
      if(PressureGradient1AFlag) then
        U2dot  = accelerationdot
      end if
    end if

    !thickness:
    qT1 = ((((Udotn+Undot))/(r*(1-Mr)**2)) +&
           ((Un*(r*Mdotr+c*(Mr-M2)))/(r**2*(1-Mr)**3)))*rho
    
    ! Changing-area term: added by Chris Hennes, May 2005
    qT2 = ((rho * Un) / (r * (1-Mr)**2)) * dAdtau
    qT = qT1*S + qT2
  end subroutine GenerateThicknessNoise

  !***************************************************************************
  !SUBROUTINE GenerateLoadingNoise(P,i,time,rhat,nhat,M,r,Mr,Mdotr,M2,S,qL)
  !  The subroutine computes the integral contribution qL corresponding to 
  ! Loading noise at point (i(P, i, time, rhat, M, r, Mr, Mdotr, &
  !                                  M2, S, dAdtau, qL)
  !on the patch P, at the source time 'time'. 
  !ARGUMENTS:
  ! - P: patch
  ! - i: integer, position on the patch
  ! - time: real, source time tau
  ! - rhat: vector, unit radiation vector
  ! - nhat: vector, unit normal vector
  ! - M: vector: Mach number vector
  ! - r: real, magnitude of the radiation vector 
  ! - Mr: real, dot product of Mach number vector and unit radiation vector
  ! - Mdotr: real: dot product of the derivative of the Mach number vector 
  !          and unit radiation vector
  ! - M2: real, dot product of Mach number vector with itself
  ! - S: real, area
  ! - qL: real: integral value
  ! - Index : integer value to determine the time step for pressure gradient formulation
  !REMARKS:
  ! - with P=rho*c**2, loadin(P, i, time, rhat, M, r, Mr, Mdotr, &
  !                                  M2, S, dAdtau, qL)
  !g and thickness are expected to have the same
  !   contribution (Isom Thickness)
  ! - all the vectors must be expressed in the same frame of reference
  !   (observer frame). The function vectbc and vectderiv performed the change of
  !   base (along others things)
  !***************************************************************************
  subroutine GenerateLoadingNoise(P, i, time, rhat, M, r, Mr, Mdotr, &
                                  M2, S, dAdtau, qL, L, Ldot, L2dot, Lr, &
                                  CBList, nSourceBase, sourcetime, keyIndexLoad)
    type (patch),intent(in)::P
    integer, intent(in)::i,time, keyIndexLoad
    type (vector), intent(in)::rhat,M
    real,intent(in)::r,Mr,Mdotr,M2,S,dAdtau,sourcetime
    real, intent(out)::qL
    type(vector), intent(out) :: Ldot,L2dot
    real, intent(out) :: Lr
    type(SingleTimeLoad), pointer:: loadDistribution
    type (vector), pointer::L
    type (vector) :: Ldotnext, Ldotprevious
    type (vector), target :: Lb0
    real::Lm,Ldotr,qL1,qL2,timedeno

    type(CBStructure), pointer :: CBList
    integer :: nSourceBase, loadSpecies
    if (associated(P%loadObject)) then
      loadSpecies = GetLoadSpecies(P%loadObject)
      select case (loadSpecies)
        case (LOAD_DATATYPE_FLOW) 
          L     => mg_flowParametersC%L(i)
          Ldot  = (mg_flowParametersR%L(i) - mg_flowParametersL%L(i)) / &
                  (mg_flowParametersR%time - mg_flowParametersL%time)

          if(PressureGradient1AFlag) then
            if(time==1) then
              timedeno= mg_flowParametersR%time - mg_flowParametersC%time
            else 
              timedeno= mg_flowParametersC%time - mg_flowParametersL%time
            end if
            L2dot=(mg_flowParametersR%L(i) -2.0*mg_flowParametersC%L(i) + &
                     mg_flowParametersL%L(i)) / timedeno**2
          end if 

          Lm    = L*M
          Lr    = L*rhat
          Ldotr = Ldot*rhat
      
        case (LOAD_DATATYPE_PRESSURE)
          L     => mg_flowParametersC%L(i)
          Ldot  = (mg_flowParametersR%L(i) - mg_flowParametersL%L(i)) / &
                  (mg_flowParametersR%time - mg_flowParametersL%time)      

          if(PressureGradient1AFlag) then
            if(time==1) then
              timedeno= mg_flowParametersR%time - mg_flowParametersC%time
            else 
              timedeno= mg_flowParametersC%time - mg_flowParametersL%time
            end if 
            L2dot=(mg_flowParametersR%L(i) -2.0*mg_flowParametersC%L(i) + &
                       mg_flowParametersL%L(i)) / timedeno**2
          end if

          Lm = L*M
          Lr = L*rhat
          Ldotr = Ldot*rhat
        case (LOAD_DATATYPE_LOADING)
          loadDistribution => GetLoading(P%loadObject,keyIndexLoad+P%iTauMin-1)
          L   => loadDistribution%loadVectors(i)
          ! HENNES NOTE: Added reference frame check on 10/17/2010: it was assumed
          ! prior to that date that loading vectors were always in the blade-fixed
          ! frame.
          if (P%refFrame == LOAD_REFFRAME_GROUND) then
            Lb0 = L
          else
            Lb0 = VectorBaseChange(L)
          end if
          Lm  = Lb0*M
          Lr  = Lb0*rhat
          nullify(L)
          L => Lb0
          call DerivativeL(P,i,keyIndexLoad+P%iTauMin-1,Ldot)

          Ldotr=Ldot*rhat   

          if(PressureGradient1AFlag) then
            if(time==1) then
              if(GetLoadType(P%loadObject)==CONSTANT_LOADING) then
                L2dot=VectorDerivative(Ldot)
     !      call SecondDerivativeL(P,i,time+P%iTauMin-1,Ldot,L2dot)
              else
                call CalculateResultantChangeOfBase(CBList, nSourceBase, &
                                                    GetLoadingSourceTime(P%loadObject,keyIndexLoad+P%iTauMin)) 
                call DerivativeL(P,i,keyIndexLoad+P%iTauMin,Ldotnext)
                call CalculateResultantChangeOfBase(CBList, nSourceBase, &
                                                    GetLoadingSourceTime(P%loadObject,keyIndexLoad+P%iTauMin-1)) 
                call DerivativeL(P,i,time+P%iTauMin-1,Ldotprevious)
                L2dot = (Ldotnext-Ldotprevious)/(GetLoadingSourceTime(P%loadObject,keyIndexLoad+P%iTauMin)- &
                                            GetLoadingSourceTime(P%loadObject,keyIndexLoad+P%iTauMin-1))
              end if
            else if(time==P%nTau) then
              if(GetLoadType(P%loadObject)==CONSTANT_LOADING) then
                 L2dot=VectorDerivative(Ldot)
      !          call SecondDerivativeL(P,i,time+P%iTauMin-1,Ldot,L2dot)
              else
                call CalculateResultantChangeOfBase(CBList, nSourceBase, & 
                                                    GetLoadingSourceTime(P%loadObject,keyIndexLoad+P%iTauMin-1))  
                call DerivativeL(P,i,keyIndexLoad+P%iTauMin-1,Ldotnext)
                call CalculateResultantChangeOfBase(CBList, nSourceBase, & 
                                                    GetLoadingSourceTime(P%loadObject,keyIndexLoad+P%iTauMin-2)) 
                call DerivativeL(P,i,keyIndexLoad+P%iTauMin-2,Ldotprevious)
                L2dot = (Ldotnext-Ldotprevious)/(GetLoadingSourceTime(P%loadObject,keyIndexLoad+P%iTauMin-1)- &
                        GetLoadingSourceTime(P%loadObject,keyIndexLoad+P%iTauMin-2))
              end if
            else
              if(GetLoadType(P%loadObject)==CONSTANT_LOADING) then
                L2dot = VectorDerivative(Ldot)
                call SecondDerivativeL(P,i,Ldot,L2dot)
              else
                call CalculateResultantChangeOfBase(CBList, nSourceBase, &
                        GetLoadingSourceTime(P%loadObject,keyIndexLoad+P%iTauMin))
                call DerivativeL(P,i,keyIndexLoad+P%iTauMin,Ldotnext) 
                call CalculateResultantChangeOfBase(CBList, nSourceBase, &
                        GetLoadingSourceTime(P%loadObject,keyIndexLoad+P%iTauMin-2))      
                call DerivativeL(P,i,keyIndexLoad+P%iTauMin-2,Ldotprevious)
                L2dot = (Ldotnext-Ldotprevious)/(GetLoadingSourceTime(P%loadObject,keyIndexLoad+P%iTauMin)- &
                        GetLoadingSourceTime(P%loadObject,keyIndexLoad+P%iTauMin-2))
              end if
            end if
            call CalculateResultantChangeOfBase(CBList, nSourceBase, sourcetime)
          end if
      end select 
    end if
    !loading:
    qL1=((Ldotr/(r*(1-Mr)**2))+Lr*(r*Mdotr+c*(Mr-M2))/(r**2*(1-Mr)**3))/c+&
         (Lr-Lm)/(r**2*(1-Mr)**2)
    ! Changing-area term: added by Chris Hennes, May 2005
    qL2= (1.0/c)*(Lr/(r*(1-Mr)**2)) *dAdTau

    qL=qL1*S + qL2

  end subroutine GenerateLoadingNoise
  
  
  function GenerateTrailingEdgeNoise(P, i, iTime, r, S, V) result(qL)
    type(patch),  intent(in) :: P
    integer,      intent(in) :: i, iTime
    real,         intent(in) :: S, r
    type(vector), intent(in) :: V
   
    type(SingleTimeLoad), pointer:: loadDistribution
    type(vector),     pointer:: L
    type(vector)             :: inplaneR, lhat
    real                     :: theta, psi, qL, D, strength
   
    loadDistribution => GetLoading (P%loadObject, iTime)
    L                => loadDistribution%loadVectors(i)
    strength         =  L%A(1)
    lhat             =  vectorSetCoordinates(sqrt(1-L%A(2)**2-L%A(3)**2),L%A(2),L%A(3))
    inplaneR         =  mg_radiationVector-(mg_radiationVector*lhat)*lhat
    theta            =  acos((V*inplaneR)/ &
                             (vectorAbsolute(V)*vectorAbsolute(inplaneR)))
    psi              =  acos((mg_radiationVector*lhat)/r)
    D                =  (cos(0.5*theta))*sqrt(sin(psi))
    qL = D*S*strength/r
    
  end function GenerateTrailingEdgeNoise
  
  subroutine GenerateIsomsThicknessNoise(rhat, nhat, nhatdot, M,&
                                         r, Mr, Mdotr, M2, S, qL)
    type (vector), intent(in)::rhat,M,nhat,nhatdot
    real,intent(in)::r,Mr,Mdotr,M2,S
    real, intent(out)::qL

    type (vector)::L,Ldot
    real::Lm,Lr,Ldotr,qL1

    L=(rho*c**2)*nhat
    Ldot=(rho*c**2)*nhatdot
    Ldotr=Ldot*rhat
    Lm=L*M
    Lr=L*rhat

    !loading:
    qL1=((Ldotr/(r*(1-Mr)**2))+&
        Lr*(r*Mdotr+c*(Mr-M2))/(r**2*(1-Mr)**3))/c+&
        (Lr-Lm)/(r**2*(1-Mr)**2)
    !qL1=Ldotr/(c*r*(1-Mr)**2)
    !qL1=Lr*(r*Mdotr+c*(Mr-M2))/(r**2*(1-Mr)**3)/c
    !qL1=(Lr-Lm)/(r**2*(1-Mr)**2)

    qL=qL1*S
  end subroutine GenerateIsomsThicknessNoise

  subroutine GenerateThicknessPGrad1A(Un,Udotn,Undot,U2dotn,Udot,ndot,U,n2dot, &
                                      r,M,Mdotr,Mr,M2,M2dotr,Mdot,rhat,dPdTTemp,S)
   real, intent(in)  :: Un,Udotn,Undot,U2dotn,r,Mdotr,Mr,M2,M2dotr,S
   type(vector), intent(in) :: Udot,ndot,U,n2dot,M,Mdot,rhat
   type(vector), intent(out) :: dPdTTemp
   real :: Q,Qdot,Q2dot,W,Wdot
   type(vector) :: I1, I2
   Q=rho*Un
   Qdot=rho*(Udotn+Undot)
   Q2dot=rho*(U2dotn+2.0*Udot*ndot+U*n2dot)
   W=Wfunction(r,Mdotr,Mr,M2)
   Wdot=Wdotfunction(r,M2dotr,M,Mdot,Mdotr,Mr,M2)
   
!   I1%A(:)=(c*(Mr*rhat%A(:)-M%A(:))/r*(Qdot*Ufunction(1,2,r,Mr)+Q*W*Ufunction(2,3,r,Mr)) )                    &
!            +rhat%A(:)*(Q2dot*Ufunction(1,2,r,Mr)+Qdot*(Vfunction(1,2,r,Mdotr,Mr,M2)+W*Ufunction(2,3,r,Mr))   &
!                  +Q*(Wdot*Ufunction(2,3,r,Mr)+W*Vfunction(2,3,r,Mdotr,Mr,M2)))

!   I1%A(:)=-1/c*I1%A(:)/(1-Mr)*S
 
  
!   I2%A(:)=((r*Mdot%A(:)-c*Mr*rhat%A(:)+c*M%A(:))/r*Q*Ufunction(2,2,r,Mr) &
!           -(rhat%A(:)-M%A(:))*(Qdot*Ufunction(2,2,r,Mr)+Q*Vfunction(2,2,r,Mdotr,Mr,M2)))/(1-Mr)*S


!------
    I1%A(:)=-1/c*(rhat%A(:)*(Q2dot/(r*(1-Mr)**3)+(3*Qdot*W+Q*Wdot)/(r**2*(1-Mr)**4) &
            +3*Q*W**2/(r**3*(1-Mr)**5)) &
            -c*M%A(:)*(Qdot/(r**2*(1-Mr)**3)+Q*W/(r**3*(1-Mr)**4)) )*S
    
    I2%A(:)=((M%A(:)-rhat%A(:))*Qdot/(r**2*(1-Mr)**3)+(-c*Mr*rhat%A(:)+c*M%A(:)+r*Mdot%A(:))*Q/(r**3*(1-Mr)**3)&
            +2*(M%A(:)-rhat%A(:))*Q*W/(r**3*(1-Mr)**4))*S

   dPdTTemp=I1 +I2
  end subroutine GenerateThicknessPGrad1A

  subroutine GenerateLoadingPGrad1A(L,Lr,Ldot,L2dot,r,M,Mdotr,Mr,M2,M2dotr,Mdot,rhat,dPdLTemp,S)
   real, intent(in)          :: r,Mdotr,Mr,M2,M2dotr,S, Lr
   type(vector), intent(in)  :: L,Ldot,L2dot,M,Mdot,rhat
   type(vector), intent(out) :: dPdLTemp
   type(vector)              :: drhatdtau
   real :: W,Wdot,dLrdtau,Ldotr,L2dotr,Ldotrdot,Lm,Ldotm,Lmdot
   type(vector) :: I3,I4,I5,I6

   W=Wfunction(r,Mdotr,Mr,M2)
   Wdot=Wdotfunction(r,M2dotr,M,Mdot,Mdotr,Mr,M2)
   drhatdtau%A(:)=c/r*(Mr*rhat%A(:)-M%A(:))
   dLrdtau=Ldot*rhat+L*drhatdtau
   Ldotr=Ldot*rhat
   L2dotr=L2dot*rhat
   Ldotrdot=Ldot*drhatdtau
   Lm=L*M
   Ldotm=Ldot*M
   Lmdot=L*Mdot

!   I3%A(:)=drhatdtau%A(:)*(Ldotr*Ufunction(1,2,r,Mr)+c*(Lr-Lm)*Ufunction(2,2,r,Mr)+Lr*W*Ufunction(2,3,r,Mr)) &
!          +rhat%A(:)*(L2dotr*Ufunction(1,2,r,Mr)+Ldotrdot*Ufunction(1,2,r,Mr)+Ldotr*Vfunction(1,2,r,Mdotr,Mr,M2) &
!                 +c*(dLrdtau-Ldotm-Lmdot)*Ufunction(2,2,r,Mr)+c*(Lr-Lm)*Vfunction(2,2,r,Mdotr,Mr,M2) &
!           +dLrdtau*W*Ufunction(2,3,r,Mr)+Lr*Wdot*Ufunction(2,3,r,Mr)+Lr*W*Vfunction(2,3,r,Mdotr,Mr,M2)) 

!   I3%A(:)=-1/c**2*I3%A(:)/(1-Mr)*S


!   I4%A(:)=(Ldot%A(:)-dLrdtau*rhat%A(:)-Lr*drhatdtau%A(:))*Ufunction(2,1,r,Mr) +&
!           (L%A(:)-Lr*rhat%A(:))*Vfunction(2,1,r,Mdotr,Mr,M2)

!   I4%A(:)=1/c*I4%A(:)/(1-Mr)*S

!   I5%A(:)=dLrdtau*(rhat%A(:)-M%A(:))*Ufunction(2,2,r,Mr)+Lr*(drhatdtau%A(:)-Mdot%A(:))*Ufunction(2,2,r,Mr) &
!           +Lr*(rhat%A(:)-M%A(:))*Vfunction(2,2,r,Mdotr,Mr,M2)

!   I5%A(:)=-1/c*I5%A(:)/(1-Mr)*S

!   I6%A(:)=(L%A(:)-3.0*Lr*rhat%A(:))/(r**3*(1-Mr))*S

!-----------------

    I3%A(:)=-1/c**2*(rhat%A(:)*(L2dotr+Ldotrdot)/(r*(1-Mr)**3)&
            +c*(-M%A(:)*Ldotr-(-Ldotr+LdotM+Lmdot)*rhat%A(:))/(r**2*(1-Mr)**3)&
            +rhat%A(:)*(3*Ldotr*W+Lr*Wdot)/(r**2*(1-Mr)**4)&
            +c**2*((2*Lr*Mr-Lm*(1+Mr))*rhat%A(:)-(Lr-Lm)*M%A(:))/(r**3*(1-Mr)**3)&
            +c*((Lr*(Mr+2)-3*Lm)*W*rhat%A(:)-Lr*W*M%A(:))/(r**3*(1-Mr)**4)&
            +3*Lr*W**2/(r**3*(1-Mr)**5)*rhat%A(:)   )*S

    I4%A(:)=1/c*((Ldot%A(:)-Ldotr*rhat%A(:))/(r**2*(1-Mr)**2)-c*((3*Lr*Mr-Lm)*rhat%A(:)-Lr*M%A(:)-Mr*L%A(:))/(r**3*(1-Mr)**2)&
            +(L%A(:)-Lr*rhat%A(:))*W/(r**3*(1-Mr)**3))*S

    I5%A(:)=-1/c*( (Ldotr*(rhat%A(:)-M%A(:))-Lr*Mdot%A(:))/r**2/(1-Mr)**3 &
            +(rhat%A(:)*(2*c*Lr*Mr-c*Lm)-M%A(:)*(c*Lr*Mr-c*Lm+c*Lr))/r**3/(1-Mr)**3 &
            +2*Lr*(rhat%A(:)-M%A(:))*W/r**3/(1-Mr)**4)*S


    I6%A(:)=(L%A(:)-3.0*Lr*rhat%A(:))/(r**3*(1-Mr))*S
   
    dPdLTemp=I3+I4+I5+I6
  end subroutine GenerateLoadingPGrad1A
 
  
  subroutine ResetGlobals()
    integer::i
    
    mg_srcVector%A        = 0
    mg_srcIndex           = 0
    mg_radiationVector%A  = 0
    mg_obsVector%A        = 0
    mg_tauLocal           = 0
    mg_tLocal             = 0
    call DestroyDataParameters (mg_flowParametersL)
    call DestroyDataParameters (mg_flowParametersC)
    call DestroyDataParameters (mg_flowParametersR)
    nullify(mg_obsCBList, mg_srcCBList, mg_flowParametersR, &
            mg_flowParametersC, mg_flowParametersL, mg_srcSurface)    
  end subroutine ResetGlobals
  
  
  subroutine CreatePatchObsTimeArrays(P, nbObs)
    implicit none
    type(patch), intent(inout):: P
    integer:: nbObs
    
    allocate(P%ObsTimeInRange(nbObs))
    
  end subroutine CreatePatchObsTimeArrays
  
  ! To save us from having to move the information about the total number of observers
  ! all through the code, the arrays were initialized early on so that the size of the 
  ! array carried the information implicitly.  Now that nTau is available we allocate the fully array.
  subroutine CreatePatchTimeRangeArrays(P, obsCount)
    implicit none
    type(patch)::P
    integer:: obsCount
    real(kind=8):: temp
    
    allocate(P%ObsTimeInRange(obsCount)%tMin(P%nTau), P%ObsTimeInRange(obsCount)%tMax(P%nTau))
    P%ObsTimeInRange(obsCount)%tMin = huge(temp)
    P%ObsTimeInRange(obsCount)%tMax = -huge(temp)

  end subroutine CreatePatchTimeRangeArrays
  

  subroutine CreateSourceTimeArray (P, CBList, obsLocation, observerCBList, &
                                    tMin, tMax, nt)
    implicit none
    type(patch), intent(inout):: P
    type(CBStructure), pointer::CBList, observerCBList
    type(Vector), intent(in):: obsLocation
    real, intent(in) :: tMin, tMax
    integer, intent(in) :: nt

    integer :: i, datatype, it, tauMinIndex, tauMaxIndex, nbNodes, iwall
    real :: minTime, dt, temptime, maxtime, tauMin, tauMax, itauMin, itauMax
    type(Matrix4) :: srcT, obsT
    type(vector) :: obsLocationAtTime, pos, etaMin, etaMax

    integer, parameter :: INDEX_SKIP = 1 !5  !ksb - I set this to 1 for debugging
    real, parameter :: PADDING_RATIO = 0.02 ! 2% padding on either side

    ! Find the source time range associated with the observer time range.
    ! Come up with some sort of dt:
    if(debugLevel.gt.13) then
      write(*,*) 'Creating source time array from patch and observer information.'
      write(*,*) 'Observer defined from ', tMin, ' to ', tMax, ' with ', nt, ' timesteps.'
      write(*,*) 'Starting observer position is ', obsLocation, ' with ', &
                GetSize(observerCBList), ' number of base chages.'
      write(*,*) 'Patch has ', GetSize(CBList), ' number of base chages.'
    end if
    dt = (tMax-tMin) / real(nt-1)
    ! Find the source time associated with the first observer time:
    obsT = IMatrix4
    if (associated (observerCBList)) then
      obsT = observerTransformationMatrix (observerCBList, tMin)
    end if
    minTime = HUGE(minTime)
    if(debugLevel.gt.13) then
      write(*,*) 'dt calculated as ', dt
      write(*,*) 'Obs location at the starting time is ', obsLocationAtTime
    end if


    mg_srcSurface => P%surfaceObject
    tauMin = GetSurfaceMinTime(mg_srcSurface)
    if(debugLevel.gt.13) then
      write(*,*) 'Minimum source time is ', tauMin
    end if
    if( getSurfaceType(mg_srcSurface) == APERIODIC_SURFACE ) then
      srcT = observerTransformationMatrix (CBList,tauMin)
    else
    ! tauMin not defined for CONSTANT or PERIODIC surfaces !ksb: should it be?
      srcT = Imatrix4
    end if
    nbNodes = GetSurfaceNbNodes(mg_srcSurface)
    !
    ! The wall loop is only needed for the starting time because the 
    ! reflected signal always must originate at an earlier time than the 
    ! direct signal if they both ARRIVE at the same time
    !
    do iwall=0,mg_NbWall
      obsLocationAtTime = obsT * obsLocation 
      if( iwall>0) then
        obsLocationAtTime = GetImagePosition(mg_wallArray(iwall), obsLocationAtTime)    
      end if
      do i=1, nbNodes, INDEX_SKIP
        if(debugLevel.gt.13) then
          write(*,*) 'Testing ', i, ' of ', GetSurfaceNbNodes(mg_srcSurface), &
                     ' skipping every ', INDEX_SKIP
        end if
        mg_srcIndex = i
        etaMin = GetPositionAtTime(mg_srcSurface,tauMin,mg_srcIndex)
        pos = srcT*etaMin
        if(debugLevel.gt.13) then
          write(*,*) 'Source position at minimum time is ', pos 
        end if
        mg_radiationVector=obsLocationAtTime- pos 
        if(debugLevel.gt.13) then
          write(*,*) 'Initial guess at radiation vector is ', mg_radiationVector
        end if
        call SourceTime (obsLocationAtTime, tMin, dt, tempTime, CBList)
        if(debugLevel.gt.13) then
          write(*,*) 'Source time calculated as ', tempTime
        end if
        minTime = MIN(minTime, tempTime)
      end do
    end do
    if(debugLevel.gt.13) then
      write(*,*) 'Minimum time calculated as ', minTime
    end if
    ! Find the source time associated with the last observer time:
    obsT = IMatrix4
    if (associated (observerCBList)) then
      obsT = observerTransformationMatrix (observerCBList, tMax)
    end if
    obsLocationAtTime = obsT * obsLocation
    maxTime = -HUGE(maxTime)
    if(debugLevel.gt.13) then
      write(*,*) 'Obs location at the last time is ', obsLocationAtTime
    end if

    tauMax = GetSurfaceMaxTime(mg_srcSurface)
    if(debugLevel.gt.13) then
      write(*,*) 'Maximum source time is ', tauMax
    end if
    if( getSurfaceType(mg_srcSurface) == APERIODIC_SURFACE ) then
      srcT = observerTransformationMatrix (CBList,tauMax)
    else
    ! tauMax not defined for CONSTANT or PERIODIC surfaces 
      srcT = Imatrix4
    end if
    do i=1, GetSurfaceNbNodes(mg_srcSurface), INDEX_SKIP
      if(debugLevel.gt.13) then
        write(*,*) 'Testing ', i, ' of ', GetSurfaceNbNodes(P%surfaceObject), &
                    ' skipping every ', INDEX_SKIP
      end if
      mg_srcIndex = i
      etaMax = GetPositionAtTime(mg_srcSurface,tauMax,mg_srcIndex)
      pos = srcT*etaMax
      if(debugLevel.gt.13) then
        write(*,*) 'Source position at maximum time is ', pos 
      end if
      mg_radiationVector=obsLocationAtTime-pos 
      if(debugLevel.gt.13) then
        write(*,*) 'Initial guess at radiation vector is ', mg_radiationVector
      end if
      call SourceTime (obsLocationAtTime, tMax, dt, tempTime, CBList)
      if(debugLevel.gt.13) then
        write(*,*) 'Source time calculated as ', tempTime
      end if
      maxTime = MAX(maxTime, tempTime)
    end do
    if(debugLevel.gt.13) then
     write(*,*) 'Maximum time calculated as ', maxTime
    end if
    !ksb - I don't think this padding is correct, but it doesn't pass the current
    !checksuite, so I don't know what all eventually get affected.  Need to try 
    !to clean this up sometime.  Setting tmin and tmax in observerIn properly
    !avoids the problem here.
    !
    ! Pad the times to make sure we got everything:
    ! (changed to use P%dTau 1/26/2011 by KSB)
    ! padding = (maxTime - minTime) * PADDING_RATIO
    minTime = minTime - 12.*P%dTau !2.*P%dTau ! padding  !ksb - minTime and maxTime are really min & max tau
    maxTime = maxTime + 12.*P%dTau !ksb deubg 2.*P%dTau 
    ! padding  - I would really like to get this fixed so no padding is needed.- added 12 to match was Ben
    ! is using in his version at Bell (which this will replace).  I get identical answers now. 7/3/2013 KSB
    if(debugLevel.gt.13) then
      write(*,*) 'With padding time range is from ', minTime, ' to ', maxTime
    end if
    ! We now have the minimum and maximum required source times, as REALs. If
    ! the position and loading data are continuous functions (i.e - constant
    ! in time except for changes of base) then we can just create the
    ! neccessary array based on this data and dTau. Alternately, if the data
    ! is periodic, then we can always find a tau range that will cover the
    ! entire t range, but we must use the discrete tau given by the input
    ! data. The last case is if we have aperiodic data: then, not only do we
    ! have to use the discrete tau given in the data, but we may not actually
    ! have enough data to do the calculation.

    ! First, figure out what type of case we have. Check the surface first:
    dataType = GetSurfaceType (P%surfaceObject)
    if (dataType == CONSTANT_SURFACE) then
      ! We have to check the other data if the surface is constant
      if (associated (P%loadObject)) then
        dataType = GetLoadType (P%loadObject)
      end if
    end if
    if(debugLevel.gt.13) then
      write(*,*) 'Data type of surface is ', dataType
    end if

    ! If the data is all constant, then we can just make up the array based on
    ! minTime, maxTime and dTau
    if (dataType == CONSTANT_SURFACE .or. dataType == CONSTANT_LOADING) then
      if(debugLevel.gt.13) then
        write(*,*) 'Data type is constant, creating tau array' 
      end if
      P%nTau = ceiling( abs(maxTime-minTime)/P%dTau ) !ksb debug: floor(abs(maxTime-minTime) / P%dTau) + 1
      allocate (P%tau(P%nTau))
      P%tau(P%nTau) = minTime + (P%nTau-1)*P%dTau
      P%tau = CreateEvenlySpacedArray(minTime, P%tau(P%nTau), P%ntau)
      if (debugLevel >= 5.and.IsMaster()) then
        write(*,*) "Created tau array from dTau = ",P%dtau 
      end if
    else if (GetSurfaceType (P%surfaceObject) /= CONSTANT_SURFACE) then
      if(debugLevel.gt.13) then
        write(*,*) 'Surface type is not constant, creating tau array from surface object' 
      end if
      call CreateSurfaceTauArray (P%surfaceObject, minTime, maxTime, P%tau, P%iTauMin)
      if (debugLevel >= 5.and.IsMaster()) then
        write (*,*) "Created tau array from surface."
      end if
    else if (associated (P%loadObject)) then
      if(debugLevel.gt.13) then
        write(*,*) 'Loading data is not constant, creating tau array from loading data object' 
      end if
      call CreateLoadingTauArray (P%loadObject, minTime, maxTime, P%tau, P%iTauMin)
      if (debugLevel >= 5.and.IsMaster()) then
        write (*,*) "Created tau array from loading."
      end if
    end if
    P%nTau = size(P%tau)
    if (P%nTau > 1) then
      P%dTau = P%tau(2) - P%tau(1)  !ksb - not sure we ever check to see tau data is evenly spaced for aperiodic input
    end if

  end subroutine CreateSourceTimeArray
 

  !****************************************************************************
  !SUBROUTINE DerivativeL(P,i,time,Ldot)
  !  The subroutine computes the derivative of the loading at the
  ! point (i,j) on the patch, at time "time".
  !ARGUMENTS:
  ! - P: patch
  ! - i,j: integer, position on the patch
  ! - time: real, source time tau
  ! - Ldot: vector, derivative of the loading
  !REMARKS:
  ! - The central difference approximation is used: f'(t)=(f(t-dt)-f(t+dt))/2dt
  ! - for the last derivative, the following formula is used:
  !   f'(t)=(f(t)-f(t-dt))/dt
  ! - for the first derivative, the following formula is used:
  !   f'(t)=(f(t+dt)-f(t))/dt//
  ! - The function "vectbc" is used to expressed to loading vectors 
  !   in the observer frame of reference
  !********************************************************************** 
  subroutine DerivativeL(P, i, itime, Ldot)
    implicit none
    type(patch), intent(in)::P
    type(vector), intent(out)::Ldot
    integer,intent(in)::i,itime

    type(vector)::ihat, jhat, khat, ihatdot, jhatdot, khatdot
    type(SingleTimeLoad), pointer:: load
    type(vector)::lDeriv

    if (P%refFrame == LOAD_REFFRAME_GROUND) then
      lDeriv = GetLoadingDerivative (P%loadObject, i, itime)
      Ldot = lDeriv
    else
      ihat = vectorSetCoordinates (1.0, 0.0, 0.0)
      jhat = vectorSetCoordinates (0.0, 1.0, 0.0)
      khat = vectorSetCoordinates (0.0, 0.0, 1.0)
      ihatdot = vectorDerivative (ihat)
      jhatdot = vectorDerivative (jhat)
      khatdot = vectorDerivative (khat)
      ihat = vectorBaseChange (ihat)
      jhat = vectorBaseChange (jhat)
      khat = vectorBaseChange (khat)
      load => GetLoading (P%loadObject, itime)
      lDeriv = GetLoadingDerivative (P%loadObject, i, itime)
      Ldot = lDeriv%A(1)*ihat + load%loadVectors(i)%A(1)*ihatdot + &
             lDeriv%A(2)*jhat + load%loadVectors(i)%A(2)*jhatdot + &
             lDeriv%A(3)*khat + load%loadVectors(i)%A(3)*khatdot
    end if
    
  end subroutine derivativeL

  !****************************************************************************
  !SUBROUTINE DerivativeL(P,i,time,Ldot)
  !  The subroutine computes the derivative of the loading at the
  ! point (i,j) on the patch, at time "time".
  !ARGUMENTS:
  ! - P: patch
  ! - i,j: integer, position on the patch
  ! - time: real, source time tau
  ! - Ldot: vector, derivative of the loading
  !REMARKS:
  ! - The central difference approximation is used: f'(t)=(f(t-dt)-f(t+dt))/2dt
  ! - for the last derivative, the following formula is used:
  !   f'(t)=(f(t)-f(t-dt))/dt
  ! - for the first derivative, the following formula is used:
  !   f'(t)=(f(t+dt)-f(t))/dt//
  ! - The function "vectbc" is used to expressed to loading vectors 
  !   in the observer frame of reference
  !********************************************************************** 
  subroutine SecondDerivativeL(P,i,Ldot, L2dot)
    implicit none
    type(patch), intent(in)::P
    type(vector), intent(in)::Ldot
    type(vector), intent(out)::L2dot
    type(vector)::ihat, jhat, khat, ihatdot, jhatdot, khatdot
    integer,intent(in)::i

    ihat = vectorSetCoordinates (1.0, 0.0, 0.0)
    jhat = vectorSetCoordinates (0.0, 1.0, 0.0)
    khat = vectorSetCoordinates (0.0, 0.0, 1.0)
    ihatdot = vectorDerivative (ihat)
    jhatdot = vectorDerivative (jhat)
    khatdot = vectorDerivative (khat)
    ihat = vectorBaseChange (ihat)
    jhat = vectorBaseChange (jhat)
    khat = vectorBaseChange (khat)
    L2dot=Ldot%A(1)*ihatdot+Ldot%A(2)*jhatdot+Ldot%A(3)*khatdot

  end subroutine SecondderivativeL

  !***********************************************************************
  !SUBROUTINE RecordSigma(P,i,iTau,qLdS,qTdS,velocity,acceleration) result(numberVariablesCounter)
  !  A patch needs to record it's information, and that is done here. After
  ! the integration is done for the time step iTau, this function is called so 
  ! that all the information is recorded.
  !AGRUEMENTS:
  ! - P: The patch that is being recorded.
  ! - i, j: The span and width location of the point we are on.
  ! - qLdS: Loading noise created by this point.
  ! - qTdS: thickness noise created by this point.
  ! - velocity: Vecocity vector at this point.
  ! - acceleration: Acceleration of this point.
  ! - numberVariablesCounter: The counter of the number of variables.
  !************************************************************************
  subroutine RecordSigma(P, i, iTau, qLdS, qTdS, velocity, acceleration, obsTime, press,&
                         keyIndexSurf, keyIndexLoad)
    implicit none
    type(patch), intent(inout)::P
    !i is the node, iTau is the source time designator (both integers)
    integer, intent(in)::i,iTau, keyIndexSurf, keyIndexLoad
    real(kind=8), dimension(:), intent(in)::qTdS,qLdS
    real(kind=8)::obsTime
    real:: vinf
    type(vector), intent(in)::velocity,acceleration
    type(SingleTimeLoad), pointer :: load, singPressureParameters, flow 
    type(DataParameters), pointer:: press 
    type(SingleSurface), pointer :: surf 

    type(Vector)::tempVector, rhat
    integer::counter, loadSpecies
    real::r,Mdotr    
    
    if (associated(P%loadObject)) then
      loadSpecies = GetLoadSpecies(P%loadObject)
    endif
    surf => GetSurface (P%surfaceObject, keyIndexSurf+P%iTauMin-1)

    counter = P%numberVariables
!
! Not sure this code section is correct.  I'm not sure that face centered data will
! always be ready in surf%nodecoord when it gets here.  KSB
!
    if(associated(surf%vectors) )then ! node centered data
      tempVector = position(surf%vectors(i))
    else if( associated(surf%nodecoord) ) then
      tempVector = position(surf%nodecoord(i)) !face centered data
    else 
      call Error('Error reading position vectors in RecordSigma')
    end if

    P%sigmaSurface(1:3,i,iTau)=tempVector%A(:)
    if( GetSurfaceIblanking(P%surfaceObject) ) then
      P%sigmaSurface(4,i,iTau)=surf%iblanks(i)
      P%sigmaSurface(5,i,iTau)=P%tau(iTau)
      P%sigmaSurface(6,i,iTau)=obsTime
    else
      P%sigmaSurface(4,i,iTau)=P%tau(iTau)
      P%sigmaSurface(5,i,iTau)=obsTime
    end if
    if(machSigmaFlag) then
       counter=counter+1
       P%sigmaSurface(counter,i,iTau)=vectorAbsolute((1/c)*velocity)
    end if
    if(loadingNoiseSigmaFlag) then
       counter=counter+1
       if(P%computeLoadingFlag.or.isomsThicknessNoiseFlag) then 
          P%sigmaSurface(counter,i,iTau)=qLdS(i)/surf%areas(i)
       else 
          P%sigmaSurface(counter,i,iTau)=0
       end if
    end if
    if(thicknessNoiseSigmaFlag) then
       counter=counter+1
       if(P%computeThicknessFlag) then 
          P%sigmaSurface(counter,i,iTau)=qTdS(i)/surf%areas(i)
       else 
          P%sigmaSurface(counter,i,iTau)=0
       end if
    end if
    if(totalNoiseSigmaFlag) then
       counter=counter+1
       if(P%computeThicknessFlag) then
          if(P%computeLoadingFlag.or.isomsThicknessNoiseFlag) then 
             P%sigmaSurface(counter,i,iTau)=(qTdS(i)+qLdS(i))/surf%areas(i)
          else 
             P%sigmaSurface(counter,i,iTau)=qTdS(i)/surf%areas(i)
          end if
       else
          if(P%computeLoadingFlag.or.isomsThicknessNoiseFlag) then 
             P%sigmaSurface(counter,i,iTau)=qLdS(i)/surf%areas(i)
          else 
             P%sigmaSurface(counter,i,iTau)=0
          end if
       end if
    end if
    if(normalSigmaFlag) then
       tempVector=VectorBaseChange(surf%normals(i))
       P%sigmaSurface((counter+1):(counter+3),i,iTau)=tempVector%A(:)
       counter=counter+3
    end if
    if(velocitySigmaFlag) then
       P%sigmaSurface((counter+1):(counter+3),i,iTau)=velocity%A(:)
       counter=counter+3
    end if
    if(accelerationSigmaFlag) then
       P%sigmaSurface((counter+1):(counter+3),i,iTau)=acceleration%A(:)
       counter=counter+3
    end if  
    if(loadingSigmaFlag) then
       if(P%computeLoadingFlag .and. associated(P%loadObject) .and. loadSpecies == LOAD_DATATYPE_LOADING) then
          load => GetLoading(P%loadObject,keyIndexLoad+P%iTauMin-1)
          if (P%refFrame /= LOAD_REFFRAME_GROUND) then
            tempVector=VectorBaseChange(load%loadVectors(i))
          else
            tempVector = load%loadVectors(i)
          end if
       elseif (P%computeLoadingFlag .and. associated(P%loadObject) .and. loadSpecies == LOAD_DATATYPE_PRESSURE) then
         tempVector=press%L(i)  
        else
          tempVector%A(:) = 0.0
       end if
       P%sigmaSurface((counter+1):(counter+3),i,iTau)=tempVector%A(:)
       counter=counter+3
    end if
    if(densitySigmaFlag) then
       counter=counter+1
       if (associated(P%loadObject) .and. loadSpecies == LOAD_DATATYPE_FLOW) then
         flow => GetLoading(P%loadObject,keyIndexLoad+P%iTauMin-1)
         P%sigmaSurface(counter,i,iTau)=flow%densVals(i)
       else
         P%sigmaSurface(counter,i,iTau)=rho
       end if
    end if
    if(momentumSigmaFlag) then
      if (associated(P%loadObject) .and. loadSpecies == LOAD_DATATYPE_FLOW) then
        flow => GetLoading(P%loadObject,keyIndexLoad+P%iTauMin-1)
        if (P%refFrame == LOAD_REFFRAME_GROUND) then
          tempVector = flow%momVectors(i)
        else if (P%refFrame == LOAD_REFFRAME_MIXED) then
          tempVector = vectorBaseChange(flow%momVectors(i))
        else if (P%refFrame == LOAD_REFFRAME_BLADE) then
          tempVector = vectorBaseChange(flow%momVectors(i)) + flow%densVals(i)*velocity      
        else if (P%refFrame == LOAD_REFFRAME_BLADE_SCALEV) then
          vinf = sqrt(velocity*velocity)
          tempVector = vinf*vectorBaseChange(flow%momVectors(i)) + flow%densVals(i)*velocity
        end if
        P%sigmaSurface((counter+1):(counter+3),i,iTau)=tempVector%A(:)
      else
        P%sigmaSurface((counter+1):(counter+3),i,iTau) = 0
      end if
      counter=counter+3
    end if
    if(pressureSigmaFlag) then
      counter=counter+1
      if (associated(P%loadObject) .and. loadSpecies == LOAD_DATATYPE_PRESSURE) then
        singPressureParameters => GetLoading(P%loadObject,keyIndexLoad+P%iTauMin-1)
        P%sigmaSurface(counter,i,iTau)=singPressureParameters%pressVals(i)
      else if (associated(P%loadObject) .and. loadSpecies == LOAD_DATATYPE_FLOW) then
        flow => GetLoading(P%loadObject,keyIndexLoad+P%iTauMin-1)
        P%sigmaSurface(counter,i,iTau)=flow%pressVals(i)
      else
        P%sigmaSurface(counter,i,iTau)=0
      end if
    end if
    if(areaSigmaFlag) then
      counter=counter+1
      P%sigmaSurface(counter,i,iTau)=surf%areas(i)
    end if  
    if(MdotrSigmaFlag) then
      counter=counter+1
      r=vectorAbsolute(mg_radiationVector)
      rhat=(1/r)*mg_radiationVector
      Mdotr=(1/c)*(velocity*rhat)
      P%sigmaSurface(counter,i,iTau)=Mdotr
    end if
    if(iblankSigmaFlag) then
      counter=counter+1
      P%sigmaSurface(counter,i,iTau)=surf%iblanks(i)
    end if
      
  end subroutine recordSigma 

  
  subroutine WritePatchPlot3DStructuredSigma (P,gridUnit, dataUnit, prefix)
    implicit none
    type(patch)::P
    integer, intent(in):: gridUnit, dataUnit
    character(len=*), intent(in) :: prefix
    integer::ivar,stat
    integer:: itau, j, istart
    integer, dimension(2):: dims
    real(kind=4):: temp
    real(kind=4), dimension(:),   allocatable:: temp1DArray
    real(kind=4), dimension(:,:), allocatable:: temp2DArray
    integer, dimension(:,:), allocatable:: itemp2DArray
    
    if (GetSurfaceGridType(P%surfaceObject) /= STRUCTURED_GRID) then
      ! This is not a structured grid: return now
      return
    end if
    
    if (.not. associated(P%sigmaSurface)) then
       call Warning('No sigma data stored for this patch.')
      return
    end if
    if (P%nTau < 1 .or. &
        P%numberVariables+P%numberFunctions < 1) then
       call Warning('No sigma data stored for this patch.')
      return
    end if
    
    call Message (prefix//trim(P%geoTitle))
!
! compact patches do not work as FieldView as iso surfaces, because the 
! marching cubes algorithm doesn't get a cube, only on side of a cube.
! The work around is to output compact patches with the minimum dimension of
! 2 (typically ni=2 in current practice) rather than 1.  This makes the 
! compact patch a volume grid rather than just a surface grid.  To do this
! I just copy the surface out twice - it doesn't seem to matter that there
! are two lines actually on top of each other.  Note:  the iso will only display
! if you use a MESH surface in FieldView - but then it works.  KSB 1/11/06
!
! ksb debug:  need to add iblanks for special cases (compact patch, etc.)
!
    call getSurfaceDimensions(P%surfaceObject,dims)
    if( dims(1)/=1 .and. dims(2)/=1 ) then   ! this is not a compact patch
      allocate (temp2DArray(size(P%sigmaSurface,2),size(P%sigmaSurface,3)))
      do ivar=1,3
        temp2DArray = P%sigmaSurface(ivar,:,:)
        call WriteBinaryReal2DArray(gridUnit, temp2DArray, stat)
        if (stat /=0) then
          call Error ("Failed to write sigma data: out of space.")
          stop
        end if
      end do
      if( GetSurfaceIblanking(P%surfaceObject) ) then
        allocate (temp2DArray(size(P%sigmaSurface,2),size(P%sigmaSurface,3)))
        itemp2DArray = P%sigmaSurface(4,:,:)
        call WriteBinaryInteger2DArray(gridUnit, itemp2DArray, stat)
        deallocate(itemp2DArray)
        istart = 5
      else
        istart = 4
      end if
      do ivar=istart, P%numberVariables+P%numberFunctions
        temp2DArray = P%sigmaSurface(ivar,:,:)
        call WriteBinaryReal2DArray(dataUnit, temp2DArray, stat)
        if (stat /=0) then
          call Error ("Failed to write sigma data: out of space.")
          stop
        end if
      end do
      deallocate(temp2DArray)
    else if( dims(1)==1 ) then  ! this case is for a compact patch ni=1
      do ivar=1,3
        do itau=1,P%nTau
          do j=1,dims(2)
            temp = P%sigmaSurface(ivar,j,itau)
            call WriteBinaryReal(gridUnit, temp, stat)
            call WriteBinaryReal(gridUnit, temp, stat)
            if (stat /=0) then
              call Error ("Failed to write sigma data: out of space.")
              stop
            end if        
          end do
        end do
      end do
      do ivar=4, P%numberVariables+P%numberFunctions
        do itau=1,P%nTau
          do j=1,dims(2)
            temp = P%sigmaSurface(ivar,j,itau)
            call WriteBinaryReal(dataUnit, temp, stat)
            call WriteBinaryReal(dataUnit, temp, stat)
            if (stat /=0) then
              call Error ("Failed to write sigma data: out of space.")
              stop
            end if
          end do
        end do
      end do
    else  ! this case is for compact patch nj=1  (not currently used, but included for completeness)
      allocate (temp1DArray(size(P%sigmaSurface,2)))
      do ivar=1,3
        do itau=1,P%nTau
          temp1DArray = P%sigmaSurface(ivar,:,itau)
          call WriteBinaryReal1DArray(gridUnit, temp1DArray, stat)
          call WriteBinaryReal1DArray(gridUnit, temp1DArray, stat)        
          if (stat /=0) then
            call Error ("Failed to write sigma data: out of space.")
          end if
        end do
      end do
      do ivar=4, P%numberVariables+P%numberFunctions
        do itau=1,P%nTau
          temp1DArray = P%sigmaSurface(ivar,:,itau)
          call WriteBinaryReal1DArray(dataUnit, temp1DArray, stat)
          call WriteBinaryReal1DArray(dataUnit, temp1DArray, stat)
          if (stat /=0) then
            call Error ("Failed to write sigma data: out of space.")
          end if 
        end do
      end do
      deallocate(temp1DArray)
    end if
  
  end subroutine WritePatchPlot3DStructuredSigma

  
  subroutine WritePatchFVUNSHeader (P, dataUnit, names)
    implicit none
    type(patch)::P
    integer, intent(in) :: dataUnit
    character (len=*), intent(in) :: names

    character (len=80) :: fullName
    integer :: stat

    ! Create the boundary name
    write (fullName, '(A80)') trim(names)//":"//trim(P%geoTitle)
    fullname = adjustl (fullName)

    ! Write the boundary section
    call WriteBinaryInteger (dataUnit, 1, stat)
    call WriteBinaryInteger (dataUnit, 1, stat)
    call WriteBinaryString (dataUnit, fullName, stat)
    
  end subroutine WritePatchFVUNSHeader

  
  subroutine WritePatchFVUNSData (P, dataUnit, patchNumber, prefix)
    implicit none
    type(patch)::P
    integer, intent(inout) :: dataUnit, patchNumber
    character(len=*), intent(in) :: prefix
    integer :: stat, t, n, ivar, nbNodes, nbAllNodes
    integer, parameter :: FV_NODES=     1001, FV_FACES=    1002, &
                          FV_ELEMENTS=  1003, FV_VARIABLES=1004, &
                          FV_BNDRY_VARS=1006

    integer, dimension(:,:), pointer :: faces =>NULL()
    integer, dimension(:,:), allocatable:: elements 
    real(kind=4), dimension(:,:), allocatable:: temp2DArray
                          
    if (GetSurfaceGridType(P%surfaceObject) /= UNSTRUCTURED_GRID) then
      ! This is not an unstructured grid: return now
      return
    end if

    if (.not. associated(P%sigmaSurface)) then
      call Warning('No sigma data stored for this patch.')
      return
    end if
    if (P%nTau < 1 .or. &
        P%numberVariables+P%numberFunctions < 1) then
       call Warning('No sigma data stored for this patch.')
      return
    end if
    
    call Message (prefix//trim(P%geoTitle))
    
    patchNumber = patchNumber + 1
    
    ! Write the header
    call WriteBinaryInteger (dataUnit, FV_NODES, stat)
    nbNodes = GetSurfaceNbNodes (P%surfaceObject)
    nbAllNodes = size(P%sigmaSurface,2)*size(P%sigmaSurface,3)
    call WriteBinaryInteger (dataUnit,nbAllNodes,stat)
    
    ! Write the nodes 
    ! P%sigmaSurface(1 for x, 2 for y, 3 for z,nbNodes,nTau)
    allocate(temp2DArray(size(P%sigmaSurface,2), size(P%sigmaSurface,3)))
    do ivar=1,3
      temp2DArray = P%sigmaSurface(ivar,:,:)
      call WriteBinaryReal2DArray(dataUnit, temp2DArray, stat)
      if (stat /=0) then
        call Error ("Failed to write sigma data: out of space.")
      end if
    end do

    ! Create the face data: we read in the connectivity at one time: we 
    ! then need to add a third "volume" dimension corresponding to the source
    ! time. We create faces for each surface at each time.
    call GetSurfaceUnstructuredFaces (P%surfaceObject, faces)
    
    ! Write the faces
    !call WriteBinaryInteger (dataUnit, FV_FACES, stat)
    !call WriteBinaryInteger (dataUnit, patchNumber, stat)
    !call WriteBinaryInteger (dataUnit, 0, stat)
    !call WriteBinaryInteger2DArray (dataUnit, sigmaFaces, stat)

    if (all(faces(:,1) == 3)) then
      ! We have triangles on the surface, or prisms in spacetime
      ! Write the element header:
      call WriteBinaryInteger (dataUnit, FV_ELEMENTS, stat)
      call WriteBinaryInteger (dataUnit, 0, stat) ! Tetrahedra
      n = size(faces,1)*(P%nTau-1)
      call WriteBinaryInteger (dataUnit, 0, stat) ! Hexahedra
      call WriteBinaryInteger (dataUnit, n, stat) ! Prisms
      call WriteBinaryInteger (dataUnit, 0, stat) ! Pyramids

      ! Create the elements:
      allocate (elements(n,6))
      do t=1,P%nTau-1
        do n=1,size(faces,1)
          elements((t-1)*size(faces,1)+n,1) = faces(n,1+1) + (t-1)*nbNodes
          elements((t-1)*size(faces,1)+n,4) = faces(n,2+1) + (t-1)*nbNodes
          elements((t-1)*size(faces,1)+n,6) = faces(n,3+1) + (t-1)*nbNodes
          elements((t-1)*size(faces,1)+n,2) = faces(n,1+1) + (t)*nbNodes
          elements((t-1)*size(faces,1)+n,3) = faces(n,2+1) + (t)*nbNodes
          elements((t-1)*size(faces,1)+n,5) = faces(n,3+1) + (t)*nbNodes
        end do
      end do
      
      ! Write out the elements:
      do n=1,size(elements,1)
        ! Write the element header. We don't care about this, so rather than dealing
        ! with Fieldview's "fv_encode_elem_header" function we just output a 1.
        call WriteBinaryInteger (dataUnit, 786432, stat)
        call WriteBinaryInteger1DArray (dataUnit, elements(n,1:6), stat)
      end do

    else if (all(faces(:,1) == 4)) then
      ! We have quadrilaterals, or hexahedra in spacetime    
      
      ! Write the element header:
      call WriteBinaryInteger (dataUnit, FV_ELEMENTS, stat)
      call WriteBinaryInteger (dataUnit, 0, stat) ! Tetrahedra
      n = size(faces,1)*(P%nTau-1)
      call WriteBinaryInteger (dataUnit, n, stat) ! Hexahedra
      call WriteBinaryInteger (dataUnit, 0, stat) ! Prisms
      call WriteBinaryInteger (dataUnit, 0, stat) ! Pyramids

      ! Create the elements:
      allocate (elements(n,8))
      do t=1,P%nTau-1
        do n=1,size(faces,1)
          elements((t-1)*size(faces,1)+n,1) = faces(n,1+1) + (t-1)*nbNodes
          elements((t-1)*size(faces,1)+n,2) = faces(n,2+1) + (t-1)*nbNodes
          elements((t-1)*size(faces,1)+n,3) = faces(n,4+1) + (t-1)*nbNodes
          elements((t-1)*size(faces,1)+n,4) = faces(n,3+1) + (t-1)*nbNodes
          elements((t-1)*size(faces,1)+n,5) = faces(n,1+1) + (t)*nbNodes
          elements((t-1)*size(faces,1)+n,6) = faces(n,2+1) + (t)*nbNodes
          elements((t-1)*size(faces,1)+n,7) = faces(n,4+1) + (t)*nbNodes
          elements((t-1)*size(faces,1)+n,8) = faces(n,3+1) + (t)*nbNodes
        end do
      end do
      
      ! Write out the elements:
      do n=1,size(elements,1)
        ! Write the element header. We don't care about this, so rather than dealing
        ! with Fieldview's "fv_encode_elem_header" function we just output a 1.
        call WriteBinaryInteger (dataUnit, 1306687, stat)
        call WriteBinaryInteger1DArray (dataUnit, elements(n,1:8), stat)
      end do

    else
      ! We have a mix, or some other type, which FV can't handle
      call Error ("Fieldview only supports hex and prism.",&
                  "All cells must be the same.")
    end if
    
    ! Write out the variables
    call WriteBinaryInteger (dataUnit, FV_VARIABLES, stat)
    do ivar=4, P%numberVariables+P%numberFunctions
      temp2DArray = P%sigmaSurface(ivar,:,:)
      call WriteBinaryReal2DArray(dataUnit, temp2DArray, stat)
      if (stat /=0) then
        call Error ("Failed to write sigma data: out of space.")
      end if
    end do

    ! Write out the boundary variables
    call WriteBinaryInteger (dataUnit, FV_BNDRY_VARS, stat)
     
    deallocate(elements, temp2DArray)

  end subroutine WritePatchFVUNSData
  
  subroutine GetPatchDimensions(P, dimensionsList)
    implicit none
    type(patch)::P
    integer, dimension(:)::dimensionsList

    call GetSurfaceDimensions(P%surfaceObject, dimensionsList(1:2))
    dimensionsList(3)=P%nTau
    dimensionsList(4)=P%numberVariables+P%numberFunctions
  
  end subroutine GetPatchDimensions
 

  subroutine GetPatchImplicitTimeRange (P, patchCBList, tauMin, tauMax, &
                                        dTau, containerNTau, obsCoordinates, &
                                        observerCBList, obsTimes)
    implicit none
    type(patch)::P
    real, intent(inout)::tauMin, tauMax, dTau
    integer, intent(in)::containerNTau
    type(vector), dimension(:,:), intent(in)::obsCoordinates
    type(CBStructure), pointer::patchCBList , observerCBList 
    real, dimension (:,:,:) :: obsTimes

    logical :: hasImplicitRange
    integer :: i, j, n, nTau, obsIMax, obsJMax, numCB, loadSpecies
    real :: minSourceTime, maxSourceTime, oTime, initialr
    type(SingleSurface), pointer :: surf =>NULL()

    ! A parameter to control how many points from the surface are evaluated. The
    ! code will check 1 out of every SURF_SKIP points.
    integer, parameter :: SURF_SKIP=1

    if (associated(P%loadObject)) then
      loadSpecies = GetLoadSpecies(P%loadObject)
    endif
    
    ! This process has two main parts: determining the lower bound on the time
    ! range, and determining the upper bound. It works by doing a retarded time
    ! calculation at each observer position for the first and last source times
    ! at every source position. We want to start recording the signal when the
    ! full signal starts reaching any observer, and to stop recording when the
    ! full signal no longer reaches any observer.

    ! Start with some initialization, so see if we can even do this:
    hasImplicitRange = .false.
    if (associated(P%loadObject)) then
      if (LoadingHasImplicitTau(P%loadObject)) then
        hasImplicitRange = .true.
        nTau = GetLoadingNTau (P%loadObject)
        !ksb debug:  This maxsourceTime line is going to fail because the loadObject only holds 4 time steps
        ! in Ben's modified version.  Ben used the 2 following lines - but that breaks the functionality
        ! i.e., the ability to let the code determine the time range implicityly.
        !
        !minSourceTime = GetLoadingSourceTime(P%loadObject, 1)     !ksb debug uncomment  %BAG 8/10/2011
        !maxSourceTime = GetLoadingSourceTime(P%loadObject, nTau)  !ksb debug 9/25/2013 - nTau will be bigger than 4
        minSourceTime = tauMin
        maxSourceTime = tauMax
      end if
    else if (SurfaceHasImplicitTau(P%surfaceObject)) then
      hasImplicitRange = .true.
      nTau = GetSurfaceNTau (P%surfaceObject)
      !minSourceTime = GetSurfaceSourceTime(P%surfaceObject, 1) !ksb debug - check these, they probably should be fixed to work.
      !maxSourceTime = GetSurfaceSourceTime(P%surfaceObject, nTau)
      minSourceTime = tauMin
      maxSourceTime = tauMax
    end if
    if ((tauMin.ne.tauMax).and.(.not.hasImplicitRange)) then
      if (tauMin.gt.tauMax) then
        call Message ("tauMin greater than tauMax", "This is not recommended")
      end if
      if (dTau.gt.0) then
        hasImplicitRange = .true.
        if (containerNTau.eq.0) then
          nTau = int(abs(tauMax-tauMin)/dTau)
        else 
          if (containerNTau.gt.0) then
            nTau = containerNTau
          else
            call Error ("nTau specified in namelist less than zero")
          end if
        end if
      else if (dTau.eq.0.0) then
        dTau = (tauMax-tauMin)/real(ntau)
      else
        call Error ("dTau specified in namelist less than zero")
      end if
      minSourceTime = tauMin
      maxSourceTime = tauMax
    end if

    if (.not. hasImplicitRange) then
      call Message ("NOTE: Patch '"//trim(P%geoTitle)//&
                    "' does not have an implicit tau range.")
      return
    end if

    ! More initialization:
    call initializeInterpolateIndex ()
    obsIMax=size(obsCoordinates,1)
    obsJMax=size(obsCoordinates,2)
    numCB=GetSize(patchCBList)

    ! Loop over the observers
    do i=1, obsIMax
      do j=1, obsJMax
        ! At a given observer location, calculate the minimum and maximum time
        ! the signal will arrive

        ! Do the minimum time first:
        surf => GetSurface (P%surfaceObject, 1)
        call CalculateResultantChangeOfBase(patchCBList, numCB, minSourceTime)
        mg_obsInitialVector=obsCoordinates(i, j)
        initialr=vectorAbsolute(mg_obsInitialVector)
        do n=1, size(surf%vectors), SURF_SKIP
          ! Estimate the observer time:
          call ObserverTime (surf%vectors(n), minSourceTime, oTime, P%dTau, &
                             P%dTau, initialr, observerCBList)
          ! If oTime is larger than the previously calculated times, then save
          ! it off as the new minTime, since we are looking for the time at
          ! which the signal from ALL sources has been captured.
          obsTimes(i,j,1) = MAX (obsTimes(i,j,1), oTime)
        end do

        ! Do the maximum time next:
        surf => GetSurface (P%surfaceObject, nTau)
        call CalculateResultantChangeOfBase(patchCBList, numCB, maxSourceTime)
        initialr=vectorAbsolute(mg_obsInitialVector)
        do n=1, size(surf%vectors), SURF_SKIP
          ! Estimate the observer time:
          call ObserverTime (surf%vectors(n), maxSourceTime, oTime, P%dTau, &
                             P%dTau, initialr, observerCBList)
          ! If oTime is smaller than the previously calculated times, then save
          ! it off as the new maxTime, since we are looking for the time at
          ! which the signal from ANY source has stopped arriving at this
          ! observer.
          obsTimes(i,j,2) = MIN (obsTimes(i,j,2), oTime)
        end do
      end do
    end do
    
  end subroutine GetPatchImplicitTimeRange
  
  subroutine RotatePatchTimeData(P)
    implicit none
    type(patch):: P
    
    if (getSurfaceType(P%surfaceObject) == APERIODIC_SURFACE) then
        call RotateSurfaceArray(P%surfaceObject)
    end if
    if (associated(P%loadObject)) then
      if(getLoadType(P%loadObject) == APERIODIC_LOADING) then
        call RotateLoadingArray(P%loadObject)
      end if
    end if
    
  end subroutine RotatePatchTimeData
  
  subroutine GetPatchTimeArrays(P, CBList, obsInitVector, ObsCBList, &
                                observerTimeArray, nTau, tauMin, geoInfo, loadInfo)
    type(patch)::P
    type(Vector), intent(in)::obsInitVector
    type(CBStructure),pointer::CBList, obsCBList    
    integer,dimension(:):: geoInfo,loadInfo
    real(kind=8), dimension(:), intent(in)::observerTimeArray
    integer:: nTau, i, patchNum
    real:: tauMin, obsTMinKind4, obsTMaxKind4, patchMinTau
    logical:: callCreateArray
    
    callCreateArray = .false.
    if (.not.associated(P%tau)) then
      callCreateArray = .true.
      if (GetAssociatedLoad(P)) then
        if (GetLoadType(P%loadObject)==APERIODIC_LOADING) then
          patchMinTau = GetLoadingSourceTime(P%loadObject,1)
          callCreateArray = .false.
        end if
      else if (GetSurfaceType(P%surfaceObject)==APERIODIC_SURFACE) then
        patchMinTau = GetSurfaceSourceTime(P%surfaceObject,1)
        callCreateArray = .false.
      end if
      if (callCreateArray) then
        ! Do a retarded time calculation to find tauMin, tauMax and nTau.'
        obsTMinKind4 = observerTimeArray(1)
        obsTMaxKind4 = observerTimeArray(size(observerTimeArray))
        call CreateSourceTimeArray (P, CBList, obsInitVector, obsCBList,&
                                    obsTMinKind4, &
                                    obsTMaxKind4, &
                                    size(observerTimeArray))
      else
!        call SourceTime (obsLocationAtTime, tMin, dt, tempTime, CBList)
        ! We already have our tau range as given from the aperiodic data file.
        P%nTau = nTau
        allocate (P%tau(P%nTau))   
        P%tau(P%nTau) = patchMinTau + (P%nTau-1)*P%dTau
        P%tau = CreateEvenlySpacedArray(patchMinTau, P%tau(P%nTau), P%ntau)
      end if
    end if
    
    P%tauMax = P%tau(P%nTau)
   
  end subroutine GetPatchTimeArrays                                


  function Ufunction(m,n,r,Mr) result(res)
   implicit none
   integer, intent(in) :: m,n
   real,    intent(in) :: r,Mr
   real                :: res
  
   res=1/(r**m*(1.0-Mr)**n)
  end function Ufunction

  function Vfunction(m,n,r,Mdotr,Mr,M2) result(res)
   implicit none
   integer, intent(in) :: m,n
   real,    intent(in) :: r,Mdotr,Mr,M2
   real                :: res

   res=(n*r*Mdotr+(n-m)*c*Mr**2+m*c*Mr-n*c*M2)/(r**(m+1)*(1-Mr)**(n+1))
  end function Vfunction

  function Wfunction(r,Mdotr,Mr,M2) result(res)
   implicit none
   real,    intent(in) :: r,Mdotr,Mr,M2
   real                :: res

   res=r*Mdotr+c*Mr-c*M2
  end function Wfunction

  function Wdotfunction(r,M2dotr,M,Mdot,Mdotr,Mr,M2) result(res)
   implicit none
   real,    intent(in) :: r,M2dotr,Mdotr,Mr,M2
   type(vector), intent(in) :: M,Mdot
   real                :: res

   res=(r**2*M2dotr-3*c*r*M*Mdot+c*(-c*M2+c*Mr**2+r*Mdotr))/r
  end function Wdotfunction
  
  function GetAssociatedLoad(P) result(load)
    implicit none
    type(patch)::P
    logical:: load
    load = .false.
    if (associated(P%loadObject)) load=.true.
  end function GetAssociatedLoad
  
  function GetPatchLoadnKey(P) result(nKey)
    implicit none
    type(patch)::P
    integer::nKey
    nKey = GetLoadnKey(P%loadObject)
  end function GetPatchLoadnKey
  
  function GetPatchSurfacenKey(P) result(nKey)
    implicit none
    type(patch)::P
    integer::nKey
    nKey = GetSurfacenKey(P%surfaceObject)
  end function GetPatchSurfacenKey  
  
  function GetPatchNTau(P) result(nTau)
    implicit none
    type(patch)::P
    integer::nTau
    nTau = P%nTau
  end function GetPatchNTau
  
  subroutine SetPatchNTau(P,nTau)
    implicit none
    type(patch), intent(inout)::P
    integer, intent(in)::nTau
    P%nTau = nTau
  end subroutine SetPatchNTau
  
  function GetPatchGeoNKey(P) result(nTau)
    implicit none
    type(patch)::P
    integer::nTau
    ntau = GetSurfacenKey(P%SurfaceObject)
  end function GetPatchGeoNKey
    
  function GetPatchDTau(P) result(dTau)
    implicit none
    type(patch)::P
    real::dTau
    dTau = P%dTau
  end function GetPatchDTau
  
  subroutine mg2Patch(P, Z, Initializing)
    implicit none
    type(patch):: P
    type(storedData):: Z
    logical:: Initializing
    
    P%obsCBList        => mg_obsCBList
    P%srcVector        =  mg_srcVector
    P%srcCBList        => mg_srcCBList
    P%srcIndex         =  mg_srcIndex
    P%tauLocal         =  mg_tauLocal
    Z%radiationVector  =  mg_radiationVector
    Z%obsVector        =  mg_obsVector
    Z%tLocal           =  mg_tLocal
    Z%radIndex	       =  mg_radIndex
    if (.not. Initializing) then
      if (P%loadSpecies.eq.LOAD_DATATYPE_PRESSURE .or. P%loadSpecies.eq.LOAD_DATATYPE_FLOW) then
        call transferDataParameters(Z%parameterL, mg_flowParametersL)
        call transferDataParameters(Z%parameterC, mg_flowParametersC)
        call transferDataParameters(Z%parameterR, mg_flowParametersR)
      end if
    end if
    
  end subroutine mg2Patch

  subroutine transferDataParameters(storedParameter, globalParameter)
    type(dataParameters), pointer:: storedParameter, globalParameter, temp
    
    temp            => storedParameter
    storedParameter => globalParameter
    globalParameter => temp
  
  end subroutine transferDataParameters  
  
  subroutine patch2MG(P, Z, iTau)
    type(patch):: P
    type(storedData):: Z
    integer:: iTau

    if (P%loadSpecies.eq.LOAD_DATATYPE_PRESSURE .or. P%loadSpecies.eq.LOAD_DATATYPE_FLOW) then
      mg_flowParametersL  => Z%parameterL
      mg_flowParametersC  => Z%parameterC
      mg_flowParametersR  => Z%parameterR
    end if
    if (iTau > 1) then
      mg_obsCBList             => P%obsCBList
      mg_srcVector             =  P%srcVector
      mg_radiationVector       =  Z%radiationVector
      mg_obsVector             =  Z%obsVector
      mg_srcCBList             => P%srcCBList
      mg_srcIndex              =  P%srcIndex
      mg_tauLocal              =  P%tauLocal
      mg_tLocal                =  Z%tLocal
      mg_radIndex	       =  Z%radIndex
    end if

  end subroutine patch2MG
  
  subroutine InitializeStoredData(P)
    implicit none
    type(patch):: P
    integer:: i
    logical:: Initializing = .TRUE.

    allocate(P%storedDataArray(0:mg_NbWall))
    do i=0,mg_NbWall
      call mg2Patch(P, P%StoredDataArray(i), Initializing)
    end do
  end subroutine InitializeStoredData
  
  function GetGeoTitle(P)  result(geoTitle)
    implicit none
    type(patch):: P
    character(len=32):: geoTitle
    
    geoTitle = P%geoTitle
  end function GetGeoTitle

  function GetLoadTitle(P) result(loadTitle)
    implicit none
    type(patch):: P
    character(len=32):: loadTitle
    
    loadTitle = P%loadTitle
  end function GetLoadTitle


  subroutine GetDimensions(P, geoStreamNumber, loadStreamNumber, loadInfo)
    implicit none
    type(patch)::P
    integer,dimension(:)::loadInfo
    integer:: geoStreamNumber, loadStreamNumber
    
    if (getSurfaceType(P%surfaceObject) == APERIODIC_SURFACE) then
    !if (isAperiodic(SURFACEOBJ)) then
      call ReadDimensions (P%surfaceObject, geoStreamNumber)
      call ReadConnectivity (P%surfaceObject, geoStreamNumber)
    end if
    if (GetAssociatedLoad(P) ) then
      if (getLoadType(P%loadObject)==APERIODIC_LOADING) then
        call ReadLoadingDimensions(P%loadObject, loadStreamNumber, loadInfo)
      end if
    end if
      
  end subroutine GetDimensions
  
  function GetPatchNbNodes (P) result (nbNodes)
    type(patch), intent(in)::P
    integer :: nbNodes
    nbNodes = P%nbNodes
  end function GetPatchNbNodes
  
  function GetPatchMinObsTimeInData(P,obsCount) result(minObsTimeInData)
    implicit none
    type(patch):: P
    real(kind=8):: minObsTimeInData
    integer:: obsCount, i
    ! We need this check to ensure we don't just grab a value
    ! that hasn't been changed from its initial huge() value
    do i=1,size(P%obsTimeInRange(obsCount)%tMin)
      minObsTimeInData = P%ObsTimeInRange(obsCount)%tMin(i)
      if (minObsTimeInData.ne.huge(minObsTimeInData)) exit
      minObsTimeInData = -minObsTimeInData  ! this is work around for compact 
      ! patch thickness noise, which should not be computed but it is if there
      ! is no loading noise either.
    end do
  end function GetPatchMinObsTimeInData
  
  function GetPatchMaxObsTimeInData(P,obsCount) result(maxObsTimeInData)
    implicit none
    type(patch):: P
    real(kind=8):: maxObsTimeInData
    integer:: obsCount, i
    ! We need this check to ensure we don't just grab a value
    ! that hasn't been changed from its initial -huge() value    
    do i=size(P%obsTimeInRange(obsCount)%tMax),1,-1
      maxObsTimeInData = P%ObsTimeInRange(obsCount)%tMax(i)
      if (maxObsTimeInData.ne.-huge(maxObsTimeInData)) exit
    end do
    maxObsTimeInData = -maxObsTimeInData  ! this is work around for compact 
      ! patch thickness noise, which should not be computed but it is if there
      ! is no loading noise either.
  end function GetPatchMaxObsTimeInData
  
  function GetPatchMaxSourceTime(P) result(maxSourceTime)
    implicit none
    type(patch):: P
    real::maxSourceTime
    maxSourceTime = P%tauMax
  end function GetPatchMaxSourceTime
  
  function GetPatchTauArray(P) result(tau)
    implicit none
    type(patch):: P
    real, dimension(:), pointer::tau
    tau => P%tau
  end function GetPatchTauArray
  
  function GetPatchITauMin(P) result(iTauMin)
    implicit none
    type(patch):: P
    integer:: iTauMin
    
    iTauMin = P%iTauMin
  end function GetPatchITauMin
  
  function PatchITauInRange(P) result(flag)
    implicit none
    type(patch):: P
    logical:: flag
    
    if( P%iTau <= P%nTau ) then
      flag = .true.
    else
      flag = .false.
    end if  
  end function PatchITauInRange
  
  function PatchLastITau(P) result(flag)
    implicit none
    type(patch):: P
    logical:: flag
    
    if( P%iTau == P%nTau ) then
      flag = .true.
    else
      flag = .false.
    end if  

  end function PatchLastITau
  
  subroutine CompareDataTimeInfo(BBnTau, BBperiod, P)
    implicit none
    integer:: BBnTau
    real:: BBperiod
    type(patch):: P

    if (BBnTau .ne. GetLoadnKey(P%loadObject)) then
      call Error(' Loading data and associated broadband data',&
                 ' must have the same number of timesteps.')
    else if (BBperiod .ne. GetLoadPeriod(P%loadObject)) then
      call Error(' Loading data and associated broadband data',&
                 ' must have the same period.')
    end if
    
  end subroutine CompareDataTimeInfo
    
  subroutine CheckPatchNSect(P, BPMData, isCompact)
    implicit none
    type(patch):: P
    type(BPM), pointer:: BPMData
    type(surface), pointer:: surf
    integer:: nSect, isCompact
    integer, dimension(2):: dims

    if (IsSurfaceCompact(P%surfaceObject)) then
      isCompact = isCompact + 1
      surf => P%surfaceObject
      call getSurfaceDimensions(surf,dims)
      nSect = dims(1)*dims(2)
      if (any(BPMData%inFile) .and. BPMData%nSect.ne.nSect) then
         call Error ('The number of segments specified in the BPM file is not',&
                    'the same as the number of points on the compact patch.')
      end if
    end if
  end subroutine CheckPatchNSect
  
  subroutine ComputePatchThrust(PeggData, P, CBList, loadTimeType, nTau, keyIndexSurf, keyIndexLoad)
    implicit none
    type(PEGG):: PeggData
    type(patch):: P
    type(CBStructure), pointer:: CBList
    type(baseChange):: BladeFrame, HubFrame
    type(SingleSurface), pointer:: currentSurface=>NULL()    
    type(singleTimeLoad), pointer:: loads=>NULL()
    type(vector):: temp, hubAxis
    integer:: nTau, i, j, k, nbTau, keyIndexSurf, keyIndexLoad, loadTimeType, BBkey
    integer:: loadSpecies, nbData, index
    
    if (associated(P%loadObject)) then
      loadSpecies = GetLoadSpecies(P%loadObject)
    end if
    loads => GetLoading(P%loadObject, 1)
    select case (loadSpecies)
      case (LOAD_DATATYPE_PRESSURE)
        nbData = size(loads%pressVals)
      case (LOAD_DATATYPE_LOADING) 
        nbData = size(loads%loadVectors)
      case (LOAD_DATATYPE_FLOW)
        call Error ("Thrust cannot be calculated from flow data.")
    end select
    if (loadTimeType.eq.LOAD_TIMETYPE_PERIODIC) nbTau = nTau-1
    do j=1,nbTau
      index =  keyIndexLoad+P%iTauMin-1+(j-1)
      currentSurface => GetSurface (P%surfaceObject, keyIndexSurf+index)
      loads => GetLoading(P%loadObject, keyIndexLoad+index)
      
      call CalculateResultantChangeOfBase(CBList, GetSize(CBList), P%tau(j))
      BladeFrame = resultantBase
      
      call CalculateResultantChangeOfBase(PeggData%CBList, GetSize(PeggData%CBList), P%tau(j))
      HubFrame = resultantBase
      if (PeggData%computeTerm(PEGG_HUBAXIS).eq.'COMPUTE') then
        hubAxis = GetHubAxis(PeggData%CBList)
      else
        BBkey = 1
        if (PeggData%timeType.eq.LOAD_TIMETYPE_PERIODIC) then
          BBkey = getBBoffset(size(PeggData%tau), PeggData%timeOffset, PeggData%period, keyIndexLoad+index)
        end if
        hubAxis = PeggData%HubAxis(GetPeggTermKey(PeggData, PEGG_HUBAXIS, BBkey))
      end if
      hubAxis = hubAxis/vectorAbsolute(hubAxis)
      hubAxis = matrixVectorMultiply(HubFrame%rotation, hubAxis)
      do i=1,nbData
        if (loadSpecies.eq.LOAD_DATATYPE_PRESSURE) then
          temp = -loads%pressVals(i)*currentSurface%areas(i)*currentSurface%normals(i)
        else
          ! At this point we can be sure that if we aren't working with pressure data
          ! we're working with loading data.
          if(P%refFrame == LOAD_REFFRAME_GROUND) then
            temp = -loads%loadVectors(i)
          else
            ! We have the hub axis the in ground frame, we need the load vector in the ground frame 
            ! V_h = g_T_b*V_b
            temp = -matrixVectorMultiply(BladeFrame%rotation,loads%loadVectors(i))
          end if
        end if
        PeggData%ShaftAxisThrust = PeggData%ShaftAxisThrust + vectorDotProduct(temp, hubAxis)
      end do
    end do
    
  end subroutine ComputePatchThrust
  
  subroutine ComputePatchArea(PeggData, P)
    implicit none
    type(PEGG):: PeggData
    type(Patch), intent(inout)::P
    type(SingleSurface), pointer:: currentSurface=>NULL()
    real:: totalBladeArea
    
    if (.not. isSurfaceCompact(P%surfaceObject)) then
      currentSurface => GetSurface (P%surfaceObject, 1)
      ! This is only an approximation since we can't assume anything
      ! about the shape of the blade.  A_b = N_b*c*R ~ (1/2)*A_surf
      PeggData%totalBladeArea = PeggData%totalBladeArea + sum(currentSurface%areas)/2.0
    end if
  end subroutine ComputePatchArea
  
  subroutine ComputePatchSpan(PeggData, P, CBList, hubPosition)
    implicit none
    type(PEGG):: PeggData
    type(Patch), intent(inout):: P
    type(CBStructure), pointer:: CBList
    type(vector):: HubPosition
    type(SingleSurface), pointer:: currentSurface=>NULL()    
    integer:: i
    
    currentSurface => GetSurface (P%surfaceObject, 1)
    call CalculateResultantChangeOfBase(CBList, GetSize(CBList), 0.0)
    do i=1,size(currentSurface%vectors)
      PeggData%BladeRadius = max(PeggData%BladeRadius, vectorAbsolute(resultantBase%position &
          + currentSurface%vectors(i)-HubPosition))
    end do
  end subroutine ComputePatchSpan  
  
  subroutine ComputePatchSectLength(BPMData, P)
    implicit none
    type(BPM):: BPMdata
    type(Patch):: P
    type(singleSurface), pointer:: currentSurface
    currentSurface => GetSurface(P%surfaceObject, 1)
    BPMData%sectLength = currentSurface%areas
       
  end subroutine ComputePatchSectLength     
  
end module patchObject
