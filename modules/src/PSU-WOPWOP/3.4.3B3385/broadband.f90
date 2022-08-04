module broadbandModule
! Written by Ben Goldman 
! Faculty advisor Dr. Kenneth S. Brentner
   
!**********************************************************************
!MODULE Broadband
!  This module contains all the Broadband subroutines.  
!  Currently only two semi-empirical methods are available:
!  The first by R.J. Pegg and the second by Brooks,Pope and Marcolini.
!**********************************************************************

  use constantsModule
  use mathModule
  use COBObject
  use IOModule
  use FrequencyDomainObject, only: GetThirdOctFreqArray
  use FrequencyAnalysisModule
  use strings
  implicit none

  !**************************************************************************
  !TYPE: Broadband
  !  The type Broadband contains any and all types associated with predicting
  !  broadband noise.  Currently only two methods are available.
  !***************************************************************************  
  type Broadband
    private
    type (PEGG), pointer:: PeggData=>null()
    type (BPM), pointer :: BPMdata=>null()
  end type Broadband

  !**************************************************************************
  !TYPE: PEGG
  !  The type PEGG contains all the data necessary to predict 
  !  broadband noise based upon the approach by R.J. Pegg
  !***************************************************************************  
  type PEGG
    integer:: nTau, nbRevs, timeType, iTau, iTauMin
    integer, dimension(:), pointer:: nbContributions=>null()
    real(kind=8), dimension(:,:), pointer:: Intensity=>null()
    real:: dTau, period, timeOffset, deltadB !maxSrcTime
    logical:: spreadBands
    real, dimension(:), pointer:: tau=>null(), obsTimeArray=>null()
    type(CBStructure), pointer::CBlist=>null()
    ! InFile just tells us whether the data is in the file, not whether
    ! the user wants it to be used for computation or not.    
    logical, dimension(PEGG_INPUT_DATA):: inFile
    character(len=32), dimension(PEGG_INPUT_DATA):: computeTerm
    ! This information is needed for Pegg's Method
    real:: totalBladeArea, BladeRadius, RotSpeed, R
    ! These terms may change over time
    type(vector), dimension(:), pointer:: hubAxis=>null()
    real, dimension(:), pointer:: ClBar=>null()
    real:: TipSpeed, theta1, shaftAxisThrust	
  end type PEGG
  
  !**************************************************************************
  !TYPE: BPM
  !  The type BPM contains all the data necessary to predict 
  !  broadband noise based upon the approach by Brooks, Pope and Marcolini
  !***************************************************************************
  type BPM
    integer:: nTau, nbRevs, timeType, nSect, iTau
    integer, dimension(:), pointer:: nbContributions=>null()
    real(kind=8), dimension(:,:,:), pointer:: Intensity=>null()
    real:: dTau, period, timeOffset !, maxSrcTime
    real, dimension(:), pointer:: tau=>null(), obsTimeArray=>null(), DopplerShift=>null()
    logical:: uniformBlade
    ! InFile just tells us whether the data is in the file, not whether
    ! the user wants it to be used for computation or not.
    logical, dimension(BPM_INPUT_DATA):: inFile
    character(len=32), dimension(BPM_INPUT_DATA):: computeTerm
    ! These terms are constant over time
    real, dimension(:), pointer:: SectChord=>null(), SectLength=>null(), TEthickness=>null(),TEflowAngle=>null(), R=>null()
    ! These terms may change over time
    real, dimension(:), pointer:: Theta=>null(), Phi=>null(), TipLCS=>null()
    real, dimension(:,:), pointer:: U=>null(), SectAOA=>null()
    integer:: BLtrip
    logical:: LBLVSnoise, TBLTEnoise, BluntNoise, BladeTipNoise, RoundBladeTip
  end type BPM
  
  type PeggNamelist
    character(len=32):: totalBladeAreaFlag, bladeRadiusFlag, rotSpeedFlag, CLBarFlag, HubAxisFlag, totalThrustFlag
    real:: totalBladeArea, bladeRadius, rotSpeed, ClBar, TotalThrust, deltadB
    logical:: spreadBands
    type(vector):: hubAxis
  end type PeggNamelist
  
  type BPMNamelist
    integer:: nSect, uniformBlade, BLTrip
    character(len=32):: sectChordFlag, sectLengthFlag, TEthicknessFlag, TEflowAngleFlag, TipLCSFlag, SectAOAFlag, UFlag
    real, dimension(:), allocatable:: sectChord, sectLength, TEThickness, TEflowAngle, SectAOA, U
    real:: TipLCS
    logical:: LBLVSnoise, TBLTEnoise, bluntNoise, bladeTipNoise, roundBladeTip
  end type BPMNamelist 
  
  private
  public:: Broadband, PEGG, BPM, PeggNamelist, BPMNamelist
  public:: readBBHeader,  readPeggNamelist, readBPMNamelist, readPeggData, readBPMdata
  public:: CreatePeggArrays, CreateBPMArrays, ComputePeggBroadbandTerms, ComputeBPMBroadbandTerms
  public:: ComputePeggIntensity, ComputeBPMintensity, ComputePeggTerm, getHubAxis, GetOmega
  public:: getBBoffset, GetPeggTermKey, getBPMtermKey, getPeggTimeType, getBPMTimeType, GetPeggData, GetBPMData
  public:: DestroyBPMNamelistData, destroyBBdata, destroyPeggData, destroyBPMData

  ! This should be plenty for BPM's method.
  integer, parameter:: mg_maxSegs = 200

contains

  subroutine ReadPeggNamelist(unitNumber, PeggNamelistData, NoiseFile)
    implicit none
    integer:: unitNumber
    type(PeggNamelist):: PeggNamelistdata
    character(len=4096):: NoiseFile, PeggNoiseFile
    character(len=32):: totalBladeAreaFlag, bladeRadiusFlag, rotSpeedFlag, CLBarFlag, hubAxisFlag, totalThrustFlag
    real:: totalBladeArea, bladeRadius, rotSpeed, CLBar, totalThrust, temp, deltadB
    type(vector):: hubAxis
    logical:: spreadBands, keepEnergyConst
    
    namelist / PeggIn /  PeggNoiseFile, totalBladeAreaFlag, bladeRadiusFlag, rotSpeedFlag, CLBarFlag, &
        hubAxisFlag, totalThrustFlag, totalBladeArea, bladeRadius, rotSpeed, CLBar, hubAxis, totalThrust,&
        deltadB, spreadBands, keepEnergyConst, EPNLPrime
                         
    temp = -huge(temp)
    ! Default values
    PeggNoiseFile      = ''
    totalBladeAreaFlag = ''
    totalBladeArea     = temp
    bladeRadiusFlag    = ''
    bladeRadius        = temp
    rotSpeedFlag       = ''
    rotSpeed           = temp
    ClBarFlag          = ''
    CLBar              = temp
    hubAxisFlag        = ''
    hubAxis            = vectorSetCoordinates(temp, temp, temp)
    totalThrustFlag    = ''
    totalThrust        = temp
    deltadB	       = 0.0
    spreadBands	       = .false.
    keepEnergyConst    = .false.
    
    read(unitNumber, nml=PeggIn)
   
    call strToUpper (totalBladeAreaFlag, len_trim(totalBladeAreaFlag))
    call strToUpper (bladeRadiusFlag,      len_trim(bladeRadiusFlag))
    call strToUpper (rotSpeedFlag,       len_trim(rotSpeedFlag))
    call strToUpper (ClBarFlag,          len_trim(ClBarFlag))
    call strToUpper (hubAxisFlag,        len_trim(hubAxisFlag))
    call strToUpper (totalThrustFlag,    len_trim(totalThrustFlag))
    
    if (totalBladeAreaFlag.eq.'USERVALUE' .and. totalBladeArea.eq.temp) then
      call Error('totalBladeAreaFlag is set to `USERVALUE` but no data is associated with totalBladeArea')
    end if
    if (bladeRadiusFlag.eq.'USERVALUE' .and. bladeRadius.eq.temp) then
      call Error('bladeRadiusFlag is set to `USERVALUE` but no data is associated with bladeRadius')                 
    end if
    if (RotSpeedFlag.eq.'USERVALUE' .and. RotSpeed.eq.temp) then
      call Error('RotSpeedFlag is set to `USERVALUE` but no data is associated with RotSpeed')
    end if
    if (CLBarFlag.eq.'USERVALUE' .and. CLBar.eq.temp) then
      call Error('CLBarFlag is set to `USERVALUE` but no data is associated with CLBar')
    end if
    if (HubAxisFlag.eq.'USERVALUE' .and. any(HubAxis%A.eq.temp)) then
      call Error('HubAxisFlag is set to `USERVALUE` but no data is associated with HubAxis')
    end if
    if (TotalThrustFlag.eq.'USERVALUE' .and. TotalThrust.eq.temp) then
      call Error('TotalThrustFlag is set to `USERVALUE` but no data is associated with TotalThrust')
    end if

    NoiseFile                           = trim(PeggNoiseFile)
    PeggNamelistdata%totalBladeAreaFlag = totalBladeAreaFlag
    PeggNamelistdata%totalBladeArea     = totalBladeArea
    PeggNamelistdata%bladeRadiusFlag    = bladeRadiusFlag
    PeggNamelistdata%bladeRadius        = bladeRadius
    PeggNamelistdata%rotSpeedFlag       = rotSpeedFlag
    PeggNamelistdata%rotSpeed           = rotSpeed
    PeggNamelistdata%CLbarFlag          = CLbarFlag
    PeggNamelistdata%ClBar              = ClBar
    PeggNamelistdata%hubAxisFlag        = hubAxisFlag
    PeggNamelistdata%hubAxis            = hubAxis
    PeggNamelistdata%TotalThrustFlag    = TotalThrustFlag
    PeggNamelistdata%TotalThrust        = TotalThrust
    PeggNamelistdata%spreadBands	= spreadBands
    
    PeggNamelistdata%deltadB = 0.
    if (spreadBands) then
      if (keepEnergyConst) then
	! This is the shift in sound level required
        ! to apply Pegg's spectral shape over all 3rd-octave
        ! bands while keeping the same amount of energy in
        ! the spectrum.
        PeggNamelistdata%deltadB	= -4.5458
      else
        PeggNamelistData%deltadB	= deltadB
      end if
    end if

  end subroutine ReadPeggNamelist

  subroutine ReadBPMNamelist(unitNumber, BPMNamelistData, NoiseFile)
    implicit none
    integer:: unitNumber
    type(BPMNamelist):: BPMNamelistdata
    character(len=4096):: NoiseFile, BPMNoiseFile  
    character(len=32):: sectChordFlag, SectLengthFlag, TEthicknessFlag, TEflowAngleFlag, TipLCSFlag, SectAOAFlag, UFlag
    integer:: nSect, uniformBlade, BLTrip, nSectTemp
    real, dimension(mg_maxSegs):: sectChord, SectLength, TEthickness, TEflowAngle, SectAOA, U
    real:: TipLCS, temp
    logical:: LBLVSnoise, TBLTEnoise, bluntNoise, bladeTipNoise, roundBladeTip
  
    namelist / BPMIn /  BPMNoiseFile, nSect, uniformBlade, BLTrip, &
                        sectChordFlag, SectLengthFlag, TEthicknessFlag, TEflowAngleFlag, TipLCSFlag, SectAOAFlag, UFlag, &
                        sectChord, SectLength, TEthickness, TEflowAngle, TipLCS, SectAOA, U, &
                        LBLVSnoise, TBLTEnoise, bluntNoise, bladeTipNoise, roundBladeTip
    
    temp = -huge(temp)
    ! Default values
    NoiseFile        = ''
    BPMNoiseFile     = ''
    nSect            = 0
    uniformBlade     = 0
    BLTrip           = 0
    sectChordFlag    = ''
    sectChord        = temp
    sectLengthFlag   = ''
    sectLength       = temp
    TEthicknessFlag  = ''
    TEthickness      = temp
    TEflowAngleFlag  = ''
    TEflowAngle      = temp
    TipLCSFlag       = ''
    TipLCS           = temp
    SectAOAFlag      = ''
    SectAOA          = temp
    UFlag            = ''
    U                = temp
    LBLVSnoise       = .false.
    TBLTEnoise       = .false.
    bluntNoise       = .false.
    bladeTipNoise    = .false.
    roundBladeTip    = .false.
    
    read(unitNumber, nml=BPMIn)
   
    call strToUpper (sectChordFlag,   len_trim(sectChordFlag))
    call strToUpper (sectLengthFlag,  len_trim(sectLengthFlag))
    call strToUpper (TEthicknessFlag, len_trim(TEthicknessFlag))
    call strToUpper (TEflowAngleFlag, len_trim(TEflowAngleFlag))
    call strToUpper (TipLCSFlag,      len_trim(TipLCSFlag))
    call strToUpper (SectAOAFlag,     len_trim(SectAOAFlag))
    call strToUpper (UFlag,           len_trim(UFlag))
    
    if (uniformBlade.ne.0 .and. uniformBlade.ne.1) then
      call Error('UniformBlade mu1st have a value of 0 or 1.')
    end if
    
    if (BLTrip.ne.0 .and. BLTrip.ne.1 .and. BLTrip.ne.2) then
      call Error('BLTrip mu1st have a value of 0, 1 or 2.')
    end if
    
    if (uniformBlade.eq.1) then
      nSectTemp = 1
    else
      nSectTemp = nSect
    end if
    
    if (SectChordFlag.eq.'USERVALUE' .and. any(SectChord(1:nSectTemp).eq.temp)) then
      call Error('SectChordFlag is set to `USERVALUE` but no data is associated with SectChord')
    end if
    if (SectLengthFlag.eq.'USERVALUE' .and. any(SectLength(1:nSectTemp).eq.temp)) then
      call Error('SectLengthFlag is set to `USERVALUE` but no data is associated with SectLength')
    end if
    if (TEThicknessFlag.eq.'USERVALUE' .and. any(TEThickness(1:nSectTemp).eq.temp)) then
      call Error('TEThicknessFlag is set to `USERVALUE` but no data is associated with TEThickness')
    end if
    if (TEflowAngleFlag.eq.'USERVALUE' .and. any(TEflowAngle(1:nSectTemp).eq.temp)) then
      call Error('TEflowAngleFlag is set to `USERVALUE` but no data is associated with TEflowAngle')
    end if
    if (TipLCSFlag.eq.'USERVALUE' .and. TipLCS.eq.temp) then
      call Error('TipLCSFlag is set to `USERVALUE` but no data is associated with TipLCS')
    end if
    if (SectAOAFlag.eq.'USERVALUE' .and. any(SectAOA(1:nSectTemp).eq.temp)) then
      call Error('SectAOAFlag is set to `USERVALUE` but no data is associated with SectAOA')
    end if
    if (UFlag.eq.'USERVALUE' .and. any(U(1:nSectTemp).eq.temp)) then
      call Error('UFlag is set to `USERVALUE` but no data is associated with U')
    end if
    
    NoiseFile                       = trim(BPMNoiseFile)
    BPMNamelistData%nSect           = nSect
    BPMNamelistData%uniformBlade    = uniformBlade
    BPMNamelistData%BLTrip          = BLTrip
    
    BPMNamelistData%sectChordFlag   = sectChordFlag
    if (sectChordFlag.eq.'USERVALUE') then
      allocate(BPMNamelistData%sectChord(nSect))   
      if (uniformBlade == 1) then
        BPMNamelistData%sectChord = sectChord(1)
      else
        BPMNamelistData%sectChord = sectChord
      end if
    end if
    BPMNamelistData%SectLengthFlag  = SectLengthFlag
    if (SectLengthFlag.eq.'USERVALUE') then
      allocate(BPMNamelistData%SectLength(nSect))
      if (uniformBlade == 1) then
        BPMNamelistData%SectLength  = SectLength(1)
      else
        BPMNamelistData%SectLength  = SectLength
      end if
    end if
    BPMNamelistData%TEThicknessFlag = TEThicknessFlag
    if (TEThicknessFlag.eq.'USERVALUE') then
      allocate(BPMNamelistData%TEThickness(nSect))
      if (uniformBlade == 1) then
        BPMNamelistData%TEThickness = TEThickness(1)
      else
        BPMNamelistData%TEThickness = TEThickness
      end if
    end if
    BPMNamelistData%TEflowAngleFlag = TEflowAngleFlag
    if (TEflowAngleFlag.eq.'USERVALUE') then
      allocate(BPMNamelistData%TEflowAngle(nSect))
      if (uniformBlade == 1) then
        BPMNamelistData%TEflowAngle = TEflowAngle(1)
      else
        BPMNamelistData%TEflowAngle = TEflowAngle
      end if
    end if
    BPMNamelistData%SectAOAFlag = SectAOAFlag
    if (SectAOAFlag.eq.'USERVALUE') then
      allocate(BPMNamelistData%SectAOA(nSect))
      if (uniformBlade == 1) then
        BPMNamelistData%SectAOA = SectAOA(1)
      else
        BPMNamelistData%SectAOA = SectAOA
      end if
    end if
    BPMNamelistData%TipLCSFlag  = TipLCSFlag
    if (TipLCSFlag.eq.'USERVALUE') then
      BPMNamelistData%TipLCS    = TipLCS
    end if
    BPMNamelistData%UFlag = UFlag
    if (UFlag.eq.'USERVALUE') then
      allocate(BPMNamelistData%U(nSect))
      if (uniformBlade == 1) then
        BPMNamelistData%U = U(1)
      else
        BPMNamelistData%U = U
      end if
    end if
    
    BPMNamelistData%LBLVSnoise      = LBLVSnoise
    BPMNamelistData%TBLTEnoise      = TBLTEnoise
    BPMNamelistData%bluntNoise      = bluntNoise
    BPMNamelistData%bladeTipNoise   = bladeTipNoise
    BPMNamelistData%roundBladeTip   = roundBladeTip
  end subroutine ReadBPMNamelist
  

  subroutine readBBHeader(NoiseFlag, NoiseFile, UnitNumber)
    implicit none
    logical:: NoiseFlag
    character(len=4096):: NoiseFile
    integer:: UnitNumber, stat, magicNum

    if (NoiseFlag .and. len_trim(NoiseFile)>0) then
      UnitNumber = GetStreamNumber()
     ! call OpenBinaryFile (UnitNumber, trim(globalFolderName)//'./'//NoiseFile, .false., stat)
     call OpenBinaryFile (UnitNumber, trim(globalFolderName)//NoiseFile, .false., stat)
      if (stat/=0) then
        call Error('Could not open the file '//trim(NoiseFile))
        stop
      end if
      call ReadBinaryInteger(UnitNumber,magicNum,stat)
      !determine endian correctness
      if (debuglevel>10) write (*,*) "Magic number is read in as ", magicNum
      if (magicNum == 704643072) then
        !file needs endian converted
        call CloseBinaryFile(UnitNumber)
       ! call OpenBinaryFile (UnitNumber, trim(globalFolderName)//'./'//NoiseFile, .true., stat)
        call OpenBinaryFile (UnitNumber, trim(globalFolderName)//NoiseFile, .true., stat)
        call ReadBinaryInteger(UnitNumber,magicNum,stat)
      end if
      if (magicNum /= 42) then
        call Error('Could not recognize file format in'//trim(NoiseFile))
        stop
      end if
    else
      UnitNumber=0
    end if

  end subroutine readBBHeader
  
  subroutine CreatePeggArrays(BBData, UnitNumber, NamelistData)
    implicit none
    type(Broadband), pointer:: BBData
    integer:: unitNumber
    type(PeggNamelist):: NamelistData
    
    type(PEGG), pointer:: PeggData
    integer:: totalBladeAreaFlag, bladeRadiusFlag, rotSpeedFlag, ClBarFlag, hubAxisFlag, totalThrustFlag
    real:: tempTotalBladeArea, tempbladeRadius, tempRotSpeed, tempClBar, tempTotalThrust
    type(vector):: tempHubAxis
    integer:: nbAllocate, nTau, i, stat

    allocate(BBdata%PeggData)
    PeggData => BBdata%PeggData
    
    PeggData%computeTerm     = ''
    PeggData%inFile          = .false.
    PeggData%nTau            = 0
    PeggData%timeType        = LOAD_TIMETYPE_CONSTANT  
    
    PeggData%computeTerm(PEGG_TOTALBLADEAREA) = NamelistData%totalBladeAreaFlag
    PeggData%computeTerm(PEGG_BLADERADIUS)    = NamelistData%bladeRadiusFlag
    PeggData%computeTerm(PEGG_ROTSPEED)       = NamelistData%RotSpeedFlag
    PeggData%computeTerm(PEGG_CLBAR)          = NamelistData%CLBarFlag
    PeggData%computeTerm(PEGG_HUBAXIS)        = NamelistData%HubAxisFlag
    PeggData%computeTerm(PEGG_TOTALTHRUST)    = NamelistData%TotalThrustFlag
    PeggData%deltadB			      = NamelistData%deltadB
    PeggData%spreadBands		      = NamelistData%spreadBands
   
    PeggData%ShaftAxisThrust = 0.
    PeggData%theta1          = 0.
    PeggData%TipSpeed        = 0.
    
    
    if (UnitNumber.gt.0) then
      call ReadBinaryInteger(unitNumber,TotalBladeAreaFlag,stat)
      call ReadBinaryInteger(unitNumber,bladeRadiusFlag,stat)
      call ReadBinaryInteger(unitNumber,rotSpeedFlag,stat)
      call ReadBinaryInteger(unitNumber,ClBarFlag,stat)
      call ReadBinaryInteger(unitNumber,hubAxisFlag,stat)
      call ReadBinaryInteger(unitNumber,totalThrustFlag,stat)
      call ReadBinaryInteger(unitNumber,PeggData%timeType,stat)
      select case (PeggData%timeType)
        case (LOAD_TIMETYPE_CONSTANT)
          nTau = 1
        case (LOAD_TIMETYPE_PERIODIC)
          call ReadBinaryInteger(unitNumber,PeggData%nbRevs,stat)
          call ReadBinaryReal(unitNumber,PeggData%period,stat)
          call ReadBinaryInteger(unitNumber,PeggData%nTau,stat)
          nTau = PeggData%nTau
        case (LOAD_TIMETYPE_APERIODIC)
          call ReadBinaryInteger(unitNumber,PeggData%nTau,stat)
          nTau = 1
        case default
          call Error('Unknown timeType used in Pegg data file.')  
      end select
    else 
      TotalBladeAreaFlag = 0
      bladeRadiusFlag = 0
      rotSpeedFlag = 0
      ClBarFlag = 0
      hubAxisFlag = 0
      totalThrustFlag = 0
    end if
    allocate(PeggData%tau(PeggData%nTau))    
      
    ! Give the 'computeTerm' variable the proper value if there is nothing present in either
    ! the namelist or data file
    call SetComputableTerm(PeggData%computeTerm, PeggData%InFile, PEGG_TOTALBLADEAREA, TotalBladeAreaFlag, nTau, nbAllocate)
    if (PeggData%computeTerm(PEGG_TOTALBLADEAREA).eq.'USERVALUE') PeggData%TotalBladeArea = NamelistData%TotalBladeArea

    call SetComputableTerm(PeggData%computeTerm, PeggData%InFile, PEGG_BLADERADIUS, BladeRadiusFlag, nTau, nbAllocate)
    if (PeggData%computeTerm(PEGG_BLADERADIUS).eq.'USERVALUE') PeggData%BladeRadius = NamelistData%BladeRadius
    
    call SetComputableTerm(PeggData%computeTerm, PeggData%InFile, PEGG_ROTSPEED, RotSpeedFlag, nTau, nbAllocate)
    if (PeggData%computeTerm(PEGG_ROTSPEED).eq.'USERVALUE') PeggData%RotSpeed = NamelistData%RotSpeed
    
    call SetComputableTerm(PeggData%computeTerm, PeggData%InFile, PEGG_CLBAR, CLBarFlag, nTau, nbAllocate)
    allocate(PeggData%CLBar(nbAllocate))
    if (PeggData%computeTerm(PEGG_CLBAR).eq.'USERVALUE') PeggData%CLBar = NamelistData%CLBar

    call SetComputableTerm(PeggData%computeTerm, PeggData%InFile, PEGG_HUBAXIS, HubAxisFlag, nTau, nbAllocate)
    allocate(PeggData%HubAxis(nbAllocate))
    if (PeggData%computeTerm(PEGG_HUBAXIS).eq.'USERVALUE') PeggData%HubAxis = NamelistData%HubAxis
    
    call SetComputableTerm(PeggData%computeTerm, PeggData%InFile, PEGG_TOTALTHRUST, TotalThrustFlag, nTau, nbAllocate)
    if (PeggData%computeTerm(PEGG_TOTALTHRUST).eq.'USERVALUE') PeggData%ShaftAxisThrust = NamelistData%TotalThrust
      
    if (UnitNumber.gt.0) then
      if (TotalBladeAreaFlag.eq.1) then
        call ReadBinaryReal(unitNumber,tempTotalBladeArea,stat)
        if (PeggData%computeTerm(PEGG_TOTALBLADEAREA).eq.'FILEVALUE') PeggData%totalBladeArea = tempTotalBladeArea
      end if
      if (bladeRadiusFlag.eq.1) then
        call ReadBinaryReal(unitNumber,tempbladeRadius,stat)
        if (PeggData%computeTerm(PEGG_bladeRadius).eq.'FILEVALUE') PeggData%bladeRadius = tempbladeRadius
      end if
      if (rotSpeedFlag.eq.1) then
        call ReadBinaryReal(unitNumber,tempRotSpeed,stat)
        if (PeggData%computeTerm(PEGG_ROTSPEED).eq.'FILEVALUE') PeggData%rotSpeed = temprotSpeed
      end if
      if (PeggData%timeType.ne.LOAD_TIMETYPE_APERIODIC) then
        call ReadPeggData(Peggdata, unitNumber, nTau)
      end if
      ! For periodic data we want to get the average over the period,
      ! but since PSU-WOPWOP periods repeat the first step we need to ignore it here.
      ! tempTotalThrust is the repeated step, so subtract it off and then divide by one
      ! less than the number of steps in file.  If the thrust is being read in
      ! then it is assumed to be along the hub axis.  If not, we'll calculate it
      ! and dot it with the hub axis.
      if (PeggData%timeType.eq.LOAD_TIMETYPE_PERIODIC) then
        PeggData%ShaftAxisThrust = (PeggData%shaftAxisThrust - tempTotalThrust)/(1.0*nTau-1.)
      end if
    end if
        
  end subroutine CreatePeggArrays
  
  subroutine CreateBPMArrays(BBdata, unitNumber, NamelistData)
    implicit none
    type(broadband), pointer:: BBdata
    integer:: unitNumber
    type(BPMNamelist):: NamelistData
    
    type(BPM), pointer:: BPMData=>null()
    integer:: sectChordFlag, sectLengthFlag, TEthicknessFlag, TEflowAngleFlag, tipLCSFlag, SectAOAFlag, Uflag, uniformBlade
    real:: tempSectChord, tempsectLength, tempTEthickness, tempTEflowAngle, tempTipLCS, tempSectAOA, tempU
    integer:: nbAllocate, i, j, stat, nTau, nSect

    allocate(BBdata%BPMData)
    BPMData => BBdata%BPMData
    
    BPMData%computeTerm     = ''
    BPMData%inFile          = .false.
    BPMData%uniformBlade    = .false.
    BPMData%nTau            = 0
    BPMData%timeType        = LOAD_TIMETYPE_CONSTANT
    
    BPMData%nSect           = NamelistData%nSect
    BPMdata%BLtrip          = NamelistData%BLTrip
    BPMData%LBLVSnoise      = NamelistData%LBLVSnoise
    BPMData%TBLTEnoise      = NamelistData%TBLTEnoise
    BPMData%bluntNoise      = NamelistData%bluntNoise
    BPMData%bladeTipNoise   = NamelistData%bladeTipNoise
    BPMData%roundBladeTip   = NamelistData%roundBladeTip
    
    BPMData%computeTerm(BPM_SECTCHORD)    = NamelistData%sectChordFlag
    BPMData%computeTerm(BPM_SECTLENGTH)   = NamelistData%sectLengthFlag
    BPMData%computeTerm(BPM_TETHICKNESS)  = NamelistData%TEthicknessFlag
    BPMData%computeTerm(BPM_TEFLOWANGLE)  = NamelistData%TEflowAngleFlag
    BPMData%computeTerm(BPM_TIPLCS)       = NamelistData%TipLCSFlag
    BPMData%computeTerm(BPM_SECTAOA)      = NamelistData%SectAOAFlag
    BPMData%computeTerm(BPM_U)            = NamelistData%UFlag    

   if (UnitNumber.gt.0) then
      call ReadBinaryInteger(unitNumber,BPMData%nSect,stat)
      call ReadBinaryInteger(unitNumber,uniformBlade,stat)
      call ReadBinaryInteger(unitNumber,sectChordFlag,stat)
      call ReadBinaryInteger(unitNumber,sectLengthFlag,stat)
      call ReadBinaryInteger(unitNumber,TEthicknessFlag,stat)
      call ReadBinaryInteger(unitNumber,TEflowAngleFlag,stat)
      call ReadBinaryInteger(unitNumber,TipLCSFlag,stat)
      call ReadBinaryInteger(unitNumber,SectAOAFlag,stat)
      call ReadBinaryInteger(unitNumber,UFlag,stat)
      call ReadBinaryInteger(unitNumber,BPMData%timeType,stat)
      select case (BPMData%timeType)
        case (LOAD_TIMETYPE_CONSTANT)
          nTau = 1
        case (LOAD_TIMETYPE_PERIODIC)
          call ReadBinaryInteger(unitNumber,BPMData%nbRevs,stat)
          call ReadBinaryReal(unitNumber,BPMData%period,stat)
          call ReadBinaryInteger(unitNumber,BPMData%nTau,stat)
          nTau = BPMData%nTau
        case (LOAD_TIMETYPE_APERIODIC)
          call ReadBinaryInteger(unitNumber,BPMData%nTau,stat)
          nTau = 1
        case default
          call Error('Unknown timeType used in BPM data file.')  
      end select
      if (BPMData%nSect.ne.NamelistData%nSect) then
        call Error(' The BPM data file must have the same number of',&
                   ' blade segments as the BPM namelist.')
      end if
    else
      sectChordFlag = 0
      sectLengthFlag = 0
      TEthicknessFlag = 0
      TEflowAngleFlag = 0
      TipLCSFlag = 0
      SectAOAFlag = 0
      UFlag = 0
      uniformBlade = NamelistData%uniformBlade
    end if
    allocate(BPMData%tau(BPMData%nTau))    
    allocate(BPMData%DopplerShift(BPMData%nSect))
    
    if (uniformBlade.eq.1) BPMData%uniformBlade = .true.

    call SetNonComputableTerm(BPMData%computeTerm, BPMData%InFile, BPM_SECTCHORD, SectAOAFlag, nTau, nbAllocate)
    allocate(BPMdata%SectChord(BPMData%nSect));
    select case(BPMData%computeTerm(BPM_SECTCHORD))
      case ('DEFAULT')
        BPMdata%SectChord = 1.0
      case ('USERVALUE')  
        BPMdata%SectChord = NamelistData%SectChord
    end select
    
    call SetComputableTerm(BPMData%computeTerm, BPMData%InFile, BPM_SECTLENGTH, sectLengthFlag, nTau, nbAllocate)
    allocate(BPMdata%sectLength(BPMData%nSect))  
    if (BPMData%computeTerm(BPM_SECTLENGTH).eq.'USERVALUE') BPMdata%sectLength = NamelistData%sectLength
    
    call SetNonComputableTerm(BPMData%computeTerm, BPMData%InFile, BPM_TETHICKNESS, TEthicknessFlag, nTau, nbAllocate)
    allocate(BPMdata%TEthickness(BPMData%nSect))
    select case(BPMData%computeTerm(BPM_TETHICKNESS))
      case ('DEFAULT')
        BPMdata%TEThickness = 0.0005
      case ('USERVALUE')  
        BPMdata%TEThickness = NamelistData%TEThickness
    end select 
    
    call SetNonComputableTerm(BPMData%computeTerm, BPMData%InFile, BPM_TEFLOWANGLE, TEflowAngleFlag, nTau, nbAllocate)
    allocate(BPMdata%TEflowAngle(BPMData%nSect));   
    select case(BPMData%computeTerm(BPM_TEFLOWANGLE))
      case ('DEFAULT')
        BPMdata%TEflowAngle = 14.0/rad2deg
      case ('USERVALUE')  
        BPMdata%TEflowAngle = NamelistData%TEflowAngle
    end select    
     
    call SetNonComputableTerm(BPMData%computeTerm, BPMData%InFile, BPM_SECTAOA, SectAOAFlag, nTau, nbAllocate)
    allocate(BPMdata%SectAOA(BPMData%nSect, nbAllocate))
    select case(BPMData%computeTerm(BPM_SECTAOA))
      case ('DEFAULT')
        BPMdata%SectAOA = 0.0
      case ('USERVALUE')  
        BPMdata%SectAOA(:,1) = NamelistData%SectAOA
    end select

    call SetNonComputableTerm(BPMData%computeTerm, BPMData%InFile, BPM_TIPLCS, TipLCSFlag, nTau, nbAllocate)
    allocate(BPMdata%tipLCS(nbAllocate))
    select case(BPMData%computeTerm(BPM_TIPLCS))
      case ('DEFAULT')
        BPMdata%TipLCS = 1.0  
      case ('USERVALUE')  
        BPMdata%TipLCS = NamelistData%TipLCS
    end select
    
    call SetComputableTerm(BPMData%computeTerm, BPMData%InFile, BPM_U, UFlag, nTau, nbAllocate)
    allocate(BPMdata%U(BPMData%nSect, nbAllocate))
    if (BPMData%computeTerm(BPM_U).eq.'USERVALUE') BPMdata%U(:,1) = NamelistData%U
        
    allocate(BPMdata%R(BPMData%nSect));      BPMData%R      = 0.
    allocate(BPMdata%Theta(BPMData%nSect)); BPMData%Theta = 0.
    allocate(BPMdata%Phi(BPMData%nSect));   BPMData%Phi   = 0.
    
    if (UnitNumber.gt.0) then
      if (BPMData%uniformBlade) then
        nSect = 1
      else
        nSect = BPMData%nSect
      end if
      do i=1,nSect
        if (SectChordFlag.eq.1) then
          call ReadBinaryReal(unitNumber,tempSectChord,stat)
          if (BPMData%computeTerm(BPM_SECTCHORD).eq.'FILEVALUE') BPMData%SectChord(i) = tempSectChord
        end if
        if (SectLengthFlag.eq.1) then
          call ReadBinaryReal(unitNumber,tempSectLength,stat)
          if (BPMData%computeTerm(BPM_SECTLENGTH).eq.'FILEVALUE') BPMData%SectLength(i) = tempSectLength
        end if
        if (TEThicknessFlag.eq.1) then
          call ReadBinaryReal(unitNumber,tempTEThickness,stat)
          if (BPMData%computeTerm(BPM_TETHICKNESS).eq.'FILEVALUE') BPMData%TEThickness(i) = tempTEThickness
        end if
        if (TEflowAngleFlag.eq.1) then
          call ReadBinaryReal(unitNumber,tempTEflowAngle,stat)
          if (BPMData%computeTerm(BPM_TEFLOWANGLE).eq.'FILEVALUE') BPMData%TEflowAngle(i) = tempTEflowAngle
        end if
      end do
      if (BPMData%timeType.ne.LOAD_TIMETYPE_APERIODIC) then
        call ReadBPMData(BPMdata, unitNumber, nTau)
      end if
      if (BPMData%uniformBlade) then
        BPMData%SectChord(:)   = BPMData%SectChord(1)
        BPMData%SectLength(:)  = BPMData%SectLength(1)
        BPMData%TEThickness(:) = BPMData%TEThickness(1)
        BPMData%TEflowAngle(:) = BPMData%TEflowAngle(1)
        BPMData%SectAOA(:,1)   = BPMData%SectAOA(:,1)
        BPMData%U(:,1)         = BPMData%U(:,1)
      end if
    end if

    BPMdata%TEflowAngle  = BPMdata%TEflowAngle*rad2deg   ! units of rad
    BPMdata%SectAOA      = BPMdata%SectAOA    *rad2deg   ! units of rad
    !ksb debug:  Baofent doesn't think this is correct tipLCS is tip lift curve slope (not an angle)
    !BPMData%tipLCS       = BPMdata%tipLCS     /rad2deg   ! units of (1/rad)
  end subroutine CreateBPMArrays
     
  subroutine SetComputableTerm(computeTerm, InFile, Term, Flag, nTau, nbAllocate)
    implicit none
    character(len=*), dimension(:):: computeTerm
    logical, dimension(:):: InFile
    integer:: Term, Flag, nTau, nbAllocate   
    
    select case (trim(computeTerm(Term)))
      case ('')
        select case(Flag)
          case(0)
            computeTerm(Term) = 'COMPUTE'
            nbAllocate = 1
          case (1)
            InFile(Term) = .true.
            computeTerm(Term) = 'FILEVALUE'
            nbAllocate = nTau
          case default
            call Error('Unknown value given for '//trim(computeTerm(Term))//' in broadband data file.')
        end select
      case ('COMPUTE')
        nbAllocate = 1
        computeTerm(Term) = 'COMPUTE'        
        select case(Flag)
          case (0)
          case (1)
            InFile(Term) = .true.
          case default
            call Error('Unknown value given for '//trim(computeTerm(Term))//' in broadband data file.')
        end select 
      case ('USERVALUE')
        nbAllocate = 1
        computeTerm(Term) = 'USERVALUE'        
        select case(Flag)
          case (0)
          case (1)
            InFile(Term) = .true.
          case default
            call Error('Unknown value given for '//trim(computeTerm(Term))//' in broadband data file.')
        end select
      case ('FILEVALUE')
        select case(Flag)
          case (0)
            computeTerm(Term) = 'COMPUTE'
            nbAllocate = 1            
          case (1)           
            InFile(Term) = .true.
            nbAllocate = nTau             
          case default
            call Error('Unknown value given for '//trim(computeTerm(Term))//' in broadband data file.')
        end select
      case default
        call Error(''//trim(computeTerm(Term))//' must be blank or have one of the',&
                   'following values: `COMPUTE`, `USERVALUE` or `FILEVALUE`.')
    end select
  end subroutine SetComputableTerm
  
  
 subroutine SetNonComputableTerm(computeTerm, InFile, Term, Flag, nTau, nbAllocate)
    implicit none
    character(len=*), dimension(:):: computeTerm
    logical, dimension(:):: InFile
    integer:: Term, Flag, nTau, nbAllocate   
    
    select case (trim(computeTerm(Term)))
      case ('')
        select case(Flag)
          case(0)
            computeTerm(Term) = 'DEFAULT'
            nbAllocate = 1
          case (1)
            InFile(Term) = .true.
            computeTerm(Term) = 'FILEVALUE'
            nbAllocate = nTau
          case default
            call Error('Unknown value given for '//trim(computeTerm(Term))//' in broadband data file.')
        end select
      case ('USERVALUE')
        nbAllocate = 1
        computeTerm(Term) = 'USERVALUE'        
        select case(Flag)
          case (0)
          case (1)
            InFile(Term) = .true.
          case default
            call Error('Unknown value given for '//trim(computeTerm(Term))//' in broadband data file.')
        end select
      case ('FILEVALUE')
        select case(Flag)
          case (0)
            computeTerm(Term) = 'DEFAULT'
            nbAllocate = 1            
          case (1)           
            InFile(Term) = .true.
            nbAllocate = nTau             
          case default
            call Error('Unknown value given for '//trim(computeTerm(Term))//' in broadband data file.')
        end select
      case default
        call Error(''//trim(computeTerm(term))//' must be blank or have one of the',&
                   'following values: `COMPUTE`, `USERVALUE` or `FILEVALUE`.')
    end select
  end subroutine SetNonComputableTerm 
  
  subroutine ReadPeggData(PeggData, UnitNumber, iTau)
    implicit none
    type(PEGG):: PeggData
    integer:: unitNumber, stat, i, iTau
    real:: tempClbar, tempTotalThrust
    type(vector):: tempHubAxis
    
    do i=1,iTau
      if (PeggData%timeType.ne.LOAD_TIMETYPE_CONSTANT) then
        call ReadBinaryReal(unitNumber,PeggData%tau(i),stat)
      end if
      if (PeggData%InFile(PEGG_CLBAR)) then
        call ReadBinaryReal(unitNumber,tempClBar,stat)
        if (PeggData%computeTerm(PEGG_CLBAR).eq.'FILEVALUE') PeggData%ClBar(i) = tempClBar
      end if
      if (PeggData%InFile(PEGG_HUBAXIS)) then
        call ReadBinaryReal1DArray(unitNumber,temphubAxis%A,stat)
        if (PeggData%computeTerm(PEGG_HUBAXIS).eq.'FILEVALUE') PeggData%hubAxis(i) = temphubAxis
      end if
      if (PeggData%InFile(PEGG_TOTALTHRUST)) then
        call ReadBinaryReal(unitNumber,tempTotalThrust,stat)
        if (PeggData%computeTerm(PEGG_TOTALTHRUST).eq.'FILEVALUE') then
          PeggData%shaftAxisThrust = PeggData%shaftAxisThrust + tempTotalThrust
        end if
      end if
    end do
  end subroutine ReadPeggData
  
  subroutine ReadBPMData(BPMData, UnitNumber, iTau)
    implicit none
    type(BPM):: BPMData
    integer:: unitNumber, nSect, iTau, i, j, stat
    real:: tempTipLCS
    real, dimension(:), allocatable:: tempSectAOA, tempU
    
    if (BPMdata%uniformBlade) then
      nSect = 1
    else
      nSect = BPMData%nSect
    end if
  
    allocate(tempSectAOA(nSect), tempU(nSect))
    do j=1,iTau
      if (BPMData%timeType.ne.LOAD_TIMETYPE_CONSTANT) then
        call ReadBinaryReal(unitNumber,BPMData%tau(iTau),stat)
      end if
      if (BPMData%InFile(BPM_SECTAOA)) then
        call ReadBinaryReal1DArray(unitNumber,tempSectAOA,stat)
        if (BPMData%computeTerm(BPM_SECTAOA).eq.'FILEVALUE') BPMData%SectAOA(:,j) = tempSectAOA
      end if
      if (BPMData%InFile(BPM_TIPLCS)) then
        call ReadBinaryReal(unitNumber,tempTipLCS,stat)
        if (BPMData%computeTerm(BPM_TIPLCS).eq.'FILEVALUE') BPMData%TipLCS(j) = tempTipLCS
      end if
      if (BPMData%InFile(BPM_U)) then
        call ReadBinaryReal1DArray(unitNumber,tempU,stat)
        if (BPMData%computeTerm(BPM_U).eq.'FILEVALUE') BPMData%U(:,j) = tempU
      end if
  
!      if (BPMdata%uniformBlade) then
!        if (BPMData%computeTerm(BPM_SectAOA).eq.'FILEVALUE') BPMData%SectAOA(:,j) = BPMData%SectAOA(1,j)
!        if (BPMData%computeTerm(BPM_U).eq.'FILEVALUE') BPMData%U(:,j) = BPMData%U(1,j)
!      end if
    end do

    deallocate(tempSectAOA, tempU)
  end subroutine ReadBPMData 
  
     
  !This routine calculates the values needed for the semi-empirical
  !approachs to predicting broadband noise.  
  subroutine ComputePeggBroadbandTerms(PeggData, radiationVector, BBkey)
    implicit none
    type(PEGG), pointer:: PeggData
    type(vector):: radiationVector, hubAxis
    real:: temp
    integer:: BBkey

    PeggData%R = vectorAbsolute(radiationVector)

    if (PeggData%computeTerm(PEGG_HUBAXIS).eq.'COMPUTE') then
      hubAxis = resultantBase%rotation*GetHubAxis(PeggData%CBList)
    else
      hubAxis = PeggData%HubAxis(GetPeggTermKey(PeggData, PEGG_HUBAXIS, BBkey))
    end if
    hubAxis = -hubAxis/vectorAbsolute(hubAxis)

    ! Angle between negative thrust axis and observer (deg)
    temp = vectorDotProduct(hubAxis,radiationVector)/vectorAbsolute(radiationVector)
    if (abs(temp).lt.(1.-epsilon) .or. abs(temp).gt.(1.+epsilon)) then
      PeggData%theta1 = acos(temp)
    else
      ! This is either 0 or pi, but Pegg uses this term with 
      ! cos(theta)^2 so it doesn't matter.
      PeggData%theta1 = 0.
    end if

    PeggData%TipSpeed = PeggData%bladeRadius*abs(PeggData%RotSpeed)
    if (PeggData%computeTerm(PEGG_CLBAR).eq.'COMPUTE') then
      !ClBar = 6*Ct/Sigma = 6*(T/(rho*A*(TipSpeed)**2))/(totalBladeArea/A) = 6*(T/(rho*totalBladeArea*(TipSpeed)**2))
      temp = 6*PeggData%ShaftAxisThrust/(rho*PeggData%totalBladeArea*(PeggData%TipSpeed)**2.)
      PeggData%ClBar(1) = 6*PeggData%ShaftAxisThrust/(rho*PeggData%totalBladeArea*(PeggData%TipSpeed)**2.)
    end if
  end subroutine ComputePeggBroadbandTerms
 
 subroutine ComputeBPMBroadbandTerms(BPMdata, surfNormal, radiationVector, srcVelocity, obsCBList, i, AOAkey, Ukey, sourceTime)
    implicit none
    type(BPM), pointer:: BPMData
    type(vector), intent(in):: SurfNormal
    type(vector):: radiationVector, srcVelocity, obsVelocity, M_src, M_obs, radUnitVector
    type(CBStructure), pointer::obsCBList
    real:: sourceTime
    integer:: i, AOAkey, Ukey

    !type(matrix):: TwistMatrix
    type(vector):: yAxisLocal, BladeXAxisGlobal, BladeYAxisGlobal, BladeZAxisGlobal, crossProd
    real:: tempReal

    BPMData%R(i) = vectorAbsolute(radiationVector)
    
    yAxisLocal       = vectorSetCoordinates(0.0, 1.0, 0.0)
    BladeZAxisGlobal = vectorBaseChange(SurfNormal)
    BladeYAxisGlobal = VectorBaseChange(yAxisLocal)
    BladeXAxisGlobal = BladeYAxisGlobal.cross.BladeZAxisGlobal

    ! Theta -> Find the angle off of the blade-fixed +x axis in the XY plane
    tempReal    = vectorDotProduct(BladeXAxisGlobal, radiationVector)/vectorAbsolute(BladeXAxisGlobal)/BPMData%R(i)
    ! This compiler has an issue with taking inverse cosine very close to +/-1
    if (abs(tempReal).lt.(1.-epsilon) .or. abs(tempReal).gt.(1.+epsilon)) then
      BPMData%Theta(i) = acos(tempReal)
    else
      if (tempReal.gt.(1.-epsilon)) then
        ! If the observer is directly beneath the rotor
        
        ! BPM base their directivity formulas on the assumption that they are dealing with 
        ! a semi-infinite plate.  Because of this, the method is faulty when the directivity angle
        ! goes to zero.  Simply move it off by a very small amount.        
        BPMData%Theta(i) = epsilon
      else
        ! If the observer is directly overhead
        BPMData%Theta(i) = pi-epsilon
      end if
    end if
    
    ! Phi -> Find the angle off of the blade-fixed +Y axis in the YZ plane
    crossProd = BladeXAxisGlobal.cross.radiationVector
    tempReal = vectorDotProduct(crossProd,BladeZAxisGlobal)/vectorAbsolute(crossProd)/vectorAbsolute(BladeZAxisGlobal)
    if (abs(tempReal).lt.(1.-epsilon) .or. abs(tempReal).gt.(1.+epsilon)) then
      BPMData%Phi(i) = acos(tempReal)
    else
      if (tempReal.gt.(1.-epsilon)) then
        ! If the observer is directly beneath the rotor
        
        ! BPM base their directivity formulas on the assumption that they are dealing with 
        ! a semi-infinite plate.  Because of this, the method is faulty when the directivity angle
        ! goes to zero.  Simply move it off by a very small amount.
        BPMData%Phi(i) = epsilon
      else
        ! If the observer is directly overhead
        BPMData%Phi(i) = pi-epsilon
      end if
    end if

    if (BPMdata%computeTerm(BPM_U).eq.'COMPUTE') then
      BPMdata%U(i,Ukey) = vectorAbsolute(srcVelocity)
    end if
    
    M_src = vectorRealDivide(srcVelocity,c)  !ksb debug: Check these - they don't seem general.
    M_obs = vectorSetCoordinates(0.0, 0.0, 0.0)
    radUnitVector = vectorRealDivide(radiationVector,vectorAbsolute(radiationVector))
    if (associated(ObsCBList)) then
      call CalculateResultantChangeOfBase(obsCBList, GetSize(obsCBList), sourceTime)
      M_obs = vectorRealDivide(resultantBase%velocity,c)
    end if
    BPMData%DopplerShift(i) = (1-vectorDotProduct(M_src,radUnitVector))/&
                              (1-vectorDotProduct(M_obs,radUnitVector))
  end subroutine ComputeBPMBroadbandTerms

  recursive function GetHubAxis(CBList) result(HubAxis)
    type(CBStructure), pointer::CBList
    type(vector)::HubAxis
    HubAxis = vectorSetCoordinates(0., 0., 0.)
    if (CBList%rotation) then
      HubAxis = CBList%AxisValue
    else
      HubAxis = GetHubAxis(CBList%next)
    end if
  end function GetHubAxis
  
  recursive function GetOmega(CBList) result(Omega)
    type(CBStructure), pointer::CBList
    real:: Omega
    Omega = 0.
    if (CBList%rotation) then
      omega = CBList%funcData%omega
    else
      Omega = GetOmega(CBList%next)
    end if
  end function GetOmega


  subroutine ComputePeggIntensity(PeggData, Intensity, BBkey)
    type(PEGG), pointer:: PeggData
    real(kind=8), dimension(:), pointer:: Intensity
    integer:: BBkey
    
    integer::i, CenterBand, fp
    real:: ClFcn, SPLConst, ThirdOctaveLimit, SPL
    real, dimension(37):: BandLevel
    !real, dimension(BB_THIRD_OCTAVE_BANDS):: frCen
    real, dimension(BB_THIRD_OCTAVE_BANDS):: frCen


    if (PeggData%spreadBands) then
      BandLevel = GetPeggInterpSpectrum()
    else
      BandLevel = GetPeggBaselineSpectrum()
    end if

    frCen = GetThirdOctFreqArray()

    ! These first two steps don't need to be done if we're using non-aperiodic data
    ! ------------------------------------------------------------------------------
    !Step 1: Compute the peak broadband noise frequency
    fp = -240.*log10(PeggData%shaftAxisThrust+epsilon**3.0)+2.448*PeggData%TipSpeed+942.       !metric

    ThirdOctaveLimit = 2.**(1./6.)
    !Step 2: Identify the peak frequency in a standard 1/3 octave band
    if (fp<frCen(1)/ThirdOctaveLimit .or. fp>frCen(size(frCen))*ThirdOctaveLimit) then
      call error("The peak broadband noise frequency is outside the limit", &
                 "of what can be handled for Pegg's approach.")
    else
      do i=1,size(frCen)
        if (fp < frCen(i)*ThirdOctaveLimit) then
	  ! Identify the center band and offset it so that
          ! it corresponds to the location of the peak 
          ! frequency correction in the BandLevel array (16)
          CenterBand = 16 - i
          exit
        end if
      end do
    end if

    !Step 3: Compute the constant term within SPL_1/3
    if (PeggData%ClBar(BBkey).le.0.48) then
      ClFcn = 10.*log10(PeggData%ClBar(BBkey)/0.4)
    else
      ClFcn = 0.8 + 80.*log10(PeggData%ClBar(BBkey)/0.48)
    end if
    SPLConst = 60.*log10(PeggData%TipSpeed/c)+10.*log10(PeggData%totalBladeArea/(PeggData%R**2.)* &
        (cos(PeggData%theta1)**2.+0.1))+ClFcn+130. + PeggData%deltadB

    ! Step 4: Compute the SPL_1/3 for all octave bands
    ! Don't calculate broadband levels for anything past the 10 kHz band
    do i=1,BB_THIRD_OCTAVE_BANDS
      if ((CenterBand+i).gt.0 .and. (Centerband+i).le.size(BandLevel) .and. BandLevel(CenterBand+i).lt.0.0) then
        SPL = BandLevel(CenterBand+i) + SPLConst
        Intensity(i) = Intensity(i) + ConvertdBtoMSP(SPL)
      end if
    end do
  end subroutine ComputePeggIntensity

  function GetPeggBaselineSpectrum() result(FreqArray)
    implicit none
    real, dimension(37):: BaseBandLevel, FreqArray

    data BaseBandLevel     / -29.0 ,     0.0  ,     0.0  ,   -24.5  ,     0.0  , &
                               0.0 ,   -19.5  ,     0.0  ,     0.0  ,   -15.3  , &
              	               0.0 ,     0.0  ,   -11.7  ,     0.0  ,     0.0  , &
                              -7.5 ,     0.0  ,     0.0  ,   -11.5  ,     0.0  , &
                               0.0 ,   -12.1  ,     0.0  ,     0.0  ,   -16.5  , &
                               0.0 ,     0.0  ,   -17.0  ,     0.0  ,     0.0  , &
                             -21.8 ,     0.0  ,     0.0  ,   -26.4  ,     0.0  , &
                               0.0 ,   -30.0  /
    FreqArray = BaseBandLevel
  end function GetPeggBaselineSpectrum

  function GetPeggInterpSpectrum() result(FreqArray)
    implicit none
    real, dimension(37):: InterpBandLevel, FreqArray

    data InterpBandLevel  / -29.0 ,  -27.83035528,   -26.35669527,  -24.5,  -23.20039475, &
			    -21.56299474, -19.5, -18.40833159, -17.03291558, -15.3, &
     			    -14.36428422, -13.18535621, -11.7, -10.60833159, -9.232915582, &
			    -7.5, -8.5396842, -9.849604208, -11.5, -11.65595263, &
			    -11.85244063, -12.1, -13.24365262, -14.68456463, -16.5, &
			    -16.62996052, -16.79370053, -17.0, -18.24762104, -19.81952505, &
		            -21.8, -22.99563683, -24.50204484, -26.4, -27.33571578, &
			    -28.51464379, -30.0  /

    FreqArray = InterpBandLevel
  end function GetPeggInterpSpectrum
  
  subroutine ComputeBPMintensity(BPMdata, Intensity, BBkey, i)
!               ----------------------------------
!               ****** Variable Definitions ******
!               ----------------------------------

!    Variable Name               Definition                     Units
!    -------------               ----------                     -----

!       TipAOA              tip angle of attack                 Degrees
!       SectAOA             segment angle of attack             Degrees
!       TipLCS              tip lift curve slope                ---
!       chord               segment chord length                meters
!       c                   speed of sound                      meters/sec
!       frCen               1/3 octave centered frequencies     hertz
!       TEthickness         segment trailing edge thickness     meters
!       BluntNoise          flag to compute bluntness noise     --- 
!       LBLNoise            flag to compute lbl noise           ---
!       TipNoise            flag to compute tip noise           ---
!       BLTrip               flag to trip boundary layer         ---
!       TBLTEnoise          flag to compute tblte noise         ---
!       L                   segment span length                 meters
!       maxFreq             maximu1m number of frequencies       ---
!       maxSeg              maximu1m number of segments          ---
!       nFreq               number of 1/3 octave frequencies    ---
!       nSeg                number of segments                  ---
!       p1                  pressure associated with            
!                             TBLTE prediction                  nt/m2
!       p2                  pressure associated with            
!                             TBLTE prediction                  nt/m2
!       p3                  pressure associated with            
!                             TBLTE prediction                  nt/m2
!       p4                  pressure associated with            
!                             TOTAL prediction                  nt/m2
!       p5                  pressure associated with            
!                             LBLVS prediction                  nt/m2
!       p6                  pressure associated with            
!                             bluntness prediction              nt/m2
!       p7                  pressure associated with            
!                             tip noise prediction              nt/m2
!       Phi                directivity angle                   degrees
!       TEflowAngle                 bluntness angle                     degrees
!       R                   segment to observer distance        meters
!       RoundTip            logical indicated rounded tip       ---
!       SPL                 total sound pressure level          dB
!       SPLAplh             sound pressure level associated     
!                             with TBLTE prediction             dB
!       SPLblnt             sound pressure level associated
!                             with bluntness prediction         dB
!       SPLlbl              sound pressure level associated
!                             with LBL prediction               dB
!       SPLp                sound pressure level associated
!                             with TBLTE prediction             dB
!       SPLs                sound pressure level associated
!                             with TBLTE prediction             dB
!       SPLtbl              sound pressure level associated
!                             with TBLTE prediction             dB
!       SPLtip              sound pressure level associated
!                             with tip noise prediction         dB
!       St                  strouhal number                     ---
!       Theta              directivity angle                   degrees
!       U                   segment freestream velocity         meters/sec
!       nu                  kinematic viscosity                 m2/sec
    type(BPM), pointer:: BPMdata
    real(kind=8), dimension(:,:):: Intensity
    real, dimension(BB_THIRD_OCTAVE_BANDS)::  frCen, shiftFrCen, St, SPLlbl, SPLtbl, &
                                SPLp, SPLs, SPLalpha, SPLblnt, &
                                SPLtip
    real:: M, RC, AOA
    integer:: i, j, BBkey, Ukey, AOAkey, LCSKey
    
  !  Set up values of 1/3 octave centered frequencies
  !  ------------------------------------------------

    frCen = GetThirdOctFreqArray()

    UKey   = getBPMTermKey(BPMData, BPM_U, BBkey)
    AOAkey = getBPMTermKey(BPMData, BPM_SECTAOA, BBkey)
    LCSKey = getBPMTermKey(BPMData, BPM_TipLCS, BBkey)

    M = -huge(c)
    SPlp     = M
    SPLs     = M
    SPLalpha = M
    SPLlbl   = M
    SPLblnt  = M
    SPLtip   = M

    !    Compute Reynolds number and mach number
    !    ---------------------------------------
    M = BPMData%U(i,UKey)/c
    RC = BPMData%U(i,UKey)*BPMData%SectChord(i)/nu
    ShiftFrCen = frCen*BPMData%DopplerShift(i)

    AOA = abs(BPMData%SectAOA(i,AOAkey))  ! This is a temporary fix until a more permanent
                                          ! method is developed to switch the pressure and 
                                          ! suction sides of the blade when AOA<0.  BAG 12/19/2012

    if (BPMdata%TBLTEnoise) then
      call TBLTE(AOA, BPMData%SectChord(i), BPMData%U(i,UKey), ShiftFrCen, BPMData%BLTrip, &
        SPLp, SPLs, SPLalpha, SPLtbl, BPMData%Theta(i), BPMData%Phi(i),  &
        BPMData%SectLength(i), BPMData%R(i), M, RC)
    end if

    if (BPMdata%LBLVSnoise) then
      call LBLVS(AOA, BPMData%SectChord(i), BPMData%U(i,UKey), ShiftFrCen, SPLlbl, &
        BPMData%Theta(i), BPMData%Phi(i), BPMData%SectLength(i), BPMData%R(i), M, RC)
    end if
      
    if (BPMdata%BluntNoise) then
      call Blunt(AOA, BPMData%SectChord(i), BPMData%U(i,UKey), ShiftFrCen, BPMData%BLtrip, &
        SPLblnt, BPMData%Theta(i), BPMData%Phi(i), BPMData%SectLength(i), BPMData%R(i), &
        BPMData%TEthickness(i), BPMData%TEflowAngle(i), M, RC)
    end if
      
    if (BPMdata%BladeTipNoise .and. (i .eq. BPMData%nSect)) then
      call TipNoise(AOA, BPMData%TipLCS(LCSKey), BPMData%SectChord(i), &
        BPMData%U(i,UKey), ShiftFrCen, SPLtip, BPMData%Theta(i), BPMData%Phi(i), &
        BPMData%R(i), BPMData%RoundBladeTip, M)
    end if

  !   Add in this segment's contribution on a mean-square pressure basis
  !   ------------------------------------------------------------------
    do j=1,BB_THIRD_OCTAVE_BANDS
      if (BPMdata%TBLTEnoise) then
        Intensity(1,j) = Intensity(1,j) + ConvertdBtoMSP(SPLp(j))
        Intensity(2,j) = Intensity(2,j) + ConvertdBtoMSP(SPLs(j))
        Intensity(3,j) = Intensity(3,j) + ConvertdBtoMSP(SPLalpha(j))
      end if

      if (BPMdata%LBLVSnoise) then
        Intensity(4,j) = Intensity(4,j) + ConvertdBtoMSP(SPLlbl(j))
      end if        
        
      if (BPMdata%BluntNoise) then
        Intensity(5,j) = Intensity(5,j) + ConvertdBtoMSP(SPLblnt(j))
      end if
        
      if (BPMdata%BladeTipNoise .and. (i .eq. BPMdata%nSect)) then
        Intensity(6,j) = Intensity(6,j) + ConvertdBtoMSP(SPLtip(j))
      end if
    end do

  ! Stopping the routine here leaves us with Intensity 
    ! --------------------------------------------------
    
    
  !  Contributions from all segments are now accounted for
  !  Compute sound pressure levels for each mechanism and
  !  for the total.
  !  -----------------------------------------------------

!    do i=1,nFreq
!      if (p1(i) .ne. 0.) SPL(1,i) = 10.*log10(p1(i))
!      if (p2(i) .ne. 0.) SPL(2,i) = 10.*log10(p2(i))
!      if (p3(i) .ne. 0.) SPL(3,i) = 10.*log10(p3(i))
!      if (p4(i) .ne. 0.) SPL(4,i) = 10.*log10(p4(i))
!      if (p5(i) .ne. 0.) SPL(5,i) = 10.*log10(p5(i))
!      if (p6(i) .ne. 0.) SPL(6,i) = 10.*log10(p6(i))
!      if (p7(i) .ne. 0.) SPL(7,i) = 10.*log10(p7(i))
!    end do
    
  !  Write output file
  !  -----------------

!    write (5,7000)
!    
!    do i=1,nFreq
!      write(5,7100) frCen(i), (SPL(j,i), j=1,3), (SPL(j,i), j=4,6), SPL(7,i)
!      
!      if (mod(i,5) .eq. 0) write(5,7200)
!    end do
!    
!    7000 format (1H1, 52X, 'One-Third Octave',/,50x,'Sound Pressure Levels', &
!                  ////,5x,'                ','    Pressure    ', &
!                           '    Suction     ','   Separation  '/, &
!                        5x,' Frequency (Hz) ','    Side TBL   ', &
!                           '   Side TBL     ','    Side TBL   ', &
!                           '    Laminar     ','    Bluntness  ', &
!                           '      Tip       ','     Total     ', &
!                           /,5x,8('---------------- '),/)
!                           
!    7100 format (8F15.3)
!    7200 format ('  ')
!    8000 format (I3)
!    8002 format(4I10)
  end subroutine ComputeBPMintensity
  
  subroutine LBLVS (SectAOA, chord, U, frCen, SPLlam, Theta, Phi, L, R, M, RC)          
                
!               ----------------------------------
!               ****** Variable Definitions ******
!               ----------------------------------

!    Variable Name               Definition                     Units
!    -------------               ----------                     -----

!       SectAOA           segment angle of attack             Degrees
!       chord               segment chord length                meters
!       c                   speed of sound                      meters/sec
!       D                   Reynolds number ratio               ---
!       DbarH               high frequency directivity          ---
!       deltaP              pressure side boundary layer
!                             thickness                         meters
!       dstrP               pressure side boundary layer
!                             displacement thickness            meters
!       dstrS               suction side boundary layer
!                             displacement thickness            meters
!       E                   Strouhal number ratio               ---
!       frCen               1/3 octave centered frequencies     hertz
!       G1                  sound pressure level function       dB
!       G2                  overall sound pressure level
!                             function                          dB
!       G3                  overall sound pressure level 
!                             function                          dB
!       BLTrip               flag to trip boundary layer         ---
!       L                   segment span length                 meters
!       M                   mach number                         ---
!       nFreq               number of 1/3 octave frequencies    ---
!       OASPL               overall sound pressure level        dB
!       Phi                directivity angle                   degrees
!       R                   segment to observer distance        meters
!       Rc                  reynolds number based on chord      ---
!       Rc0                 reference reynolds number           ---
!       scale               geometric scaling term              
!       SPLlam              sound pressure level due to 
!                             laminar mechanism                 dB
!       StPrim              strouhal number based on pressure
!                             side boundary layer thickness     ---
!       St1Prim             reference strouhal number           ---
!       StPkPrm             peak strouhal number                ---
!       Theta              directivity angle                   degrees
!       U                   segment freestream velocity         meters/sec
!       nu                  kinematic viscosity                 m2/sec
    implicit none
    real, dimension(BB_THIRD_OCTAVE_BANDS):: StPrim, SPLlam, frCen
    real:: SectAOA, chord, U, L, M, RC, Theta, Phi, R, term, deltaP, dstrS, dstrP, DbarH
    real:: St1Prim, StPkPrm, Rc0, D, G2, G3, scale, E, G1
    integer:: i, BLTrip
    
    term = M**5.*L/R**2.
    
!    Compute boundary layer thickness
!    --------------------------------
    call Thick(chord, SectAOA, BLTrip, deltaP, dstrS, dstrP, M, RC)
    
!    Compute directivity function
!    ----------------------------
    call DirectH(M, Theta, Phi, DbarH)
    
!    Compute reference Strouhal Number
!    ---------------------------------
    if  (RC .le. 1.3E+05) then
      St1Prim = 0.18
    else if (RC .le. 4.0E+05) then
      St1Prim=0.001756*RC**0.3931
    else 
      St1Prim = 0.28
    end if
    
    StPkPrm = 10.**(-0.04*SectAOA)*St1Prim
    
!    Compute Reference Reynolds Number
!    ---------------------------------
    if (SectAOA .le. 3.0) Rc0 = 10.**(0.215*SectAOA+4.978)
    if (SectAOA .gt. 3.0) Rc0 = 10.**(0.120*SectAOA+5.263)
    
!    Compute peak scaled spectrum level
!    ----------------------------------
    D = Rc/Rc0
    if (D .le. 0.3237) then
      G2 = 77.852*log10(D)+15.328
    else if (D .le. 0.5689) then
      G2 = 65.188*log10(D)+9.125
    else if (D .le. 1.7579) then
      G2 = -114.052*log10(D)**2.
    else if (D .le. 3.0889) then
      G2 = -65.188*log10(D)+9.125
    else
      G2 = -77.852*log10(D)+15.328
    end if
    
    G3 = 171.04-3.03*SectAOA
    scale = 10.*log10(deltaP*DbarH*term)
    
    SPLlam = G2+G3+scale

!   Compute scaled sound pressure levels for each strouhal number
!   -------------------------------------------------------------
    do i=1,BB_THIRD_OCTAVE_BANDS
      StPrim(i) = frCen(i)*deltaP/U
      E = StPrim(i)/StPkPrm
      
      if (E .lt. 0.5974) then
        G1 = 39.8*log10(E)-11.12
      else if (E .le. 0.8545) then
        G1 = 98.409*log10(E)+2.0
      else if (E .lt. 1.17) then
        G1 = -5.076+sqrt(2.484-506.25*(log10(E))**2.)
      else if (E .lt. 1.674) then
        G1 = -98.409*log10(E)+2.0
      else 
        G1 = -39.8*log10(E)-11.12
      end if
      
      SPLlam(i) = SPLlam(i) + G1
    end do            
    return      
  end subroutine LBLVS
  
  subroutine TBLTE(SectAOA, chord, U, frCen, BLTrip, SPLp, SPLs, &
                    SPLalpha, SPLtbl, Theta, Phi, L, R, M, RC)
                    
!               ----------------------------------
!               ****** Variable Definitions ******
!               ----------------------------------

!    Variable Name               Definition                     Units
!    -------------               ----------                     -----

!       A                   strouhal number ratio               ---
!       A0                  function used in 'A' calculation    ---
!       A02                 function used in 'A' calculation    ---
!       AA                  'A' spectrum shape evaluated at
!                             strouhal number ratio             dB
!       SectAOA           segment angle of attack             Degrees
!       AmaxA               maximu1m 'A' curve evaluated at
!                             strouhal number ratio             dB
!       AmaxA0              maximu1m 'A' curve evaulated at A0   dB
!       AmaxA02             maximu1m 'A' curve evaluated at A02  dB
!       AmaxB               maximu1m 'A' curve evaluated at B    dB
!       AminA               minimu1m 'A' curve evaluated at 
!                             strouhal number ratio             dB
!       AminA0              minimu1m 'A' curve evaluated at A0   dB
!       AminA02             minimu1m 'A' curve evaluated at A02  dB
!       AminB               minimu1m 'A' curve evaluated at B    dB
!       ArA0                interpolation factor                ---
!       ArA02               interpolation factor                ---
!       B                   strouhal number ratio               ---
!       B0                  function used in 'B' calculation    ---
!       BB                  'B' spectrum shape evaluated at
!                             strouhal number ratio             dB
!       beta                used in 'B' computation             ---
!       beta0               used in 'B' computation             ---
!       BmaxB               maximu1m 'B' evaluated at B          dB
!       BmaxB0              maximu1m 'B' evaluated at B0         dB
!       BminB               minimu1m 'B' evaluated at B          dB
!       BminB0              minimu1m 'B' evaluated at B0         dB
!       BrB0                interpolation factor                ---
!       chord               segment chord length                meters
!       c                   speed of sound                      meters/sec
!       DbarH               high frequency directivity          ---
!       DbarL               low frequency directivity           ---
!       delK1               correction to amplitude function    dB
!       deltaP              pressure side boundary layer
!                             thickness                         meters
!       dstrP               pressure side displacement thickness meters
!       dstrS               suction side displacement thickness meters
!                             layer thickness
!       frCen               1/3 octave centered frequencies     hertz
!       gamma               used in 'B' computation             ---
!       gamma0              used in 'B' computation             ---
!       BLTrip               flag to trip boundary layer         ---
!       K1                  amplitude function                  dB
!       K2                  amplitude function                  dB
!       L                   span                                meters
!       M                   mach number                         ---
!       nFreq               number of 1/3 octave frequencies    ---
!       Phi                directivity angle                   degrees
!       p1                  pressure side pressure              nt/m2
!       p2                  suction side pressure               nt/m2
!       p4                  pressure from angle of attack 
!                             contribution                      nt/m2
!       R                   segment to observer distance        meters
!       Rc                  reynolds number based on chord      ---
!       RdstrP              reynolds number based on pressure
!                             side displacement thickness       ---
!       RdstrS              reynolds number based on suction
!                             side displacement thickness       ---
!       SPLAplh             sound pressure level associated     
!                             with TBLTE prediction             dB
!       SPLp                sound pressure level associated
!                             with TBLTE prediction             dB
!       SPLs                sound pressure level associated
!                             with TBLTE prediction             dB
!       SPLtbl              sound pressure level associated
!                             with TBLTE prediction             dB
!       StP                 pressure side strouhal number       ---
!       StS                 suction side strouhal number        ---
!       St1                 peak strouhal number                ---
!       St1Prim             peak strouhal number                ---
!       St2                 peak strouhal number                ---
!       StPeak              peak strouhal number                ---
!       SWITCH              logical for computation of angle
!                             of attack contribution            ---
!       Theta              directivity angle                   degrees
!       U                   segment freestream velocity         meters/sec
!       nu                  kinematic viscosity                 m2/sec                    
!       xCheck              used to check for angle of attack
!                             contribution                      --
    implicit none    
    real, dimension(BB_THIRD_OCTAVE_BANDS):: SPLtbl, SPLp, SPLs, SPLalpha, StP, StS, frCen
    logical:: SWITCH
    real:: SectAOA, chord, L, M, K1, K2, RC, term, U, Theta, Phi, R
    real:: deltaP, dstrS, dstrP, DbarL, DbarH, RdstrS, RdstrP, St1, St2, St1Prim, StPeak
    real:: delK1, gamma, gamma0, beta, beta0, xCheck
    real:: A, AA, AminA, AmaxA, A0, A02, AminA0, AmaxA0, AminA02, AmaxA02, Ara0, Ara02
    real:: AminB, AmaxB
    real:: B, BB, BminB, BmaxB, B0, BminB0, BmaxB0, Brb0
    real:: p1, p2, p3
    integer:: BLTrip
    
    integer:: i
    
    term = M**5.*L/R**2.
!   Compute Boundary Layer Thickness
!   --------------------------------
    call Thick(chord, SectAOA, BLTrip, deltaP, dstrS, dstrP, M, RC)
    
!   Compute Directivity Function
!   ----------------------------
    call directL(M, Theta, Phi, DbarL)
    call directH(M, Theta, Phi, DbarH)
    
!   Calculate the reynolds numbers based on pressure 
!   and suction displacement thickness        
!   ------------------------------------------------
    RdstrS = dstrS*U/nu
    RdstrP = dstrP*U/nu

!   Determine peak strouhal nubmers to be used for
!   'A' and 'B' curve calculations
!   ----------------------------------------------
    St1 = 0.02*M**(-0.6)
    
    if (SectAOA .le. 1.333) then
      St2 = St1
    else if (SectAOA .le. 12.5) then
      St2 = St1*10.**(0.0054*(SectAOA-1.333)**2.)
    else
      St2 = 4.72*St1
    end if
    
    St1Prim = (St1+St2)/2.
    
    call A0Comp(Rc, A0)
    call A0Comp(3.*Rc, A02)
    
!   Evaluate minimu1m and maximu1m 'A' curves at A0
!   ---------------------------------------------
    call Amin(A0, AminA0)
    call Amax(A0, AmaxA0)
    call Amin(A02, AminA02)
    call Amax(A02, AmaxA02)
    
!   Compute 'A' max/min ratio
!   -------------------------
    Ara0  = (20.+AminA0)/(AminA0-AmaxA0)
    Ara02 = (20.+AminA02)/(AminA02-AmaxA02)
    
!   Compute B0 to be used in 'B' curve calculations    
!   -----------------------------------------------
    if (Rc .lt. 9.52E+04) then 
      B0 = 0.3
    else if (Rc .lt. 8.57E+05) then
      B0 = (-4.48E-13)*(Rc-8.57E+05)**2.+0.56
    else
      B0 = 0.56
    end if
      
!   Evaluate minimu1m and maximu1m 'B' curves at B0
!   ---------------------------------------------
    call Bmin(B0, BminB0)
    call Bmax(B0, BmaxB0)
    
!   Compute 'B' max/min ratio
!   -------------------------
    BrB0 = (20.+BminB0)/(BminB0-BmaxB0)
    
!   For each center frequency, compute an
!   'A' prediction for the pressure side
!   -------------------------------------
    StPeak = St1
    
    do i=1,BB_THIRD_OCTAVE_BANDS
      StP(i) = frCen(i)*dstrP/U
      A      = log10(StP(i)/StPeak)
      call Amin(A, AminA)
      call Amax(A, AmaxA)
      AA     = AminA + ArA0*(AmaxA-AminA)

      if (Rc .lt. 2.47E+05) then
        K1 = -4.31*log10(Rc) + 156.3
      else if (Rc .lt. 8.0E+05) then
        K1 = -9.0*log10(Rc) + 181.6
      else
        K1 = 128.5
      end if
      
      if (RdstrP .le. 5000.) then
        delK1 = -SectAOA*(5.29 - 1.43*log10(RdstrP))
      else
        delK1 = 0.0
      end if
      
      SPLp(i) = AA + K1-3. + 10.*log10(dstrP*DbarH*term)+delK1

      gamma  =  27.094*M + 3.31
      beta   =  72.650*M + 10.74
      gamma0 =  23.430*M + 4.651
      beta0  = -34.190*M - 13.820
      
      if (SectAOA .le. (gamma0-gamma)) then
        K2 = -1000.0
      else if (SectAOA .le. (gamma0+gamma)) then
        K2 = sqrt(beta**2.-(beta/gamma)**2.*(SectAOA-gamma0)**2.)+beta0
      else
        K2 = -12.0
      end if
      
      K2 = K2 + K1
      
      StS(i) = frCen(i)*dstrS/U
      
!   Check for 'A' computation for suction side
!   ------------------------------------------
      xCheck = gamma0
      SWITCH = .false.
      if ((SectAOA .ge. xCheck).or.(SectAOA .gt. 12.5)) SWITCH=.TRUE.
      if (.not. SWITCH) then
        A = log10(StS(i)/St1Prim)
        call Amin(A, AminA)
        call Amax(A, AmaxA)
        AA = AminA + ArA0*(AmaxA-AminA)
        
        SPLs(i) = AA + K1 - 3. + 10.*log10(dstrS*DbarH*term)
        
!       'B' curve computation
!      ---------------------
        B = abs(log10(StS(i)/St2))
        call Bmin(B, BminB)
        call Bmax(B, BmaxB)
        BB = BminB + BrB0*(BmaxB-BminB)
        SPLalpha(i) = BB + K2 + 10.*log10(dstrS*DbarH*term)        
      else
!      The 'A' computation is dropped if 'SWITCH' is true
!      --------------------------------------------------
        SPLs(i) = 10.*log10(dstrS*DbarL*term)
        SPLp(i) = 10.*log10(dstrS*DbarL*term) 
        
        B = abs(log10(StS(i)/St2))
        call Amin(B, AminB)
        call Amax(B, AmaxB)
        BB = AminB + ArA02*(AmaxB-AminB)
        SPLalpha(i) = BB + K2 + 10.*log10(dstrS*DbarL*term)
      end if
      
!   Sum all contributions from 'A' and 'B' on both
!   pressure and suction side on a mean-square pressure basis
!   ---------------------------------------------------------
      if (SPLp(i)    .lt. -100.) SPLp(i)    = -100.
      if (SPLs(i)    .lt. -100.) SPLs(i)    = -100.
      if (SPLalpha(i).lt. -100.) SPLalpha(i) = -100.

!     determine whether the segment is in the inverse flow region     
!     For this case, the TBLTE noise is set to be zero.
!     Because the AOA is a very large negative value, .GT. -50 degree
!     Also the noise is not TBLTE noise at all. 
!     Modifided by Baofeng Cheng, 10/18/2012
      if (abs(sectAOA) >= 50.) then ! Not sure 50 deg is a good value - perhaps it should be smaller. 
          SPLp(i) = 0.0
          SPLs(i) = 0.0
          SPLalpha(i) = 0.0
      end if   

      p1 = 10.**(SPLp(i)   /10.)
      p2 = 10.**(SPLs(i)   /10.)
      p3 = 10.**(SPLalpha(i)/10.)
      
      SPLtbl(i) = 10.*log10(p1 + p2 + p3)
    end do
  
  end subroutine TBLTE

!  This  subroutine defines the curve fit corresponding
!  to the A-curve for the minimu1m allowed reynolds number.  
  subroutine Amin(A, AminA)
    real:: A, AminA, x1

    x1 = abs(A)
    if (x1 .le. 0.204) then
      AminA = sqrt(67.552 - 886.788*x1**2.) - 8.219
    else if (x1 .le. 0.244) then
      AminA = -32.665*x1 + 3.981
    else
      AminA = -142.795*x1**3. + 103.656*x1**2. - 57.757*x1 + 6.006
    end if
  end subroutine Amin

!  This  subroutine defines the curve fit corresponding
!  to the A-curve for the maximu1m allowed reynolds number.   
  subroutine Amax(A, AmaxA)
    real:: A, AmaxA, x1
    
    x1 = abs(A)
    if (x1 .le. 0.13) then
      AmaxA = sqrt(67.552 - 886.788*x1**2.) - 8.219
    else if (x1 .le. 0.321) then
      AmaxA = -15.901*x1 + 1.098
    else
      AmaxA = -4.669*x1**3. + 3.491*x1**2. - 16.669*x1 + 1.149
    end if
  end subroutine Amax
  
!  This  subroutine defines the curve fit corresponding
!  to the B-curve for the minimu1m allowed reynolds number.
  subroutine Bmin(B, BminB)
    real:: B, BminB, x1

    x1 = abs(B)
    if (x1 .le. 0.13) then
      BminB = sqrt(16.888 - 886.788*x1**2.) - 4.109
    else if (x1 .le. 0.145) then
      BminB = -83.607*x1 + 8.138
    else
      BminB = -817.81*x1**3. + 355.21*x1**2. - 135.024*x1 + 10.619
    end if
  end subroutine Bmin
  
!  This  subroutine defines the curve fit corresponding
!  to the B-curve for the maximu1m allowed reynolds number.
  subroutine Bmax(B, BmaxB)
    real:: B, BmaxB, x1  
  
    x1 = abs(B)
    if (x1 .le. 0.1) then
      BmaxB = sqrt(16.888 - 886.788*x1**2.) - 4.109
    else if (x1 .le. 0.187) then
      BmaxB = -31.313*x1 + 1.854
    else
      BmaxB = -80.541*x1**3. + 44.174*x1**2. - 39.381*x1 + 2.344
    end if
  end subroutine Bmax
  
!  This  subroutine determines where the A-curve
!  takes on a value of -20dB
  subroutine A0Comp(Rc, A0)
    real:: Rc, A0
    
    if (Rc .lt. 9.52E+04) then
      A0 = 0.57
    else if (Rc .lt. 8.57E+05) then
      A0 = (-9.57E-13)*(Rc-8.57E+05)**2. + 1.13
    else
      A0 = 1.13
    end if
  end subroutine A0Comp
  
!  This subroutine computes the high frequency 
!  directivity function for the input observer location  
  subroutine directH(M, Theta, Phi, Dbar)
    real:: M, Theta, Phi, Dbar
    real:: Mc

    Mc = 0.8*M
    
    Dbar = 2.*sin(Theta/2.)**2.*sin(Phi)**2./((1.+M*cos(Theta))*(1.+(M-Mc)*cos(Theta))**2.)
  end subroutine directH
  
!  This subroutine computes the low frequency 
!  directivity function for the input observer location  
  subroutine directL(M, Theta, Phi, Dbar)
    real:: M, Theta, Phi, Dbar
    
    Dbar = (sin(Theta)*sin(Phi))**2./(1.+M*cos(Theta))**4.
  end subroutine directL
  
  subroutine Blunt(SectAOA, chord, U, frCen, BLTrip, SPLblnt, Theta, Phi, &
                    L, R, TEthickness, TEflowAngle, M, RC)
                    
!               ----------------------------------
!               ****** Variable Definitions ******
!               ----------------------------------

!    Variable Name               Definition                     Units
!    -------------               ----------                     -----

!       SectAOA           segment angle of attack             Degrees
!       aTerm               used to compute peak strouhal no.   ---
!       chord               segment chord length                meters
!       c                   speed of sound                      meters/sec
!       DbarH               high frequency directivity          ---
!       deltaP              pressure side boundary layer
!                             thickness                         meters
!       DstarH              average displacement thickness
!                             over trailing edge bluntness      ---
!       dstrP               pressure side displacement thickness meters
!       dstrS               suction side displacement thickness meters
!                             layer thickness
!       eta                 ratio of strouhal numbers           ---
!       frCen               1/3 octave centered frequencies     hertz
!       F4temp              G5 evaluated at minimu1m HDstarP     dB
!       G4                  scaled spectrum level               dB
!       G5                  spectrum shape function             dB
!       G50                 G5 evaluated at TEflowAngle=0.0             dB
!       G514                G5 evaluated at TEflowAngle=14.0            dB
!       TEthickness         trailing edge bluntness             meters
!       HDstar              bluntness over average displacement
!                             thickness                         ---
!       HDstarL             minimu1m allowed value of HDstar     ---
!       HDstarP             modified value of HDstar            ---
!       BLTrip               trigger for boundary layer tripping ---
!       L                   span                                meters
!       M                   mach number                         ---
!       nFreq               number of 1/3 octave frequencies    ---
!       Phi                directivity angle                   degrees
!       TEflowAngle                 trailing edge angle                 degrees
!       R                   segment to observer distance        meters
!       Rc                  reynolds number based on chord      ---
!       scale               scaling factor                      ---
!       SPLblnt             sound pressure levels due to
!                           bluntness                         ---
!       StPeak              peak strouhal number                ---
!       StPPP               strouhal number                     ---
!       Theta              directivity angle                   degrees
!       U                   segment freestream velocity         meters/sec
!       nu                  kinematic viscosity                 m2/sec
    implicit none
    real, dimension(BB_THIRD_OCTAVE_BANDS):: SPLblnt, frCen, StPPP
    real:: SectAOA, chord, U, Theta, Phi, L, R, TEthickness, TEflowAngle, M, RC
    real:: deltaP, dstrS, dstrP, dstrAvg, DbarH, HDstar, DstarH, aTerm, StPeak 
    real:: eta, HDstarL, G514, HDstarP, G4, G5, G50, F4temp, scale
    integer:: BLTrip, i
    
!   Compute boundary layer thickness
!   --------------------------------
    call Thick(chord, SectAOA, BLTrip, deltaP, dstrS, dstrP, M, RC)    
    
!   Compute average displacement thickness
!   --------------------------------------
    dstrAvg = (dstrS + dstrP)/2.
    
    HDstar = TEthickness/dstrAvg
    DstarH = 1./HDstar
    
!   Compute directivity function
!   ----------------------------
    call directH(M, Theta, Phi, DbarH)
    
!   Compute peak strouhal number
!   ----------------------------
    aTerm = 0.212 - 0.0045*TEflowAngle
    
    if (HDstar .ge. 0.2) then
      StPeak = aTerm/(1.0 + 0.235*DstarH - 0.0132*DstarH**2.0)
    else
      StPeak = 0.1*HDstar + 0.095 - 0.00243*TEflowAngle
    end if
      
!   Compute scaled spectrum level
!   -----------------------------
    if (HDstar .le. 5.0) then
      G4 = 17.5*log10(HDstar) + 157.5 - 1.114*TEflowAngle
    else
      G4 = 169.7 - 1.114*TEflowAngle
    end if        
    
!   For each frequency, compute spectrum shape referenced to 0 dB
!   -------------------------------------------------------------
    do i=1,BB_THIRD_OCTAVE_BANDS
      StPPP(i) = frCen(i)*TEthickness/U
      eta      = log10(StPPP(i)/StPeak)
      
      HDstarL = HDstar      
      call G5comp(HDstarL, eta, G514)      
      HDstarP = 6.724*HDstar**2. - 4.019*HDstar + 1.107      
      call G5comp(HDstarP, eta, G50)
      
      G5 = G50 + 0.0714*TEflowAngle*(G514-G50)
      if (G5 .gt. 0.) G5 = 0.
      call G5comp(0.25, eta, F4temp)
      if (G5 .gt. F4temp) G5 = F4temp
      
      scale = 10.*log10(M**5.5*TEthickness*DbarH*L/R**2.0)
      SPLblnt(i) = G4 + G5 + scale
    end do
                    
  end subroutine Blunt
  
  subroutine G5comp(HDstar, eta, G5)
    implicit none
    real:: HDstar, eta, G5, M, K, mu1, eta0
    
    if (HDstar .lt. 0.25) then
      mu1 = 0.1211
    else if (HDstar .le. 0.62) then
      mu1 = -0.2175*HDstar + 0.1755
    else if (HDstar .lt. 1.15) then
      mu1 = -0.0308*HDstar + 0.0596
    else
      mu1 = 0.0242
    end if
    
    if (HDstar .le. 0.02) then
      M = 0.0
    else if (HDstar .lt. 0.5) then
      M = 68.724*HDstar - 1.35
    else if (HDstar .le. 0.62) then
      M = 308.475*HDstar - 121.23
    else if (HDstar .le. 1.15) then
      M = 224.811*HDstar - 69.354
    else if (HDstar .lt. 1.2) then
      M = 1583.28*HDstar - 1631.592
    else
      M = 268.344
    end if
    
    if (M .lt. 0.) M = 0.
    
    eta0 = -sqrt((M**2.)*(mu1**4.)/(6.25+(M*mu1)**2.))
    
    K = 2.5*sqrt(1.0 - (eta0/mu1)**2.) - 2.5 - M*eta0
    
    if (eta .le. eta0) then
      G5 = M*eta + K
    else if (eta .le. 0) then
      G5 = 2.5*sqrt(1.0 - (eta/mu1)**2.) - 2.5
    else if (eta .le. 0.03615995) then
      G5 = sqrt(1.5625 - 1194.99*eta**2.) - 1.25
    else
      G5 = -155.543*eta + 4.375
    end if
  end subroutine G5comp
  
  subroutine tipNoise(TipAOA, TipLCS, chord, U, frCen, SPLtip, Theta, Phi, &
                        R, RoundTip, M)

!               ----------------------------------
!               ****** Variable Definitions ******
!               ----------------------------------

!    Variable Name               Definition                     Units
!    -------------               ----------                     -----

!       TipAOA              tip angle of attack                 Degrees
!       TipLCS              tip lift curve slope                ---
!       alpTipP             corrected tip angle of attack       degrees
!       chord               segment chord length                meters
!       c                   speed of sound                      meters/sec
!       DbarH               directivity                         ---
!       frCen               1/3 octave centered frequencies     hertz
!       L                   segment span length                 meters
!       M                   mach number                         ---
!       Mm                  maximu1m mach number                 ---
!       nFreq               number of centered frequencies      ---
!       Phi                directivity angle                   degrees
!       R                   segment to observer distance        meters
!       RoundTip            logical indicated rounded tip       ---
!       scale               scaling term                        ---
!       SPLtip              sound pressure level associated
!                             with tip noise prediction         dB
!       StPP                strouhal number                     ---
!       term                scaling term                        ---
!       Theta              directivity angle                   degrees
!       U                   segment freestream velocity         meters/sec
!       Um                  maximu1m velocity                    meter/sec
!       nu                  kinematic viscosity                 m2/sec
    implicit none
    real, dimension(BB_THIRD_OCTAVE_BANDS):: SPLtip, frCen
    real:: TipAOA, TipLCS, Chord, U, L, M, Mm, Theta, Phi, R
    real:: alpTipP, DbarH, Um, term, scale, StPP, term2
    logical:: RoundTip
    
    integer:: i
    
    ! This is step is a very big problem
    ! if either the TipAOA or TipLCS is negative 
    ! since it directly leads to the log of a negative number.
    alpTipP = TipAOA*TipLCS
    
    call directH(M, Theta, Phi, DbarH)
    
    if (RoundTip) then
      L = 0.008*alpTipP*chord
    else
      if (abs(alpTipP) .le. 2.0) then
        L = (0.023 + 0.0169*alpTipP)*chord
      else
        L = (0.0378 + 0.0095*alpTipP)*chord
      end if
    end if
    
    Mm = (1.0 + 0.036*alpTipP)*M
    Um = Mm*c
    
    term = ((M*L/R)**2.)*(Mm**3.)*DbarH
    if (term .ne. 0.) then
      scale = 10.*log10(term)
    else
      scale = 0.
    end if
    
    do i=1,BB_THIRD_OCTAVE_BANDS
      StPP = frCen(i)*L/Um
      SPLtip(i) = 126.0 - 30.5*(log10(StPP+epsilon**3.) + 0.3)**2.0 + scale
    end do
  end subroutine tipNoise
  
  subroutine Thick(chord, SectAOA, BLTrip, deltaP, dstrS, dstrP, M, RC)
  
!               ----------------------------------
!               ****** Variable Definitions ******
!               ----------------------------------

!    Variable Name               Definition                     Units
!    -------------               ----------                     -----

!       SectAOA           segment angle of attack             Degrees
!       chord               segment chord length                meters
!       c                   speed of sound                      meters/sec
!       delta0              boundary layer thickness at
!                             zero angle of attack              meters
!       deltaP              pressure side boundary layer
!                             thickness                         meters
!       dstr0               displacement thickness at zero
!                             angle of attack                   meters
!       dstrP               pressure side displacement
!                             thickness                         meters
!       dstrS               suction side displacement 
!                             thickness                         meters                    
!       BLTrip               flag to trip boundary layer         ---
!       M                   mach number                         ---
!       Rc                  reynolds number based on chord      ---
!       U                   segment freestream velocity         meters/sec
!       nu                  kinematic viscosity                 m2/sec
    implicit none
    real:: chord, SectAOA, deltaP, dstrS, dstrP, M, RC
    real:: dstr0
    integer:: BLTrip
    real:: delta0

    delta0 = 10.**(1.6569 - 0.9045*log10(Rc) + 0.0596*log10(Rc)**2.)*chord
    if (BLTrip .eq. 2) delta0 = 0.6*delta0
    
!   Compute pressure side boundary layer thickness
!   ----------------------------------------------
    deltaP = 10.**(-0.04175*sectAOA + 0.00106*sectAOA**2.)*delta0
    
!   Compute zero angle of attack displacement thickness
!   ---------------------------------------------------
    if ((BLTrip .eq. 1) .or. (BLTrip .eq. 2)) then
      if (Rc .le. 0.3E+06) then
        dstr0 = 0.0601*Rc**(-0.114)*chord
      else
        dstr0 = 10.**(3.411 - 1.5397*log10(Rc) + 0.1059*log10(Rc)**2.)*chord
      end if
      if (BLTrip .eq. 2) dstr0 = dstr0*0.6
    else
      dstr0 = 10.**(3.0187 - 1.5397*log10(Rc) + 0.1059*log10(Rc)**2.)*chord
    end if
    
!   pressure side displacement thickness
!   ------------------------------------
    dstrP = 10.**(-0.0432*sectAOA + 0.00113*sectAOA**2.)*dstr0
    if (BLTrip .eq. 2) dstrP = dstrP*1.48
    
!   suction side displacement thickness
!   -----------------------------------
    if (BLTrip .eq. 1) then
      if (sectAOA .le. 5.) then
        dstrS = 10.**(0.0679*sectAOA)*dstr0        
      else if (sectAOA .le. 12.5) then
        dstrS = 0.381*10.**(0.1516*sectAOA)*dstr0
      else
        dstrS = 14.296*10.**(0.0258*sectAOA)*dstr0
      end if
    else
      if (sectAOA .le. 7.5) then
        dstrS = 10.**(0.0679*sectAOA)*dstr0
      else if (sectAOA .le. 12.5) then
        dstrS = 0.0162*10.**(0.3066*sectAOA)*dstr0
      else
        dstrS = 52.42*10.**(0.0258*sectAOA)*dstr0
      end if
    end if
  end subroutine Thick
  
  function ComputePeggTerm(PeggData, Term) result(computeTerm)
    implicit none
    type(PEGG):: PeggData
    integer:: term
    logical:: computeTerm
    computeTerm = .false.
    if (PeggData%computeTerm(term).eq.'COMPUTE') computeTerm = .true.
  
  end function ComputePeggTerm
  
  function getBBoffset(nTau, DataTimeOffset, period, time) result(realTime)
    implicit none
    integer, intent(in):: nTau, time
    real:: period, DataTimeOffset
    
    real:: p
    integer:: timeOffset, offsetTime, realTime
    
    timeOffset = INT((DataTimeOffset/period)*(nTau-1))
    offsetTime = time + timeOffset
    if (offsetTime <= 0) then
      p = (offsetTime-1)/(nTau-1)
      realTime = offsetTime + (abs(p)+1)*(nTau-1)
    else
      realTime = modulo (offsetTime, nTau-1)
    end if
    if (realTime == 0) realTime = nTau-1
    if (realTime == nTau) realTime = 1
  
  end function getBBoffset
  
  function getPeggTermKey(PeggData, Term, BBkey) result(key)
    implicit none
    type(PEGG):: PeggData
    integer:: Term, BBkey, key 
    if (PeggData%computeTerm(Term).ne.'FILEVALUE') then
      key = 1
    else
      key = BBkey
    end if  
  
  end function getPeggTermKey
  
  
  function getBPMTermKey(BPMData, Term, BBkey) result(key)
    implicit none
    type(BPM):: BPMData
    integer:: Term, BBkey, key 
    if (BPMData%computeTerm(Term).ne.'FILEVALUE') then
      key = 1
    else
      key = BBkey
    end if
  
  end function getBPMTermKey
  
  function getPeggTimetype(PeggData) result(timeType)
    implicit none
    type(PEGG):: PeggData
    integer:: timeType

    timeType = PeggData%timetype
    
  end function getPeggTimetype
  
  function getBPMTimetype(BPMData) result(timeType)
    implicit none
    type(BPM):: BPMData
    integer:: timeType

    timeType = BPMData%timetype
    
  end function getBPMTimetype 
  
  function GetPeggData(BBdata) result(PeggData)
    implicit none
    type(broadband), target:: BBdata
    type(PEGG), pointer:: PeggData
    
    PeggData => BBData%PeggData
  end function GetPeggData
  
  function GetBPMData(BBdata) result(BPMData)
    implicit none
    type(broadband):: BBdata
    type(BPM), pointer:: BPMData
    
    BPMData => BBData%BPMData
  end function GetBPMData
  
  subroutine DestroyBPMNamelistData(NamelistData)
    implicit none
    type(BPMNamelist):: NamelistData
    if (allocated(NamelistData%sectChord))   deallocate(NamelistData%sectChord)
    if (allocated(NamelistData%sectLength))  deallocate(NamelistData%sectLength)
    if (allocated(NamelistData%TEThickness)) deallocate(NamelistData%TEThickness)
    if (allocated(NamelistData%TEflowAngle)) deallocate(NamelistData%TEflowAngle)
    if (allocated(NamelistData%SectAOA))     deallocate(NamelistData%SectAOA)
    if (allocated(NamelistData%U))           deallocate(NamelistData%U)
  end subroutine DestroyBPMNamelistData
  
  subroutine DestroyBBdata(BBdata)
    implicit none
    type(Broadband):: BBdata
    if (associated(BBdata%PeggData)) then
      call DestroyPeggData(BBData%PeggData)
      deallocate(BBdata%PeggData)
    end if
    if (associated(BBdata%BPMData)) then
      call DestroyBPMData(BBData%BPMData)
      deallocate(BBdata%BPMData)
    end if
  end subroutine DestroyBBdata
  
  subroutine DestroyPeggData(PeggData)
    implicit none
    type(PEGG), pointer:: PeggData
    if (associated(PeggData%Tau))             	   deallocate(PeggData%tau)
    if (associated(PeggData%obsTimeArray))    	   deallocate(PeggData%obsTimeArray)
    if (associated(PeggData%hubAxis))         	   deallocate(PeggData%hubAxis)
    if (associated(PeggData%ClBar))                deallocate(PeggData%ClBar)
    if (associated(PeggData%nbContributions)) 	   deallocate(PeggData%nbContributions)
    if (associated(PeggData%Intensity))		   deallocate(PeggData%intensity)
  end subroutine DestroyPeggData
  
  subroutine DestroyBPMData(BPMData)
    implicit none
    type(BPM), pointer:: BPMdata
    
    if (associated(BPMData%tau))             	  deallocate(BPMData%tau)
    if (associated(BPMData%obsTimeArray))    	  deallocate(BPMData%obsTimeArray)
    if (associated(BPMData%nbContributions)) 	  deallocate(BPMData%nbContributions)
    if (associated(BPMdata%Intensity))		  deallocate(BPMData%intensity)
    if (associated(BPMData%SectChord))      	  deallocate(BPMData%SectChord)
    if (associated(BPMData%SectLength))     	  deallocate(BPMData%SectLength)
    if (associated(BPMData%TEThickness))     	  deallocate(BPMData%TEthickness)
    if (associated(BPMData%TEflowAngle))    	  deallocate(BPMData%TEflowAngle)
    if (associated(BPMData%Theta))         	  deallocate(BPMData%Theta)
    if (associated(BPMData%Phi))            	  deallocate(BPMData%Phi)
    if (associated(BPMData%U))               	  deallocate(BPMData%U)
    if (associated(BPMData%SectAOA))         	  deallocate(BPMData%SectAOA)
    if (associated(BPMData%TipLCS))          	  deallocate(BPMData%TipLCS)
  end subroutine
  
end module broadbandModule
