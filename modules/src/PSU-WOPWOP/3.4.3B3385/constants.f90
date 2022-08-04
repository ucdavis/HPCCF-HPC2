module constantsModule
! PSU-WOPWOP
! $Id: constants.f90 3371 2017-08-14 13:29:16Z brentner $
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


!*****************************************************
!TYPE: Constants module
!  This module encapsulates all environmental constants that
! need to be used on all levels of the program.  There are two 
! functions here for setting constants. The first is 
! setEnvironmentalConstants which depends on user input, the second
! is setConstants which does not depend on any outside information.
! This was necessary because the grid generation program needed
! some of these constants but would have no user input.
!******************************************************

  use MPIModule
  use debugModule
  implicit none

  integer::integrationType, sigmaStructuredFileType, sigmaUnstructuredFileType
  integer,parameter:: aperiodicArraySize=4
  !ksb debug real::rho, c, Re, M, nu, pi, gamma, deltaT, mh, epsilon, P_ref, rad2deg
  real::rho, c, nu, pi, gamma, deltaT, epsilon, P_ref, rad2deg
  real::gasConstant, Temperature, RelHumidity
  character(len=4096)::globalFolderName
  
  integer            :: THICK_APTH,      &
                        LOAD_APTH,       &
                        BB_PEGG_APTH,    &
                        BB_BPM_APTH,     &
                        TOTAL_APTH,      &
                        THICK_PGX,       &
                        THICK_PGY,       &
                        THICK_PGZ,       &
                        LOAD_PGX,        &
                        LOAD_PGY,        &
                        LOAD_PGZ,        &
                        TOTAL_PGX,       &
                        TOTAL_PGY,       &
                        TOTAL_PGZ
                        
  character(len=4096),dimension(:),allocatable::PPRIMETITLEARRAY                        
                        
  ! File format enumeration
  integer, parameter :: FORMAT_PLOT3D    = 100, &
                        FORMAT_FIELDVIEW = 101, &
                        FORMAT_TECPLOT   = 102
  
  !File Version 1.0 enumerations
  integer, parameter:: GEO_CHECK        = 1, &
                       GEO_NBZONES      = 2, &
                       GEO_GRIDTYPE     = 3, &
                       GEO_TIMETYPE     = 4, &
                       GEO_NORMTYPE     = 5, &
                       GEO_PRECISION    = 6, &
                       GEO_IBLANK       = 7, &
                       GEO_FUTURE2      = 8, &
                       LOAD_CHECK       = 1, &
                       LOAD_NBZONES     = 2, &
                       LOAD_GRIDTYPE    = 3, &
                       LOAD_TIMETYPE    = 4, &
                       LOAD_NORMTYPE    = 5, &
                       LOAD_DATATYPE    = 6, &
                       LOAD_REFFRAME    = 7, &
                       LOAD_PRECISION   = 8, &
                       LOAD_FUTURE1     = 9, &
                       LOAD_FUTURE2     = 10,&
                       GEO_SIZE         = 8, &
                       LOAD_SIZE        = 10

  integer, parameter:: GEO_GRIDTYPE_STRUCTURED   = 1, &
                       GEO_GRIDTYPE_UNSTRUCTURED = 2, &
                       GEO_TIMETYPE_CONSTANT     = 1, &
                       GEO_TIMETYPE_PERIODIC     = 2, &
                       GEO_TIMETYPE_APERIODIC    = 3, &
                       GEO_TIMETYPE_MTFAPERIODIC = 4, &
                       GEO_TIMETYPE_QPERIODIC    = 5, &
                       GEO_TIMETYPE_MTFQPERIODIC = 6, &
                       GEO_PRECISION_SINGLE      = 1, &
                       GEO_PRECISION_DOUBLE      = 2, &
                       GEO_NORMTYPE_NODES        = 1, &
                       GEO_NORMTYPE_CELLS        = 2

  integer, parameter:: LOAD_GRIDTYPE_STRUCTURED   = 1, &
                       LOAD_GRIDTYPE_UNSTRUCTURED = 2, &
                       LOAD_TIMETYPE_CONSTANT     = 1, &
                       LOAD_TIMETYPE_PERIODIC     = 2, &
                       LOAD_TIMETYPE_APERIODIC    = 3, &
                       LOAD_TIMETYPE_MTFAPERIODIC = 4, &
                       LOAD_TIMETYPE_QPERIODIC    = 5, &
                       LOAD_TIMETYPE_MTFQPERIODIC = 6, &                       
                       LOAD_PRECISION_SINGLE      = 1, &
                       LOAD_PRECISION_DOUBLE      = 2, &
                       LOAD_DATATYPE_PRESSURE     = 1, &
                       LOAD_DATATYPE_LOADING      = 2, &
                       LOAD_DATATYPE_FLOW         = 3, &
                       LOAD_REFFRAME_GROUND       = 1, &
                       LOAD_REFFRAME_MIXED        = 2, &
                       LOAD_REFFRAME_BLADE        = 3, &
                       LOAD_REFFRAME_BLADE_SCALEV = 4, &
                       LOAD_NORMTYPE_NODES        = 1, &
                       LOAD_NORMTYPE_CELLS        = 2
  integer, parameter:: VERSION_0                = 0, &
                       VERSION_1                = 1

  !Integration type enumerations
  integer, parameter:: TIME_DOMAIN                = 0,&
                       FREQUENCY_DOMAIN           = 1
  !integer, parameter:: SURFACE                    = 1,&
  !                     LOADING                    = 2   
  integer, parameter:: THICKNESS_OUT             = 1, &
                       LOADING_OUT               = 2

  ! Broadband parameters
 integer, parameter::  BB_THIRD_OCTAVE_BANDS = 27, &
                       
                       PEGG_OCTAVE_BANDS     = 9, &
                       PEGG_INPUT_DATA       = 6, &
                       PEGG_TOTALBLADEAREA   = 1, &
                       PEGG_BLADERADIUS      = 2, &
                       PEGG_ROTSPEED         = 3, &
                       PEGG_CLBAR            = 4, &
                       PEGG_HUBAXIS          = 5, &
                       PEGG_TOTALTHRUST      = 6, &
                       
                       BPM_INPUT_DATA        = 8, &
                       BPM_SECTCHORD         = 1, &
                       BPM_SECTLENGTH        = 2, &
                       BPM_TETHICKNESS       = 3, &
                       BPM_TEFLOWANGLE       = 4, &
                       BPM_TIPLCS            = 5, &                       
                       BPM_SECTAOA           = 6, &
                       BPM_TIPAOA            = 7, &
                       BPM_U                 = 8, &
                       BPM_COMPUTE_FLAGS     = 5, &
                       BPM_LBLVSNOISE        = 1, &
                       BPM_TBLTENOISE        = 2, &
                       BPM_BLUNTNOISE        = 3, &
                       BPM_BLADETIPNOISE     = 4, &
                       BPM_ROUNDBLADETIP     = 5, &
                       
                       BPM_TERMS             = 7, &
                       BPM_SPLP              = 1, &
                       BPM_SPLS              = 2, &
                       BPM_SPLALPHA          = 3, &
                       BPM_SPLLBLVS          = 4, &
                       BPM_SPLBLUNT          = 5, &
                       BPM_SPLTIP            = 6, &
                       BPM_SPLTOTAL          = 7
                           
  ! Calculation flags\
  logical::thicknessNoiseFlag,      & 
           loadingNoiseFlag,        &
           totalNoiseFlag,          &
           PressureGradientFlag,    &
           PressureGradient1AFlag,  &           
           segmentFlag,             &
           readInPressureFlag,      &
           globalPeggNoiseFlag,     &
           globalBPMNoiseFlag

  ! Sigma surface flags
  logical::sigmaFlag,               &
           loadingNoiseSigmaFlag,   &
           thicknessNoiseSigmaFlag, &
           isomsThicknessNoiseFlag, &
           normalSigmaFlag,         &
           machSigmaFlag,           &
           totalNoiseSigmaFlag,     &
           observerSigmaFlag,       &
           observerContainerSigmaFlag,&
           velocitySigmaFlag,       &
           accelerationSigmaFlag,   &
           loadingSigmaFlag,        &
           densitySigmaFlag,        &
           momentumSigmaFlag,       &
           pressureSigmaFlag,       &
           areaSigmaFlag,           &
           MdotrSigmaFlag,          &
           iblankSigmaFlag

  ! Other input/output flags
  logical::OASPLdBFlag,             &
           OASPLdBAFlag,            &
           OASPLFlag,               &
           SPLdBFlag,               &
           SPLdBAFlag,              &
           SPLFlag,                 &
           iBlankFlag,              &
           acousticPressureFlag,    &
           broadbandFlag,           &
           narrowbandFlag,          &
           audioFlag,               &
           spectrumFlag,            &
           phaseFlag,               &
           octaveApproxFlag,        &
           PNLFlag,                 &
           PNLTFlag,                &
           PartH36Flag,             &
           EPNLFlag,                &
           SELFlag,                 &
           octaveFlag,              &         
           complexPressureFlag,     &
           FSCOutputFlag,           &
           FSCThicknessFlag,        &
           FSCLoadingFlag,          &
           FAACertFlag,             &
           newSegLHS,               &
           newSegRHS,               &
           dataLimitLHS,            &
           dataLimitRHS,            &
           forceEPNL,               &
           forceSEL,                &
           HitMaxExpansions,        &
           EPNLPrime,               &
           AtmAtten
           
  logical::ASCIIOutputFlag
  ! Check for whether the case meets
  ! the PNLTM-10dB and LAmax-10dB criteria for FAA cert.
  ! It will default to true until 
  ! we've determined that it HASN'T been met.
  logical, dimension(2)::PNLTdBDrop, &
                         SELdBDrop
                         

  
  ! FFTW 3.3.4 constants
      INTEGER FFTW_R2HC
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_WISDOM_ONLY
      PARAMETER (FFTW_WISDOM_ONLY=2097152)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_NO_DFT_R2HC
      PARAMETER (FFTW_NO_DFT_R2HC=512)
      INTEGER FFTW_NO_NONTHREADED
      PARAMETER (FFTW_NO_NONTHREADED=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)
      INTEGER FFTW_NO_SLOW
      PARAMETER (FFTW_NO_SLOW=262144)
      INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
      PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
      INTEGER FFTW_ALLOW_PRUNING
      PARAMETER (FFTW_ALLOW_PRUNING=1048576)
contains

subroutine BroadcastEnvironmentalConstants()
!ksb debug: I don't think this subroutine is used.
 
  call BroadcastInteger(integrationType)
  call BroadcastInteger(debugLevel)
  call BroadcastLogical(sigmaFlag)  
  call BroadcastLogical(OASPLdBFlag)
  call BroadcastLogical(OASPLdBAFlag)
  call BroadcastLogical(SPLdBFlag)
  call BroadcastLogical(SPLdBAFlag)
  call BroadcastLogical(spectrumFlag)
  call BroadcastLogical(audioFlag)
  call BroadcastLogical(phaseFlag)
  call BroadcastLogical(iBlankFlag)
  call BroadcastLogical(acousticPressureFlag)
  call BroadcastLogical(broadbandFlag)
  call BroadcastLogical(narrowbandFlag)
  call BroadcastLogical(thicknessNoiseFlag)
  call BroadcastLogical(loadingNoiseFlag)
  call BroadcastLogical(totalNoiseFlag)
  call BroadcastReal   (rho)
  call BroadcastReal   (c)
  call BroadcastReal   (pi)
  call BroadcastReal   (gamma)
  call BroadcastReal   (nu)
  call BroadcastReal   (epsilon)
  call BroadcastReal   (P_ref)
  call BroadcastReal   (RelHumidity)
  call BroadcastReal   (gasConstant)
  
end subroutine BroadcastEnvironmentalConstants
     
!*****************************************************************************************
!SUBROUTINE setEvironmentalConstants(unitNumber)
!  This sets all the environemtal constants including rho and c 
! as well as time parameters
!ARGUMENTS:
!  - unitNumber:  cooresponding to the namelist
!*****************************************************************************************
subroutine SetEnvironmentalConstants(unitNumber)
  integer, intent(in)::unitNumber
  real:: nu0, mu
  namelist / EnvironmentConstants / rho, c, gamma, mu, nu, integrationType, epsilon, P_ref, &
  &  RelHumidity, gasConstant
  ! Integration Type
  ! 1 = frequency domain
  ! 0 = time domain (default)
  ! Set the default values
  gasConstant     = 287.04  ! J kg^-1 K-1
  RelHumidity     = 70.0
  rho             = 1.225
  c               = 342.0
  gamma           = 1.4
  integrationType = 0
  pi              = 4.0*atan(1.) !3.14159265359
  rad2deg         = 180.0/pi
  mu              = 0.0000151
  deltaT          = 0.5
  nu0             = mu/rho; nu = nu0
  epsilon         = 1.e-5
  P_ref           = 20.e-6
  read(unitNumber, nml=EnvironmentConstants)
  
  if ( c == 0 .or. rho == 0 .or. mu == 0 .or. nu == 0) then
    call Error('Fluid parameters are incorrect. Stopping.', &
               'c = '//RealToString(c)//' , rho  = '//RealToString(rho), &
               'mu = '//RealToString(mu)//' , nu = '//RealToString(nu)//'')
    stop
  end if
  !
  ! Consistency check for mu, nu, and rho
  !
  if( abs(nu0 - nu) > epsilon ) then
      ! nu has changed through namelist input
      if( abs(mu/rho - nu) > epsilon ) then
          call Warning('Value of mu is inconsistent with the value of nu.', &
                       'Using nu = '//RealToString(nu)//'' )
      end if
  else
      if( abs(mu/rho - nu) > epsilon ) then
          ! mu or rho has changed, but nu was not specified.  Update nu here.
          nu = mu/rho
      end if
  end if
      

  Temperature = (c**2.)/(gamma*gasConstant) !normally Kelvin

end subroutine SetEnvironmentalConstants

subroutine ResetPPrimeIndices()
  integer::i !This keeps track of the index
  if (.not.allocated(PPRIMETITLEARRAY)) then
    allocate(PPRIMETITLEARRAY(GetNewSizeOfPPrime()))
  end if
  i = 1
  if(totalNoiseFlag)then
    if((.not.thicknessNoiseFlag).and.(.not.loadingNoiseFlag))then
      !ksb debug: 12/2/16
      TOTAL_APTH = 1
      TOTAL_PGX  = 2
      TOTAL_PGY  = 3
      TOTAL_PGZ  = 4
    else if(loadingNoiseFlag.and.(.not.thicknessNoiseFlag))then
      !ksb debug: 12/2/26
      LOAD_APTH  = 1
      TOTAL_APTH = 2
      LOAD_PGX   = 3
      LOAD_PGY   = 4
      LOAD_PGZ   = 5
      TOTAL_PGX  = 6
      TOTAL_PGY  = 7
      TOTAL_PGZ  = 8
    else if(thicknessNoiseFlag.and.(.not.loadingNoiseFlag))then
      THICK_APTH = 1
      TOTAL_APTH = 2
      THICK_PGX  = 3
      THICK_PGY  = 4
      THICK_PGZ  = 5
      TOTAL_PGX  = 6
      TOTAL_PGY  = 7
      TOTAL_PGZ  = 8
     
      !ksb debug: 12/2/16 - I'm not sure why there is not branch for if both thickness AND loading are .true.
    end if
  end if
end subroutine ResetPPrimeIndices          

subroutine SetPPrimeIndices()
  integer::i !This keeps track of the index
  if (.not.allocated(PPRIMETITLEARRAY)) then
    allocate(PPRIMETITLEARRAY(GetSizeOfPPrime()))
  end if
  i=0  
  if(thicknessNoiseFlag.or.totalNoiseFlag)then
    i=i+1
    THICK_APTH=i
    PPRIMETITLEARRAY(i)='Thickness'
  end if
  if(loadingNoiseFlag.or.totalNoiseFlag)then
    i=i+1
    LOAD_APTH=i
    PPRIMETITLEARRAY(i)='Loading'
  end if
  if(totalNoiseFlag)then
    i=i+1
    TOTAL_APTH=i
    PPRIMETITLEARRAY(i)='Total'
  end if
  if(PressureGradientFlag.or.PressureGradient1AFlag)then
    if(thicknessNoiseFlag.or.totalNoiseFlag)then
      i=i+1
      THICK_PGX=i
      PPRIMETITLEARRAY(i)='ThicknessPGX'
      i=i+1
      THICK_PGY=i
      PPRIMETITLEARRAY(i)='ThicknessPGY'
      i=i+1
      THICK_PGZ=i
      PPRIMETITLEARRAY(i)='ThicknessPGZ'
    end if
    if(loadingNoiseFlag.or.totalNoiseFlag)then
      i=i+1
      LOAD_PGX=i
      PPRIMETITLEARRAY(i)='LoadingPGX'
      i=i+1
      LOAD_PGY=i
      PPRIMETITLEARRAY(i)='LoadingPGY'      
      i=i+1
      LOAD_PGZ=i
      PPRIMETITLEARRAY(i)='LoadingPGZ'
    end if
    if(totalNoiseFlag)then
      i=i+1
      TOTAL_PGX=i
      PPRIMETITLEARRAY(i)='TotalPGX'
      i=i+1
      TOTAL_PGY=i
      PPRIMETITLEARRAY(i)='TotalPGY'
      i=i+1
      TOTAL_PGZ=i
      PPRIMETITLEARRAY(i)='TotalPGZ'      
    end if
  end if  
end subroutine SetPPrimeIndices

function GetSizeOfPPrime() result(i)
  integer::i
  i=0
  if(thicknessNoiseFlag.or.totalNoiseFlag)then
    i=i+1
  end if
  if(loadingNoiseFlag.or.totalNoiseFlag)then
    i=i+1
  end if
  if(totalNoiseFlag)then
    i=i+1
  end if
  if(PressureGradientFlag.or.PressureGradient1AFlag)then
    if(thicknessNoiseFlag.or.totalNoiseFlag)then
      i=i+3
    end if
    if(loadingNoiseFlag.or.totalNoiseFlag)then
      i=i+3
    end if
    if(totalNoiseFlag)then
      i=i+3
    end if    
  end if  
end function GetSizeOfPPrime

function GetNewSizeOfPPrime() result(i)
  integer::i
  i=0
  if(thicknessNoiseFlag)then
    i=i+1
  end if
  if(loadingNoiseFlag)then
    i=i+1
  end if
  if(totalNoiseFlag)then
    i=i+1
  end if
  if(PressureGradientFlag.or.PressureGradient1AFlag)then
    if(thicknessNoiseFlag)then
      i=i+3
    end if
    if(loadingNoiseFlag)then
      i=i+3
    end if
    if(totalNoiseFlag)then
      i=i+3
    end if    
  end if  
end function GetNewSizeOfPPrime

end module ConstantsModule
