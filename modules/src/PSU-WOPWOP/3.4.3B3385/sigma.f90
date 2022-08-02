module SigmaModule
! PSU-WOPWOP 
! $Id: sigma.f90 3301 2013-10-30 01:43:19Z brentner $ 
! $LastChangedDate: 2013-10-29 21:43:19 -0400 (Tue, 29 Oct 2013) $
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
 
 
  use containerObject
  use observerContainerObject
  use observerObject
  use debugModule
  use constantsModule
  use IOModule
  implicit none

  private
  public :: WriteSigmaData

  ! File formats currently supported:
  !
  !  Structured:
  !    Plot3D
  !
  !  Unstructured:
  !    Fieldview unstructured (FV-UNS)

contains

  ! This function handles writing all of the sigma data in the appropriate file
  ! formats. It takes in the top-level container list and the base for the
  ! filenames.
  subroutine WriteSigmaData (containers, filename, obsCont)
    type(container), dimension(:), intent(in) :: containers
    character(len=*), intent(in) :: filename
    type(observerContainer), intent(inout) :: obsCont

    integer, dimension(:,:), allocatable :: dims

    integer::nbStructuredPatches, nbUnstructuredPatches,i, counter
    
    ! Calculate how many of each type of patches there are:
    nbStructuredPatches = 0
    nbUnstructuredPatches = 0
    do i=1, size(containers)
      call GetNbContainerSurfaces (containers(i), &
                                   nbStructuredPatches, &
                                   nbUnstructuredPatches)
    end do
    ! Write the structured data:
    if (nbStructuredPatches > 0) then
      if (sigmaStructuredFileType == FORMAT_PLOT3D) then
        ! We need to get the dimensions of all the structured patches:
        allocate(dims(nbStructuredPatches, 4))
        counter = 0
        do i=1, size(containers)
          call GetContainerStructuredDims(containers(i), dims, counter)
        end do
        call WritePlot3DStructured (containers, filename, nbStructuredPatches, &
                                    dims, obsCont)
        deallocate(dims)                                
      else
        call Warning ("Only Plot3D output is currently supported for "//&
                      "structured sigma data.", "No output will be created.")
      end if
    end if
    ! Write the unstructured data:
    if (nbUnstructuredPatches > 0) then
      if (sigmaUnstructuredFileType == FORMAT_FIELDVIEW) then
        call WriteFieldviewUnstructured(containers,filename,nbUnstructuredPatches)
      else
        call Warning ("Only Fieldview output is currently supported for "//&
                      "unstructured sigma data.", "No output will be created.")
      end if
    end if

  end subroutine WriteSigmaData

  
  subroutine WritePlot3DStructured (containers, filename, nbPatches, dims, obsCont)
    type(container), dimension(:) :: containers
    character(len=*), intent(in) :: filename
    integer :: nbPatches, xUnit, fnUnit, namUnit, stat
    integer, dimension(:,:) :: dims
    type(observerContainer), intent(inout) :: obsCont
    integer :: i
    ! See if there is anything to output:
    if (nbPatches == 0 .and. .not. observerSigmaFlag) then
      return
    end if    
    ! Create the observer sigma:
    if ((observerSigmaFlag.or.observerContainerSigmaFlag).and.(obsCont%nbObs==1)) then
      call EnlargeObsContainer(obsCont)
    end if
    ! Open the files
    xUnit = GetStreamNumber()
    fnUnit = GetStreamNumber()
    namUnit = GetStreamNumber()
    call CreateBinaryFile (xUnit, trim(filename)//".x", .false., stat)
    call CreateBinaryFile (fnUnit, trim(filename)//".fn", .false., stat)    
    open (unit=namUnit, file=trim(filename)//".nam", status='REPLACE')
    if (stat /= 0) then
      call Error ("Could not open sigma output files. Make sure the directory",&
                  "exists and that you can write to it and try again.")
      stop
    end if    
    ! Write the headers
    if (observerSigmaFlag.or.observerContainerSigmaFlag) nbPatches = nbPatches+1
    call WriteBinaryInteger (xUnit, nbPatches, stat)
    do i=1,size(dims,1)
     if(dims(i,1)==1 .and. dims(i,2)==1) then !This is for a point source
      call WriteBinaryInteger(xUnit,1,stat)
      call WriteBinaryInteger(xUnit,1,stat)  
     else 
      if( dims(i,1) /= 1 ) then  ! this is a mod compact sigma surface - I don't want them to be only a line
        call WriteBinaryInteger (xUnit, dims(i,1), stat)
      else
        call WriteBinaryInteger (xUnit, 2, stat)
      end if
      if( dims(i,2) /= 1 ) then  ! this is a mod compact sigma surface - I don't want them to be only a line
        call WriteBinaryInteger (xUnit, dims(i,2), stat)
      else
        call WriteBinaryInteger (xUnit, 2, stat)
      end if
     end if
      call WriteBinaryInteger (xUnit, dims(i,3), stat)
    end do
    if (observerSigmaFlag.or.observerContainerSigmaFlag) then
      call WriteBinaryInteger (xUnit, obsCont%NbObsDim1, stat)
      call WriteBinaryInteger (xUnit, obsCont%NbObsDim2, stat)
      call WriteBinaryInteger (xUnit, obsCont%obs(1)%nt, stat)
    end if
    call WriteBinaryInteger (fnUnit, nbPatches, stat)
    do i=1,size(dims,1)
      if(dims(i,1)==1 .and. dims(i,2)==1) then !This is for a point source
      call WriteBinaryInteger(fnUnit,1,stat)
      call WriteBinaryInteger(fnUnit,1,stat)  
     else  
      if( dims(i,1) /=1  ) then  ! mod to ensure compact sigma surfaces are not a line
        call WriteBinaryInteger (fnUnit, dims(i,1), stat)
      else
        call WriteBinaryInteger (fnUnit, 2, stat)
      end if
      if( dims(i,2) /=1  ) then  ! mod to ensure compact sigma surfaces are not a line
        call WriteBinaryInteger (fnUnit, dims(i,2), stat)
      else
        call WriteBinaryInteger (fnUnit, 2, stat)
      end if   
     end if    
      call WriteBinaryInteger (fnUnit, dims(i,3), stat)
      call WriteBinaryInteger (fnUnit, dims(i,4)-3, stat)
    end do
    if (observerSigmaFlag.or.observerContainerSigmaFlag) then
      call WriteBinaryInteger (fnUnit, obsCont%NbObsDim1, stat)
      call WriteBinaryInteger (fnUnit, obsCont%NbObsDim2, stat)
      call WriteBinaryInteger (fnUnit, obsCont%obs(1)%nt, stat)
      call WriteBinaryInteger (fnUnit, dims(1,4)-3, stat) ! Not quite correct
    end if

    call Message ("Writing structured sigma data:")
    ! Loop over the zones and output them
    do i=1, size(containers)
      call WriteContainerStructuredSigma (containers(i), xUnit, fnUnit,"  ")
    end do
    
    ! Write the observer sigma:
    if(observerSigmaFlag.or.observerContainerSigmaFlag) then
      call WriteObsContSigmaSurface (obsCont,dims(1,4), xUnit, fnUnit)
    end if

    ! Write the namelist:
    call WriteSigmaNameList (namUnit)
    close (namUnit)
    call CloseBinaryFile (xUnit)
    call CloseBinaryFile (fnUnit)
    
  end subroutine WritePlot3DStructured

  
  subroutine WriteSigmaNameList(namUnit) 
    integer, intent(in) :: namUnit
 
    write(namUnit,*) 'Source Time' 
    write(namUnit,*) 'Observer Time' 
    if(machSigmaFlag)           write(namUnit,*) 'Mach Number' 
    if(loadingNoiseSigmaFlag)   write(namUnit,*) 'Loading Noise' 
    if(thicknessNoiseSigmaFlag) write(namUnit,*) 'Thickness Noise' 
    if(totalNoiseSigmaFlag)     write(namUnit,*) 'Total Noise' 
    if(normalSigmaFlag) then 
      write(namUnit,*) 'Nx; Normal vector' 
      write(namUnit,*) 'Ny' 
      write(namUnit,*) 'Nz' 
    end if 
    if(velocitySigmaFlag) then 
      write(namUnit,*) 'Vx; Velocity vector' 
      write(namUnit,*) 'Vy' 
      write(namUnit,*) 'Vz' 
    end if 
    if(accelerationSigmaFlag) then 
      write(namUnit,*) 'Ax; Acceleration vector' 
      write(namUnit,*) 'Ay' 
      write(namUnit,*) 'Az' 
    end if 
    if(loadingSigmaFlag) then 
      write(namUnit,*) 'Lnx; Loading vector' 
      write(namUnit,*) 'Lny' 
      write(namUnit,*) 'Lnz' 
    end if 
    if(densitySigmaFlag) then 
      write(namUnit,*) 'Density' 
    end if 
    if(momentumSigmaFlag) then 
      write(namUnit,*) 'MomentumX; Momentum vector' 
      write(namUnit,*) 'MomentumY' 
      write(namUnit,*) 'MomentumZ' 
    end if 
    if(pressureSigmaFlag) then 
      write(namUnit,*) 'Pressure' 
    end if 
    if(areaSigmaFlag) then
      write(namUnit,*) 'Area'
    end if
    if(MdotrSigmaFlag) then
      write(namUnit,*) 'Mdotr'
    end if
    if(iblankSigmaFlag) then
      write(namUnit,*) 'iblankVal'
    end if
     
  end subroutine writeSigmaNameList

  
  subroutine WriteFieldviewUnstructured (containers, filename, nbPatches)
    type(container), dimension(:) :: containers
    character(len=*), intent(in) :: filename
    integer :: nbPatches, fnUnit, stat, i, patchNumber
    character(len=80) :: tempString
    integer, parameter:: FV_MAGIC=66051,&
                         FV_VARIABLES=1004,FV_BNDRY_VARIABLES=1006
    real(kind=4):: temp
    
    ! See if there is anything to output:
    if (nbPatches == 0) then
      return
    end if
    fnUnit = GetStreamNumber()
    ! Open the file:
    call CreateBinaryFile (fnUnit, trim(filename)//".uns", .false., stat)
    if (stat /= 0) then
      call Error ("Could not open sigma output files. Make sure the directory",&
                  "exists and that you can write to it and try again.", &
                  "File:"//trim(filename)//".uns")
      stop
    end if

    ! Write the magic number:
    call WriteBinaryInteger (fnUnit, FV_MAGIC, stat)
    
    ! Write the header string:
    tempString = "FIELDVIEW"
    call WriteBinaryString (fnUnit, tempString, stat)

    ! Write the format version number (really two integers)
    call WriteBinaryInteger (fnUnit, 2, stat)
    call WriteBinaryInteger (fnUnit, 5, stat)

    ! Write zeros for the time, fsmach, alpha and re
    temp = 0.0
    call WriteBinaryReal (fnUnit, temp, stat)
    call WriteBinaryReal (fnUnit, temp, stat)
    call WriteBinaryReal (fnUnit, temp, stat)
    call WriteBinaryReal (fnUnit, temp, stat)

    ! Write the number of grids
    call WriteBinaryInteger (fnUnit, nbPatches, stat)
    
    ! Write the number of boundaries (equal to the number of grids here)
    !call WriteBinaryInteger (fnUnit, nbPatches, stat)
    call WriteBinaryInteger (fnUnit, 0, stat)

    ! Write out the boundary names: these are basically just the zone names
    !do i=1,size(containers)
    !  call WriteContainerFVUNSHeader (containers(i), fnUnit)
    !end do

    ! Write the variables
    call WriteSigmaFVUNSVariables (fnUnit) ! Volume variables    \  same
    !call WriteSigmaFVUNSVariables (fnUnit) ! Boundary variables  /  values
    call WriteBinaryInteger (fnUnit, 0, stat)
    
    ! Finally, write out the zone data:
    patchNumber = 0
    call Message ("Writing unstructured sigma data:")
    do i=1,size(containers)
      call WriteContainerFVUNSData (containers(i), fnUnit, patchNumber,"  ")
    end do

    ! Close the file
    call CloseBinaryFile (fnUnit)

  end subroutine WriteFieldviewUnstructured 
  
  
  subroutine WriteSigmaFVUNSVariables(namUnit)
    integer, intent(in) :: namUnit
    integer :: numberFunctions, stat

    ! Count up the functions
    numberFunctions = 2 ! Source time and observer time
    if(loadingNoiseSigmaFlag)   numberFunctions = numberFunctions + 1
    if(thicknessNoiseSigmaFlag) numberFunctions = numberFunctions + 1
    if(totalNoiseSigmaFlag)     numberFunctions = numberFunctions + 1
    if(normalSigmaFlag)         numberFunctions = numberFunctions + 3
    if(machSigmaFlag)           numberFunctions = numberFunctions + 1
    if(velocitySigmaFlag)       numberFunctions = numberFunctions + 3
    if(accelerationSigmaFlag)   numberFunctions = numberFunctions + 3
    if(loadingSigmaFlag)        numberFunctions = numberFunctions + 3
    if(densitySigmaFlag)        numberFunctions = numberFunctions + 1
    if(momentumSigmaFlag)       numberFunctions = numberFunctions + 3
    if(pressureSigmaFlag)       numberFunctions = numberFunctions + 1
    if(areaSigmaFlag)           numberFunctions = numberFunctions + 1
    if(MdotrSigmaFlag)          numberFunctions = numberFunctions + 1
    if(iblankSigmaFlag)         numberFunctions = numberFunctions + 1
 
    ! Write the number of functions
    call WriteBinaryInteger (namUnit, numberFunctions, stat)
    
    call Write80CharString (namUnit, 'Source Time', stat)
    call Write80CharString (namUnit, 'Observer Time' , stat)
    if(machSigmaFlag)           call Write80CharString (namUnit, 'Mach Number' , stat)
    if(loadingNoiseSigmaFlag)   call Write80CharString (namUnit, 'Loading Noise' , stat)
    if(thicknessNoiseSigmaFlag) call Write80CharString (namUnit, 'Thickness Noise' , stat)
    if(totalNoiseSigmaFlag)     call Write80CharString (namUnit, 'Total Noise' , stat)
    if(normalSigmaFlag) then 
      call Write80CharString (namUnit, 'Nx;Normal vector' , stat)
      call Write80CharString (namUnit, 'Ny' , stat)
      call Write80CharString (namUnit, 'Nz' , stat)
    end if 
    if(velocitySigmaFlag) then 
      call Write80CharString (namUnit, 'Vx;Velocity vector' , stat)
      call Write80CharString (namUnit, 'Vy' , stat)
      call Write80CharString (namUnit, 'Vz' , stat)
    end if 
    if(accelerationSigmaFlag) then 
      call Write80CharString (namUnit, 'Ax;Acceleration vector' , stat)
      call Write80CharString (namUnit, 'Ay' , stat)
      call Write80CharString (namUnit, 'Az' , stat)
    end if 
    if(loadingSigmaFlag) then 
      call Write80CharString (namUnit, 'Lnx;Loading vector' , stat)
      call Write80CharString (namUnit, 'Lny' , stat)
      call Write80CharString (namUnit, 'Lnz' , stat)
    end if 
    if(densitySigmaFlag) then 
      call Write80CharString (namUnit, 'Density' , stat)
    end if 
    if(momentumSigmaFlag) then 
      call Write80CharString (namUnit, 'MomentumX;Momentum vector' , stat)
      call Write80CharString (namUnit, 'MomentumY' , stat)
      call Write80CharString (namUnit, 'MomentumZ' , stat)
    end if 
    if(pressureSigmaFlag) then 
      call Write80CharString (namUnit, 'Pressure' , stat)
    end if 
    if(areaSigmaFlag) then
      call Write80CharString (namUnit, 'Area', stat)
    end if
    if(MdotrSigmaFlag) then
      call Write80CharString (namUnit, 'Mdotr', stat)
    end if
    if(iblankSigmaFlag) then
       call Write80CharString (namUnit, 'IblankVal', stat)
    end if 
     
  end subroutine WriteSigmaFVUNSVariables

  
  subroutine Write80CharString (u, str, stat)
    integer, intent(in) :: u
    integer, intent(out) :: stat
    character(len=*), intent(in) :: str
    character(len=80) :: newStr
    write (newStr,'(A80)') trim(adjustl(str))
    newStr = adjustl(newStr)
    call WriteBinaryString (u, newStr, stat)
  end subroutine Write80CharString

end module SigmaModule
