!!!!!!EXTRA

   ! id =bmif%bmif_create(nxyz, nthreads)
   ! print *, "debug1"
	! status = RM_InitializeYAML(id, yaml_file)     ! Deprecated, but working
   status = bmif%bmif_initialize(yaml_file)  ! but this one not working!! Now working! 09/10/2024
   print *, "initalize"
   ! nxyz = 40
   ! print *, "debug2"


!!!   Set properties
   ! status = RM_SetErrorOn(id, 1)
	status = bmif%SetErrorOn(1)
   ! status = RM_SetErrorHandlerMode(id, 2)  ! exit on error
	status = bmif%SetErrorHandlerMode(1)
   ! status = RM_SetComponentH2O(id, 0)
   ! status = RM_SetRebalanceFraction(id, 0.5d0)
	status =bmif%SetRebalanceFraction(0.5d0)
   ! status = RM_SetRebalanceByCell(id, 1)
	status =bmif%SetRebalanceByCell(1)
   ! status = RM_UseSolutionDensityVolume(id, 0)
	status = bmif%UseSolutionDensityVolume(0)
   ! status = RM_SetPartitionUZSolids(id, 0) #not implemented

!!! Set concentration units
   ! status = RM_SetUnitsSolution(id, 2)      ! 1, mg/L; 2, mol/L; 3, kg/kgs
	status =bmif%SetUnitsSolution(2)
   ! status = RM_SetUnitsPPassemblage(id, 1)  ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status =bmif%SetUnitsSSassemblage(1)
   ! status = RM_SetUnitsExchange(id, 1)      ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
   ! status = RM_SetUnitsSurface(id, 1)       ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
   ! status = RM_SetUnitsGasPhase(id, 1)      ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
   ! status = RM_SetUnitsSSassemblage(id, 1)  ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
   ! status = RM_SetUnitsKinetics(id, 1)      ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status =bmif%SetUnitsKinetics(1)

!!! Set conversion from seconds to user units (days)
   ! status = RM_SetTimeConversion(id, dble(1.0 / 86400.0))
	t = 1.0 / 86400.0
	! status = RM_SetTimeConversion(id, t)     ! Deprecated
	status = bmif%SetTimeConversion(t)
   
!!! Set representative volume
   allocate(rv(nxyz))
   rv = 1.0
   ! status = RM_SetRepresentativeVolume(id, rv)
   ! DoubleVector = 1.0d0
	! status = RM_SetRepresentativeVolume(id, DoubleVector)     ! Deprecated
	status = bmif%SetRepresentativeVolume(rv)
   
!!! Set initial porosity ! For transport? No need right?
   allocate(por(nxyz))
   ! por = 0.2
   ! status = RM_SetPorosity(id, por)
	status = bmif%SetPorosity(rv)

!!! Set initial saturation
   allocate(sat(nxyz))
   sat = 1.0
   ! status = RM_SetSaturationUser(id, sat)
   status = bmif%SetSaturationUser(sat)

!!! --------------------------------------------------------------------------
!!! Set initial conditions
!!! --------------------------------------------------------------------------

!!! Set printing of chemistry file to false
   ! status = RM_SetPrintChemistryOn(id, 0, 1, 0)  ! workers, initial_phreeqc, utility
   status = bmif%SetPrintChemistryOn(0, 1, 0)  ! workers, initial_phreeqc, utility

!!! Load database
   ! status = RM_LoadDatabase(id, "aed/" // dbase) 
   status = bmif%LoadDatabase("aed/" // dbase)

!!! Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
  ! CALL register_basic_callback_fortran()

!!! Run file to define solutions and reactants for initial conditions, selected output
  ! There are three types of IPhreeqc instances in PhreeqcRM
  ! Argument 1 refers to the worker IPhreeqcs for doing reaction calculations for transport
  ! Argument 2 refers to the InitialPhreeqc instance for accumulating initial and boundary conditions
  ! Argument 3 refers to the Utility instance available for processing

!!! Run file to define solutions and reactants for initial conditions, selected output
!   status = RM_RunFile(id, 1, 1, 1, "advect.pqi")
!   status = bmif%RunFile(1, 1, 1, "advect.pqi")

!!! Define Phreeqc string to clear contents of workers and utility
  string = "DELETE; -all"

!!! Run string to clear contents of workers and utility
!   status = RM_RunString(id, 1, 0, 1, string)  ! workers, initial_phreeqc, utility
  status = bmif%RunString(1, 0, 0, string)

!!! Determine number of components to transport
!   ncomps = RM_FindComponents(id)
  ncomps = bmif%FindComponents()
  
!!! Print some of the reaction module information		
!   write(string1, "(A,I10)") "Number of threads:                                ", RM_GetThreadCount(id)
!   print *, string1
!   write(string1, "(A,I10)") "Number of MPI processes:                          ", RM_GetMpiTasks(id)
!   print *, string1
!   write(string1, "(A,I10)") "MPI task number:                                  ", RM_GetMpiMyself(id)
!   print *, string1
!   status = RM_GetFilePrefix(id, alloc_string)
!   write(string1, "(A,A)") "File prefix:                                      ", alloc_string
!   print *, string1
!   write(string1, "(A,I10)") "Number of grid cells in the user's model:         ", RM_GetGridCellCount(id)
!   print *, string1
!   write(string1, "(A,I10)") "Number of chemistry cells in the reaction module: ", RM_GetChemistryCellCount(id)
!   print *, string1
!   write(string1, "(A,I10)") "Number of components for transport:               ", RM_GetComponentCount(id)
!   print *, string1 
  
!!! Get component information
!   status = RM_GetGfw(id, gfw)
  status = bmif%GetGfw(gfw)
!   status = RM_GetComponents(id, components)
  status = bmif%GetComponents(components)
  do i = 1, ncomps
     write(string,"(A10, F15.4)") components(i), gfw(i)
   !   print *, string
  enddo

!!! Example Selected Output
!   call example_selected_output(id)

!!! Set initial conditions
  ! Four ways to set initial conditions
  ! 1. Mixing
  !allocate(ic1(nxyz,7), ic2(nxyz,7), f1(nxyz,7))
  !ic1 = -1
  !ic2 = -1
  !f1 = 1.0
  !do i = 1, nxyz
  !   ic1(i,1) = 1       ! Solution 1
  !   ic1(i,2) = -1      ! Equilibrium phases none
  !   ic1(i,3) = 1       ! Exchange 1
  !   ic1(i,4) = -1      ! Surface none
  !   ic1(i,5) = -1      ! Gas phase none
  !   ic1(i,6) = -1      ! Solid solutions none
  !   ic1(i,7) = -1      ! Kinetics none
  !enddo
  !status = RM_InitialPhreeqc2Module(id, ic1, ic2, f1)
  ! 2. No mixing is defined, so the following is equivalent
  ! status = RM_InitialPhreeqc2Module(id, ic1)
  ! 3. Simplest for these conditions
!   allocate(solutions(nxyz), exchanges(nxyz))
!   solutions = 1
! !   status = RM_InitialSolutions2Module(id, solutions)
!   status = bmif%InitialSolutions2Module(solutions)
!   exchanges = 1
! !   status = RM_InitialExchanges2Module(id, exchanges)
!   status = bmif%InitialExchanges2Module(exchanges)
!   ! 4. alternative for setting initial conditions
!   ! cell number in second argument (-1 indicates last solution, 40 in this case)
!   ! in advect.pqi and any reactants with the same number--
!   ! Equilibrium phases, exchange, surface, gas phase, solid solution, and (or) kinetics--
!   ! will be written to cells 18 and 19 (0 based)
!   allocate (module_cells(2))
!   module_cells(1) = 18
!   module_cells(2) = 19
! !   status = RM_InitialPhreeqcCell2Module(id, -1, module_cells, 2)
!   status = bmif%InitialPhreeqcCell2Module(-1, module_cells, 2)


!!! Initial equilibration of cells
  print *, "pre"
  time = 0.0
  time_step = 0.0
  allocate(c(nxyz, ncomps))
!   status = RM_SetTime(id, time)
  status = bmif%SetTime(time)
!   status = RM_SetTimeStep(id, time_step)
  status = bmif%SetTimeStep(time_step)
!   status = RM_RunCells(id) 
  status = bmif%RunCells()
  print *, "post"
!   status = RM_GetConcentrations(id, c) ! Get concentrations after equilibration to array 'c'
  status = bmif%GetConcentrations(c) ! Get concentrations after equilibration to array 'c'

!   Print initial concentrations for each cell
!   print *, 'Concentrations: '
!   do i = 1, nxyz  ! Loop through cells
!       print *, 'Cell ', i
!       do j = 1, ncomps  ! Loop through components
!           print *, components(j), ':', c(i,j)
!       end do
!       print *  ! Add an empty line for better readability between cells
!   end do

!!! --------------------------------------------------------------------------
!!! #TODO: Set boundary condition ! SHOULD ATMOSPHERE COME HERE?
!!! --------------------------------------------------------------------------

  nbound = 1
  allocate(bc1(nbound), bc2(nbound), bc_f1(nbound))
  allocate(bc_conc(nbound, ncomps))  
  bc1 = 0           ! solution 0 from Initial IPhreeqc instance
  bc2 = -1          ! no bc2 solution for mixing
  bc_f1 = 1.0       ! mixing fraction for bc1 
  status = RM_InitialPhreeqc2Concentrations(id, bc_conc, nbound, bc1, bc2, bc_f1)


!!! Register state variables: Dissolved Components for libaed-water

!   status = RM_GetComponents(id, components)
!   print *,'NCOMPS:', ncomps
!   print *,'COMPONENTS:'

  ALLOCATE(data%id_phreeqcrm_dv(ncomps))
  !ALLOCATE(data%id_phreeqcrm_sv(num_sv))
!   status = RM_GetComponents(id, names)
  DO i=1,ncomps
   ! print *, names(i)
   data%id_phreeqcrm_dv(i) = aed_define_variable('CO_'//trim(components(i)) , 'mmol/m3', 'long name', dr_init(i), dr_min(i), dr_max(i), 0.)
  ENDDO


!!!! Register state variables: Equilibrium Phases, for non-cohesive? 
   ! nep = RM_GetEquilibriumPhasesCount(id)
   nep = bmif%GetEquilibriumPhasesCount()
   ALLOCATE(data%id_phreeqcrm_sv(nep))
   ! status = RM_GetEquilibriumPhasesNames(id, names)
   status = bmif%GetEquilibriumPhasesNames(names)
   DO i=1,nep
    data%id_phreeqcrm_sv(i) = aed_define_variable('EP_'//trim(names(i)) , 'mmol/m3', 'long name', dr_init(i), dr_min(i), dr_max(i), 0.)
   ENDDO