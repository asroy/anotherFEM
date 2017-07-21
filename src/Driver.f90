program PlayStation
      use KindDefinition_mod, only : DP
      use FEM_Generic_mod
      use FEM_ElementType_Lagrange_mod
      use FEM_equation_Poisson_mod
      use FEM_equation_Euler_mod
      use LinearSystem_Generic_mod
      use LinearSystem_Solver_GMRES_mod
      use LinearSystem_Solver_LUSGS_mod
      use LinearSystem_Preconditioner_ILU_mod
      use Mesh_Generic_mod
      use Mesh_Importer_Grd_mod
      use TimeSolver_Generic_mod
      use TimeSolver_Steady_Newton_mod
      use TimeSolver_Steady_TimeMarch_mod
      use TimeSolver_Steady_TimeMarch_RmsMinimization_mod
      use GlobalVariableForDebug_mod
      implicit none

      ! mesh
      type (mesh_t) :: mesh
      class (mesh_importer_t), allocatable :: mesh_importer

      ! discretization
      class (FEM_t), allocatable :: FEM

      ! equation
      class (equation_t), allocatable :: equation

      ! time solver
      class (time_solver_t), allocatable :: time_solver

      ! linear system
      type (linear_system_t) :: linear_system
      class (linear_solver_t), allocatable :: linear_solver
      class (preconditioner_t), allocatable :: preconditioner


      open(fpd1,file='output/fpd1.dat')
      open(fpd2,file='output/fpd2.dat')
      open(fpd3,file='output/fpd3.dat')
      open(fpd4,file='output/fpd4.dat')

      ! mesh
      allocate ( grd_mesh_importer_t :: mesh_importer )

      ! discretization
      allocate ( lagrange_t :: FEM )

      ! equation
     !allocate ( poisson_t        :: equation )
      allocate ( euler_equation_t :: equation )

      ! non-linear solver
     !allocate ( steady_newton_t     :: time_solver )
     !allocate ( steady_time_march_t :: time_solver )
      allocate ( steady_time_march_rms_minimization_t :: time_solver )

      ! linear solver
     !allocate ( LUSGS_t :: linear_solver )
      allocate ( GMRES_t :: linear_solver )

      ! preconditioner for linear system
      allocate ( ILU_t :: preconditioner )

      ! import mesh
     !call mesh_importer%Import_mesh ( mesh, fpIO, 'input/one_triangle.grd' )
     !call mesh_importer%Import_mesh ( mesh, fpIO, 'input/onegrid-16.grd' )
     !call mesh_importer%Import_mesh ( mesh, 1, fpIO, 'input/naca0012-coarse.grd' )
     !call mesh_importer%Import_mesh ( mesh, 1, fpIO, 'input/naca0012-fine.grd' )
      call mesh_importer%Import_mesh ( mesh, 1, fpIO, 'input/naca0012-trex.grd' )
     !call mesh_importer%Import_mesh ( mesh, 1, fpIO, 'input/v4-naca0012_tuo-p1.grd' )
      call mesh%Output_mesh_tecplot ( fpIO, 'output/mesh.plt')

      ! equation
      call equation%Reset_equation

      ! generate discretization from mesh
      call FEM%Generate_FEM ( mesh, equation )
      write(*,*) 'PlayStation:', FEM%nequation, FEM%nfbasis
      call FEM%depen_fbasis%Output_depen_tecplot ( fpIO, 'output/depen_fbasis.plt' )
      call FEM%depen_matrix%Output_depen_tecplot ( fpIO, 'output/depen_matrix.plt' )


      ! initialize solution
      call equation%Initialize_solution ( FEM )

      call FEM % Output_FEM_solution_tecplot( fpIO,'output/initial_solution.plt' )

      ! prepare linear system
      call linear_system%Reallocate_linear_system ( FEM%nequation, FEM%nfbasis, FEM%depen_matrix )
      call linear_solver%Reallocate_linear_solver ( linear_system )
     !call preconditioner%Reallocate_preconditioner ( linear_system )

      ! non-linear solver
      call time_solver%Reset_time_solver
      call time_solver%Do_time_solver ( FEM, equation, linear_system, linear_solver, preconditioner )

      call FEM % Output_FEM_solution_tecplot( fpIO,'output/solution.plt' )
end program PlayStation
