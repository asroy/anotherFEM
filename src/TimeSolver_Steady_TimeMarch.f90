module TimeSolver_Steady_TimeMarch_mod
      use KindDefinition_mod, only : DP
      use Common_Function_mod
      use IO_mod
      use TimeSolver_Generic_mod
      use GlobalVariableForDebug_mod
      save

      type, extends(time_solver_t) :: steady_time_march_t
        logical :: Use_dt_local
        real(DP) :: dt_global

        real(DP) :: cfl_min
        real(DP) :: cfl_max
        integer :: iramp

        integer :: niteration
        real(DP) :: tolerance
      contains
        procedure, pass(this) :: Reset_time_solver => Reset_steady_time_march
        procedure, pass(this) :: Do_time_solver => Do_steady_time_march
      end type steady_time_march_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reset_steady_time_march ( this )
      implicit none
      class(steady_time_march_t), intent(out) :: this

      this % Use_dt_local = .true.
      this % dt_global = 0.

      this % cfl_min = 1.
      this % cfl_max = 1000.
      this % iramp   = 100

      ! niteration, tolerance
      this % niteration = 1000
      this % tolerance = 1.e-15

end subroutine Reset_steady_time_march
!=====================================================================
!
!=====================================================================
subroutine Do_steady_time_march ( this, FEM, equation, linear_system, linear_solver, preconditioner )
      implicit none
      class(steady_time_march_t),  intent(in)    :: this
      class(FEM_t),               intent(inout) :: FEM
      class(equation_t),          intent(inout) :: equation
      class(linear_system_t),     intent(inout) :: linear_system
      class(linear_solver_t),     intent(inout) :: linear_solver
      class(preconditioner_t),    intent(inout) :: preconditioner

      type(time_control_t) :: time_control
      real(DP), dimension(nequation_MAX) :: rms
      integer :: ndim, n
      integer :: iter

      ndim = linear_system % ndim
      n    = linear_system % n

      time_control % Include_time = .true.
      time_control % ntimeStep_time_low = 0
      time_control % ntimeStep_time_up  = 1
      time_control % alph(0:1) = [ -1., 1. ]
      time_control % Use_dt_local = .true.
      time_control % dt_global = 0.
      time_control % Include_space = .true.
      time_control % itimeStep_space = 1
      time_control % ratio_space = 1.

      do iter = 1, this % niteration
        if (iter < this % iramp ) then
          time_control % cfl = this % cfl_min + ( iter - 1 )*( this % cfl_max - this % cfl_min)/( this % iramp - 1 )
        else
          time_control % cfl = this % cfl_max
        end if

       !if( iter == 3 ) call FEM % Read_q( 1, fpIO, 'output/read_q-2.bin')

       !call equation % Calculate_dt_local   ( FEM, time_control )
        call equation % Calculate_dt_local_2 ( FEM, time_control )

        call FEM % Calculate_FEM_residual( linear_system % rhs, equation, time_control )
        linear_system % rhs(:,:) = - linear_system % rhs(:,:)

        rms = vector_rms ( linear_system % rhs, ndim, n )/n
        write(*,   200)  'Do_steady_time_march:', iter, time_control % cfl, rms(1:ndim)
        write(fpd1,200)  'Do_steady_time_march:', iter, time_control % cfl, rms(1:ndim)

        if( rms(1) < this % tolerance ) exit

        call FEM % Calculate_FEM_jacobian( linear_system % A, equation, time_control )

        call linear_solver % Solve_linear_system( linear_system, preconditioner )

        call FEM % Update_FEM_solution( linear_system % x, 1. , 1 )

        call FEM % Output_FEM_solution_tecplot( fpIO,'output/Do_steady_time_march-q.plt' )
        call equation % Check_solution ( FEM, 1, fpIO, 'output/check_solution.plt' )

        call FEM % Update_FEM_solution_history ( 0, 1 )

       !if( iter == 2 ) call FEM%Write_q(1, fpIO, 'output/write_q-2.bin')
      end do


100   format(a,20(es16.5))
200   format(a,i7,20(es16.5))

end subroutine Do_steady_time_march

end module TimeSolver_Steady_TimeMarch_mod
