module TimeSolver_Steady_Newton_mod
      use KindDefinition_mod, only : DP
      use Common_Function_mod
      use IO_mod
      use TimeSolver_Generic_mod
      use GlobalVariableForDebug_mod
      save

      type, extends(time_solver_t) :: steady_newton_t
        integer :: niteration
        real(DP) :: tolerance
      contains
        procedure, pass(this) :: Reset_time_solver => Reset_steady_newton
        procedure, pass(this) :: Do_time_solver => Do_steady_newton
      end type steady_newton_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reset_steady_newton ( this )
      implicit none
      class(steady_newton_t), intent(out) :: this

      ! niteration, tolerance
      this % niteration = 1
      this % tolerance = 1.e-14

end subroutine Reset_steady_newton
!=====================================================================
!
!=====================================================================
subroutine Do_steady_newton ( this, FEM, equation, linear_system, linear_solver, preconditioner )
      implicit none
      class(steady_newton_t),  intent(in)    :: this
      class(FEM_t),            intent(inout) :: FEM
      class(equation_t),       intent(inout) :: equation
      class(linear_system_t),  intent(inout) :: linear_system
      class(linear_solver_t),  intent(inout) :: linear_solver
      class(preconditioner_t), intent(inout) :: preconditioner

      type(time_control_t) :: time_control
      real(DP) :: rms(nequation_MAX)
      integer :: ndim, n
      integer :: iter

      ndim = linear_system % ndim
      n    = linear_system % n

      time_control % Include_time = .false.
      time_control % ntimeStep_time_low = 0
      time_control % ntimeStep_time_up  = 0
      time_control % alph(:) = 0.
      time_control % Use_dt_local = .false.
      time_control % dt_global = 0.
      time_control % Include_space = .true.
      time_control % itimeStep_space = 1
      time_control % ratio_space = 1.

      do iter = 1, this % niteration
        call FEM % Calculate_FEM_residual( linear_system % rhs, equation, time_control )
        linear_system % rhs(:,:) = - linear_system % rhs(:,:)

        rms = vector_rms ( linear_system % rhs, ndim, n )/n
        write(*   ,100)  'Do_steady_newton:', iter, rms(1:ndim)
        write(fpd1,100)  'Do_steady_newton:', iter, rms(1:ndim)

        if( rms(1) < this % tolerance ) exit

        call FEM % Calculate_FEM_jacobian( linear_system % A, equation, time_control )

        call linear_solver % Solve_linear_system( linear_system, preconditioner )

        call FEM % Update_FEM_solution( linear_system % x, 1., 1 )

        call FEM % Output_FEM_solution_tecplot( fpIO,'output/Do_steady_newton-q.plt' )
      end do

100   format(a,i6,20(es16.4))

end subroutine Do_steady_newton

end module TimeSolver_Steady_Newton_mod
