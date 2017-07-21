module TimeSolver_Generic_mod
      use KindDefinition_mod, only : DP
      use FEM_Generic_mod
      use LinearSystem_Generic_mod

      save

! - time solver
      type :: time_solver_t
      contains
        procedure, pass(this) :: Reset_time_solver
        procedure, pass(this) :: Do_time_solver
      end type time_solver_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reset_time_solver ( this )
      implicit none
      class(time_solver_t), intent(out) :: this
end subroutine Reset_time_solver
!=====================================================================
!
!=====================================================================
subroutine Do_time_solver ( this, FEM, equation, linear_system, linear_solver, preconditioner )
      implicit none
      class(time_solver_t),      intent(in)    :: this
      class(FEM_t),              intent(inout) :: FEM
      class(equation_t),         intent(inout) :: equation
      class(linear_system_t),    intent(inout) :: linear_system
      class(linear_solver_t),    intent(inout) :: linear_solver
      class(preconditioner_t),   intent(inout) :: preconditioner
end subroutine Do_time_solver

end module TimeSolver_Generic_mod
