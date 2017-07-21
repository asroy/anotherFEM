module LinearSystem_Generic_mod
      use KindDefinition_mod, only : DP
      use Common_Type_CRS_mod
      use Common_Type_Dependency_mod
      implicit none
      save

      type :: linear_system_t
        integer :: ndim
        integer :: n
        type(CRS_t) :: A
        real(DP), allocatable :: rhs(:,:)
        real(DP), allocatable :: x(:,:)
      contains
        procedure, pass(this) :: Reallocate_linear_system
      end type linear_system_t

      type :: linear_solver_t
        integer :: ndim
        integer :: n
      contains
        procedure, pass(this) :: Reallocate_linear_solver
        procedure, pass(this) :: Solve_linear_system
      end type linear_solver_t

      type :: preconditioner_t
        integer :: ndim
        integer :: n
      contains
        procedure, pass(this) :: Reallocate_preconditioner
        procedure, pass(this) :: Calculate_preconditioner
        procedure, pass(this) :: Apply_preconditioner
      end type preconditioner_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reallocate_linear_system ( this, ndim, n, depen )
      implicit none
      class(linear_system_t), intent(inout) :: this
      integer,                intent(in)    :: ndim, n
      class(depen_t),         intent(in)    :: depen

      this%ndim = ndim
      this%n = n

      if( allocated (this%x) ) then
        deallocate( this%x )
        deallocate( this%rhs )
        call this%A%Deallocate_CRS
      end if

      allocate( this%x(ndim,n) )
      allocate( this%rhs(ndim,n) )

      call this%A%Reallocate_CRS ( ndim, n, depen%nnz, depen%ia, depen%ja )

end subroutine Reallocate_linear_system
!=====================================================================
!
!=====================================================================
subroutine Reallocate_linear_solver ( this, LS )
      implicit none
      class(linear_solver_t), intent(out) :: this
      class(linear_system_t), intent(in) :: LS
end subroutine Reallocate_linear_solver
!=====================================================================
!
!=====================================================================
subroutine Solve_linear_system ( this, LS, preconditioner )
      implicit none
      class(linear_solver_t), intent(inout) :: this
      class(linear_system_t), intent(inout) :: LS
      class(preconditioner_t), intent(inout) :: preconditioner
end subroutine Solve_linear_system
!=====================================================================
!
!=====================================================================
subroutine Reallocate_preconditioner ( this, LS )
      implicit none
      class(preconditioner_t), intent(out) :: this
      class(linear_system_t),  intent(in)  :: LS
end subroutine Reallocate_preconditioner
!=====================================================================
!
!=====================================================================
subroutine Calculate_preconditioner ( this, LS )
      implicit none
      class(preconditioner_t), intent(inout) :: this
      class(linear_system_t),  intent(in)    :: LS
end subroutine Calculate_preconditioner
!=====================================================================
!
!=====================================================================
subroutine Apply_preconditioner ( this, x, y )
      implicit none
      class(preconditioner_t),                 intent(in)  :: this
      real(DP),                dimension(:,:), intent(out) :: x
      real(DP),                dimension(:),   intent(in)  :: y
end subroutine Apply_preconditioner

end module LinearSystem_Generic_mod
