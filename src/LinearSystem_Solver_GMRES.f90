module LinearSystem_Solver_GMRES_mod
      use KindDefinition_mod, only : DP
      use LinearSystem_Generic_mod
      use IO_mod
      use GlobalVariableForDebug_mod
      use GMRES_mod
      use LinearSystem_Preconditioner_ILU_mod      ! this is tucktaping, should not use non-generic preconditioner here. Change the code later
      implicit none
      save

      type, extends(linear_solver_t) :: GMRES_t
        private
        integer :: nrestart, nsearch
        real(DP) :: tolerance

        real(DP), dimension(:,:), allocatable :: F
        real(DP), dimension(:),   allocatable :: AP
        real(DP), dimension(:),   allocatable :: X
        real(DP), dimension(:,:), allocatable :: phi
        real(DP), dimension(:,:), allocatable :: tres

        type(ILU_t) :: ILU                          ! this is tucktaping, should not use non-generic preconditioner here, the code should be changed later
      contains
        procedure, pass(this) :: Reallocate_linear_solver => Reallocate_GMRES
        procedure, pass(this) :: Solve_linear_system => GMRES
      end type GMRES_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reallocate_GMRES( this, LS )
      implicit none
      class(GMRES_t), intent(out) :: this
      class(linear_system_t), intent(in) :: LS
      integer :: ndim
      integer :: n
      integer :: nsearch

      if( allocated (this%F) )     deallocate( this%F )
      if( allocated (this%AP) )    deallocate( this%AP )
      if( allocated (this%X) )     deallocate( this%X )
      if( allocated (this%phi) )   deallocate( this%phi )
      if( allocated (this%tres) )  deallocate( this%tres )

      this%nrestart = 3
      this%nsearch = 50
      this%tolerance = 1.d-14

      nsearch = this%nsearch

      ndim = LS%ndim
      n = LS%n

      allocate( this%F   (ndim,nsearch+1) )
      allocate( this%AP  (ndim*n) )
      allocate( this%X   (ndim*n) )
      allocate( this%phi (ndim,n) )
      allocate( this%tres(ndim,n) )

      ! this is tucktaping, should not use non-generic preconditioner here
      call this % ILU % Reallocate_preconditioner ( LS )

end subroutine Reallocate_GMRES
!=====================================================================
!
!=====================================================================
subroutine GMRES ( this, LS, preconditioner )
      implicit none
      class(GMRES_t),          intent(inout) :: this
      class(linear_system_t),  intent(inout) :: LS
      class(preconditioner_t), intent(inout) :: preconditioner

      integer :: n, ndim

      n = LS%n
      ndim = LS%ndim

!     initialize solution
      LS%x(:,:) = 0.

      LS%rhs(:,:) = -LS%rhs(:,:)

      ! this is tucktaping, should not use non-generic preconditioner here. Change the code later
      call gmres_ ( LS%x, LS%A%e, LS%rhs, LS%A%ia, LS%A%iau, LS%A%ja, this % ILU % ALU%e, this % ILU % ALU%ia, this % ILU % ALU%iau, this % ILU % ALU%ja, n ,LS%A%nnz, this % ILU % ALU%nnz, ndim, this%tolerance, this%nrestart, this%nsearch)


end subroutine GMRES

end module LinearSystem_Solver_GMRES_mod
