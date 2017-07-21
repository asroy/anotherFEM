module FEMEntry_LocalNode_Default_mod
      use KindDefinition_mod, only : DP
      use FEMEntry_Generic_mod
      implicit none
      save

      type, extends (local_node_t) :: default_local_node_t
      contains
        procedure, pass(this) :: Reset_local_node => Reset_default_local_node
      end type default_local_node_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reset_default_local_node( this, basis_order )
      implicit none
      class(default_local_node_t), intent(inout) :: this
      integer,                     intent(in)    :: basis_order

      this%basis_order = basis_order

      select case( basis_order )
      case(1)
        this%nnd = 3
        this%rnd( 1) = 0.    ;  this%snd( 1) = 0.
        this%rnd( 2) = 1.    ;  this%snd( 2) = 0.
        this%rnd( 3) = 0.    ;  this%snd( 3) = 1.
      case(2)
        this%nnd = 6
        this%rnd( 1) = 0.    ;  this%snd( 1) = 0.
        this%rnd( 2) = 1.    ;  this%snd( 2) = 0.
        this%rnd( 3) = 0.    ;  this%snd( 3) = 1.
        this%rnd( 4) = 1./2. ;  this%snd( 4) = 0.
        this%rnd( 5) = 1./2. ;  this%snd( 5) = 1./2.
        this%rnd( 6) = 0.    ;  this%snd( 6) = 1./2.
      case(3)
        this%nnd = 10
        this%rnd( 1) = 0.    ;  this%snd( 1) = 0.
        this%rnd( 2) = 1.    ;  this%snd( 2) = 0.
        this%rnd( 3) = 0.    ;  this%snd( 3) = 1.
        this%rnd( 4) = 1./3. ;  this%snd( 4) = 0.
        this%rnd( 5) = 2./3. ;  this%snd( 5) = 1./3.
        this%rnd( 6) = 0.    ;  this%snd( 6) = 2./3.
        this%rnd( 7) = 2./3. ;  this%snd( 7) = 0.
        this%rnd( 8) = 1./3. ;  this%snd( 8) = 2./3.
        this%rnd( 9) = 0.    ;  this%snd( 9) = 1./3.
        this%rnd(10) = 1./3. ;  this%snd(10) = 1./3.
      end select

end subroutine Reset_default_local_node

end module FEMEntry_LocalNode_Default_mod
