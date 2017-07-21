module FEMEntry_ElementType_Lagrange_mod
      use KindDefinition_mod, only : DP
      use FEMEntry_Generic_mod
      use GlobalVariableForDebug_mod
      implicit none
      save

! - Lagranage mapping of element and bedge
      type, extends (element_t) :: lagrange_element_t
      contains
        procedure, pass(this) :: Initialize_element => Initialize_lagrange_element
        procedure, pass(this) :: Get_element_N => Get_lagrange_element_N
        procedure, pass(this) :: Get_element_Nr_Ns => Get_lagrange_element_Nr_Ns
        procedure, pass(this) :: Get_element_Nrr_Nrs_Nss => Get_lagrange_element_Nrr_Nrs_Nss
      end type lagrange_element_t

      type, extends (bedge_t) :: lagrange_bedge_t
      contains
        procedure, pass(this) :: Initialize_bedge => Initialize_lagrange_bedge
        procedure, pass(this) :: Get_bedge_to_element_rs => Get_lagrange_bedge_to_element_rs
        procedure, pass(this) :: Get_bedge_N_nz => Get_lagrange_bedge_N_nz
        procedure, pass(this) :: Get_bedge_Nr_nz => Get_lagrange_bedge_Nr_nz
      end type lagrange_bedge_t

contains

!=====================================================================
!
!=====================================================================
subroutine Initialize_lagrange_element ( this, basis_order )
      implicit none
      class(lagrange_element_t), intent(inout) :: this
      integer,                   intent(in)    :: basis_order

      this % basis_order = basis_order

      select case ( basis_order )
      case(1)
        this % nfb = 3
        this % nnd = 3
      case(2)
        this % nfb = 6
        this % nnd = 6
      case(3)
        this % nfb = 10
        this % nnd = 10
      case default
        write(*,*) 'wrong! Initialize_lagrange_element: basis_order', basis_order
        stop
      end select

      this % area_rs = 0.5

end subroutine Initialize_lagrange_element
!=====================================================================
!
!=====================================================================
subroutine Get_lagrange_element_N ( this, N, r, s )
      implicit none
      class(lagrange_element_t),               intent(in)  :: this
      real(DP),                  dimension(:), intent(out) :: N
      real(DP),                                intent(in)  :: r, s

      real(DP) :: r2, rs, s2, r3, r2s, rs2, s3

      select case( this%basis_order )
      case(1)
        N(1) = 1. - r - s
        N(2) =      r
        N(3) =          s
      case(2)
        r2 = r*r
        rs = r*s
        s2 = s*s

        N(1) = 1. - 3.*r - 3.*s + 2.*r2 + 4.*rs + 2.*s2
        N(2) =    -    r        + 2.*r2
        N(3) =           -    s                 + 2.*s2
        N(4) =      4.*r        - 4.*r2 - 4.*rs
        N(5) =                            4.*rs
        N(6) =             4.*s         - 4.*rs - 4.*s2
      case(3)
        r2 = r*r
        rs = r*s
        s2 = s*s
        r3  = r2*r
        r2s = r2*s
        rs2 = r*s2
        s3  = s2*s

        N(1)  = 1. - 5.5*r - 5.5*s +  9. *r2 + 18. *rs +  9. *s2 -  4.5*r3 - 13.5*r2s - 13.5*rs2 -  4.5*s3
        N(2)  =          r         -  4.5*r2                     +  4.5*r3
        N(3)  =                  s                     -  4.5*s2                                 +  4.5*s3
        N(4)  =      9. *r         - 22.5*r2 - 22.5*rs           + 13.5*r3 + 27. *r2s + 13.5*rs2
        N(5)  =                              -  4.5*rs                     + 13.5*r2s
        N(6)  =            - 4.5*s           +  4.5*rs + 18. *s2                      - 13.5*rs2 - 13.5*s3
        N(7)  =    - 4.5*r         + 18. *r2 +  4.5*rs           - 13.5*r3 - 13.5*r2s
        N(8)  =                              -  4.5*rs                                + 13.5*rs2
        N(9)  =              9. *s           - 22.5*rs - 22.5*s2           + 13.5*r2s + 27. *rs2 + 13.5*s3
        N(10) =                                27. *rs                     - 27. *r2s - 27. *rs2
      case default
        write(*,*) 'wrong! Get_lagrange_element_N: not defined'
        stop
      end select

end subroutine Get_lagrange_element_N
!=====================================================================
!
!=====================================================================
subroutine Get_lagrange_element_Nr_Ns ( this, Nr, Ns, r, s )
      implicit none
      class(lagrange_element_t),               intent(in)  :: this
      real(DP),                  dimension(:), intent(out) :: Nr, Ns
      real(DP),                                intent(in)  :: r, s

      real(DP) :: r2, rs, s2

      select case( this%basis_order )
      case(1)
        Nr(1) = - 1.
        Nr(2) =   1.
        Nr(3) =   0.

        Ns(1) = - 1.
        Ns(2) =   0.
        Ns(3) =   1.
      case(2)
        Nr(1) = - 3. + 4.*r + 4.*s
        Nr(2) = - 1. + 4.*r
        Nr(3) =   0.
        Nr(4) =   4. - 8.*r - 4.*s
        Nr(5) =               4.*s
        Nr(6) =             - 4.*s

        Ns(1) = - 3. + 4.*r + 4.*s
        Ns(2) =   0.
        Ns(3) = - 1         + 4.*s
        Ns(4) =      - 4.*r
        Ns(5) =        4.*r
        Ns(6) =   4. - 4.*r - 8.*s
      case(3)
        r2 = r*r
        rs = r*s
        s2 = s*s

        Nr(1)  = - 5.5 + 18. *r + 18. *s - 13.5*r2 - 27. *rs - 13.5*s2
        Nr(2)  =   1.  -  9. *r          + 13.5*r2
        Nr(3)  =   0.
        Nr(4)  =   9.  - 45. *r - 22.5*s + 40.5*r2 + 54. *rs + 13.5*s2
        Nr(5)  =                -  4.5*s           + 27. *rs
        Nr(6)  =                   4.5*s                     - 13.5*s2
        Nr(7)  = - 4.5 + 36. *r +  4.5*s - 40.5*r2 - 27. *rs
        Nr(8)  =                -  4.5*s                     + 13.5*s2
        Nr(9)  =                - 22.5*s           + 27. *rs + 27. *s2
        Nr(10) =                  27. *s           - 54. *rs - 27. *s2

        Ns(1)  = - 5.5 + 18. *r + 18. *s - 13.5*r2 - 27. *rs - 13.5*s2
        Ns(2)  =   0.
        Ns(3)  =   1.           -  9. *s                     + 13.5*s2
        Ns(4)  =       - 22.5*r          + 27. *r2 + 27. *rs
        Ns(5)  =       -  4.5*r          + 13.5*r2
        Ns(6)  = - 4.5 +  4.5*r + 36. *s           - 27. *rs - 40.5*s2
        Ns(7)  =          4.5*r          - 13.5*r2
        Ns(8)  =       -  4.5*r                    + 27. *rs
        Ns(9)  =   9.  - 22.5*r - 45. *s + 13.5*r2 + 54. *rs + 40.5*s2
        Ns(10) =         27. *r          - 27. *r2 - 54. *rs
      case default
        write(*,*) 'wrong! Get_lagrange_element_Nr_Ns: not defined'
        stop
      end select

end subroutine Get_lagrange_element_Nr_Ns
!=====================================================================
!
!=====================================================================
subroutine Get_lagrange_element_Nrr_Nrs_Nss ( this, Nrr, Nrs, Nss, r, s )
      implicit none
      class(lagrange_element_t),               intent(in)  :: this
      real(DP),                  dimension(:), intent(out) :: Nrr, Nrs, Nss
      real(DP),                                intent(in)  :: r, s

      select case( this % basis_order )
      case(1)
        Nrr(1:3) = 0.
        Nrs(1:3) = 0.
        Nss(1:3) = 0.
      case(2)
        Nrr(1) =   4.
        Nrr(2) =   4.
        Nrr(3) =   0.
        Nrr(4) = - 8.
        Nrr(5) =   0.
        Nrr(6) =   0.

        Nrs(1) =   4.
        Nrs(2) =   0.
        Nrs(3) =   0.
        Nrs(4) = - 4.
        Nrs(5) =   4.
        Nrs(6) = - 4.

        Nss(1) =   4.
        Nss(2) =   0.
        Nss(3) =   4.
        Nss(4) =   0.
        Nss(5) =   0.
        Nss(6) = - 8.
      case(3)
        Nrr(1)  =   18.  - 27. *r - 27. *s
        Nrr(2)  = -  9.  + 27. *r
        Nrr(3)  =    0.
        Nrr(4)  = - 45.  + 81. *r + 54. *s
        Nrr(5)  =                   27. *s
        Nrr(6)  =    0.
        Nrr(7)  =   36.  - 81. *r - 27. *s
        Nrr(8)  =    0.
        Nrr(9)  =                   27. *s
        Nrr(10) =                 - 54. *s

        Nrs(1)  =   18.  - 27. *r - 27. *s
        Nrs(2)  =    0.
        Nrs(3)  =    0.
        Nrs(4)  = - 22.5 + 54. *r + 27. *s
        Nrs(5)  = -  4.5 + 27. *r
        Nrs(6)  =    4.5          - 27. *s
        Nrs(7)  =    4.5 - 27. *r
        Nrs(8)  = -  4.5          + 27. *s
        Nrs(9)  = - 22.5 + 27. *r + 54. *s
        Nrs(10) =   27.  - 54. *r - 54. *s

        Nss(1)  =   18.  - 27. *r - 27. *s
        Nss(2)  =    0.
        Nss(3)  = -  9.           + 27. *s
        Nss(4)  =          27. *r
        Nss(5)  =    0.
        Nss(6)  =   36.  - 27. *r - 81. *s
        Nss(7)  =    0.
        Nss(8)  =          27. *r
        Nss(9)  = - 45.  + 54. *r + 81. *s
        Nss(10) =        - 54. *r
      case default
        write(*,*) 'wrong! Get_lagrange_element_Nrr_Nrs_Nss: not defined'
        stop
      end select

end subroutine Get_lagrange_element_Nrr_Nrs_Nss
!=====================================================================
!
!=====================================================================
subroutine Initialize_lagrange_bedge ( this, basis_order )
      implicit none
      class(lagrange_bedge_t), intent(inout) :: this
      integer,                 intent(in)    :: basis_order

      this % basis_order = basis_order

      select case( basis_order )
      case(1)
        this % nfb   = 3
        this % nfbnz = 2
      case(2)
        this % nfb   = 6
        this % nfbnz = 3
      case(3)
        this % nfb   = 10
        this % nfbnz = 4
      case default
        write(*,*) 'wrong! Initialize_lagrange_bedge: basis_order', basis_order
        stop
      end select

      this % length_r = 1.

end subroutine Initialize_lagrange_bedge
!=====================================================================
!
!=====================================================================
subroutine Get_lagrange_bedge_to_element_rs ( this, r, s, rb )
      implicit none
      class(lagrange_bedge_t), intent(in)  :: this
      real(DP),                intent(out) :: r, s
      real(DP),                intent(in)  :: rb

      select case( this % elementSide )
      case(1)
        r = rb
        s = 0.
      case(2)
        r = 1. - rb
        s = rb
      case(3)
        r = 0.
        s = 1. - rb
      case default
        write(*,*) 'wrong! Get_lagrange_bedge_to_element_rs'
        stop
      end select

end subroutine Get_lagrange_bedge_to_element_rs
!=====================================================================
!
!=====================================================================
subroutine Get_lagrange_bedge_N_nz ( this, N_nz, rb )
      implicit none
      class(lagrange_bedge_t),               intent(in)  :: this
      real(DP),                dimension(:), intent(out) :: N_nz
      real(DP),                              intent(in)  :: rb

      real(DP) :: rb2, rb3

      select case( this % basis_order )
      case(1)
        N_nz(1) = 1. - rb
        N_nz(2) =      rb
      case(2)
        rb2 = rb*rb

        N_nz(1) = 1. - 3.*rb + 2.*rb2
        N_nz(2) =    - 1.*rb + 2.*rb2
        N_nz(3) =      4.*rb - 4.*rb2
      case(3)
        rb2 = rb*rb
        rb3 = rb2*rb

        N_nz(1) = 1. - 5.5*rb +  9. *rb2 -  4.5*rb3
        N_nz(2) =          rb -  4.5*rb2 +  4.5*rb3
        N_nz(3) =      9. *rb - 22.5*rb2 + 13.5*rb3
        N_nz(4) =    - 4.5*rb + 18. *rb2 - 13.5*rb3
      case default
        write(*,*) 'wrong! Get_lagrange_bedge_N_nz'
        stop
      end select

end subroutine Get_lagrange_bedge_N_nz
!=====================================================================
!
!=====================================================================
subroutine Get_lagrange_bedge_Nr_nz ( this, Nr_nz, rb )
      implicit none
      class(lagrange_bedge_t),               intent(in)  :: this
      real(DP),                dimension(:), intent(out) :: Nr_nz
      real(DP),                              intent(in)  :: rb

      real(DP) :: rb2

      select case( this%basis_order )
      case(1)
        Nr_nz(1) = - 1.
        Nr_nz(2) =   1.
      case(2)
        rb2 = rb*rb

        Nr_nz(1) = - 3. + 4.*rb
        Nr_nz(2) = - 1. + 4.*rb
        Nr_nz(3) =   4. - 8.*rb
      case(3)
        rb2 = rb*rb

        Nr_nz(1) = - 5.5 + 18. *rb - 13.5*rb2
        Nr_nz(2) =   1.  -  9. *rb + 13.5*rb2
        Nr_nz(3) =   9.  - 45. *rb + 40.5*rb2
        Nr_nz(4) = - 4.5 + 36. *rb - 40.5*rb2
      case default
        write(*,*) 'wrong! Get_lagrange_bedge_Nr_nz'
        stop
      end select

end subroutine Get_lagrange_bedge_Nr_nz

end module FEMEntry_ElementType_Lagrange_mod
