module FEMEntry_Generic_mod
      use KindDefinition_mod, only : DP
      use Constant_mod
      use Common_Function_mod
      implicit none
      save

! - element, bedge, basis function, node
      ! element
      type :: element_t
        integer :: basis_order                                       ! basis functions order
        integer :: nfb                                               ! # of basis functions
        integer :: ifb(nfb_MAX)                                      ! ID of basis functions in element
        integer :: nnd                                               ! # of nodes in the element
        integer :: ind(nnd_MAX)                                      ! ID of nodes
        real(DP) :: area_rs                                          ! area of element in isoparametric coordinate
      contains
        procedure, pass(this) :: Initialize_element
        procedure, pass(this) :: Get_element_N
        procedure, pass(this) :: Get_element_Nr_Ns
        procedure, pass(this) :: Get_element_Nrr_Nrs_Nss
        procedure, pass(this) :: Calculate_element_Jel_area_Nx_Ny
        procedure, pass(this) :: Calculate_element_Nxx_Nxy_Nyy
      end type element_t

      ! bedge
      type :: bedge_t
        integer :: basis_order                                         ! basis functions order
        integer :: nfbnz                                               ! # of non-zero basis functions on bedge
        integer :: ifbnz(nfbnz_MAX)
        integer :: nfb                                                 ! # of basis functions of the element
        integer :: ifb(nfb_MAX)
        integer :: iel                                                 ! belonging element
        integer :: elementSide                                         ! which side of belonging element
        integer :: bcType                                              ! boundary condition type
        real(DP) :: length_r                                           ! length of bedge in isoparametric coordinate
      contains
        procedure, pass(this) :: Initialize_bedge
        procedure, pass(this) :: Get_bedge_to_element_rs
        procedure, pass(this) :: Get_bedge_N_nz
        procedure, pass(this) :: Get_bedge_Nr_nz
        procedure, pass(this) :: Calculate_bedge_xynorm_Jsf_length
      end type bedge_t

      ! basis function
      type :: fbasis_t
        real(DP) :: x, y                                                     ! coordinate
        real(DP) :: q(nequation_MAX,ntimeStep_LOWEST:ntimeStep_UPPEST)       ! q
        real(DP) :: dt_local                                                 ! local time step of time marching for steady-state problem
        real(DP) :: e                                                        ! for debug
        real(DP) :: evector(nequation_MAX)                                   ! for debug
        real(DP) :: ematrix(nequation_MAX,nequation_MAX)                     ! for debug
      end type fbasis_t

      ! node
      type :: node_t
        real(DP) :: x, y                                             ! node coordinate
      end type node_t

! - generic mapping of 2D and 1D quadrature
      type :: quadrature_2d_t
        integer :: polynomial_order                                  ! order of polynomial to be integrated
        integer :: nqd                                               ! # of quadrature
        real(DP) :: rqd(nqd_2d_MAX), sqd(nqd_2d_MAX)                 ! (r,s) of quadrature
        real(DP) :: wqd(nqd_2d_MAX)                                  ! weight of quadrature
      contains
        procedure, pass(this) :: Reset_quadrature_2d
      end type quadrature_2d_t

      type :: quadrature_1d_t
        integer :: polynomial_order                                  ! order of polynomial to be integrated
        integer :: nqd                                               ! # of quadrature
        real(DP) :: rqd(nqd_1d_MAX)                                  ! (r) of quadrature
        real(DP) :: wqd(nqd_1d_MAX)                                  ! weight of quadrature
      contains
        procedure, pass(this) :: Reset_quadrature_1d
      end type quadrature_1d_t

! - generic mapping of local nodes
      type :: local_node_t
        integer :: basis_order                                       ! basis functions order
        integer :: nnd                                               ! # of nodes in the element
        real(DP) :: rnd(nnd_MAX), snd(nnd_MAX)                       ! (r,s) of nodes
      contains
        procedure, pass(this) :: Reset_local_node
      end type local_node_t

contains

!=====================================================================
!
!=====================================================================
subroutine Initialize_element ( this, basis_order )
      implicit none
      class(element_t), intent(inout) :: this
      integer,          intent(in)    :: basis_order
end subroutine Initialize_element
!=====================================================================
!
!=====================================================================
subroutine Get_element_N( this, N, r, s )
      implicit none
      class(element_t),               intent(in)  :: this
      real(DP),         dimension(:), intent(out) :: N
      real(DP),                       intent(in)  :: r, s
end subroutine Get_element_N
!=====================================================================
!
!=====================================================================
subroutine Get_element_Nr_Ns( this, Nr, Ns, r, s )
      implicit none
      class(element_t),               intent(in)  :: this
      real(DP),         dimension(:), intent(out) :: Nr, Ns
      real(DP),                       intent(in)  :: r, s
end subroutine Get_element_Nr_Ns
!=====================================================================
!
!=====================================================================
subroutine Get_element_Nrr_Nrs_Nss ( this, Nrr, Nrs, Nss, r, s )
      implicit none
      class(element_t),               intent(in)  :: this
      real(DP),         dimension(:), intent(out) :: Nrr, Nrs, Nss
      real(DP),                       intent(in)  :: r, s
end subroutine Get_element_Nrr_Nrs_Nss
!=====================================================================
!
!=====================================================================
subroutine Calculate_element_Jel_area_Nx_Ny ( this, Jel_area, Nx, Ny, xfb, yfb, r, s )
      implicit none
      class(element_t),               intent(in)  :: this
      real(DP),                       intent(out) :: Jel_area
      real(DP),         dimension(:), intent(out) :: Nx, Ny
      real(DP),         dimension(:), intent(in)  :: xfb, yfb
      real(DP),                       intent(in)  :: r, s

      integer :: nfb
      real(DP), dimension(nfb_MAX) :: Nr, Ns
      real(DP) :: xr, xs, yr, ys
      real(DP) :: rx, ry, sx, sy
      real(DP) :: Jel, Jinv

      nfb = this % nfb

      call this % Get_element_Nr_Ns( Nr, Ns, r, s )

      xr = dot_product_my ( xfb, Nr, nfb )
      xs = dot_product_my ( xfb, Ns, nfb )
      yr = dot_product_my ( yfb, Nr, nfb )
      ys = dot_product_my ( yfb, Ns, nfb )

      Jel = xr*ys - xs*yr
      Jinv = 1./Jel

      Jel_area = Jel * this % area_rs

      rx =   ys*Jinv
      ry = - xs*Jinv
      sx = - yr*Jinv
      sy =   xr*Jinv

      Nx(1:nfb) = Nr(1:nfb)*rx + Ns(1:nfb)*sx
      Ny(1:nfb) = Nr(1:nfb)*ry + Ns(1:nfb)*sy

end subroutine Calculate_element_Jel_area_Nx_Ny
!=====================================================================
!
!=====================================================================
subroutine Calculate_element_Nxx_Nxy_Nyy ( this, Nxx, Nxy, Nyy, xfb, yfb, r, s )
      implicit none
      class(element_t),               intent(in)  :: this
      real(DP),         dimension(:), intent(out) :: Nxx, Nxy, Nyy
      real(DP),         dimension(:), intent(in)  :: xfb, yfb
      real(DP),                       intent(in)  :: r, s

      integer :: nfb
      real(DP), dimension(nfb_MAX) :: Nr, Ns
      real(DP), dimension(nfb_MAX) :: Nrr, Nrs, Nss
      real(DP), dimension(nfb_MAX) :: Nxr, Nxs, Nyr, Nys
      real(DP) :: xr, xs, yr, ys
      real(DP) :: rx, ry, sx, sy
      real(DP) :: xrr, xrs, xss, yrr, yrs, yss
      real(DP) :: rxr, rxs, ryr, rys, sxr, sxs, syr, sys
      real(DP) :: Jel, Jinv, J2inv
      real(DP) :: Jr, Js

      nfb = this % nfb

      ! 1st derivatives
      call this % Get_element_Nr_Ns( Nr, Ns, r, s )

      xr = dot_product_my ( xfb, Nr, nfb )
      xs = dot_product_my ( xfb, Ns, nfb )
      yr = dot_product_my ( yfb, Nr, nfb )
      ys = dot_product_my ( yfb, Ns, nfb )

      Jel = xr*ys - xs*yr
      Jinv = 1./Jel

      rx =   ys*Jinv
      ry = - xs*Jinv
      sx = - yr*Jinv
      sy =   xr*Jinv

      ! 2nd derivatives
      call this % Get_element_Nrr_Nrs_Nss( Nrr, Nrs, Nss, r, s )

      xrr = dot_product_my ( xfb, Nrr, nfb )
      xrs = dot_product_my ( xfb, Nrs, nfb )
      xss = dot_product_my ( xfb, Nss, nfb )

      yrr = dot_product_my ( yfb, Nrr, nfb )
      yrs = dot_product_my ( yfb, Nrs, nfb )
      yss = dot_product_my ( yfb, Nss, nfb )

      Jr = xrr*ys + xr*yrs - xrs*yr - xs*yrr
      Js = xrs*ys + xr*yss - xss*yr - xs*yrs

      J2inv = Jinv*Jinv

      rxr =   Jinv*yrs - J2inv*Jr*ys
      rxs =   Jinv*yss - J2inv*Js*ys
      ryr = - Jinv*xrs + J2inv*Jr*xs
      rys = - Jinv*xss + J2inv*Js*xs
      sxr = - Jinv*yrr + J2inv*Jr*yr
      sxs = - Jinv*yrs + J2inv*Js*yr
      syr =   Jinv*xrr - J2inv*Jr*xr
      sys =   Jinv*xrs - J2inv*Js*xr

      Nxr(1:nfb) = Nrr(1:nfb)*rx + Nr(1:nfb)*rxr + Nrs(1:nfb)*sx + Ns(1:nfb)*sxr
      Nxs(1:nfb) = Nrs(1:nfb)*rx + Nr(1:nfb)*rxs + Nss(1:nfb)*sx + Ns(1:nfb)*sxs
      Nyr(1:nfb) = Nrr(1:nfb)*ry + Nr(1:nfb)*ryr + Nrs(1:nfb)*sy + Ns(1:nfb)*syr
      Nys(1:nfb) = Nrs(1:nfb)*ry + Nr(1:nfb)*rys + Nss(1:nfb)*sy + Ns(1:nfb)*sys

      Nxx(1:nfb) = Nxr(1:nfb)*rx + Nxs(1:nfb)*sx
      Nxy(1:nfb) = Nxr(1:nfb)*ry + Nxs(1:nfb)*sy
      Nyy(1:nfb) = Nyr(1:nfb)*ry + Nys(1:nfb)*sy

end subroutine Calculate_element_Nxx_Nxy_Nyy
!=====================================================================
!
!=====================================================================
subroutine Initialize_bedge ( this, basis_order )
      implicit none
      class(bedge_t), intent(inout) :: this
      integer,        intent(in)    :: basis_order
end subroutine Initialize_bedge
!=====================================================================
!
!=====================================================================
subroutine Get_bedge_to_element_rs ( this, r, s, rb )
      implicit none
      class(bedge_t), intent(in)  :: this
      real(DP),       intent(out) :: r, s
      real(DP),       intent(in)  :: rb
end subroutine Get_bedge_to_element_rs
!=====================================================================
!
!=====================================================================
subroutine Get_bedge_N_nz ( this, N_nz, rb )
      implicit none
      class(bedge_t),               intent(in)  :: this
      real(DP),       dimension(:), intent(out) :: N_nz
      real(DP),                     intent(in)  :: rb
end subroutine Get_bedge_N_nz
!=====================================================================
!
!=====================================================================
subroutine Get_bedge_Nr_nz ( this, Nr_nz, rb )
      implicit none
      class(bedge_t),               intent(in)  :: this
      real(DP),       dimension(:), intent(out) :: Nr_nz
      real(DP),                     intent(in)  :: rb
end subroutine Get_bedge_Nr_nz
!=====================================================================
!
!=====================================================================
subroutine Calculate_bedge_xynorm_Jsf_length ( this, xnorm, ynorm, Jsf_length, xfbnz, yfbnz, rb )
      implicit none
      class(bedge_t),               intent(in)  :: this
      real(DP),                     intent(out) :: xnorm, ynorm
      real(DP),                     intent(out) :: Jsf_length
      real(DP),       dimension(:), intent(in)  :: xfbnz, yfbnz
      real(DP),                     intent(in)  :: rb

      integer :: nfbnz
      real(DP), dimension(nfbnz_MAX) :: Nr_nz
      real(DP) :: xr, yr
      real(DP) :: Jsf

      nfbnz = this % nfbnz

      call this % Get_bedge_Nr_nz ( Nr_nz, rb )

      xr = dot_product_my ( xfbnz, Nr_nz, nfbnz )
      yr = dot_product_my ( yfbnz, Nr_nz, nfbnz )

      Jsf = sqrt ( xr*xr + yr*yr )
      Jsf_length = Jsf * this % length_r

      xnorm =  yr/Jsf
      ynorm = -xr/Jsf

end subroutine Calculate_bedge_xynorm_Jsf_length
!=====================================================================
!
!=====================================================================
subroutine Reset_quadrature_2d( this, polynomial_order )
      implicit none
      class(quadrature_2d_t), intent(inout) :: this
      integer, intent(in) :: polynomial_order
end subroutine Reset_quadrature_2d
!=====================================================================
!
!=====================================================================
subroutine Reset_quadrature_1d( this, polynomial_order )
      implicit none
      class(quadrature_1d_t), intent(inout) :: this
      integer, intent(in) :: polynomial_order
end subroutine Reset_quadrature_1d
!=====================================================================
!
!=====================================================================
subroutine Reset_local_node( this, basis_order )
      implicit none
      class(local_node_t), intent(inout) :: this
      integer, intent(in) :: basis_order
end subroutine Reset_local_node

end module FEMEntry_Generic_mod
