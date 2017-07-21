module FEM_Equation_Poisson_mod
      use KindDefinition_mod, only : DP
      use Constant_mod
      use FEMEntry_Generic_mod
      use FEM_Generic_mod
      use Common_Overload_mod
      use Common_Function_mod
      use GlobalVariableForDebug_mod
      implicit none
      save

      type, extends (equation_t) :: poisson_t
      contains
        procedure, pass(this) :: Reset_equation => Reset_poisson_equation
        procedure, nopass     :: equation_element_polynomial_order => poisson_equation_element_polynomial_order
        procedure, nopass     :: equation_bedge_polynomial_order => poisson_equation_bedge_polynomial_order
        procedure, pass(this) :: Initialize_solution => Initialize_poisson_solution
        procedure, pass(this) :: Calculate_dt_local => Calculate_poisson_dt_local
        procedure, pass(this) :: Calculate_interior_term => Calculate_interior_poisson_term
        procedure, pass(this) :: Calculate_interior_term_linearization => Calculate_interior_poisson_term_linearization
        procedure, pass(this) :: Calculate_interior_term_linearization_hand => Calculate_interior_poisson_term_linearization_hand
      end type poisson_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reset_poisson_equation( this )
      implicit none
      class(poisson_t), intent(out) :: this

      this%nequation = 1

end subroutine Reset_poisson_equation
!=====================================================================
!
!=====================================================================
function poisson_equation_element_polynomial_order( basis_order )  result ( polynomial_order )
      implicit none
      integer, intent(in) :: basis_order
      integer             :: polynomial_order

      polynomial_order = 2*( basis_order - 1 )

end function poisson_equation_element_polynomial_order
!=====================================================================
!
!=====================================================================
function poisson_equation_bedge_polynomial_order( basis_order )  result ( polynomial_order )
      implicit none
      integer, intent(in) :: basis_order
      integer             :: polynomial_order

      polynomial_order = basis_order

end function poisson_equation_bedge_polynomial_order
!=====================================================================
!
!=====================================================================
subroutine Initialize_poisson_solution ( this, FEM )
      implicit none
      class(poisson_t), intent(in)    :: this
      class(FEM_t),     intent(inout) :: FEM

      integer :: i

      do i = 1, FEM%nfbasis
        FEM%fbasis(i)%q(1,:) = 0.
      end do

end subroutine Initialize_poisson_solution
!=====================================================================
!
!=====================================================================
subroutine Calculate_poisson_dt_local ( this, FEM, time_control )
      implicit none
      class(poisson_t),     intent(in)    :: this
      class(FEM_t),         intent(inout) :: FEM
      type(time_control_t), intent(in)    :: time_control

      FEM % fbasis(:) % dt_local = 100000.

end subroutine Calculate_poisson_dt_local
!=====================================================================
!
!=====================================================================
subroutine Calculate_interior_poisson_term ( this, res, element, fbasis_nfb, time_control, r, s )
      implicit none
      class(poisson_t),                         intent(in)  :: this
      real(DP),             dimension(:,:),     intent(out) :: res
      class(element_t),                         intent(in)  :: element
      type(fbasis_t),       dimension(:),       intent(in)  :: fbasis_nfb
      type(time_control_t),                     intent(in)  :: time_control
      real(DP),                                 intent(in)  :: r, s

      integer :: nfb
      integer :: ntime_low, ntime_up, itime_space
      logical :: Include_time
      logical :: Use_dt_local
      real(DP) :: cfl
      real(DP) :: dt_global
      real(DP) :: ratio_space
      real(DP) :: Jel_area
      real(DP) :: dt
      real(DP), dimension(nfb_MAX) :: N, Nx, Ny
      real(DP), dimension(1,ntimeStep_LOWEST:ntimeStep_UPPEST) :: q_ntime
      real(DP), dimension(1) :: q, qx, qy
      real(DP), dimension(1) :: qt
      real(DP) :: source
      integer :: i, itime

      Include_time  = time_control % Include_time
      ntime_low     = time_control % ntimeStep_time_low
      ntime_up      = time_control % ntimeStep_time_up
      Use_dt_local  = time_control % Use_dt_local
      cfl           = time_control % cfl
      dt_global     = time_control % dt_global
      itime_space   = time_control % itimeStep_space
      ratio_space   = time_control % ratio_space

      nfb = element%nfb

      call element%Get_element_N ( N, r, s )
      call element%Calculate_element_Jel_area_Nx_Ny ( Jel_area, Nx, Ny, fbasis_nfb(:) % x, fbasis_nfb(:) % y, r, s )

      ! {q}, qx := {dq/dx}, qy := {dq/dy}
      q (1) = dot_product_my ( N , fbasis_nfb(:) % q(1,itime_space), nfb )
      qx(1) = dot_product_my ( Nx, fbasis_nfb(:) % q(1,itime_space), nfb )
      qy(1) = dot_product_my ( Ny, fbasis_nfb(:) % q(1,itime_space), nfb )

      ! set {res} to 0
      res(1,1:nfb) = 0.

      ! qt := {dq/dt}
      if( Include_time ) then
        ! q_ntime
        q_ntime(1,itime_space) = q(1)

        do itime = ntime_low, ntime_up
          if( itime == itime_space )  cycle

          q_ntime (1,itime) = dot_product_my ( N , fbasis_nfb(:) % q(1,itime), nfb )
        end do

        ! dt
        if ( Use_dt_local ) then
          dt = cfl * dot_product_my ( N, fbasis_nfb(:) % dt_local, nfb )
        else
          dt = dt_global
        end if

        ! qt
        qt(:) = time_control % time_discretization ( q_ntime, dt, 1 )
      else
        qt(:) = 0.
      end if

      ! term 1: N*qt
      if ( Include_time ) then
        do i = 1, nfb
          res(1,i) = res(1,i) + N(i)*qt(1)
        end do
      end if

      ! term 2: Nx*qx + Ny*qy
      res(1,1:nfb) = res(1,1:nfb) + ( ratio_space * Nx(1:nfb) ) * qx(1) + ( ratio_space * Ny(1:nfb) ) * qy(1)

      ! term 3: N*source
      source = 1.
      res(1,1:nfb) = res(1,1:nfb) + ratio_space*N(1:nfb)*source

      ! multiply Jacobian
      res(1,1:nfb) = Jel_area*res(1,1:nfb)

end subroutine Calculate_interior_poisson_term
!=====================================================================
!
!=====================================================================
subroutine Calculate_interior_poisson_term_linearization ( this, aa, element, fbasis_nfb, time_control, r, s )
      implicit none
      class(poisson_t),                         intent(in)  :: this
      real(DP),             dimension(:,:,:,:), intent(out) :: aa
      class(element_t),                         intent(in)  :: element
      type(fbasis_t),       dimension(:),       intent(in)  :: fbasis_nfb
      type(time_control_t),                     intent(in)  :: time_control
      real(DP),                                 intent(in)  :: r, s

      integer :: nfb
      integer :: ntime_low, ntime_up, itime_space
      logical :: Include_time, Include_space
      logical :: Use_dt_local
      real(DP) :: cfl
      real(DP) :: dt_global
      real(DP) :: ratio_space
      real(DP) :: Jel_area
      real(DP) :: dt
      real(DP), dimension(nfb_MAX) :: N, Nx, Ny
      type(dType), dimension(1,ntimeStep_LOWEST:ntimeStep_UPPEST) :: q_ntime
      type(dType), dimension(1) :: q, qx, qy
      type(dType), dimension(1) :: qt
      type(dType) :: source
      type(dType), dimension(1,nfb_MAX) :: res
      real(DP), dimension(1,1,nfb_MAX) :: dres_q, dres_qx, dres_qy
      integer :: i, j, itime

      Include_time  = time_control % Include_time
      ntime_low     = time_control % ntimeStep_time_low
      ntime_up      = time_control % ntimeStep_time_up
      Use_dt_local  = time_control % Use_dt_local
      cfl           = time_control % cfl
      dt_global     = time_control % dt_global
      Include_space = time_control % Include_space
      itime_space   = time_control % itimeStep_space
      ratio_space   = time_control % ratio_space

      nfb = element%nfb

      call element%Get_element_N ( N, r, s )
      call element%Calculate_element_Jel_area_Nx_Ny ( Jel_area, Nx, Ny, fbasis_nfb(:) % x, fbasis_nfb(:) % y, r, s )

! - get res' derivative w.r.t. to (q, qx, qy) using overload
      ! initialize {q}, qx := {dq/dx}, qy := {dq/dy}
      overloadLength = 3

      call initialize( q (1), 1, dot_product_my ( N , fbasis_nfb(:) % q(1,itime_space), nfb ) )
      call initialize( qx(1), 2, dot_product_my ( Nx, fbasis_nfb(:) % q(1,itime_space), nfb ) )
      call initialize( qy(1), 3, dot_product_my ( Ny, fbasis_nfb(:) % q(1,itime_space), nfb ) )

      ! set {res} to 0
      res(1,1:nfb) = 0.

      ! qt := {dq/dt}
      if( Include_time ) then
        ! q_ntime
        q_ntime(1,itime_space) = q(1)

        do itime = ntime_low, ntime_up
          if( itime == itime_space )  cycle

          q_ntime (1,itime) = dot_product_my ( N , fbasis_nfb(:) % q(1,itime), nfb )
        end do

        ! dt
        if ( Use_dt_local ) then
          dt = cfl * dot_product_my ( N, fbasis_nfb(:) % dt_local, nfb )
        else
          dt = dt_global
        end if

        ! qt
        qt(:) = time_control % time_discretization ( q_ntime, dt, 1 )
      else
        qt(:) = 0.
      end if

      ! term 1: N*qt
      if ( Include_time ) then
        do i = 1, nfb
          res(1,i) = res(1,i) + N(i)*qt(1)
        end do
      end if

      ! term 2: Nx*qx + Ny*qy
      res(1,1:nfb) = res(1,1:nfb) + ( ratio_space * Nx(1:nfb) ) * qx(1) + ( ratio_space * Ny(1:nfb) ) * qy(1)

      ! term 3: N*source
      source = 1.
      res(1,1:nfb) = res(1,1:nfb) + ratio_space*N(1:nfb)*source

      ! multiply Jacobian
      res(1,1:nfb) = Jel_area*res(1,1:nfb)

! - take answer out
      do i = 1, nfb
        dres_q (1,1,i) = res(1,i)%deriv(1)
        dres_qx(1,1,i) = res(1,i)%deriv(2)
        dres_qy(1,1,i) = res(1,i)%deriv(3)
      end do

! - get res' derivative w.r.t qfb using chain rule
      ! (q, qx, qy)'s derivative w.r.t. qfb is (N, Nx, Ny)
      do i = 1, nfb
      do j = 1, nfb
        aa(1,1,i,j) = dres_q(1,1,i)*N(j) + dres_qx(1,1,i)*Nx(j) + dres_qy(1,1,i)*Ny(j)
      end do
      end do


      ! >>> debug
     !write(*,*) 'aaaaaaaaa, r, s', r, s
     !do i = 1, nfb
     !  write(*,100)'aa', aa(1,1,i,1:nfb) 
     !end do
      ! <<< debug

100   format(a,20(es16.5))      

end subroutine Calculate_interior_poisson_term_linearization
!=====================================================================
!
!=====================================================================
subroutine Calculate_interior_poisson_term_linearization_hand ( this, aa, element, fbasis_nfb, time_control, r, s )
      implicit none
      class(poisson_t),                         intent(in)  :: this
      real(DP),             dimension(:,:,:,:), intent(out) :: aa
      class(element_t),                         intent(in)  :: element
      type(fbasis_t),       dimension(:),       intent(in)  :: fbasis_nfb
      type(time_control_t),                     intent(in)  :: time_control
      real(DP),                                 intent(in)  :: r, s

      real(DP), dimension(nfb_MAX) :: Nx, Ny
      real(DP) :: Jel_area
      integer :: nfb
      integer :: i, j

      nfb = element%nfb

      aa(1,1,1:nfb,1:nfb) = 0.

      call element%Calculate_element_Jel_area_Nx_Ny ( Jel_area, Nx, Ny, fbasis_nfb(:) % x, fbasis_nfb(:) % y, r, s )

      ! term 1
      do j = 1, nfb
      do i = 1, nfb
        aa(1,1,i,j) = aa(1,1,i,j) + Nx(i)*Nx(j) + Ny(i)*Ny(j)
      end do
      end do

      ! multiply Jacobian
      aa(1,1,1:nfb,1:nfb) = Jel_area*aa(1,1,1:nfb,1:nfb)


      ! >>> debug
     !write(*,*) 'bbbbbbb, r, s', r, s
     !do i = 1, nfb
     !  write(*,100)'aa', aa(1,1,i,1:nfb) 
     !end do
      ! <<< debug


100   format(a,20(es16.5))      


end subroutine Calculate_interior_poisson_term_linearization_hand

end module FEM_Equation_Poisson_mod
