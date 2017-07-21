module FEM_Equation_Euler_mod
      use KindDefinition_mod, only : DP
      use Constant_mod
      use FEMEntry_Generic_mod
      use FEM_Generic_mod
      use Euler_mod
      use Common_Overload_mod
      use Common_Function_mod
      use GlobalVariableForDebug_mod
      implicit none
      save

      type, extends (equation_t) :: euler_equation_t
        real(DP) :: q_farfield(4)

        logical :: Do_shock_capture
        real(DP) :: k_shock = 1.              ! for weak shock
       !real(DP) :: k_shock = 0.1             ! for strong shock
        real(DP) :: epsilon_shock = 0.1       ! epislon for smoothing shock indicator

        logical :: Do_2nd_derivative
      contains
        procedure, pass(this) :: Reset_equation => Reset_euler_equation
        procedure, nopass     :: equation_element_polynomial_order => euler_equation_element_polynomial_order
        procedure, nopass     :: equation_bedge_polynomial_order => euler_equation_bedge_polynomial_order
        procedure, pass(this) :: Initialize_solution => Initialize_euler_solution
        procedure, pass(this) :: Enforce_essential_boundary_condition => Enforce_euler_essential_boundary_condition
        procedure, pass(this) :: Calculate_dt_local => Calculate_euler_dt_local
        procedure, pass(this) :: Calculate_dt_local_2 => Calculate_euler_dt_local_2
        procedure, pass(this) :: Calculate_interior_term => Calculate_euler_interior_term
        procedure, pass(this) :: Calculate_interior_term_linearization => Calculate_euler_interior_term_linearization
        procedure, pass(this) :: Calculate_natural_boundary_term => Calculate_euler_natural_boundary_term
        procedure, pass(this) :: Calculate_natural_boundary_term_linearization => Calculate_euler_natural_boundary_term_linearization
        procedure, pass(this) :: max_omega_for_dq => euler_max_omega_for_dq
        procedure, pass(this) :: Check_solution => Check_euler_solution
      end type euler_equation_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reset_euler_equation( this )
      implicit none
      class(euler_equation_t), intent(out) :: this

      real(DP) :: onedeg = acos(0.)/90.

      this % nequation = 4

     !this % q_farfield(1:4) = [ 1., 0.6, 0.  , 1. ]
     !this % q_farfield(1:4) = [ 1., 0.6, 0.05, 1. ]
     !this % q_farfield(1:4) = [ 1., 0.8, 0.  , 1. ]
     !this % q_farfield(1:4) = [ 1., 0.8, 0.05, 1. ]
     !this % q_farfield(1:4) = [ 1., 2. , 0.  , 1. ]
     !this % q_farfield(1:4) = [ 1., 3. , 0.  , 1. ]
     !this % q_farfield(1:4) = [ 1., 4. , 0.  , 1. ]
      this % q_farfield(1:4) = [ 1., 0.5*cos(1.25*onedeg), 0.5*sin(1.25*onedeg)  , 1. ]
     !this % q_farfield(1:4) = [ 1., 0.8*cos(1.25*onedeg), 0.8*sin(1.25*onedeg)  , 1. ]

      this % Do_shock_capture = .false.


      if ( this % Do_shock_capture ) then
        this % Do_2nd_derivative = .true.
      else
        this % Do_2nd_derivative = .false.
      end if

      ! >>> debug
     !write(*,*) 'debug! Reset_euler_equation'
     !this % Do_shock_capture = .true.
     !this % Do_2nd_derivative = .true.
      ! <<< debug


      write(*,100) 'Reset_euler_equation: q_farfield:', this % q_farfield(1:4)
      write(*,*)   'Reset_euler_equation: shock capture', this % Do_shock_capture
      write(*,*)   'Reset_euler_equation: k_shock', this % k_shock
      write(*,*)   'Reset_euler_equation: epsilon_shock', this % epsilon_shock
      write(*,*)   'Reset_euler_equation: 2nd derivative', this % Do_2nd_derivative

100   format(a,20(es16.5))

end subroutine Reset_euler_equation
!=====================================================================
!
!=====================================================================
function euler_equation_element_polynomial_order( basis_order )  result ( polynomial_order )
      implicit none
      integer, intent(in) :: basis_order
      integer             :: polynomial_order

      polynomial_order = 4*basis_order - 1

end function euler_equation_element_polynomial_order
!=====================================================================
!
!=====================================================================
function euler_equation_bedge_polynomial_order( basis_order )  result ( polynomial_order )
      implicit none
      integer, intent(in) :: basis_order
      integer             :: polynomial_order

      polynomial_order = 4*basis_order

end function euler_equation_bedge_polynomial_order
!=====================================================================
!
!=====================================================================
subroutine Initialize_euler_solution ( this, FEM )
      implicit none
      class(euler_equation_t), intent(in)    :: this
      class(FEM_t),            intent(inout) :: FEM

      integer :: i, ii

      do i = 1, FEM % nfbasis
        do ii = 1, 4
          FEM % fbasis(i) % q(ii,:) = this % q_farfield(ii)
        end do

        ! >>> debug
       !write(*,*) 'debug! Initialize_euler_solution'
       !FEM % fbasis(i) % q(:,1) = manufactured_solution ( FEM % fbasis(i) % x, FEM % fbasis(i) % y )
       !FEM % fbasis(i) % q(:,0) = FEM % fbasis(i) % q(:,1)
        ! <<< debug
      end do


end subroutine Initialize_euler_solution
!=====================================================================
!
!=====================================================================
subroutine Enforce_euler_essential_boundary_condition ( this, FEM, itimeStep )
      implicit none
      class(euler_equation_t), intent(in)    :: this
      class(FEM_t),            intent(inout) :: FEM
      integer,                 intent(in)    :: itimeStep

      integer :: i, ibedge, ifb
      integer :: bcType

      do ibedge = 1, FEM % nbedge
        bcType = FEM % bedge(ibedge) % bcType
        if ( bcType /= 0 )  cycle

        do i = 1, FEM % bedge(ibedge) % nfbnz
          ifb = FEM % bedge(ibedge) % ifbnz(i)

          FEM % fbasis(ifb) % q(1:4,itimeStep) = this % q_farfield(1:4)
        end do
      end do

end subroutine Enforce_euler_essential_boundary_condition
!=====================================================================
!
!=====================================================================
subroutine Calculate_euler_dt_local ( this, FEM, time_control )
      implicit none
      class(euler_equation_t), intent(in)    :: this
      class(FEM_t),            intent(inout) :: FEM
      type(time_control_t),    intent(in)    :: time_control

      integer :: nq, nfbasis
      integer :: itime_space
      logical, save :: Is_1st_time = .true.
      real(DP), dimension(:), allocatable, save :: v1, v2
      integer :: basis_order, nfb
      real(DP) :: Jel_area
      type(fbasis_t), dimension(nfb_MAX) :: fbasis_nfb
      real(DP),       dimension(nfb_MAX) :: N, Nx, Ny
      integer,        dimension(nfb_MAX) :: ifb_nfb
      real(DP),       dimension(nfb_MAX) :: dt_local
      real(DP),       dimension(4) :: q
      real(DP) :: r, s, w
      real(DP) :: tmp
      integer :: i, k, ii, iel, ifb


      nq = this%nequation

      itime_space = time_control % itimeStep_space

      if ( Is_1st_time ) then
        allocate ( v1( FEM % nfbasis ) )
        allocate ( v2( FEM % nfbasis ) )
      end if

      Is_1st_time = .false.

      v1(:) = 0.
      v2(:) = 0.

      do iel = 1, FEM % nelement
        basis_order = FEM % element(iel) % basis_order

        call FEM % Get_element_fbasis_ifb_nfb ( fbasis_nfb, ifb_nfb, nfb, iel )

        do k = 1, FEM % quadrature_2d (basis_order) % nqd
          r = FEM % quadrature_2d (basis_order) % rqd(k)
          s = FEM % quadrature_2d (basis_order) % sqd(k)
          w = FEM % quadrature_2d (basis_order) % wqd(k)

          call FEM % element(iel) % Get_element_N ( N, r, s )
          call FEM % element(iel) % Calculate_element_Jel_area_Nx_Ny ( Jel_area, Nx, Ny, fbasis_nfb(:) % x, fbasis_nfb(:) % y, r, s )

          do ii = 1, nq
            q(ii) = dot_product_my ( N, fbasis_nfb(:) % q(ii,itime_space), nfb )
          end do

          do i = 1, nfb
           !tmp = abs( Nx(i)*q(2) + Ny(i)*q(3) ) + sqrt( q(4) )*sqrt( Nx(i)**2 + Ny(i)**2 )
            tmp = sqrt( q(2)**2 + q(3)**2 ) + sqrt( q(4) )*sqrt( Nx(i)**2 + Ny(i)**2 )

            ifb = ifb_nfb(i)
            v1(ifb) = v1(ifb) + w*Jel_area*tmp
            v2(ifb) = v2(ifb) + w*Jel_area
          end do
        end do
      end do

      FEM % fbasis(:) % dt_local = 0.5*v2(:)/v1(:)


      ! >>> debug
      FEM % fbasis(:) % e = FEM % fbasis(:) % dt_local

      call FEM % Output_FEM_solution_tecplot ( fpIO, 'output/dt_local.plt')
      ! <<< debug

100   format(a,20(es16.5))


end subroutine Calculate_euler_dt_local
!=======================================================================
!
!=======================================================================
subroutine Calculate_euler_dt_local_2 ( this, FEM, time_control )
      implicit none
      class(euler_equation_t), intent(in)    :: this
      class(FEM_t),            intent(inout) :: FEM
      type(time_control_t),    intent(in)    :: time_control

      integer :: nq
      integer :: itime_space
      logical, save :: Is_1st_time = .true.
      logical,  dimension(:), allocatable, save :: Is_p1_node
      real(DP), dimension(:), allocatable, save :: v1, v2
      integer :: basis_order, nfb, nfbnz
      type(fbasis_t), dimension(nfb_MAX)   :: fbasis_nfb
      type(fbasis_t), dimension(nfbnz_MAX) :: fbasis_nfbnz
      integer,        dimension(nfb_MAX)   :: ifb_nfb
      integer,        dimension(nfbnz_MAX) :: ifb_nfbnz
      integer :: ifb_p1(4)
      real(DP) :: dt1, dt2, dt3
      real(DP) :: q(4)
      real(DP) :: tmp
      real(dp) :: x(4), y(4)
      real(dp) :: xave,yave,xmid,ymid
      real(dp) :: xnorm,ynorm,rlen, area
      integer :: i, iel, ibedge, l, ifb


      nq = this%nequation

      itime_space = time_control % itimeStep_space

      if ( Is_1st_time ) then
        allocate ( v1( FEM % nfbasis ) )
        allocate ( v2( FEM % nfbasis ) )
        allocate ( Is_p1_node( FEM % nfbasis ) )

        v2(:) = 0.
        Is_p1_node(:) = .false.

        do iel = 1, FEM % nelement
          call FEM % Get_element_fbasis_ifb_nfb ( fbasis_nfb, ifb_nfb, nfb, iel )

          x(1:3) = fbasis_nfb(1:3) % x
          y(1:3) = fbasis_nfb(1:3) % y

          area = 0.5*( ( x(2) - x(1) )*( y(3) - y(1) ) - ( x(3) - x(1) )*( y(2) - y(1) ) )

          do i = 1, 3
            ifb = ifb_nfb(i)
            Is_p1_node(ifb) = .true.
            v2(ifb) = v2(ifb) + area/3.
          end do

        end do
      end if

      Is_1st_time = .false.

      ! dt_local for p1 node
      v1(:) = 0.

      do iel = 1, FEM % nelement
        call FEM % Get_element_fbasis_ifb_nfb ( fbasis_nfb, ifb_nfb, nfb, iel )

        ifb_p1(1:4) = [ 1, 2, 3, 1 ]

        x(1:4) = fbasis_nfb(ifb_p1(1:4)) % x
        y(1:4) = fbasis_nfb(ifb_p1(1:4)) % y

        xave = ( x(1) + x(2) + x(3) )/3.
        yave = ( y(1) + y(2) + y(3) )/3.

        do l = 1, 3
          xmid = 0.5*(x(l) + x(l+1))
          ymid = 0.5*(y(l) + y(l+1))

          xnorm =   yave - ymid
          ynorm = -(xave - xmid)
          rlen  = sqrt(xnorm*xnorm + ynorm*ynorm)
          xnorm = xnorm/rlen
          ynorm = ynorm/rlen

          q(1:nq) = 0.5*( fbasis_nfb(ifb_p1(l)) % q(1:nq,itime_space) + fbasis_nfb(ifb_p1(l+1)) % q(1:nq,itime_space)  )
          tmp = ( abs( xnorm*q(2) + ynorm*q(3) ) + sqrt(q(4)) )*rlen

          v1(ifb_nfb(ifb_p1(l)  )) = v1(ifb_nfb(ifb_p1(l)  )) + tmp
          v1(ifb_nfb(ifb_p1(l+1))) = v1(ifb_nfb(ifb_p1(l+1))) + tmp
        end do
      end do

      ! bedge
      do ibedge = 1, FEM % nbedge
        call FEM % Get_bedge_fbasis_ifb_nfbnz ( fbasis_nfbnz, ifb_nfbnz, nfbnz, ibedge )

        x(1:2) = fbasis_nfbnz(1:2) % x
        y(1:2) = fbasis_nfbnz(1:2) % y

        xnorm =   y(2) - y(1)
        ynorm = -(x(2) - x(1))
        rlen  = sqrt(xnorm*xnorm + ynorm*ynorm)
        xnorm = xnorm/rlen
        ynorm = ynorm/rlen
        rlen  = 0.5*rlen

        q(1:nq) = 0.5*( fbasis_nfbnz(1) % q(1:nq,itime_space) + fbasis_nfb(2) % q(1:nq,itime_space) )
        tmp = ( abs( xnorm*q(2) + ynorm*q(3) ) + sqrt(q(4)) )*rlen

        v1(ifb_nfbnz(1)) = v1(ifb_nfbnz(1)) + tmp
        v1(ifb_nfbnz(2)) = v1(ifb_nfbnz(2)) + tmp
      end do

      !
      do i = 1, FEM % nfbasis
        if( Is_p1_node(i) )  FEM % fbasis(i) % dt_local = v2(i)/v1(i)
      end do


      ! for non-p1 node
      do iel = 1, FEM % nelement
        call FEM % Get_element_fbasis_ifb_nfb ( fbasis_nfb, ifb_nfb, nfb, iel )

        select case ( FEM % element(iel) % basis_order )
        case(1)
          continue
        case(2)
          dt1 = FEM % fbasis( ifb_nfb(1) ) % dt_local
          dt2 = FEM % fbasis( ifb_nfb(2) ) % dt_local
          dt3 = FEM % fbasis( ifb_nfb(3) ) % dt_local

          FEM % fbasis( ifb_nfb(4) ) % dt_local = ( dt1 + dt2 )/2.
          FEM % fbasis( ifb_nfb(5) ) % dt_local = ( dt2 + dt3 )/2.
          FEM % fbasis( ifb_nfb(6) ) % dt_local = ( dt3 + dt1 )/2.
        case(3)
          dt1 = FEM % fbasis( ifb_nfb(1) ) % dt_local
          dt2 = FEM % fbasis( ifb_nfb(2) ) % dt_local
          dt3 = FEM % fbasis( ifb_nfb(3) ) % dt_local

          FEM % fbasis( ifb_nfb(4) ) % dt_local = ( 2.*dt1 + dt2 )/3.
          FEM % fbasis( ifb_nfb(5) ) % dt_local = ( 2.*dt2 + dt3 )/3.
          FEM % fbasis( ifb_nfb(6) ) % dt_local = ( 2.*dt3 + dt1 )/3.

          FEM % fbasis( ifb_nfb(7) ) % dt_local = ( dt1 + 2.*dt2 )/3.
          FEM % fbasis( ifb_nfb(8) ) % dt_local = ( dt2 + 2.*dt3 )/3.
          FEM % fbasis( ifb_nfb(9) ) % dt_local = ( dt3 + 2.*dt1 )/3.

          FEM % fbasis( ifb_nfb(10) ) % dt_local = ( dt1 + dt2 + dt3 )/3.
        case default
          write(*,*) 'wrong! Calculate_euler_dt_local_2'
          stop
        end select
      end do


      ! >>> debug
      FEM % fbasis(:) % e = FEM % fbasis(:) % dt_local

      call FEM % Output_FEM_solution_tecplot ( fpIO, 'output/dt_local_2.plt')
      ! <<< debug


end subroutine Calculate_euler_dt_local_2
!=====================================================================
!
!=====================================================================
subroutine Calculate_euler_interior_term ( this, res, element, fbasis_nfb, time_control, r, s )
      implicit none
      class(euler_equation_t),                 intent(in)  :: this
      real(DP),                dimension(:,:), intent(out) :: res
      class(element_t),                        intent(in)  :: element
      type(fbasis_t),          dimension(:),   intent(in)  :: fbasis_nfb
      type(time_control_t),                    intent(in)  :: time_control
      real(DP),                                intent(in)  :: r, s

      integer :: nfb, nq
      integer :: ntime_low, ntime_up, itime_space
      logical :: Include_time
      logical :: Use_dt_local
      real(DP) :: cfl
      real(DP) :: dt_global
      real(DP) :: ratio_space
      real(DP), dimension(nfb_MAX) :: N, Nx, Ny, Nxx, Nxy, Nyy
      real(DP) :: Jel_area
      real(DP) :: dt
      real(DP), dimension(4,ntimeStep_LOWEST:ntimeStep_UPPEST) :: q_ntime, cq_ntime
      real(DP), dimension(4)   :: q, qx, qy, qxx, qxy, qyy
      real(DP), dimension(4)   :: cqt, cqx, cqy
      real(DP), dimension(4)   :: f, g, fx, gy
      real(DP), dimension(4)   :: f_shock, g_shock, fx_shock, gy_shock
      real(DP), dimension(4)   :: v1, v2, v3, v4
      real(DP), dimension(4,4) :: Taui, Ajac, Bjac, M
      real(DP) :: h_length
      integer :: i, ii, itime

      ! for debug
      real(DP), dimension(4,4) :: emat


      nq  = this % nequation
      nfb = element % nfb

      Include_time  = time_control % Include_time
      ntime_low     = time_control % ntimeStep_time_low
      ntime_up      = time_control % ntimeStep_time_up
      Use_dt_local  = time_control % Use_dt_local
      cfl           = time_control % cfl
      dt_global     = time_control % dt_global
      itime_space   = time_control % itimeStep_space
      ratio_space   = time_control % ratio_space

      res(1:nq,1:nfb) = 0.

      call element % Get_element_N ( N, r, s )
      call element % Calculate_element_Jel_area_Nx_Ny ( Jel_area, Nx, Ny, fbasis_nfb(:) % x, fbasis_nfb(:) % y, r, s )

      if( this % Do_2nd_derivative ) then
        call element % Calculate_element_Nxx_Nxy_Nyy ( Nxx, Nxy, Nyy, fbasis_nfb(:) % x, fbasis_nfb(:) % y, r, s )
      else
        Nxx(1:nfb) = 0.
        Nxy(1:nfb) = 0.
        Nyy(1:nfb) = 0.
      end if

      ! {q}, qx := {dq/dx}, qy := {dq/dy}
      do ii = 1, nq
        q  (ii) = dot_product_my ( N  , fbasis_nfb(:) % q(ii,itime_space), nfb )
        qx (ii) = dot_product_my ( Nx , fbasis_nfb(:) % q(ii,itime_space), nfb )
        qy (ii) = dot_product_my ( Ny , fbasis_nfb(:) % q(ii,itime_space), nfb )

        if( this % Do_2nd_derivative ) then
          qxx(ii) = dot_product_my ( Nxx, fbasis_nfb(:) % q(ii,itime_space), nfb )
          qxy(ii) = dot_product_my ( Nxy, fbasis_nfb(:) % q(ii,itime_space), nfb )
          qyy(ii) = dot_product_my ( Nyy, fbasis_nfb(:) % q(ii,itime_space), nfb )
        else
          qxx(ii) = 0.
          qxy(ii) = 0.
          qyy(ii) = 0.
        end if
      end do


      if (Do_print) then
        write(*,100) 'r,s', r,s
        write(*,100) 'N  ', N  (1:nfb)
        write(*,100) 'Nx ', Nx (1:nfb)
        write(*,100) 'Ny ', Ny (1:nfb)
        write(*,100) 'Nxx', Nxx(1:nfb)
        write(*,100) 'Nxy', Nxy(1:nfb)
        write(*,100) 'Nyy', Nyy(1:nfb)
        write(*,100) 'Jel_area', Jel_area
        write(*,100) 'q',q(1:nq)
      end if


      ! cqt := {dcq/dt}
      if( Include_time ) then
        ! q_ntime
        q_ntime(:,itime_space) = q(:)

        do itime = ntime_low, ntime_up
          if( itime == itime_space )  cycle

          do ii = 1, nq
            q_ntime (ii,itime) = dot_product_my ( N(:) , fbasis_nfb(:) % q(ii,itime), nfb )
          end do
        end do

        ! cq_ntime
        do itime = ntime_low, ntime_up
          cq_ntime(:,itime) = conservative_variable_cq ( q_ntime(:,itime) )
        end do

        ! dt
        if ( Use_dt_local ) then
          dt = cfl * dot_product_my ( N, fbasis_nfb(:) % dt_local, nfb )
         !dt = cfl * 1./sum( sqrt( Nx(1:nfb)**2 + Ny(1:nfb)**2 ) ) / ( sqrt( q(4) ) + sqrt( q(2)**2 + q(3)**2 )  )/nfb
        else
          dt = dt_global
        end if

        ! cqt
        cqt(:) = time_control % time_discretization ( cq_ntime, dt, nq )
      else
        cqt(:) = 0.
      end if

      ! term 1: N*{dcq/dt}
      if ( Include_time ) then
        do i = 1, nfb
          res(:,i) = res(:,i) + N(i)*cqt(:)
        end do
      end if

1000  continue

      ! term 2: Nx*{F} + Ny*{G}
      f(:) = convective_flux_f ( q )
      g(:) = convective_flux_g ( q )


      ! >>> debug
      if( Do_print ) then
        write(*,*) 'convective f, g'
        write(*,100) 'f', f
        write(*,100) 'g', g
      end if
      ! <<< debug


     ! extra diffusive flux for shock capturing
      if( this % Do_shock_capture ) then
        h_length = 1./sum( sqrt( Nx(1:nfb)**2 + Ny(1:nfb)**2 ) )

       !call shock_capture_flux_f_g          ( f_shock, g_shock, q, qx, qy, h_length, this % k_shock )
        call shock_capture_flux_f_g_modified ( f_shock, g_shock, q, qx, qy, h_length, this % k_shock, this % epsilon_shock )

        f(:) = f(:) - f_shock(:)
        g(:) = g(:) - g_shock(:)

        ! >>> debug
        if( Do_print ) then
          write(*,*) 'f_shock, g_shock'
          write(*,100) 'f_shock', f_shock
          write(*,100) 'g_shock', g_shock
        end if
        ! <<< debug

      end if

      do i = 1, nfb
        res(:,i) = res(:,i) - ( ratio_space*Nx(i) )*f(:) - ( ratio_space*Ny(i) )*g(:)
      end do


      ! >>> debug
      if( Do_print ) then
        write(*,100) '-Nx*f - Ny*g'
        do i = 1, nfb
          write(*,100) '-Nx*f - Ny*g', - Nx(i)*f(:) - Ny(i)*g(:)
        end do
      end if
      ! <<< debug


2000  continue

      ! term 3: ( Nx*[A] + Ny*[B] )*[tau]*( {dcq/dt} + {dF/dx} + {dG/dy} ), where {dF/dx} = [A]*{dcq/dx}, {dG/dy} = [B]*{dcq/dy}
      v1(:) = 0.

      Taui(:,:) = get_tau_inverse ( q, Nx, Ny, nfb )

      Ajac(:,:) = convective_jacobian_A ( q )
      Bjac(:,:) = convective_jacobian_B ( q )

      M(:,:) = cq_derivative_wrt_q (q)

      ! >>> debug
      if( Do_print ) then
        do ii = 1, nq
          write(*,100) 'Ajac', Ajac(ii,:)
        end do

        do ii = 1, nq
          write(*,100) 'Bjac', Bjac(ii,:)
        end do

       !do ii = 1, nq
       !  write(*,100) 'Nx(3)*Ajac + Ny(3)*Bjac', Nx(3)*Ajac(ii,:) + Ny(3)*Bjac(ii,:)
       !end do

       !emat = matmul( Ajac, M )
       !do ii = 1, nq
       !  write(*,100) 'Ajac*M', emat(ii,:)
       !end do

       !emat = matmul( Bjac, M )
       !do ii = 1, nq
       !  write(*,100) 'Bjac*M', emat(ii,:)
       !end do
      end if
      ! <<< debug


      ! dcq/dt
      if( Include_time ) then
        v1(:) = v1(:) + cqt(:)
      end if


      cqx(:) = matrix_vector_multiply ( M, qx, nq, nq )
      cqy(:) = matrix_vector_multiply ( M, qy, nq, nq )

      fx(:) = matrix_vector_multiply ( Ajac, cqx, nq, nq )
      gy(:) = matrix_vector_multiply ( Bjac, cqy, nq, nq )

      v2(:) = fx(:) + gy(:)

      ! extra flux for shock capture
      if( this % Do_shock_capture ) then
        h_length = 1./sum( sqrt( Nx(1:nfb)**2 + Ny(1:nfb)**2 ) )

        call shock_capture_flux_fx_gy_modified ( fx_shock, gy_shock, q, qx, qy, qxx, qxy, qyy, h_length, this % k_shock, this % epsilon_shock )

        v2(:) = v2(:) - fx_shock(:) - gy_shock(:)
      end if

      v1(:) = v1(:) + ratio_space * v2(:)


      if (Do_print) then
        write(*,100) 'v1: cqt + fx + gy', v1(1:nq)
      end if


      ! inv([tau])*{v1}
      v2(:) = GaussJordan ( Taui, v1, nq )  ! Taui, v1 will be modified by GaussJordan


      if (Do_print) then
        write(*,100) 'v2: taui*{cqt + fx + gy}', v2(1:nq)
      end if


      v3(:) = matrix_vector_multiply ( Ajac, v2, nq, nq )
      v4(:) = matrix_vector_multiply ( Bjac, v2, nq, nq )


      if (Do_print) then
        write(*,100) 'v3: A*tau*{cqt + fx + gy}', v3(1:nq)
        write(*,100) 'v4: B*tau*{cqt + fx + gy}', v4(1:nq)
      end if


      do i = 1, nfb
        res(:,i) = res(:,i) + Nx(i)*v3(:) + Ny(i)*v4(:)
      end do

3000  continue

      ! multiply Jacobian
      res(:,1:nfb) = Jel_area*res(:,1:nfb)

100   format(a,20(es16.5))
200   format(a,i7,20(es16.5))

end subroutine Calculate_euler_interior_term
!=====================================================================
!
!=====================================================================
subroutine Calculate_euler_interior_term_linearization ( this, aa, element, fbasis_nfb, time_control, r, s )
      implicit none
      class(euler_equation_t),                     intent(in)  :: this
      real(DP),                dimension(:,:,:,:), intent(out) :: aa
      class(element_t),                            intent(in)  :: element
      type(fbasis_t),          dimension(:),       intent(in)  :: fbasis_nfb
      type(time_control_t),                        intent(in)  :: time_control
      real(DP),                                    intent(in)  :: r, s

      integer :: nfb, nq
      integer :: ntime_low, ntime_up, itime_space
      logical :: Include_time
      logical :: Use_dt_local
      real(DP) :: cfl
      real(DP) :: dt_global
      real(DP) :: ratio_space
      real(DP), dimension(nfb_MAX) :: N, Nx, Ny, Nxx, Nxy, Nyy
      real(DP) :: Jel_area
      real(DP) :: dt
      type(dType), dimension(4,ntimeStep_LOWEST:ntimeStep_UPPEST) :: q_ntime, cq_ntime
      type(dType), dimension(4) :: q, qx, qy, qxx, qxy, qyy
      type(dType), dimension(4) :: cqt, cqx, cqy
      type(dType), dimension(4) :: f, g, fx, gy
      type(dType), dimension(4) :: f_shock, g_shock, fx_shock, gy_shock
      type(dType), dimension(4) :: v1, v2, v3, v4
      type(dType), dimension(4,4) :: Taui, Ajac, Bjac, M
      type(dType), dimension(4,nfb_MAX) :: res
      real(DP),    dimension(4,4,nfb_MAX) :: dres_q, dres_qx, dres_qy, dres_qxx, dres_qxy, dres_qyy
      real(DP) :: h_length
      integer :: i, j, ii, jj, itime


      nq  = this % nequation
      nfb = element % nfb

      Include_time  = time_control % Include_time
      ntime_low     = time_control % ntimeStep_time_low
      ntime_up      = time_control % ntimeStep_time_up
      Use_dt_local  = time_control % Use_dt_local
      cfl           = time_control % cfl
      dt_global     = time_control % dt_global
      itime_space   = time_control % itimeStep_space
      ratio_space   = time_control % ratio_space

      res(1:nq,1:nfb) = 0.

      call element % Get_element_N ( N, r, s )
      call element % Calculate_element_Jel_area_Nx_Ny ( Jel_area, Nx, Ny, fbasis_nfb(:) % x, fbasis_nfb(:) % y, r, s )

      if( this % Do_2nd_derivative ) then
        call element % Calculate_element_Nxx_Nxy_Nyy ( Nxx, Nxy, Nyy, fbasis_nfb(:) % x, fbasis_nfb(:) % y, r, s )
      else
        Nxx(1:nfb) = 0.
        Nxy(1:nfb) = 0.
        Nyy(1:nfb) = 0.
      end if

! - get res' derivative w.r.t. to (q, qx, qy) using overload
      ! initialize {q}, qx := {dq/dx}, qy := {dq/dy}
      if( .not. this % Do_2nd_derivative ) then
        overloadLength = 3*nq
      else
        overloadLength = 6*nq
      end if

      do ii = 1, nq
        call initialize( q (ii),        ii, dot_product_my ( N , fbasis_nfb(:) % q(ii,itime_space), nfb ) )
        call initialize( qx(ii),   nq + ii, dot_product_my ( Nx, fbasis_nfb(:) % q(ii,itime_space), nfb ) )
        call initialize( qy(ii), 2*nq + ii, dot_product_my ( Ny, fbasis_nfb(:) % q(ii,itime_space), nfb ) )

        if( this % Do_2nd_derivative ) then
          call initialize( qxx(ii), 3*nq + ii, dot_product_my ( Nxx, fbasis_nfb(:) % q(ii,itime_space), nfb ) )
          call initialize( qxy(ii), 4*nq + ii, dot_product_my ( Nxy, fbasis_nfb(:) % q(ii,itime_space), nfb ) )
          call initialize( qyy(ii), 5*nq + ii, dot_product_my ( Nyy, fbasis_nfb(:) % q(ii,itime_space), nfb ) )
        else
          qxx(ii) = 0.
          qxy(ii) = 0.
          qyy(ii) = 0.
        end if
      end do


      ! cqt := {dcq/dt}
      if( Include_time ) then
        ! q_ntime
        q_ntime(:,itime_space) = q(:)

        do itime = ntime_low, ntime_up
          if( itime == itime_space )  cycle

          do ii = 1, nq
            q_ntime (ii,itime) = dot_product_my ( N(:) , fbasis_nfb(:) % q(ii,itime), nfb )
          end do
        end do

        ! cq_ntime
        do itime = ntime_low, ntime_up
          cq_ntime(:,itime) = conservative_variable_cq ( q_ntime(:,itime) )
        end do

        ! dt
        if ( Use_dt_local ) then
          dt = cfl * dot_product_my ( N, fbasis_nfb(:) % dt_local, nfb )
         !dt = cfl * 1./sum( sqrt( Nx(1:nfb)**2 + Ny(1:nfb)**2 ) ) / ( sqrt( q(4) ) + sqrt( q(2)**2 + q(3)**2 )  )/nfb
        else
          dt = dt_global
        end if

        ! cqt
        cqt(:) = time_control % time_discretization ( cq_ntime, dt, nq )
      else
        cqt(:) = 0.
      end if

      ! term 1: N*{dcq/dt}
      if ( Include_time ) then
        do i = 1, nfb
          res(:,i) = res(:,i) + N(i)*cqt(:)
        end do
      end if

1000  continue

      ! term 2: Nx*{F} + Ny*{G}
      f(:) = convective_flux_f ( q )
      g(:) = convective_flux_g ( q )


     ! extra diffusive flux for shock capturing
      if( this % Do_shock_capture ) then
        h_length = 1./sum( sqrt( Nx(1:nfb)**2 + Ny(1:nfb)**2 ) )

       !call shock_capture_flux_f_g          ( f_shock, g_shock, q, qx, qy, h_length, this % k_shock )
        call shock_capture_flux_f_g_modified ( f_shock, g_shock, q, qx, qy, h_length, this % k_shock, this % epsilon_shock )

        f(:) = f(:) - f_shock(:)
        g(:) = g(:) - g_shock(:)
      end if

      do i = 1, nfb
        res(:,i) = res(:,i) - ( ratio_space*Nx(i) )*f(:) - ( ratio_space*Ny(i) )*g(:)
      end do


2000  continue


      ! term 3: ( Nx*[A] + Ny*[B] )*[tau]*( {dcq/dt} + {dF/dx} + {dG/dy} ), where {dF/dx} = [A]*{dcq/dx}, {dG/dy} = [B]*{dcq/dy}
      v1(:) = 0.

      Taui(:,:) = get_tau_inverse ( q, Nx, Ny, nfb )

      Ajac(:,:) = convective_jacobian_A ( q )
      Bjac(:,:) = convective_jacobian_B ( q )

      M(:,:) = cq_derivative_wrt_q (q)


      ! dcq/dt
      if( Include_time ) then
        v1(:) = v1(:) + cqt(:)
      end if


      cqx(:) = matrix_vector_multiply ( M, qx, nq, nq )
      cqy(:) = matrix_vector_multiply ( M, qy, nq, nq )

      fx(:) = matrix_vector_multiply ( Ajac, cqx, nq, nq )
      gy(:) = matrix_vector_multiply ( Bjac, cqy, nq, nq )

      v2(:) = fx(:) + gy(:)

      ! extra flux for shock capture
      if( this % Do_shock_capture ) then
        h_length = 1./sum( sqrt( Nx(1:nfb)**2 + Ny(1:nfb)**2 ) )

        call shock_capture_flux_fx_gy_modified ( fx_shock, gy_shock, q, qx, qy, qxx, qxy, qyy, h_length, this % k_shock, this % epsilon_shock )

        v2(:) = v2(:) - fx_shock(:) - gy_shock(:)
      end if

      v1(:) = v1(:) + ratio_space * v2(:)


      if (Do_print) then
        write(*,100) 'v1: cqt + fx + gy', v1(1:nq)
      end if


      ! inv([tau])*{v1}
      v2(:) = GaussJordan ( Taui, v1, nq )  ! Taui, v1 will be modified by GaussJordan


      if (Do_print) then
        write(*,100) 'v2: taui*{cqt + fx + gy}', v2(1:nq)
      end if


      v3(:) = matrix_vector_multiply ( Ajac, v2, nq, nq )
      v4(:) = matrix_vector_multiply ( Bjac, v2, nq, nq )


      if (Do_print) then
        write(*,100) 'v3: A*tau*{cqt + fx + gy}', v3(1:nq)
        write(*,100) 'v4: B*tau*{cqt + fx + gy}', v4(1:nq)
      end if


      do i = 1, nfb
        res(:,i) = res(:,i) + Nx(i)*v3(:) + Ny(i)*v4(:)
      end do

3000  continue

      ! multiply Jacobian
      res(:,1:nfb) = Jel_area*res(:,1:nfb)


! - take answer out
      do i = 1, nfb
        do ii = 1, nq
        do jj = 1, nq
          dres_q  (ii,jj,i) = res(ii,i) % deriv(     jj)
          dres_qx (ii,jj,i) = res(ii,i) % deriv(  nq+jj)
          dres_qy (ii,jj,i) = res(ii,i) % deriv(2*nq+jj)

          if( this % Do_2nd_derivative ) then
            dres_qxx(ii,jj,i) = res(ii,i)%deriv(3*nq+jj)
            dres_qxy(ii,jj,i) = res(ii,i)%deriv(4*nq+jj)
            dres_qyy(ii,jj,i) = res(ii,i)%deriv(5*nq+jj)
          else
            dres_qxx(ii,jj,i) = 0.
            dres_qxy(ii,jj,i) = 0.
            dres_qyy(ii,jj,i) = 0.
          end if
        end do
        end do
      end do


! - get res' derivative w.r.t qfb using chain rule
      ! (q, qx, qy)'s derivative w.r.t. qfb is (N, Nx, Ny)
      do i = 1, nfb
      do j = 1, nfb
        aa(:,:,i,j) = dres_q(:,:,i)*N(j) + dres_qx(:,:,i)*Nx(j) + dres_qy(:,:,i)*Ny(j)

        if( this % Do_2nd_derivative ) then
          aa(:,:,i,j) = aa(:,:,i,j) + dres_qxx(:,:,i)*Nxx(j) + dres_qxy(:,:,i)*Nxy(j) + dres_qyy(:,:,i)*Nyy(j)
        end if
      end do
      end do


100   format(a,20(es16.5))

end subroutine Calculate_euler_interior_term_linearization
!=====================================================================
!
!=====================================================================
subroutine Calculate_euler_natural_boundary_term ( this, res, bedge, fbasis_nfbnz,  element, fbasis_nfb, time_control, rb )
      implicit none
      class(euler_equation_t),                 intent(in)  :: this
      class(bedge_t),                          intent(in)  :: bedge
      type(fbasis_t),          dimension(:),   intent(in)  :: fbasis_nfbnz
      class(element_t),                        intent(in)  :: element
      type(fbasis_t),          dimension(:),   intent(in)  :: fbasis_nfb
      type(time_control_t),                    intent(in)  :: time_control
      real(DP),                dimension(:,:), intent(out) :: res
      real(DP),                                intent(in)  :: rb

      integer :: nq, nfb, nfbnz
      integer :: itime_space
      real(DP) :: ratio_space
      real(DP), dimension(nfbnz_MAX) :: N_nz
      real(DP), dimension(nfb_MAX  ) :: N, Nx, Ny
      real(DP) :: xnorm, ynorm, Jsf_length
      real(DP) :: r, s
      real(DP), dimension(4) :: q, qb, rr
      real(DP) :: spre
      integer :: i, ii


      nq = this % nequation
      nfb = element % nfb
      nfbnz = bedge % nfbnz

      ! time control
      itime_space = time_control % itimeStep_space
      ratio_space = time_control % ratio_space


      call bedge % Calculate_bedge_xynorm_Jsf_length ( xnorm, ynorm, Jsf_length, fbasis_nfbnz(:) % x, fbasis_nfbnz(:) % y, rb )
      call bedge % Get_bedge_N_nz ( N_nz, rb )

      select case( bedge % bcType )
      ! inviscid wall
      case(1000)
        do ii = 1, nq
          q(ii) = dot_product_my( N_nz(:), fbasis_nfbnz(:) % q(ii,itime_space), nfbnz )
        end do

        spre = gi*q(1)*q(4)

        rr(1) = 0.
        rr(2) = spre*xnorm
        rr(3) = spre*ynorm
        rr(4) = 0.

        rr(:) = (ratio_space*Jsf_length)*rr(:)

        do i = 1, nfbnz
          res(:,i) = rr(:)*N_nz(i)
        end do

       !write(*,*) '1000-bc'
       !write(*,100) 'rb', rb
       !write(*,100) 'x', fbasis_nfbnz(1:nfbnz) % x
       !write(*,100) 'y', fbasis_nfbnz(1:nfbnz) % y
       !write(*,100) 'xnorm, ynorm', xnorm, ynorm
       !write(*,100) 'N_nz', N_nz(1:nfbnz)
       !write(*,100) 'Jst_length', Jsf_length
       !write(*,100) 'spre', spre
       !write(*,100) 'rr', rr(1:4)

      case(2000)
        do ii = 1, nq
          q(ii) = dot_product_my( N_nz(:), fbasis_nfbnz(:) % q(ii,itime_space), nfbnz )
        end do

        qb(1:4) = this % q_farfield(1:4)

        rr(:) = roe_flux_face( q, qb, xnorm, ynorm )

        rr(:) = (ratio_space*Jsf_length)*rr(:)

        do i = 1, nfbnz
          res(:,i) = rr(:)*N_nz(i)
        end do
      case default
        write(*,*) 'wrong! Calculate_euler_natural_boundary_term:', bedge % bcType
        stop
      end select


100   format(a,20(es16.5))

end subroutine Calculate_euler_natural_boundary_term
!=====================================================================
!
!=====================================================================
subroutine Calculate_euler_natural_boundary_term_linearization ( this, aa, bedge, fbasis_nfbnz,  element, fbasis_nfb, time_control, rb )
      implicit none
      class(euler_equation_t),                        intent(in)  :: this
      class(bedge_t),                                 intent(in)  :: bedge
      type(fbasis_t),       dimension(:),             intent(in)  :: fbasis_nfbnz
      class(element_t),                               intent(in)  :: element
      type(fbasis_t),       dimension(:),             intent(in)  :: fbasis_nfb
      type(time_control_t),                           intent(in)  :: time_control
      real(DP),             dimension(:,:,:,:),       intent(out) :: aa
      real(DP),                                       intent(in)  :: rb

      integer :: nq, nfb, nfbnz
      integer :: itime_space
      real(DP) :: ratio_space
      real(DP), dimension(nfbnz_MAX) :: N_nz
      real(DP), dimension(nfb_MAX)   :: N, Nx, Ny
      real(DP) :: xnorm, ynorm, Jsf_length
      real(DP) :: r, s
      type(dType), dimension(4)   :: q, qb, rr
      type(dType) :: spre
      real(DP),    dimension(4,4) :: drr_q, drr_qx, drr_qy
      integer :: i, j, ii, jj

      nq = this % nequation
      nfb = element % nfb
      nfbnz = bedge % nfbnz

      ! time control
      itime_space = time_control % itimeStep_space
      ratio_space = time_control % ratio_space


      call bedge % Calculate_bedge_xynorm_Jsf_length ( xnorm, ynorm, Jsf_length, fbasis_nfbnz(:) % x, fbasis_nfbnz(:) % y, rb )
      call bedge % Get_bedge_N_nz ( N_nz, rb )

      call bedge % Get_bedge_to_element_rs ( r, s, rb )
      call element % Get_element_N ( N, r, s )

      select case( bedge % bcType )
      ! inviscid wall
      case(1000)
        ! {rr}'s derivative w.r.t. {q} on quadrature point, using overload
        overloadLength = nq

        do ii = 1, nq
          call initialize( q(ii), ii, dot_product_my ( N_nz(:), fbasis_nfbnz(:) % q(ii,itime_space), nfbnz ) )
        end do

        spre = gi*q(1)*q(4)

        rr(1) = 0.
        rr(2) = spre*xnorm
        rr(3) = spre*ynorm
        rr(4) = 0.

        rr(:) = (ratio_space*Jsf_length)*rr(:)

        ! take out {rr}'s derivative w.r.t. {q}
        do ii = 1, nq
        do jj = 1, nq
          drr_q(ii,jj) = rr(ii)%deriv(jj)
        end do
        end do

        ! use chain rule to get {res}'s derivative w.r.t. qfb
        do i = 1, nfbnz
        do j = 1, nfb
          aa(:,:,i,j) = N_nz(i)*N(j)*drr_q(:,:)
        end do
        end do

      ! far-field
      case(2000)
        ! {rr}'s derivative w.r.t. {q} on quadrature point, using overload
        overloadLength = nq

        do ii = 1, nq
          call initialize( q(ii), ii, dot_product_my ( N_nz(:), fbasis_nfbnz(:) % q(ii,itime_space), nfbnz ) )
        end do

        qb(1:4) = this % q_farfield(1:4)

        rr(:) = roe_flux_face( q, qb, xnorm, ynorm )

        rr(:) = (ratio_space*Jsf_length)*rr(:)

        ! take out {rr}'s derivative w.r.t. {q}
        do ii = 1, nq
        do jj = 1, nq
          drr_q(ii,jj) = rr(ii)%deriv(jj)
        end do
        end do

        ! use chain rule to get {res}'s derivative w.r.t. qfb
        do i = 1, nfbnz
        do j = 1, nfb
          aa(:,:,i,j) = N_nz(i)*N(j)*drr_q(:,:)
        end do
        end do
      case default
        write(*,*) 'wrong! Calculate_euler_natural_boundary_term_linearization:', bedge % bcType
        stop
      end select

100   format(a,20(es16.5))
200   format(a,i7,i7,20(es16.5))

end subroutine Calculate_euler_natural_boundary_term_linearization
!=====================================================================
!
!=====================================================================
function euler_max_omega_for_dq ( this, FEM, dq, itimeStep )  result( omega_max )
      implicit none
      class(euler_equation_t), intent(in)  :: this
      class(FEM_t),            intent(in)  :: FEM
      real(DP),                intent(in)  :: dq(FEM%nequation,FEM%nfbasis)
      integer,                 intent(in)  :: itimeStep
      real(DP)                             :: omega_max

      real(DP), parameter :: coef = 0.5
      real(DP) :: omega_r, omega_t
      integer :: i


      omega_r = huge(omega_r)
      omega_t = huge(omega_t)

      omega_r = 1.5
      omega_t = 1.5

      do i = 1, FEM % nfbasis
        if( dq(1,i) < 0. ) then
          omega_r = min( omega_r, abs( coef/( dq(1,i)/FEM % fbasis(i) % q(1,itimeStep) ) ) )
        end if

        if( dq(4,i) < 0. ) then
          omega_t = min( omega_t, abs( coef/( dq(4,i)/FEM % fbasis(i) % q(4,itimeStep) ) ) )
        end if
      end do

      omega_max = min( omega_r, omega_t )

end function euler_max_omega_for_dq
!=====================================================================
!
!=====================================================================
subroutine Check_euler_solution ( this, FEM, itimeStep, fp, filename )
      implicit none
      class(euler_equation_t), intent(in) :: this
      class(FEM_t),            intent(in) :: FEM
      integer,                 intent(in) :: itimeStep
      integer,                 intent(in) :: fp
      character(len=*),        intent(in) :: filename

      integer :: nq
      integer :: i
      logical :: Has_invalid

      nq = this % nequation

      Has_invalid = .false.

      open(fp, file=trim(filename))
      write(fp,*) 'variables="x","y","q1","q2","q3","q4"'

      do i = 1, FEM % nfbasis
        if( FEM % fbasis(i) % q(1,itimeStep) < 0. ) then
          write(*,200) 'q(1) < 0: itimeStep, x, y, q', itimeStep, FEM % fbasis(i) % x, FEM % fbasis(i) % y, FEM % fbasis(i) % q(1:nq,itimeStep)
          write(fp,100)  FEM % fbasis(i) % x, FEM % fbasis(i) % y, FEM % fbasis(i) % q(1:nq,itimeStep)

          Has_invalid = .true.
        end if

        if( FEM % fbasis(i) % q(4,itimeStep) < 0. ) then
          write(*,200) 'q(4) < 0: itimeStep, x, y, q', itimeStep, FEM % fbasis(i) % x, FEM % fbasis(i) % y, FEM % fbasis(i) % q(1:nq,itimeStep)
          write(fp,100)  FEM % fbasis(i) % x, FEM % fbasis(i) % y, FEM % fbasis(i) % q(1:nq,itimeStep)

          Has_invalid = .true.
        end if
      end do

      close(fp)

      if( Has_invalid )  stop

100   format(20(es16.5))
200   format(a,i7,20(es16.5))

end subroutine Check_euler_solution
!=============================================================================
!
!=============================================================================
function manufactured_solution ( x, y )  result(q)
      implicit none
      real(DP), intent(in)  :: x, y
      real(DP)              :: q(4)

      real(DP), parameter :: pi = 2.*acos(0.)
      real(DP) :: rho0, u0, v0, t0, rho, rhou, rhov, rhoT

      rho0 = 1.1
      u0   = 0.6
      v0   = 0.1
      t0   = 0.8

      rho  = rho0    * (1.d0 + sin(pi*x)*cos(pi*x)*sin(pi*y)*cos(pi*y))
      rhou = rho0*u0 * (1.d0 + sin(1.5*pi*x)*cos(1.5*pi*x)*sin(1.5*pi*y)*cos(1.5*pi*y))
      rhov = rho0*v0 * (1.d0 + sin(1.5*pi*x)*cos(1.5*pi*x)*sin(1.5*pi*y)*cos(1.5*pi*y))
      rhoT = rho0*T0 * (1.d0 + sin(pi*x)**2*sin(pi*y)**2)

      q(1) = rho
      q(2) = rhou/rho
      q(3) = rhov/rho
      q(4) = rhoT/rho

end function manufactured_solution

end module FEM_Equation_Euler_mod
