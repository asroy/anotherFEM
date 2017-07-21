module FEM_ElementType_Lagrange_mod
      use KindDefinition_mod, only : DP
      use FEM_Generic_mod
      use FEMEntry_Generic_mod
      use FEMEntry_ElementType_Lagrange_mod
      use Mesh_Generic_mod
      implicit none
      save

      type, extends (FEM_t) :: lagrange_t
      contains
        procedure, pass(this) :: Generate_FEM_connectivity => Generate_lagrange_FEM_connectivity
        procedure, pass(this) :: Generate_FEM_connectivity_high_order_mesh => Generate_lagrange_FEM_connectivity_high_order_mesh
        procedure, pass(this) :: Output_FEM_solution_tecplot => Output_lagrange_FEM_solution_tecplot
      end type lagrange_t

contains

!=====================================================================
!
!=====================================================================
subroutine Generate_lagrange_FEM_connectivity ( this, mesh )
      implicit none
      class(lagrange_t), intent(inout) :: this
      class(mesh_t),     intent(in)    :: mesh

      integer :: i, ifb, iel, inode1, inode2, inode3, icell1, icell2, iside1, iside2, iside

      integer :: ibedge
      real(DP) :: x, y

      write(fpd1,*) 'Generate_lagrange_FEM_connectivity: Start'

! - generate connectivty
      write(fpd1,*) 'mesh%ntriangle, mesh%nbedge, mesh%nedge, mesh%nnode'
      write(fpd1,*)  mesh%ntriangle, mesh%nbedge, mesh%nedge, mesh%nnode

      this%nelement = mesh%ntriangle
      this%nbedge   = mesh%nbedge

      select case (this%basis_order_START)
      case(1)
        this%nfbasis = mesh%nnode
        this%nnode = mesh%nnode
      case(2)
        this%nfbasis = mesh%nnode + mesh%nedge
        this%nnode = mesh%nnode + mesh%nedge
      case(3)
        this%nfbasis = mesh%nnode + 2*mesh%nedge + mesh%ntriangle
        this%nnode = mesh%nnode + 2*mesh%nedge + mesh%ntriangle
      case default
        write(*,*) 'wrong! Generate_lagrange_FEM_connectivity: basis_order_START', this%basis_order_START
        stop
      end select


      write(fpd1,*) 'this%nelement, this%nbedge, this%nnode'
      write(fpd1,*)  this%nelement, this%nbedge, this%nnode

      allocate ( lagrange_element_t :: this%element(this%nelement) )
      allocate ( lagrange_bedge_t :: this%bedge(this%nbedge) )
      allocate ( this%fbasis( this%nfbasis ) )
      allocate ( this%node(this%nnode) )

      do i = 1, this % nelement
        call this % element(i) % Initialize_element ( this % basis_order_START )
      end do

      do i = 1, this % nbedge
        call this % bedge(i) % Initialize_bedge ( this % basis_order_START )

        this % bedge(i) % iel         = mesh % bedge(i) % cell
        this % bedge(i) % elementSide = mesh % bedge(i) % cellSide
        this % bedge(i) % bcType      = mesh % bedge(i) % bcType
      end do

      ! vortex basis functions
      do i = 1, mesh%nnode
        this%fbasis(i)%x = mesh%node(i)%x
        this%fbasis(i)%y = mesh%node(i)%y
      end do

      do i = 1, mesh%ntriangle
        this%element(i)%ifb(1:3) = mesh%triangle(i)%node(1:3)
      end do

      ! edge and bubble basis functions
      select case ( this%basis_order_START )
      case(1)
        continue
      case(2)
        ifb = mesh%nnode

        ! edge basis functions
        do i = 1, mesh%nedge
          inode1 = mesh%edge(i)%node1
          inode2 = mesh%edge(i)%node2

          this%fbasis(ifb+1)%x = ( mesh%node(inode1)%x + mesh%node(inode2)%x )/2.
          this%fbasis(ifb+1)%y = ( mesh%node(inode1)%y + mesh%node(inode2)%y )/2.

          icell1 = mesh%edge(i)%cell1
          iside1 = mesh%edge(i)%cellSide1

          this%element(icell1)%ifb(3+iside1) = ifb + 1

          icell2 = mesh%edge(i)%cell2
          iside2 = mesh%edge(i)%cellSide2

          if( icell2 > 0 ) then
            this%element(icell2)%ifb(3+iside2) = ifb + 1
          end if

          ifb = ifb + 1
        end do
      case(3)
        ifb = mesh%nnode

        ! edge basis functions
        do i = 1, mesh%nedge
          inode1 = mesh%edge(i)%node1
          inode2 = mesh%edge(i)%node2

          this%fbasis(ifb+1)%x = ( 2.*mesh%node(inode1)%x +    mesh%node(inode2)%x )/3.
          this%fbasis(ifb+1)%y = ( 2.*mesh%node(inode1)%y +    mesh%node(inode2)%y )/3.

          this%fbasis(ifb+2)%x = (    mesh%node(inode1)%x + 2.*mesh%node(inode2)%x )/3.
          this%fbasis(ifb+2)%y = (    mesh%node(inode1)%y + 2.*mesh%node(inode2)%y )/3.

          icell1 = mesh%edge(i)%cell1
          iside1 = mesh%edge(i)%cellSide1

          this%element(icell1)%ifb(3+iside1) = ifb + 1
          this%element(icell1)%ifb(6+iside1) = ifb + 2

          icell2 = mesh%edge(i)%cell2
          iside2 = mesh%edge(i)%cellSide2

          if ( icell2 > 0 ) then
            this%element(icell2)%ifb(3+iside2) = ifb + 2
            this%element(icell2)%ifb(6+iside2) = ifb + 1
          end if

          ifb = ifb + 2
        end do

        ! bubble basis functions
        do i = 1, mesh%ntriangle
          inode1 = mesh%triangle(i)%node(1)
          inode2 = mesh%triangle(i)%node(2)
          inode3 = mesh%triangle(i)%node(3)

          this%fbasis(ifb+1)%x = ( mesh%node(inode1)%x + mesh%node(inode2)%x + mesh%node(inode3)%x )/3.
          this%fbasis(ifb+1)%y = ( mesh%node(inode1)%y + mesh%node(inode2)%y + mesh%node(inode3)%y )/3.

          this%element(i)%ifb(10) = ifb + 1

          ifb = ifb + 1
        end do
      end select

      ! bedge basis functions
      do i = 1, this%nbedge
        iel   = this%bedge(i)%iel
        iside = this%bedge(i)%elementSide

        select case(iside)
        case(1)
          this%bedge(i)%ifbnz(1) = this%element(iel)%ifb(1)
          this%bedge(i)%ifbnz(2) = this%element(iel)%ifb(2)
        case(2)
          this%bedge(i)%ifbnz(1) = this%element(iel)%ifb(2)
          this%bedge(i)%ifbnz(2) = this%element(iel)%ifb(3)
        case(3)
          this%bedge(i)%ifbnz(1) = this%element(iel)%ifb(3)
          this%bedge(i)%ifbnz(2) = this%element(iel)%ifb(1)
        case default
          write(*,*) 'wrong! Generate_lagrange_FEM_connectivity'
          stop
        end select

        if ( this%element(iel)%basis_order == 2 ) then
          this%bedge(i)%ifbnz(3) = this%element(iel)%ifb(3+iside)
        else if ( this%element(iel)%basis_order == 3 ) then
          this%bedge(i)%ifbnz(3) = this%element(iel)%ifb(3+iside)
          this%bedge(i)%ifbnz(4) = this%element(iel)%ifb(6+iside)
        end if
      end do

      ! >>> debug
      write(*,*) 'snap naca0012'

      do ibedge = 1, this % nbedge
        if( this % bedge(ibedge)%bcType /= 1000 )  cycle

        do i = 1, this % bedge(ibedge)%nfbnz
          ifb = this % bedge(ibedge)%ifbnz(i)

          x = this % fbasis(ifb)%x
          y = this % fbasis(ifb)%y

          if( abs(y) < 1.e-6 )  cycle
          if( x > 0.437 ) cycle

          x = x + 0.5
          if( y > 0. ) then
            y =   0.12/0.2*( 0.2969*sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
          else
            y = - 0.12/0.2*( 0.2969*sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
          end if

          this % fbasis(ifb) % y = y
        end do
      end do
      ! <<< debug


      write(*,*) 'Generate_lagrange_FEM_connectivity: End'

end subroutine Generate_lagrange_FEM_connectivity
!=====================================================================
!
!=====================================================================
subroutine Generate_lagrange_FEM_connectivity_high_order_mesh ( this, mesh )
      implicit none
      class(lagrange_t), intent(inout) :: this
      class(mesh_t),     intent(in)    :: mesh

      integer :: i, ifb, iel, inode1, inode2, inode3, icell1, icell2, iside1, iside2, iside
      integer :: nfb

      write(fpd1,*) 'Generate_lagrange_FEM_connectivity: Start'

      if( this % basis_order_START /= mesh % mesh_basis_order ) then
        write(*,*) 'wrong! Generate_lagrange_FEM_connectivity: mesh_basis_order'
        stop
      end if

! - generate connectivty
      write(fpd1,*) 'mesh%ntriangle, mesh%nbedge, mesh%nedge, mesh%nnode'
      write(fpd1,*)  mesh%ntriangle, mesh%nbedge, mesh%nedge, mesh%nnode

      this%nelement = mesh%ntriangle
      this%nbedge   = mesh%nbedge

      this%nfbasis = mesh%nnode
      this%nnode = mesh%nnode

      write(fpd1,*) 'this%nelement, this%nbedge, this%nnode'
      write(fpd1,*)  this%nelement, this%nbedge, this%nnode

      allocate ( lagrange_element_t :: this%element(this%nelement) )
      allocate ( lagrange_bedge_t :: this%bedge(this%nbedge) )
      allocate ( this%fbasis( this%nfbasis ) )
      allocate ( this%node(this%nnode) )

      do i = 1, this % nelement
        call this % element(i) % Initialize_element ( this % basis_order_START )
      end do

      do i = 1, this % nbedge
        call this % bedge(i) % Initialize_bedge ( this % basis_order_START )

        this % bedge(i) % iel         = mesh % bedge(i) % cell
        this % bedge(i) % elementSide = mesh % bedge(i) % cellSide
        this % bedge(i) % bcType      = mesh % bedge(i) % bcType
      end do

      ! vortex basis functions
      do i = 1, mesh%nnode
        this%fbasis(i)%x = mesh%node(i)%x
        this%fbasis(i)%y = mesh%node(i)%y
      end do

      select case( this % basis_order_start)
      case(1)
        nfb = 3
      case(2)
        nfb = 6
      case(3)
        nfb = 10
      case default
        write(*,*) 'wrong! Ge'
        stop
      end select

      do i = 1, mesh%ntriangle
        this%element(i)%ifb(1:nfb) = mesh%triangle(i)%node(1:nfb)
      end do


      ! bedge basis functions
      do i = 1, this%nbedge
        iel   = this%bedge(i)%iel
        iside = this%bedge(i)%elementSide

        select case(iside)
        case(1)
          this%bedge(i)%ifbnz(1) = this%element(iel)%ifb(1)
          this%bedge(i)%ifbnz(2) = this%element(iel)%ifb(2)
        case(2)
          this%bedge(i)%ifbnz(1) = this%element(iel)%ifb(2)
          this%bedge(i)%ifbnz(2) = this%element(iel)%ifb(3)
        case(3)
          this%bedge(i)%ifbnz(1) = this%element(iel)%ifb(3)
          this%bedge(i)%ifbnz(2) = this%element(iel)%ifb(1)
        case default
          write(*,*) 'wrong! Generate_lagrange_FEM_connectivity'
          stop
        end select

        if ( this%element(iel)%basis_order == 2 ) then
          this%bedge(i)%ifbnz(3) = this%element(iel)%ifb(3+iside)
        else if ( this%element(iel)%basis_order == 3 ) then
          this%bedge(i)%ifbnz(3) = this%element(iel)%ifb(3+iside)
          this%bedge(i)%ifbnz(4) = this%element(iel)%ifb(6+iside)
        end if
      end do

      write(*,*) 'Generate_lagrange_FEM_connectivity: End'

end subroutine Generate_lagrange_FEM_connectivity_high_order_mesh
!=====================================================================
!
!=====================================================================
subroutine Output_lagrange_FEM_solution_tecplot( this, fp, filename )
      implicit none
      class(lagrange_t), intent(in) :: this
      integer, intent(in) :: fp
      character(len=*), intent(in) :: filename

      integer :: nq
      integer :: ntri, ntriangle
      integer, dimension(:,:), pointer :: node
      integer :: i, j, k


      nq = this % nequation

      ntriangle = 0
      do i = 1, this%nelement
        call px_triangle_to_p1_triangles ( ntri, node, this%element(i)%basis_order  )
        ntriangle = ntriangle + ntri
      end do

      open(fp, file=trim(filename))
      write(fp,200) 'variables="x","y","q1","q2","q3","q4","e","ev1","ev2","ev3","ev4", "em11", "em12", "em13", "em14", "em21", "em22", "em23", "em24", "em31", "em32", "em33", "em34", "em41", "em42", "em43", "em44"'
     !write(fp,200) 'variables="x","y","q1","q2","q3","q4","e","ev1","ev2","ev3","ev4", "em11", "em22", "em33", "em44"'
     !write(fp,200) 'variables="x","y","q1","q2","q3","q4","e"'
     !write(fp,200) 'variables="x","y","q1","e"'
      write(fp,*) 'zone i=',this%nfbasis,', e=',ntriangle,',DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
      do i = 1, this%nfbasis
        write(fp,100) this%fbasis(i)%x, this%fbasis(i)%y, this%fbasis(i)%q(1:nq,1), this%fbasis(i)%e, this%fbasis(i)%evector(1:nq), ((this%fbasis(i)%ematrix(j,k),k=1,nq),j=1,nq)
       !write(fp,100) this%fbasis(i)%x, this%fbasis(i)%y, this%fbasis(i)%q(1:nq,1), this%fbasis(i)%e, this%fbasis(i)%evector(1:nq),  (this%fbasis(i)%ematrix(j,j),j=1,nq)
       !write(fp,100) this%fbasis(i)%x, this%fbasis(i)%y, this%fbasis(i)%q(1:nq,1), this%fbasis(i)%e
      end do
      do i = 1, this%nelement
        call px_triangle_to_p1_triangles ( ntri, node, this%element(i)%basis_order  )

        do j = 1, ntri
          write(fp,*) this%element(i)%ifb(node(1:3,j))
        end do
      end do
      close(fp)

100   format(50(es23.15))
200   format(a200)

      contains

      !=====================================================================
      !
      ! map high order element to p1 pieces
      !
      !=====================================================================
      subroutine px_triangle_to_p1_triangles ( ntriangle, node, basis_order )
            implicit none
            integer,                          intent(out) :: ntriangle
            integer, dimension(:,:), pointer, intent(out) :: node
            integer,                          intent(in)  :: basis_order

            integer :: ntriangle_p1 = 1
            integer :: ntriangle_p2 = 4
            integer :: ntriangle_p3 = 9
            integer, dimension(3,1), target :: node_p1 = reshape ( (/ 1, 2, 3 /), (/ 3, 1 /) )
            integer, dimension(3,4), target :: node_p2 = reshape ( (/ 1, 4, 6, &
                                                                      2, 5, 4, &
                                                                      3, 6, 5, &
                                                                      4, 5, 6 /), (/ 3, 4 /) )
            integer, dimension(3,9), target :: node_p3 = reshape ( (/ 1, 4, 9, &
                                                                      2, 5, 7, &
                                                                      3, 6, 8, &
                                                                      4, 7,10, &
                                                                      5, 8,10, &
                                                                      6, 9,10, &
                                                                      4,10, 9, &
                                                                      5,10, 7, &
                                                                      6,10, 8 /), (/ 3, 9 /) )

            select case( basis_order )
            case (1)
              ntriangle = ntriangle_p1
              node => node_p1
            case (2)
              ntriangle = ntriangle_p2
              node => node_p2
            case (3)
              ntriangle = ntriangle_p3
              node => node_p3
            case default
              write(*,*) 'wrong!: px_triangle_to_p1_triangles'
              stop
            end select

      end subroutine px_triangle_to_p1_triangles

end subroutine Output_lagrange_FEM_solution_tecplot

end module FEM_ElementType_Lagrange_mod
