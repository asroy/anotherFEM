module Mesh_Generic_mod
      use KindDefinition_mod, only : DP
      use Constant_mod
      implicit none
      save

      type :: mesh_triangle_t
        integer :: node(nnd_MAX)
      end type mesh_triangle_T

      type :: mesh_edge_t
        integer :: node1, node2
        integer :: cell1, cell2
        integer :: cellSide1, cellSide2
      end type mesh_edge_t

      type :: mesh_bedge_t
        integer :: cell
        integer :: cellSide
        integer :: bcType
      end type mesh_bedge_t

      type :: mesh_node_t
        real(DP) :: x, y
      end type mesh_node_T

      ! mesh
      type :: mesh_t
        integer :: mesh_basis_order

        integer :: nzone

        integer :: ntriangle
        integer :: nedge
        integer :: nbedge
        integer :: nnode

        integer, allocatable :: ntrianglez(:)
        integer, allocatable :: nedgez(:)
        integer, allocatable :: nbedgez(:)
        integer, allocatable :: nnodez(:)

        type (mesh_triangle_t), allocatable :: triangle(:)
        type (mesh_edge_t), allocatable :: edge(:)
        type (mesh_bedge_t), allocatable :: bedge(:)
        type (mesh_node_t), allocatable :: node(:)
      contains
        procedure, pass(this) :: Output_mesh_tecplot
      end type mesh_t

      ! generic mesh importer
      type :: mesh_importer_t
      contains
        procedure, nopass :: Import_mesh
      end type mesh_importer_t

contains

!=====================================================================
!
!=====================================================================
subroutine Output_mesh_tecplot ( this, fp, filename )
      implicit none
      class(mesh_t), intent(in) :: this
      integer, intent(in) :: fp
      character(len=*), intent(in) :: filename

      integer, dimension(:), allocatable :: node_bcType
      integer :: i, ibedge, icell, l
      integer :: local_node(4)

      allocate( node_bcType( this%nnode ) )

      node_bcType(:) = 0

      do ibedge = 1, this % nbedge
        icell = this % bedge(ibedge) % cell

        local_node(1:3) = this % triangle(icell) % node(1:3)
        local_node(4) = local_node(1)

        l = this % bedge(ibedge) % cellSide

        node_bcType(local_node(l  )) = this % bedge(ibedge) % bcType
        node_bcType(local_node(l+1)) = this % bedge(ibedge) % bcType
      end do

      open(fp, file=trim(filename))
      write(fp,*) 'variables="x","y", "bcType"'
      write(fp,*) 'zone i=',this%nnode,', e=',this%ntriangle,',DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
      do i = 1, this%nnode
        write(fp,*) this%node(i)%x, this%node(i)%y, node_bcType(i)
      end do
      do i = 1, this%ntriangle
        write(fp,*) this%triangle(i)%node(1:3)
      end do
      close(fp)

      deallocate( node_bcType )

end subroutine Output_mesh_tecplot
!=====================================================================
!
!=====================================================================
subroutine Import_mesh ( mesh, mesh_basis_order, fpmesh, meshfile )
      implicit none
      type(mesh_t),     intent(out) :: mesh
      integer,          intent(in)  :: mesh_basis_order
      integer,          intent(in)  :: fpmesh
      character(len=*), intent(in)  :: meshfile
end subroutine Import_mesh

end module Mesh_Generic_mod
