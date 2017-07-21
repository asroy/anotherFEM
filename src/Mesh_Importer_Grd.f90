module Mesh_Importer_Grd_mod
      use KindDefinition_mod, only : DP
      use Mesh_Generic_mod
      implicit none
      save

      ! grd mesh importer
      type, extends (mesh_importer_t) :: grd_mesh_importer_t
      contains
        procedure, nopass :: Import_mesh => Import_mesh_from_grd
      end type grd_mesh_importer_t

contains

!=====================================================================
!
!=====================================================================
subroutine Import_mesh_from_grd ( mesh, mesh_basis_order, fpmesh, meshfile )
      implicit none
      type(mesh_t),     intent(out) :: mesh
      integer,          intent(in)  :: mesh_basis_order
      integer,          intent(in)  :: fpmesh
      character(len=*), intent(in)  :: meshfile

      type :: psp_t
        integer :: npsp
        integer, allocatable, dimension(:) :: edge, dir, node2, cell1, cell2, cellSide1, cellSide2
      endtype psp_t
      type(psp_t), allocatable :: psp(:)
      integer, parameter :: mnbr = 20
      integer :: ndump
      integer :: nnd_tri
      integer :: i, j, t, izone, l, nn(4), counter, iedge, n1, n2, nl, nr

      real(DP) :: tx, ty, c, s

      mesh % mesh_basis_order = mesh_basis_order


      select case(mesh_basis_order)
      case(1)
        nnd_tri = 3
      case(2)
        nnd_tri = 6
      case(3)
        nnd_tri = 10
      case default
        write(*,*) 'wrong! Mesh_Importer_Grd, mesh_basis_order'
        stop
      end select


      ! read *.grd meshfile
      open(fpmesh, file = meshfile, status = 'old')
      read(fpmesh,*)
      read(fpmesh,*)
      read(fpmesh,*)
      read(fpmesh,*)
      read(fpmesh,*)  ndump, mesh%nzone, mesh%nnode, ndump

      allocate( mesh%ntrianglez(mesh%nzone), mesh%nbedgez(mesh%nzone) )

      read(fpmesh,*)
      do i = 1, mesh%nzone
        read(fpmesh,*)  ndump, mesh%ntrianglez(i), mesh%nbedgez(i)
      end do

      mesh%ntriangle = 0
      mesh%nbedge = 0
      do izone = 1, mesh%nzone
        mesh%ntriangle = mesh%ntriangle + mesh%ntrianglez(izone)
        mesh%nbedge = mesh%nbedge + mesh%nbedgez(izone)
      end do

      allocate( mesh%triangle(mesh%ntriangle) )
      allocate( mesh%bedge(mesh%nbedge) )
      allocate( mesh%node(mesh%nnode) )

      read(fpmesh,*)
      do i = 1, mesh%ntriangle
        read(fpmesh,*)  ndump, ndump, mesh%triangle(i)%node(1:nnd_tri)
      end do

      read(fpmesh,*)
      do i = 1, mesh%nnode
        read(fpmesh,*) mesh%node(i)%x, mesh%node(i)%y
      end do
      

      ! >>> debug
     !write(*,*) 'debug! Mesh_Importer_Grd'
     !mesh % node(:) % x = mesh % node(:) % x/20.
     !mesh % node(:) % y = mesh % node(:) % y/20.

     !c = cos(0.8)
     !s = sin(0.8)
     !do i = 1, mesh % nnode
     !  tx = mesh % node(i) % x
     !  ty = mesh % node(i) % y

     !  mesh % node(i) % x = c*tx - s*ty
     !  mesh % node(i) % y = s*tx + c*ty
     !end do
      ! <<< debug


      read(fpmesh,*)
      do i = 1, mesh%nbedge
        read(fpmesh,*)  mesh%bedge(i)%bcType, mesh%bedge(i)%cell, mesh%bedge(i)%cellSide

        ! >>> debug
       !write(*,*) 'debug! Mesh_Importer_Grd'
       !mesh%bedge(i)%bcType = 0
        if ( mesh % bedge(i) % bcType == 3000 )  mesh % bedge(i) % bcType = 0
       !if ( mesh % bedge(i) % bcType == 2000 )  mesh % bedge(i) % bcType = 0
       !mesh%bedge(i)%bcType = 1000
        if ( mesh % bedge(i) % bcType == 4000 )  mesh % bedge(i) % bcType = 1000
       !mesh%bedge(i)%bcType = 2000
        ! <<< debug
      end do
      close(fpmesh)


      ! get the rest mesh information
      allocate( mesh%nedgez(mesh%nzone), mesh%nnodez(mesh%nzone) )


      allocate( psp(mesh%nnode) )
      do i = 1, mesh%nnode
        allocate( psp(i)%edge(mnbr), psp(i)%dir(mnbr), psp(i)%node2(mnbr), psp(i)%cell1(mnbr), psp(i)%cell2(mnbr), psp(i)%cellSide1(mnbr), psp(i)%cellSide2(mnbr) )
        psp(i)%npsp = 0
        psp(i)%dir(:) = 0
        psp(i)%node2(:) = 0
        psp(i)%cell1(:) = 0
        psp(i)%cell2(:) = 0
      enddo

      t = 0
      mesh%nedge = 0
      mesh%nedgez(:) = 0
      do 1 izone = 1, mesh%nzone
      do 1 i = 1, mesh%ntrianglez(izone)
        t = t + 1

        nn(1:3) = mesh%triangle(t)%node(1:3)
        nn(4) = nn(1)

        do 11 l = 1, 3
          nl = nn(l)
          nr = nn(l+1)
          n1 = min( nl, nr )
          n2 = max( nl, nr )

          counter = 1
          do j = 1, psp(n1)%npsp
            if( n2 == psp(n1)%node2(j) )  exit
            counter = counter + 1
          enddo

          if( counter == psp(n1)%npsp+1 )  then
            mesh%nedge = mesh%nedge + 1
            mesh%nedgez(izone) = mesh%nedgez(izone) + 1

            psp(n1)%npsp = counter
            psp(n1)%edge(counter) = mesh%nedge
            psp(n1)%node2(counter) = n2
          endif

          if( n1 == nl )  then
            psp(n1)%dir(counter) = 1
            psp(n1)%cell1(counter) = t
            psp(n1)%cellSide1(counter) = l
          else
            psp(n1)%dir(counter) = 2
            psp(n1)%cell2(counter) = t
            psp(n1)%cellSide2(counter) = l
          endif

11      continue
1     continue

!     edge
      allocate( mesh%edge(mesh%nedge) )
      do 2 i = 1, mesh%nnode
      do 2 j = 1, psp(i)%npsp
        iedge = psp(i)%edge(j)

        if( psp(i)%dir(j) == 1 )  then
          mesh%edge(iedge)%node1 = i
          mesh%edge(iedge)%node2 = psp(i)%node2(j)
          mesh%edge(iedge)%cell1 = psp(i)%cell1(j)
          mesh%edge(iedge)%cell2 = psp(i)%cell2(j)
          mesh%edge(iedge)%cellSide1 = psp(i)%cellSide1(j)
          mesh%edge(iedge)%cellSide2 = psp(i)%cellSide2(j)
        else
          mesh%edge(iedge)%node1 = psp(i)%node2(j)
          mesh%edge(iedge)%node2 = i
          mesh%edge(iedge)%cell1 = psp(i)%cell2(j)
          mesh%edge(iedge)%cell2 = psp(i)%cell1(j)
          mesh%edge(iedge)%cellSide1 = psp(i)%cellSide2(j)
          mesh%edge(iedge)%cellSide2 = psp(i)%cellSide1(j)
        endif
2     continue


      do i = 1, mesh%nnode
        deallocate( psp(i)%edge, psp(i)%dir, psp(i)%node2, psp(i)%cell1, psp(i)%cell2, psp(i)%cellSide1, psp(i)%cellSide2 )
      end do
      deallocate( psp )

end subroutine Import_mesh_from_grd

end module Mesh_Importer_Grd_mod
