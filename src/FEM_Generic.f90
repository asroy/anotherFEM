module FEM_Generic_mod
      use KindDefinition_mod, only : DP
      use Constant_mod
      use Common_Type_CRS_mod
      use Common_Type_Dependency_mod
      use FEMEntry_Generic_mod
      use FEMEntry_LocalNode_Default_mod
      use FEMEntry_Quadrature_Gauss_mod
      use TimeControl_mod
      use Mesh_Generic_mod
      use GlobalVariableForDebug_mod
      implicit none
      save

! - FEM
      type :: FEM_t
        integer :: nequation

        integer :: nelement
        integer :: nbedge
        integer :: nfbasis
        integer :: nnode

        integer :: basis_order_START

        class(element_t), allocatable :: element(:)
        class(bedge_t),   allocatable :: bedge(:)
        type(fbasis_t),   allocatable :: fbasis(:)
        type(node_t),     allocatable :: node(:)

        class(quadrature_2d_t), allocatable :: quadrature_2d(:)
        class(quadrature_1d_t), allocatable :: quadrature_1d(:)
        class(local_node_t), allocatable :: local_node(:)

        type(depen_t) :: depen_fbasis
        type(depen_t) :: depen_matrix

        integer, allocatable :: ifb_to_imat(:)
        integer, allocatable :: imat_to_ifb(:)
      contains
        procedure, pass(this) :: Generate_FEM
        procedure, pass(this) :: Generate_FEM_connectivity
        procedure, pass(this) :: Generate_FEM_connectivity_high_order_mesh
        procedure, pass(this) :: Calculate_FEM_residual
        procedure, pass(this) :: Calculate_FEM_jacobian
        procedure, pass(this) :: Get_element_fbasis_ifb_nfb
        procedure, pass(this) :: Get_bedge_fbasis_ifb_nfbnz
        procedure, pass(this) :: Update_FEM_solution
        procedure, pass(this) :: Update_FEM_solution_history
        procedure, pass(this) :: Output_FEM_solution_tecplot
        procedure, pass(this) :: Write_q
        procedure, pass(this) :: Read_q
      end type FEM_t

! - equation
      type :: equation_t
        integer :: nequation
      contains
        procedure, pass(this) :: Reset_equation
        procedure, nopass     :: equation_element_polynomial_order
        procedure, nopass     :: equation_bedge_polynomial_order
        procedure, pass(this) :: Initialize_solution
        procedure, pass(this) :: Enforce_essential_boundary_condition
        procedure, pass(this) :: Calculate_dt_local
        procedure, pass(this) :: Calculate_dt_local_2
        procedure, pass(this) :: Calculate_interior_term
        procedure, pass(this) :: Calculate_interior_term_linearization
        procedure, pass(this) :: Calculate_interior_term_linearization_hand
        procedure, pass(this) :: Calculate_natural_boundary_term
        procedure, pass(this) :: Calculate_natural_boundary_term_linearization
        procedure, pass(this) :: max_omega_for_dq
        procedure, pass(this) :: Check_solution
      end type equation_t

contains

!=====================================================================
!
!=====================================================================
subroutine Generate_FEM ( this, mesh, equation )
      implicit none
      class(FEM_t),         intent(inout) :: this
      class(mesh_t),        intent(in)    :: mesh
      class(equation_t),    intent(in)    :: equation

      integer, allocatable :: npsp(:), psp(:,:)
      integer, parameter :: mnbr = 60
      logical :: Has_stored
      integer :: i, is, k, ii, nc, ifb, ifbs, imat, jfb, jmat, iel, ibasis_order, ipoly_order

! - get number of equation
      this%nequation = equation%nequation

! - define quadrature mapping, local node mapping
      allocate ( gauss_quadrature_2d_t :: this%quadrature_2d(nbasis_order_MAX) )
      allocate ( gauss_quadrature_1d_t :: this%quadrature_1d(nbasis_order_MAX) )
      allocate ( default_local_node_t  :: this%local_node   (nbasis_order_MAX) )

      do ibasis_order = 1, nbasis_order_MAX
        ipoly_order = equation%equation_element_polynomial_order( ibasis_order )
        call this%quadrature_2d(ibasis_order)%Reset_quadrature_2d( ipoly_order )

        ipoly_order = equation%equation_bedge_polynomial_order( ibasis_order )
        call this%quadrature_1d(ibasis_order)%Reset_quadrature_1d( ipoly_order )

        call this%local_node(ibasis_order)%Reset_local_node(ibasis_order)
      end do

! - basis order
      this%basis_order_START = 3

! - generate FEM connectivty
      call this%Generate_FEM_connectivity ( mesh )
     !call this%Generate_FEM_connectivity_high_order_mesh ( mesh )

! - depen_fbasis
      allocate( npsp(this%nfbasis), psp(mnbr,this%nfbasis) )

      npsp(:) = 0
      do iel = 1, this%nelement
        do i = 1, this%element(iel)%nfb
          ifb = this%element(iel)%ifb(i)
          do is = 1, this%element(iel)%nfb
            ifbs = this%element(iel)%ifb(is)
            Has_stored = .false.
            do ii = 1, npsp(ifb)
              if( ifbs == psp(ii,ifb) ) then
                Has_stored = .true.
                exit
              end if
            end do
            if( .not. Has_stored ) then
              npsp(ifb) = npsp(ifb) + 1
              if( npsp(ifb) > mnbr ) then
                write(*,*) 'wrong! Generate_FEM: mnbr too small'
                stop
              end if
              psp( npsp(ifb), ifb ) = ifbs
            end if
          end do
        end do
      end do

      nc = 0
      do i = 1,this%nfbasis
        nc = nc + npsp(i)
      end do

      call this%depen_fbasis%Reallocate_depen ( this%nfbasis, nc )

      this%depen_fbasis%ia(1) = 1
      do i = 1,this%nfbasis
        this%depen_fbasis%ia(i+1) = this%depen_fbasis%ia(i) + npsp(i)
      end do

      k = 0
      do i = 1,this%nfbasis
        do ii = 1, npsp(i)
          k = k + 1
          is = psp(ii,i)
          this%depen_fbasis%ja(k) = is
        end do
      end do

      deallocate ( npsp, psp )

! - ifb_to_imat, imat_to_ifb
      allocate ( this%ifb_to_imat(this%nfbasis) )
      allocate ( this%imat_to_ifb(this%nfbasis) )

      call cuthill_mckee( this%depen_fbasis%n, this%depen_fbasis%nnz, this%depen_fbasis%ia, this%depen_fbasis%ja, this%imat_to_ifb )

      do imat = 1, this%nfbasis
        ifb = this%imat_to_ifb(imat)
        this%ifb_to_imat(ifb) = imat
      end do


      ! >>> debug: turn off reorder
      write(*,*) 'debug! Generate_FEM: reorder turned off'
     !do i = 1, this % nfbasis
     !  this % imat_to_ifb(i) = i
     !  this % ifb_to_imat(i) = i
     !end do
      ! <<< debug


! - depen_matrix
      call this%depen_matrix%Reallocate_depen ( this%depen_fbasis%n, this%depen_fbasis%nnz )

      nc = 0
      this%depen_matrix%ia(1) = 1
      do imat = 1, this%depen_matrix%n
        ifb = this%imat_to_ifb(imat)
        this%depen_matrix%ia(imat+1) = this%depen_matrix%ia(imat) + this%depen_fbasis%ia(ifb+1) - this%depen_fbasis%ia(ifb)

        do k = this%depen_fbasis%ia(ifb), this%depen_fbasis%ia(ifb+1)-1
          jfb = this%depen_fbasis%ja(k)
          jmat = this%ifb_to_imat(jfb)

          nc = nc + 1
          this%depen_matrix%ja(nc) = jmat
        end do

        call sort ( this%depen_matrix%ja, this%depen_matrix%nnz, this%depen_matrix%ia(imat), this%depen_matrix%ia(imat+1)-1 )

      end do


      contains

      !=============================================================================
      !
      ! reorder using Cuthill-McKee algorithm
      !
      !=============================================================================
      subroutine cuthill_mckee( n, nnz, ia, ja, newtoold )
            implicit none
      ! - in
            integer, intent(in) :: n, nnz
            integer, intent(in) :: ia(n+1), ja(nnz)
      ! - out
            integer, intent(out) :: newtoold(n)
      ! - local
            integer :: mnbr = 100
            integer, allocatable :: psp_tmp(:,:), npsp_tmp(:)
            integer, allocatable :: oldtonew(:)
      !     index, counter
            integer :: i, j, k, i1, i2, i1new, i2new
            integer :: is, p, c, nei, counter_r, nnew, nold
            integer :: mnei, num_g, num_r
            integer :: nel, ned
            integer :: flag
            integer :: ntmp


      !
      ! - A. generate point-to-neigbhor_point map
            allocate( npsp_tmp(n) )
            allocate( psp_tmp(mnbr,n) )

            npsp_tmp(:) = 0
            psp_tmp(:,:) = -1


            do 1 i = 1, n
      !       add onegrid connectivity
              do 11  k = ia(i), ia(i+1)-1
                j = ja(k)
                if( j == i )  goto 11

                ntmp = npsp_tmp(i) + 1
                if( ntmp > mnbr ) then
                  write(*,*) 'increase mnbr in Cuthill-McKee!'
                  stop
                end if
                npsp_tmp(i) = ntmp
                psp_tmp(ntmp,i) = j
      11      continue
      1     continue

      !
      ! - B. bandwith reduction using Cuthill-McKee algorithm
      ! - reordering: get new sequence of nodes
            allocate( oldtonew(n) )

            counter_r = 0
            num_r = 0
            num_g = n
            do i=1,n
              oldtonew(i) = i
            end do

            do 6 while( num_g > 0 )
      !       find node p in oldtonew with the least neigbhors
              mnei = 2*n
              ntmp = -1
              do 61 i=1,n
                if( oldtonew(i) /= -1 .and. npsp_tmp(i) < mnei )  then
                  mnei = npsp_tmp(i)
                  ntmp = i
                end if
      61      continue
              p = ntmp

      !       move p from oldtonew into newtoold
              num_g = num_g - 1
              oldtonew(p) = -1

              num_r = num_r + 1
              newtoold(num_r) = p


              do 62 while( counter_r < num_r )
      !         pick out c: 1st unsorted element in newtoold
                counter_r = counter_r + 1
                c = newtoold(counter_r)


      !         move neighbour of c from oldtonew into newtoold
                do 621 is=1,npsp_tmp(c)
                  mnei = 2*n
                  ntmp = -1
                  do 6211 i=1,npsp_tmp(c)
                    j = psp_tmp(i,c)
                    if( oldtonew(j) /= -1 .and. npsp_tmp(j) < mnei )  then
                      mnei = npsp_tmp(j)
                      ntmp = j
                    end if
      6211        continue

                  if( ntmp == -1 )  goto 621
                  i = ntmp

                  num_g = num_g - 1
                  oldtonew(i) = -1

                  num_r = num_r + 1
                  newtoold(num_r) = i

      621       continue
      62      continue
      6     continue

            deallocate( oldtonew, psp_tmp, npsp_tmp )

      end subroutine cuthill_mckee
      !=====================================================================
      !
      !=====================================================================
      subroutine sort ( ja, nnz, istart, iend )
            integer, intent(in)    :: nnz, istart, iend
            integer, intent(inout) :: ja(nnz)

            integer i, j, min, minsave, jsave

            do i = istart, iend
              min = ja(i)
              minsave = ja(i)
              jsave = i
              do j = i+1, iend
                if( ja(j) < min ) then
                  min = ja(j)
                  jsave = j
                end if
              end do
              ja(i) = min
              ja(jsave) = minsave

            end do
      end subroutine sort

end subroutine Generate_FEM
!=====================================================================
!
!=====================================================================
subroutine Generate_FEM_connectivity ( this, mesh )
      implicit none
      class(FEM_t),      intent(inout) :: this
      class(mesh_t),     intent(in)    :: mesh
end subroutine Generate_FEM_connectivity
!=====================================================================
!
!=====================================================================
subroutine Generate_FEM_connectivity_high_order_mesh ( this, mesh )
      implicit none
      class(FEM_t),      intent(inout) :: this
      class(mesh_t),     intent(in)    :: mesh
end subroutine Generate_FEM_connectivity_high_order_mesh
!=====================================================================
!
!=====================================================================
subroutine Calculate_FEM_residual ( this, res, equation, time_control )
      implicit none
      class(FEM_t),                                                  intent(inout) :: this
      real(DP),              dimension(this%nequation,this%nfbasis), intent(out)   :: res
      class(equation_t),                                             intent(in)    :: equation
      type(time_control_t),                                          intent(in)    :: time_control

      type(fbasis_t), dimension(nfb_MAX)   :: fbasis_nfb
      type(fbasis_t), dimension(nfbnz_MAX) :: fbasis_nfbnz
      integer, dimension(nfb_MAX)   :: ifb_nfb, imat_nfb
      integer, dimension(nfbnz_MAX) :: ifb_nfbnz, imat_nfbnz
      real(DP), dimension(nequation_MAX,nfb_MAX) :: rri, rr
      integer :: nq, nfbasis
      integer :: nfb, nfbnz
      integer :: basis_order, bcType
      real(DP) :: r, s, w
      integer :: i, k, iel, ifb, imat, ibedge

      write(fpd1,*) 'Calculate_FEM_residual: Start'

      nq = this%nequation
      nfbasis = this%nfbasis

      res(:,:) = 0.

     !goto 1000         ! debug


! - interior integral
      do iel = 1, this % nelement
        basis_order = this % element(iel) % basis_order

        call this % Get_element_fbasis_ifb_nfb ( fbasis_nfb, ifb_nfb, nfb, iel )

        imat_nfb(1:nfb) = this % ifb_to_imat(ifb_nfb(1:nfb))

        rri(1:nq,1:nfb) = 0.

        do k = 1, this % quadrature_2d(basis_order) % nqd
          r = this % quadrature_2d(basis_order) % rqd(k)
          s = this % quadrature_2d(basis_order) % sqd(k)
          w = this % quadrature_2d(basis_order) % wqd(k)


          Do_print = .false.

          
          if( Do_print ) then
            write(*,*) 'before 1 gauss'
            write(*,*) 'iel,ifb_nfb',iel,',', ifb_nfb(1:nfb)
            write(*,200) 'k,r,s,w', k, r, s, w
          end if

          call equation % Calculate_interior_term ( rr, this % element(iel), fbasis_nfb, time_control, r, s )


          if( Do_print ) then
            write(*,*) 'after 1 gauss'
            do i = 1, nfb
              write(*,100) 'rr', rr(1:nq,i)
            end do
          end if


          Do_print = .false.

          rri(1:nq,1:nfb) = rri(1:nq,1:nfb) + w*rr(1:nq,1:nfb)
        end do


        ! >>> debug
       !if( any( ifb_nfb(1:nfb) == 333 ) ) then
       !  write(*,*) 'iel,ifb_nfb',iel,',', ifb_nfb(1:nfb)
       !  do i = 1, nfb
       !    write(*,200) 'rri', ifb_nfb(i), rri(1:nq,i)
       !  end do
       !end if
        ! <<< debug


        do i = 1,nfb
          imat = imat_nfb(i)
          res(1:nq,imat) = res(1:nq,imat) + rri(1:nq,i)
        end do
      end do


      ! >>> debug
      do i = 1, this%nfbasis
        imat = this%ifb_to_imat(i)
        this%fbasis(i)%evector(1:nq) = res(1:nq,imat)
      end do

      call this%Output_FEM_solution_tecplot ( fpIO, 'output/residual-no_boundary.plt')
      ! <<< debug


1000  continue


! - natural BCs
      do ibedge = 1, this % nbedge

        bcType = this % bedge(ibedge) % bcType

        if ( bcType >= 0 .and. bcType <= 999 )  cycle

        basis_order = this % bedge(ibedge) % basis_order
        iel         = this % bedge(ibedge) % iel

        call this % Get_bedge_fbasis_ifb_nfbnz ( fbasis_nfbnz, ifb_nfbnz, nfbnz, ibedge )
        call this % Get_element_fbasis_ifb_nfb ( fbasis_nfb,   ifb_nfb,   nfb,   iel )

        imat_nfbnz(1:nfbnz) = this % ifb_to_imat(ifb_nfbnz(1:nfbnz))
        imat_nfb  (1:nfb)   = this % ifb_to_imat(ifb_nfb  (1:nfb))

        rri(1:nq,1:nfbnz) = 0.

        do k = 1, this % quadrature_1d(basis_order) % nqd
          r = this % quadrature_1d(basis_order) % rqd(k)
          w = this % quadrature_1d(basis_order) % wqd(k)

          call equation % Calculate_natural_boundary_term ( rr, this % bedge(ibedge), fbasis_nfbnz,  this % element(iel), fbasis_nfb, time_control, r )

          rri(1:nq,1:nfbnz) = rri(1:nq,1:nfbnz) + w*rr(1:nq,1:nfbnz)
        end do

        do i = 1, nfbnz
          imat = imat_nfbnz(i)
          res(1:nq,imat) = res(1:nq,imat) + rri(1:nq,i)
        end do
      end do


2000  continue


! - essential BCs (should be applied after natural BCs
      do ibedge = 1, this % nbedge
        bcType = this % bedge(ibedge) % bcType

        if ( bcType >= 0 .and. bcType <= 999 ) then
          do i = 1, this % bedge(ibedge) % nfbnz
            ifb = this % bedge(ibedge) % ifbnz(i)
            imat = this % ifb_to_imat(ifb)
            res(1:nq,imat) = 0.
          end do
        end if
      end do


      ! >>> debug
      do i = 1, this%nfbasis
        imat = this%ifb_to_imat(i)
        this%fbasis(i)%evector(1:nq) = res(1:nq,imat)
      end do

      call this%Output_FEM_solution_tecplot ( fpIO, 'output/residual.plt')
      ! <<< debug


      write(fpd1,*) 'Calculate_FEM_residual: End'

100   format(a,20(es16.5))
200   format(a,i7,20(es16.5))

end subroutine Calculate_FEM_residual
!=====================================================================
!
!=====================================================================
subroutine Calculate_FEM_jacobian ( this, A, equation, time_control )
      implicit none
      class(FEM_t),         intent(inout) :: this
      class(CRS_t),         intent(inout) :: A
      class(equation_t),    intent(in)    :: equation
      type(time_control_t), intent(in)    :: time_control

      type(fbasis_t), dimension(nfb_MAX)   :: fbasis_nfb
      type(fbasis_t), dimension(nfbnz_MAX) :: fbasis_nfbnz
      integer, dimension(nfb_MAX)   :: ifb_nfb, imat_nfb
      integer, dimension(nfbnz_MAX) :: ifb_nfbnz, imat_nfbnz
      real(DP), dimension(nequation_MAX,nequation_MAX,nfb_MAX,nfb_MAX):: aai, aa
      integer :: nq
      integer :: nfb, nfbnz
      integer :: basis_order, bcType
      real(DP) :: r, s, w
      integer :: i, j, k, iel, ifb, imat, jmat, ibedge

      real(DP) :: emat(nequation_MAX,nequation_MAX)


      write(fpd1,*) 'Calculate_FEM_jacobian: Start'

      nq = this % nequation

      call A % Set_CRS_to_0


     !goto 1000       ! debug


! - element integral
      do iel = 1, this % nelement
        basis_order = this % element(iel) % basis_order

        call this % Get_element_fbasis_ifb_nfb ( fbasis_nfb, ifb_nfb, nfb, iel )

        do i = 1, nfb
          imat_nfb(i) = this % ifb_to_imat(ifb_nfb(i))
        end do

        aai(1:nq,1:nq,1:nfb,1:nfb) = 0.

        do k = 1, this % quadrature_2d(basis_order) % nqd
          r = this % quadrature_2d(basis_order) % rqd(k)
          s = this % quadrature_2d(basis_order) % sqd(k)
          w = this % quadrature_2d(basis_order) % wqd(k)

          call equation % Calculate_interior_term_linearization      ( aa, this % element(iel), fbasis_nfb, time_control, r, s )
         !call equation % Calculate_interior_term_linearization_hand ( aa, this % element(iel), fbasis_nfb, time_control, r, s )

          aai(1:nq,1:nq,1:nfb,1:nfb) = aai(1:nq,1:nq,1:nfb,1:nfb) + w*aa(1:nq,1:nq,1:nfb,1:nfb)
        end do

        do j = 1, nfb
          jmat = imat_nfb(j)
          do i = 1, nfb
            imat = imat_nfb(i)
            call A % Add_a_block_to_CRS( imat, jmat, aai(:,:,i,j) )
          end do
        end do
      end do


      ! >>> debug
      ! add row-wise
      do ifb = 1, this%nfbasis
        this%fbasis(ifb)%ematrix(:,:) = 0.
      end do

      do i = 1, A%n
        ifb = this%imat_to_ifb(i)
        do k = A%ia(i), A%ia(i+1)-1
          this%fbasis(ifb)%ematrix(1:nq,1:nq) = this%fbasis(ifb)%ematrix(1:nq,1:nq) + abs( A%e(1:nq,1:nq,k) )
         !this%fbasis(ifb)%ematrix(1:nq,1:nq) = this%fbasis(ifb)%ematrix(1:nq,1:nq) +      A%e(1:nq,1:nq,k)
        end do
      end do

      call this%output_fem_solution_tecplot( fpIO, 'output/jacobian_add_row-no_boundary.plt' )

      ! add column-wise
      do ifb = 1, this%nfbasis
        this%fbasis(ifb)%ematrix(:,:) = 0.
      end do

      do k = 1, A%nnz
        i = A%ja(k)
        ifb = this%imat_to_ifb(i)
        this%fbasis(ifb)%ematrix(1:nq,1:nq) = this%fbasis(ifb)%ematrix(1:nq,1:nq) + abs (A%e(1:nq,1:nq,k) )
       !this%fbasis(ifb)%ematrix(1:nq,1:nq) = this%fbasis(ifb)%ematrix(1:nq,1:nq) +      A%e(1:nq,1:nq,k)
      end do

      call this%output_fem_solution_tecplot( fpIO, 'output/jacobian_add_column-no_boundary.plt' )
      ! <<< debug


1000  continue


! - natural BCs
      do ibedge = 1, this % nbedge

        bcType = this % bedge(ibedge) % bcType

        if ( bcType >= 0 .and. bcType <= 999 )  cycle

        basis_order = this % bedge(ibedge) % basis_order
        iel         = this % bedge(ibedge) % iel

        call this % Get_bedge_fbasis_ifb_nfbnz ( fbasis_nfbnz, ifb_nfbnz, nfbnz, ibedge )
        call this % Get_element_fbasis_ifb_nfb ( fbasis_nfb,   ifb_nfb,   nfb,   iel )

        imat_nfbnz(1:nfbnz) = this % ifb_to_imat(ifb_nfbnz(1:nfbnz))
        imat_nfb  (1:nfb)   = this % ifb_to_imat(ifb_nfb  (1:nfb))

        aai(1:nq,1:nq,1:nfbnz,1:nfb) = 0.

        do k = 1, this % quadrature_1d(basis_order) % nqd
          r = this % quadrature_1d(basis_order) % rqd(k)
          w = this % quadrature_1d(basis_order) % wqd(k)

          call equation % Calculate_natural_boundary_term_linearization ( aa, this % bedge(ibedge), fbasis_nfbnz,  this % element(iel), fbasis_nfb, time_control, r )

          aai(1:nq,1:nq,1:nfbnz,1:nfb) = aai(1:nq,1:nq,1:nfbnz,1:nfb) + w*aa(1:nq,1:nq,1:nfbnz,1:nfb)
        end do

        do j = 1, nfb
          jmat = imat_nfb(j)
          do i = 1, nfbnz
            imat = imat_nfbnz(i)
            call A % Add_a_block_to_CRS( imat, jmat, aai(:,:,i,j) )
          end do
        end do

      end do

2000  continue


! - essential BCs (should be applied after natural BCs
      do ibedge = 1, this % nbedge
        bcType = this % bedge(ibedge) % bcType

        if ( bcType >= 0 .and. bcType <= 999 ) then
          do i = 1, this % bedge(ibedge) % nfbnz
            ifb = this % bedge(ibedge) % ifbnz(i)
            imat = this % ifb_to_imat(ifb)
            call A % Set_a_line_to_unity (imat)
          end do
        end if
      end do

      ! >>> debug
      ! add row-wise
      do ifb = 1, this%nfbasis
        this%fbasis(ifb)%ematrix(:,:) = 0.
      end do

      do i = 1, A%n
        ifb = this%imat_to_ifb(i)
        do k = A%ia(i), A%ia(i+1)-1
          this%fbasis(ifb)%ematrix(1:nq,1:nq) = this%fbasis(ifb)%ematrix(1:nq,1:nq) + abs( A%e(1:nq,1:nq,k) )
         !this%fbasis(ifb)%ematrix(1:nq,1:nq) = this%fbasis(ifb)%ematrix(1:nq,1:nq) +      A%e(1:nq,1:nq,k)
        end do
      end do

      call this%output_fem_solution_tecplot( fpIO, 'output/jacobian_add_row.plt' )


      ! add column-wise
      do ifb = 1, this%nfbasis
        this%fbasis(ifb)%ematrix(:,:) = 0.
      end do

      do k = 1, A%nnz
        i = A%ja(k)
        ifb = this%imat_to_ifb(i)
        this%fbasis(ifb)%ematrix(1:nq,1:nq) = this%fbasis(ifb)%ematrix(1:nq,1:nq) + abs (A%e(1:nq,1:nq,k) )
       !this%fbasis(ifb)%ematrix(1:nq,1:nq) = this%fbasis(ifb)%ematrix(1:nq,1:nq) +      A%e(1:nq,1:nq,k)
      end do

      call this%output_fem_solution_tecplot( fpIO, 'output/jacobian_add_column.plt' )
      ! <<< debug


      ! >>> debug
     !emat(1:nq,1:nq) = 0.
     !do i = 1, this % nfbasis
     !  emat(1:nq,1:nq) = emat(1:nq,1:nq) + abs( this % fbasis(i) % ematrix(1:nq,1:nq) )
     !end do
     !write(*,*) 'dddddddddddddddddddddd', emat(1:nq,1:nq)
      ! <<< debug


      write(fpd1,*) 'Calculate_FEM_jacobian: End'

100   format(a,20(es16.5))

end subroutine Calculate_FEM_jacobian
!=====================================================================
!
!=====================================================================
subroutine Get_element_fbasis_ifb_nfb ( this, fbasis_nfb, ifb_nfb, nfb, iel )
      implicit none
      class(FEM_t),                 intent(in)  :: this
      integer,                      intent(out) :: nfb
      type(fbasis_t), dimension(:), intent(out) :: fbasis_nfb
      integer,        dimension(:), intent(out) :: ifb_nfb
      integer,                      intent(in)  :: iel

      integer :: i, ifb

      nfb = this % element(iel) % nfb

      do i = 1, nfb
        ifb = this % element(iel) % ifb(i)
        fbasis_nfb(i) = this%fbasis(ifb)
        ifb_nfb(i) = ifb
      end do

end subroutine Get_element_fbasis_ifb_nfb
!=====================================================================
!
!=====================================================================
subroutine Get_bedge_fbasis_ifb_nfbnz ( this, fbasis_nfbnz, ifb_nfbnz, nfbnz, ibedge )
      implicit none
      class(FEM_t),                 intent(in)  :: this
      integer,                      intent(in)  :: ibedge
      type(fbasis_t), dimension(:), intent(out) :: fbasis_nfbnz
      integer,        dimension(:), intent(out) :: ifb_nfbnz
      integer,                      intent(out) :: nfbnz

      integer :: i, ifb

      nfbnz = this % bedge(ibedge) % nfbnz

      do i = 1, nfbnz
        ifb = this % bedge(ibedge) % ifbnz(i)
        fbasis_nfbnz(i) = this % fbasis(ifb)
        ifb_nfbnz(i) = ifb
      end do

end subroutine Get_bedge_fbasis_ifb_nfbnz
!=====================================================================
!
!=====================================================================
subroutine Update_FEM_solution ( this, dq, omega, itimeStep )
      implicit none
      class(FEM_t),                 intent(inout) :: this
      real(DP),     dimension(:,:), intent(in)    :: dq
      real(DP),                     intent(in)    :: omega
      integer,                      intent(in)    :: itimeStep

      integer :: imat, ifb
      integer :: nq

      nq = this % nequation

      do ifb = 1, this % nfbasis
        imat = this % ifb_to_imat(ifb)
        this % fbasis(ifb) % q(1:nq,itimeStep) = this % fbasis(ifb) % q(1:nq,itimeStep) + omega*dq(1:nq,imat)
      end do

end subroutine Update_FEM_solution
!=====================================================================
!
!=====================================================================
subroutine Update_FEM_solution_history ( this, ntimeStep_low, ntimeStep_up )
      implicit none
      class(FEM_t),  intent(inout) :: this
      integer,       intent(in)    :: ntimeStep_low, ntimeStep_up

      integer :: nq
      integer :: itimeStep, ifb

      nq = this%nequation

      do ifb = 1, this%nfbasis
      do itimeStep = ntimeStep_low, ntimeStep_up-1
        this%fbasis(ifb)%q(1:nq,itimeStep) = this%fbasis(ifb)%q(1:nq,itimeStep+1)
      end do
      end do

end subroutine Update_FEM_solution_history
!=====================================================================
!
!=====================================================================
subroutine Output_FEM_solution_tecplot( this, fp, filename )
      implicit none
      class(FEM_t),     intent(in) :: this
      integer,          intent(in) :: fp
      character(len=*), intent(in) :: filename
end subroutine Output_FEM_solution_tecplot
!=====================================================================
!
!=====================================================================
subroutine Write_q ( this, itimeStep, fp, filename )
      implicit none
      class(FEM_t),      intent(in) :: this
      integer,           intent(in) :: itimeStep
      integer,           intent(in) :: fp
      character(len=*),  intent(in) :: filename

      integer :: nq
      integer :: i

      nq = this % nequation

      open(fp, file=trim(filename), form='unformatted', status='unknown')
      do i = 1, this % nfbasis
        write(fp) this % fbasis(i) % q(1:nq,itimeStep)
      end do
      close(fp)

      write(*,*) 'Write_q:'//trim(filename)

end subroutine Write_q
!=====================================================================
!
!=====================================================================
subroutine Read_q ( this, itimeStep, fp, filename )
      implicit none
      class(FEM_t),      intent(inout) :: this
      integer,           intent(in)    :: itimeStep
      integer,           intent(in)    :: fp
      character(len=*),  intent(in)    :: filename

      integer :: nq
      integer :: i

      nq = this % nequation

      open(fp, file=trim(filename), form='unformatted', status='old')
      do i = 1, this % nfbasis
        read(fp) this % fbasis(i) % q(1:nq,itimeStep)
      end do
      close(fp)

      write(*,*) 'Read_q:'//trim(filename)

end subroutine Read_q
!=====================================================================
!
!=====================================================================
subroutine Reset_equation( this )
      implicit none
      class(equation_t), intent(out) :: this
end subroutine Reset_equation
!=====================================================================
!
!=====================================================================
function equation_nequation ( this )  result( nequation )
      implicit none
      class(equation_t), intent(out) :: this
      integer                        :: nequation
end function equation_nequation
!=====================================================================
!
!=====================================================================
function equation_element_polynomial_order( basis_order )  result ( polynomial_order )
      implicit none
      integer, intent(in) :: basis_order
      integer             :: polynomial_order
end function equation_element_polynomial_order
!=====================================================================
!
!=====================================================================
function equation_bedge_polynomial_order( basis_order )  result ( polynomial_order )
      implicit none
      integer, intent(in) :: basis_order
      integer             :: polynomial_order
end function equation_bedge_polynomial_order
!=====================================================================
!
!=====================================================================
subroutine Initialize_solution ( this, FEM )
      implicit none
      class(equation_t), intent(in)    :: this
      class(FEM_t),      intent(inout) :: FEM
end subroutine Initialize_solution
!=====================================================================
!
!=====================================================================
subroutine Enforce_essential_boundary_condition ( this, FEM, itimeStep )
      implicit none
      class(equation_t), intent(in)    :: this
      class(FEM_t),      intent(inout) :: FEM
      integer,           intent(in)    :: itimeStep
end subroutine Enforce_essential_boundary_condition
!=====================================================================
!
!=====================================================================
subroutine Calculate_dt_local ( this, FEM, time_control )
      implicit none
      class(equation_t),    intent(in)    :: this
      class(FEM_t),         intent(inout) :: FEM
      type(time_control_t), intent(in)    :: time_control
end subroutine Calculate_dt_local
!=====================================================================
!
!=====================================================================
subroutine Calculate_dt_local_2 ( this, FEM, time_control )
      implicit none
      class(equation_t),    intent(in)    :: this
      class(FEM_t),         intent(inout) :: FEM
      type(time_control_t), intent(in)    :: time_control
end subroutine Calculate_dt_local_2
!=====================================================================
!
!=====================================================================
subroutine Calculate_interior_term ( this, res, element, fbasis_nfb, time_control, r, s )
      implicit none
      class(equation_t),                    intent(in)  :: this
      real(DP),             dimension(:,:), intent(out) :: res
      class(element_t),                     intent(in)  :: element
      type(fbasis_t),       dimension(:),   intent(in)  :: fbasis_nfb
      type(time_control_t),                 intent(in)  :: time_control
      real(DP),                             intent(in)  :: r, s
end subroutine Calculate_interior_term
!=====================================================================
!
!=====================================================================
subroutine Calculate_interior_term_linearization ( this, aa, element, fbasis_nfb, time_control, r, s )
      implicit none
      class(equation_t),                        intent(in)  :: this
      real(DP),             dimension(:,:,:,:), intent(out) :: aa
      class(element_t),                         intent(in)  :: element
      type(fbasis_t),       dimension(:),       intent(in)  :: fbasis_nfb
      type(time_control_t),                     intent(in)  :: time_control
      real(DP),                                 intent(in)  :: r, s
end subroutine Calculate_interior_term_linearization
!=====================================================================
!
!=====================================================================
subroutine Calculate_interior_term_linearization_hand ( this, aa, element, fbasis_nfb, time_control, r, s )
      implicit none
      class(equation_t),                        intent(in)  :: this
      real(DP),             dimension(:,:,:,:), intent(out) :: aa
      class(element_t),                         intent(in)  :: element
      type(fbasis_t),       dimension(:),       intent(in)  :: fbasis_nfb
      type(time_control_t),                     intent(in)  :: time_control
      real(DP),                                 intent(in)  :: r, s

      write(*,*) 'wrong! Caclulate_interior_term_linearization_hand: undefined'
      stop
end subroutine Calculate_interior_term_linearization_hand
!=====================================================================
!
!=====================================================================
subroutine Calculate_natural_boundary_term ( this, res, bedge, fbasis_nfbnz,  element, fbasis_nfb, time_control, rb )
      implicit none
      class(equation_t),                    intent(in)  :: this
      class(bedge_t),                       intent(in)  :: bedge
      type(fbasis_t),       dimension(:),   intent(in)  :: fbasis_nfbnz
      class(element_t),                     intent(in)  :: element
      type(fbasis_t),       dimension(:),   intent(in)  :: fbasis_nfb
      type(time_control_t),                 intent(in)  :: time_control
      real(DP),             dimension(:,:), intent(out) :: res
      real(DP),                             intent(in)  :: rb
end subroutine Calculate_natural_boundary_term
!=====================================================================
!
!=====================================================================
subroutine Calculate_natural_boundary_term_linearization ( this, aa, bedge, fbasis_nfbnz,  element, fbasis_nfb, time_control, rb )
      implicit none
      class(equation_t),                        intent(in)  :: this
      class(bedge_t),                           intent(in)  :: bedge
      type(fbasis_t),       dimension(:),       intent(in)  :: fbasis_nfbnz
      class(element_t),                         intent(in)  :: element
      type(fbasis_t),       dimension(:),       intent(in)  :: fbasis_nfb
      type(time_control_t),                     intent(in)  :: time_control
      real(DP),             dimension(:,:,:,:), intent(out) :: aa
      real(DP),                                 intent(in)  :: rb
end subroutine Calculate_natural_boundary_term_linearization
!=====================================================================
!
!=====================================================================
function max_omega_for_dq ( this, FEM, dq, itimeStep )  result( omega_max )
      implicit none
      class(equation_t), intent(in)  :: this
      class(FEM_t),      intent(in)  :: FEM
      real(DP),          intent(in)  :: dq(FEM%nequation,FEM%nfbasis)
      integer,           intent(in)  :: itimeStep
      real(DP)                       :: omega_max
end function max_omega_for_dq
!=====================================================================
!
!=====================================================================
subroutine Check_solution ( this, FEM, itimeStep, fp, filename )
      implicit none
      class(equation_t), intent(in) :: this
      class(FEM_t),      intent(in) :: FEM
      integer,           intent(in) :: itimeStep
      integer,           intent(in) :: fp
      character(len=*),  intent(in) :: filename
end subroutine Check_solution

end module FEM_Generic_mod
