module Common_Type_CRS_mod
      use KindDefinition_mod, only : DP
      implicit none
      save

      type :: CRS_t
        integer :: n
        integer :: ndim
        integer :: nnz
        integer, allocatable :: ia(:)
        integer, allocatable :: iau(:)
        integer, allocatable :: ja(:)
        real(DP), allocatable :: e(:,:,:)
      contains
        procedure, pass(A) :: Reallocate_CRS
        procedure, pass(A) :: Deallocate_CRS
        procedure, pass(A) :: Set_CRS_to_0
        procedure, pass(A) :: Set_CRS_to_unity
        procedure, pass(A) :: Add_a_block_to_CRS
        procedure, pass(A) :: Set_a_line_to_unity
        procedure, pass(A) :: Output_CRS
        procedure, pass(A) :: Output_CRS_tecplot
      end type CRS_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reallocate_CRS( A, ndim, n, nnz, ia, ja )
      implicit none
      class(CRS_t), intent(inout) :: A
      integer,      intent(in) :: ndim
      integer,      intent(in) :: n
      integer,      intent(in) :: nnz
      integer,      intent(in) :: ia(:)
      integer,      intent(in) :: ja(:)

      integer :: i, k

      if( allocated(A % ia) ) then
        deallocate( A % ia )
        deallocate( A % iau )
        deallocate( A % ja )
        deallocate( A % e )
      end if

      A % n = n
      A % ndim = ndim
      A % nnz = nnz

      allocate( A % ia(n+1) )
      allocate( A % iau(n) )
      allocate( A % ja(nnz) )
      allocate( A % e(ndim,ndim,nnz) )

      A % ia(:) = ia(:)
      A % ja(:) = ja(:)

      A % iau(:) = 0
      do i = 1, n
        do k = A % ia(i), A % ia(i+1)-1
          if( A % ja(k) == i ) then
            A % iau(i) = k
            exit
          end if
        end do
      end do

end subroutine Reallocate_CRS
!=====================================================================
!
!=====================================================================
subroutine Deallocate_CRS( A )
      implicit none
      class(CRS_t), intent(inout) :: A

      A % ndim = 0
      A % n = 0
      A % nnz = 0

      if( allocated(A % ia) ) then
        deallocate( A % ia )
        deallocate( A % iau )
        deallocate( A % ja )
        deallocate( A % e )
      end if

end subroutine Deallocate_CRS
!=====================================================================
!
!=====================================================================
subroutine Set_CRS_to_0( A )
      implicit none
      class(CRS_t), intent(inout) :: A

      A % e(:,:,:) = 0.

end subroutine Set_CRS_to_0
!=====================================================================
!
!=====================================================================
subroutine Set_CRS_to_unity( A )
      implicit none
      class(CRS_t), intent(inout) :: A

      integer :: i, k, ii

      A % e(:,:,:) = 0.

      do i = 1, A % n
        k = A % iau(i)
        do ii = 1, A % ndim
          A % e(ii,ii,k) = 1.
        end do
      end do

end subroutine Set_CRS_to_unity
!=====================================================================
!
!=====================================================================
subroutine Add_a_block_to_CRS( A, i, j, aa )
      implicit none
      class(CRS_t), intent(inout) :: A
      integer,      intent(in)    :: i, j
      real(DP),     intent(in)    :: aa(:,:)

      integer :: ndim
      logical :: Is_found
      integer :: k

      ndim = A % ndim

      Is_found = .false.
      do k = A % ia(i),A % ia(i+1)-1
        if( A % ja(k) == j ) then
          Is_found = .true.
          A % e(1:ndim,1:ndim,k) = A % e(1:ndim,1:ndim,k) + aa(1:ndim,1:ndim)
          exit
        end if
      end do

      if( .not. Is_found ) then
        write(*,*) 'wrong! Add_a_block_to_CRS', i, j
        write(*,*) A % ia(i), A % ia(i+1)-1
        write(*,*) A % ja( A % ia(i) : A % ia(i+1) - 1 )
        stop
      end if

end subroutine Add_a_block_to_CRS
!=====================================================================
!
!=====================================================================
subroutine Set_a_line_to_unity( A, i )
      implicit none
      class(CRS_t), intent(inout) :: A
      integer,      intent(in)    :: i

      integer :: k, ii

      A % e(:,:,A%ia(i):A%ia(i+1)-1) = 0.

      k = A % iau(i)
      do ii = 1, A % ndim
        A % e(ii,ii,k) = 1.
      end do

end subroutine Set_a_line_to_unity
!=====================================================================
!
!=====================================================================
subroutine Output_CRS( A, fp )
      implicit none
      class(CRS_t), intent(in) :: A
      integer,      intent(in) :: fp

      integer :: n, ndim
      integer :: i, j, k, ii, jj
      real(DP) :: aa(A%ndim)
      real(DP), parameter :: TOL = 1d-15
      logical :: Is_found

      write(fp,*) 'Output_CRS: Start'

      n = A%n
      ndim = A%ndim

      do 1 i = 1,n
        write(fp,110,advance='no') i,A%ia(i+1)-A%ia(i),':'
        do 11 j = 1,n
          Is_found = .false.
          do k = A%ia(i),A%ia(i+1)-1
            if( A%ja(k) == j ) then
              Is_found = .true.
              write(fp,120,advance='no') '*'
              exit
            end if
          end do
          if( .not. Is_found ) write(fp,120,advance='no') ' '
11      end do
        write(fp,*)
1     end do

110   format(i6,i3,a1)
120   format(a2)

      write(fp,*)

      do 2 i = 1,n
        do 21 ii = 1,ndim
          do 211 j = 1,n
            aa(:) = 0.
            Is_found = .false.
            do 2111 k = A%ia(i),A%ia(i+1)-1
              if( A%ja(k) == j )  then
                aa(:) = A%e(ii,:,k)
                Is_found = .true.
                exit
              end if
2111        end do

            do 2112 jj = 1,ndim
              if( Is_found ) then
                if( abs(aa(jj)) < TOL ) then
                  write(fp,220,advance='no') '         *     '
                else
                  write(fp,210,advance='no') aa(jj)
                end if
              else
                write(fp,220,advance='no') '               '
              end if
2112        end do

            write(fp,230,advance='no') '|'
211       end do
          write(fp,*)
21      end do

        do 22 j = 1,n
        do 22 jj = 1,ndim
          write(fp,240,advance='no') '___________________'
22      end do

        write(fp,*)
2     end do

210   format(ES16.6)
220   format(a16)
230   format(a2)
240   format(a18)

      write(fp,*) 'Output_CRS: End'
      write(fp,*)

end subroutine Output_CRS
!=====================================================================
!
!=====================================================================
subroutine Output_CRS_tecplot ( A, fp, prefix )
      implicit none
      class(CRS_t),     intent(in) :: A
      integer,          intent(in) :: fp
      character(len=*), intent(in) :: prefix

      character(len=100) :: filename
      integer :: i, j, k

      write(filename, '(a)') trim(prefix)
      open(fp, file=filename, form='formatted', status='replace')

      write(fp,*) 'variables ="j","-i"'
      write(fp,'(a,a)') 'zone T=',trim(prefix)

      do i = 1, A % n
        do k = A % ia(i), A % ia(i+1)-1
          j = A % ja(k)
          write(fp,*) j, -i
        end do
      end do

      close(fp)

end subroutine Output_CRS_tecplot

end module Common_Type_CRS_mod
