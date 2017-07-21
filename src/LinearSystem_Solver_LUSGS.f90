module LinearSystem_Solver_LUSGS_mod
      use KindDefinition_mod, only : DP
      use LinearSystem_Generic_mod
      use IO_mod
      use GlobalVariableForDebug_mod
      implicit none
      save

      type, extends(linear_solver_t) :: LUSGS_t
        private
        integer :: niter
        real(DP) :: etol
        integer,  allocatable :: P(:,:)
        real(DP), allocatable :: y(:)
        real(DP), allocatable :: z(:)
        real(DP), allocatable :: Ax(:)
        real(DP), allocatable :: Uxold(:,:)
        real(DP), allocatable :: Uxolder(:)
      contains
        procedure, pass(this) :: Reallocate_linear_solver => Reallocate_LUSGS
        procedure, pass(this) :: Solve_linear_system => LUSGS
      end type LUSGS_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reallocate_LUSGS( this, LS )
      implicit none
      class(LUSGS_t), intent(out) :: this
      class(linear_system_t), intent(in) :: LS
      integer :: ndim
      integer :: n


      if( allocated(this%P) )  then
        deallocate( this%P )
        deallocate( this%y )
        deallocate( this%z )
        deallocate( this%Ax )
        deallocate( this%Uxold )
        deallocate( this%Uxolder )
      endif

      this%niter = 10000
      this%etol  = 1.d-14

      ndim = LS%ndim
      n = LS%n

      allocate( this%P(ndim,n) )
      allocate( this%y(ndim) )
      allocate( this%z(ndim) )
      allocate( this%Ax(ndim) )
      allocate( this%Uxold(ndim,n) )
      allocate( this%Uxolder(ndim) )

end subroutine Reallocate_LUSGS
!=====================================================================
!
!=====================================================================
subroutine LUSGS( this, LS, preconditioner )
      implicit none
      class(LUSGS_t),          intent(inout) :: this
      class(linear_system_t),  intent(inout) :: LS
      class(preconditioner_t), intent(inout) :: preconditioner

      real(DP) :: err0, edrop
      integer :: n, ndim
      integer :: i, j, k, kau, ii, iis, iip, iip_mat, ii_mat, iis_mat, iter
      real(DP) :: pvt, tmp, err

      n = LS%n
      ndim = LS%ndim

!     initialize solution
      LS%x(:,:) = 0.

!     LU decomposition of diaganal block matirx
      do 1 i = 1,n
        k = LS%A%iau(i)

        do ii_mat = 1,ndim
          this%P(ii_mat,i) = ii_mat
        enddo

        do 11 ii_mat = 1,ndim-1
          ii = this%P(ii_mat,i)

          iip = ii
          tmp = abs( LS%A%e(ii,ii_mat,k) )
          do iis_mat = ii_mat+1,ndim
            iis = this%P(iis_mat,i)

            if( tmp < abs( LS%A%e(iis,ii_mat,k) ) )  then
              iip = iis
              iip_mat = iis_mat
              tmp = abs( LS%A%e(iis,ii_mat,k) )
            endif
          enddo

          if( iip /= ii )  then
            this%P(ii_mat,i) = iip
            this%P(iip_mat,i) = ii
          endif

          iip = this%P(ii_mat,i)
          pvt = LS%A%e(iip,ii_mat,k)
          do 111 iis_mat = ii_mat+1,ndim
            iis = this%P(iis_mat,i)
            tmp = LS%A%e(iis,ii_mat,k)/pvt
            LS%A%e(iis,ii_mat,k) = tmp
            LS%A%e(iis,ii_mat+1:ndim,k) = LS%A%e(iis,ii_mat+1:ndim,k) - tmp*LS%A%e(iip,ii_mat+1:ndim,k)
111       enddo
11      enddo
1     enddo

!     initialize for error calculation
      this%Uxold(:,:) = 0.

!     this
!     calculated error is the error of previous iteration
      iter = 0
      edrop = 1d20
      do 2 while ( iter == 1 .or. ( edrop > this%etol .and. iter < this%niter ) )
        iter = iter + 1

!       forward block Gauss-Seidel
        err = 0.
        do 21 i = 1,n
          kau = LS%A%iau(i)

!         for error calculation
          this%Uxolder(:) = this%Uxold(:,i)
          this%Uxold(:,i) = 0.

          this%y(:) = LS%rhs(:,i)
!         LS%A
          do k = LS%A%ia(i),LS%A%ia(i+1)-1
            if( k == kau ) cycle
            j = LS%A%ja(k)
            this%Ax(:) = matmul(LS%A%e(:,:,k),LS%x(:,j))
            this%y(:) = this%y(:) - this%Ax(:)
            if( j > i ) this%Uxold(:,i) = this%Uxold(:,i) + this%Ax(:)  ! for error calculation
          enddo

!         pivoted-LU solver inside [ndim*ndim] block
!         forward substution
          do ii_mat = 1,ndim
            ii = this%P(ii_mat,i)
            this%z(ii_mat) = this%y(ii)
            do iis_mat = 1,ii_mat-1
              this%z(ii_mat) = this%z(ii_mat) - LS%A%e(ii,iis_mat,kau)*this%z(iis_mat)
            enddo
          enddo
!         backward substution
          do ii_mat = ndim,1,-1
            ii = this%P(ii_mat,i)
            LS%x(ii_mat,i) = this%z(ii_mat)
            do iis_mat = ii_mat+1,ndim
              LS%x(ii_mat,i) = LS%x(ii_mat,i) - LS%A%e(ii,iis_mat,kau)*LS%x(iis_mat,i)
            enddo
            LS%x(ii_mat,i) = LS%x(ii_mat,i)/LS%A%e(ii,ii_mat,kau)
          enddo

!         acculumate error
          do ii_mat = 1,ndim
            err = err + abs( this%Uxold(ii_mat,i) - this%Uxolder(ii_mat) )
          enddo
21      enddo

        err = err/(n*ndim)
        if( iter == 2 ) then
          err0 = err
          write(*,*) 'err0:', err0
        end if

        if( iter == 2  .or. mod(iter,100) == 0 ) then
          edrop = err/err0
          write(*,*) iter, edrop, err
        end if
2     enddo
      write(*,*) iter, edrop, err

end subroutine LUSGS

end module LinearSystem_Solver_LUSGS_mod
