module Common_Function_mod
      use KindDefinition_mod, only : DP
      use Common_overload_mod
      implicit none
      save

      private

      public :: vector_rms
      public :: dot_product_my, matrix_vector_multiply
      public :: GaussJordan

      interface dot_product_my
        module procedure dot_product_real_real
        module procedure dot_product_dtype_real
        module procedure dot_product_real_dtype
        module procedure dot_product_dtype_dtype
      end interface

      interface matrix_vector_multiply
        module procedure matrix_vector_multiply_real_real
        module procedure matrix_vector_multiply_dtype_real
        module procedure matrix_vector_multiply_real_dtype
        module procedure matrix_vector_multiply_dtype_dtype
      end interface

      interface GaussJordan
        module procedure GaussJordan_real
        module procedure GaussJordan_dtype
      end interface

contains

!=====================================================================
!
!=====================================================================
function vector_rms ( u, ndim, n )  result ( rms )
      implicit none
      integer,                  intent(in)  :: ndim, n
      real(DP), dimension(:,:), intent(in)  :: u
      real(DP), dimension(ndim)             :: rms

      integer :: i, ii

      rms(1:ndim) = 0.
      do i = 1, n
      do ii = 1,ndim
        rms(ii) = rms(ii) + u(ii,i)**2
      end do
      end do

      rms(1:ndim) = sqrt ( rms(1:ndim) )

end function vector_rms
!=======================================================================
!
!=======================================================================
function dot_product_real_real ( a, b, n )  result ( f )
      implicit none
      integer,                   intent(in) :: n
      real(DP),    dimension(:), intent(in) :: a
      real(DP),    dimension(:), intent(in) :: b
      real(DP)                              :: f

      integer :: i

      f = 0.
      do i = 1, n
        f = f + a(i)*b(i)
      end do

end function dot_product_real_real
!=======================================================================
!
!=======================================================================
function dot_product_dtype_real ( a, b, n )  result ( f )
      implicit none
      integer,                   intent(in) :: n
      type(dType), dimension(:), intent(in) :: a
      real(DP),    dimension(:), intent(in) :: b
      type(dType)                           :: f

      integer :: i

      f = 0.
      do i = 1, n
        f = f + a(i)*b(i)
      end do

end function dot_product_dtype_real
!=======================================================================
!
!=======================================================================
function dot_product_real_dtype ( a, b, n )  result ( f )
      implicit none
      integer,                   intent(in) :: n
      real,        dimension(:), intent(in) :: a
      type(dType), dimension(:), intent(in) :: b
      type(dType)                           :: f

      integer :: i

      f = 0.
      do i = 1, n
        f = f + a(i)*b(i)
      end do

end function dot_product_real_dtype
!=======================================================================
!
!=======================================================================
function dot_product_dtype_dtype ( a, b, n )  result ( f )
      implicit none
      integer,                   intent(in) :: n
      type(dType), dimension(:), intent(in) :: a
      type(dType), dimension(:), intent(in) :: b
      type(dType)                           :: f

      integer :: i

      f = 0.
      do i = 1, n
        f = f + a(i)*b(i)
      end do

end function dot_product_dtype_dtype
!=======================================================================
!
!=======================================================================
function matrix_vector_multiply_real_real ( A, b, m, n )  result ( v )
      implicit none
      integer,                     intent(in) :: m, n
      real(DP),    dimension(:,:), intent(in) :: A
      real(DP),    dimension(:),   intent(in) :: b
      real(DP),    dimension(m)               :: v

      integer :: i, j

      v = 0.
      do i = 1, m
      do j = 1, n
        v(i) = v(i) + A(i,j)*b(j)
      end do
      end do

end function matrix_vector_multiply_real_real
!=======================================================================
!
!=======================================================================
function matrix_vector_multiply_dtype_real ( A, b, m, n )  result ( v )
      implicit none
      integer,                     intent(in) :: m, n
      type(dType), dimension(:,:), intent(in) :: A
      real(DP),    dimension(:),   intent(in) :: b
      type(dType), dimension(m)               :: v

      integer :: i, j

      v = 0.
      do i = 1, m
      do j = 1, n
        v(i) = v(i) + A(i,j)*b(j)
      end do
      end do

end function matrix_vector_multiply_dtype_real
!=======================================================================
!
!=======================================================================
function matrix_vector_multiply_real_dtype ( A, b, m, n )  result ( v )
      implicit none
      integer,                     intent(in) :: m, n
      real(DP),    dimension(:,:), intent(in) :: A
      type(dType), dimension(:),   intent(in) :: b
      type(dType), dimension(m)               :: v

      integer :: i, j

      v = 0.
      do i = 1, m
      do j = 1, n
        v(i) = v(i) + A(i,j)*b(j)
      end do
      end do

end function matrix_vector_multiply_real_dtype
!=======================================================================
!
!=======================================================================
function matrix_vector_multiply_dtype_dtype ( A, b, m, n )  result ( v )
      implicit none
      integer,                     intent(in) :: m, n
      type(dType), dimension(:,:), intent(in) :: A
      type(dType), dimension(:),   intent(in) :: b
      type(dType), dimension(m)               :: v

      integer :: i, j

      v = 0.
      do i = 1, m
      do j = 1, n
        v(i) = v(i) + A(i,j)*b(j)
      end do
      end do

end function matrix_vector_multiply_dtype_dtype
!=======================================================================
!
! Solver linear equations using Gaussian elimination with row pivoting
!   A - LHS (will be modified)
!   b - RHS (will be modified)
!   x - solution
!
!=======================================================================
function GaussJordan_real ( A, b, n )  result( x )
      integer,                  intent(in)    :: n
      real(DP), dimension(:,:), intent(inout) :: A
      real(DP), dimension(:),   intent(inout) :: b
      real(DP), dimension(n)                  :: x

      real(DP), parameter :: epsilon = 1.e-13
      integer :: i,j,k
      integer :: error
      integer :: ipeak
      real(DP) :: factor
      real(DP) :: temp

      do i = 1,n
!
! Find pivot row
!
        ipeak = i
        do j = i+1,n
          if(abs(A(j,i)).gt.abs(A(ipeak,i)))then
            ipeak = j
          end if
        end do
!
! Check for singular
!
        if(abs(A(ipeak,i)).lt.epsilon)then
          error = 1
          write(*,'("Singular matrix in gaussJordan: stopping")')
          stop
        end if
        if(ipeak.ne.i)then
          do k = 1,n
            temp = A(ipeak,k)
            A(ipeak,k) = A(i,k)
            A(i,k) = temp
          end do
          temp = b(ipeak)
          b(ipeak) = b(i)
          b(i) = temp
        end if
!
! Zero out appropriate columns
!
        do j = 1,n
          if(j.ne.i)then
            factor = -A(j,i)/A(i,i)
            do k = 1,n
              A(j,k) = A(i,k)*factor + A(j,k)
            end do
            b(j) = b(i)*factor + b(j)
          end if
        end do
      end do
!
! Now obtain the solution by simply dividing b(:) by the diagonal
!
      do i = 1,n
        x(i) = b(i)/A(i,i)
      end do
      return

end function GaussJordan_real
!=======================================================================
!
! Solver linear equations using Gaussian elimination with row pivoting
!   A - LHS (will be modified)
!   b - RHS (will be modified)
!   x - solution
!
!=======================================================================
function GaussJordan_dtype ( A, b, n )  result( x )
      integer,                     intent(in)    :: n
      type(dType), dimension(:,:), intent(inout) :: A
      type(dType), dimension(:),   intent(inout) :: b
      type(dType), dimension(n)                  :: x

      real(DP), parameter :: epsilon = 1.e-13
      integer :: i,j,k
      integer :: error
      integer :: ipeak
      type(dType) :: factor
      type(dType) :: temp

      do i = 1,n
!
! Find pivot row
!
        ipeak = i
        do j = i+1,n
          if(abs(A(j,i)).gt.abs(A(ipeak,i)))then
            ipeak = j
          end if
        end do
!
! Check for singular
!
        if(abs(A(ipeak,i)).lt.epsilon)then
          error = 1
          write(*,'("Singular matrix in gaussJordan: stopping")')
          stop
        end if
        if(ipeak.ne.i)then
          do k = 1,n
            temp = A(ipeak,k)
            A(ipeak,k) = A(i,k)
            A(i,k) = temp
          end do
          temp = b(ipeak)
          b(ipeak) = b(i)
          b(i) = temp
        end if
!
! Zero out appropriate columns
!
        do j = 1,n
          if(j.ne.i)then
            factor = -A(j,i)/A(i,i)
            do k = 1,n
              A(j,k) = A(i,k)*factor + A(j,k)
            end do
            b(j) = b(i)*factor + b(j)
          end if
        end do
      end do
!
! Now obtain the solution by simply dividing b(:) by the diagonal
!
      do i = 1,n
        x(i) = b(i)/A(i,i)
      end do
      return

end function GaussJordan_dtype

end module Common_Function_mod
