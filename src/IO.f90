module IO_mod
      use KindDefinition_mod, only : DP
      implicit none
      save

contains

!=====================================================================
!
!=====================================================================
subroutine Output_vector( x, n, ndim, fp )
      implicit none
      integer, intent(in) :: n, ndim
      real(DP), intent(in) :: x(ndim,n)
      integer, intent(in) :: fp
      integer :: i, ii
      real(DP), parameter :: TOL = 1d-15

      write(fp,*)  'Output_vector: Start'

      do i = 1,n
        do ii = 1,ndim
          if( abs(x(ii,i)) < TOL ) then
            write(fp,20,advance='no') '         *     '
          else
            write(fp,10,advance='no') x(ii,i)
          end if
        end do
        write(fp,30,advance='no') '|'
      end do
      write(fp,*)

10    format(ES16.6)
20    format(a16)
30    format(a2)

      write(fp,*)  'Output_vector: End'

end subroutine Output_vector

end module IO_mod
