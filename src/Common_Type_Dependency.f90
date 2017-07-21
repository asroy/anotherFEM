module Common_Type_Dependency_mod
      use KindDefinition_mod, only : DP
      implicit none
      save

      type :: depen_t
        integer :: n
        integer :: nnz
        integer, allocatable :: ia(:)
        integer, allocatable :: ja(:)
      contains
        procedure, pass(this) :: Reallocate_depen
        procedure, pass(this) :: Deallocate_depen
        procedure, pass(this) :: Output_depen
        procedure, pass(this) :: Output_depen_tecplot
      end type depen_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reallocate_depen( this, n, nnz )
      implicit none
      class(depen_t), intent(inout) :: this
      integer,        intent(in)    :: n
      integer,        intent(in)    :: nnz

      if( allocated(this%ia) ) then
        deallocate( this%ia )
        deallocate( this%ja )
      end if

      this%n = n
      this%nnz = nnz

      allocate (this%ia(n+1) )
      allocate (this%ja(nnz) )

end subroutine Reallocate_depen
!=====================================================================
!
!=====================================================================
subroutine Deallocate_depen( this )
      implicit none
      class(depen_t), intent(inout) :: this

      this%n = 0
      this%nnz = 0

      if( allocated(this%ia) ) then
        deallocate( this%ia )
        deallocate( this%ja )
      end if

end subroutine Deallocate_depen
!=====================================================================
!
!=====================================================================
subroutine Output_depen ( this, fp )
      implicit none
      class(depen_t), intent(in) :: this
      integer,        intent(in) :: fp

      integer :: i, j, k
      logical :: Is_found

      write(fp,*) 'Output_depen: Start'

      do 1 i = 1,this%n
        write(fp,10,advance='no') i,this%ia(i+1)-this%ia(i),':'
        do 11 j = 1,this%n
          do k = this%ia(i),this%ia(i+1)-1
            if( this%ja(k) == j ) then
              write(fp,20,advance='no') j
              exit
            end if
          end do
11      end do
        write(fp,*)
1     end do

      write(fp,*)

      do 2 i = 1,this%n
        write(fp,10,advance='no') i,this%ia(i+1)-this%ia(i),':'
        do 21 j = 1,this%n
          Is_found = .false.
          do k = this%ia(i),this%ia(i+1)-1
            if( this%ja(k) == j ) then
              Is_found = .true.
              write(fp,30,advance='no') '*'
              exit
            end if
          end do
          if( .not. Is_found ) write(fp,30,advance='no') ' '
21      end do
        write(fp,*)
2     end do

      write(fp,*) 'Output_depen: End'
      write(fp,*)

10    format(i6,i3,a1)
20    format(i6)
30    format(a2)

end subroutine Output_depen
!=====================================================================
!
!=====================================================================
subroutine Output_depen_tecplot ( this, fp, prefix )
      implicit none
      class(depen_t),   intent(in) :: this
      integer,          intent(in) :: fp
      character(len=*), intent(in) :: prefix

      character(len=100) :: filename
      integer :: i, j, k

      write(filename, '(a)') trim(prefix)
      open(fp, file=filename, form='formatted', status='replace')

      write(fp,*) 'variables ="j","-i"'
      write(fp,'(a,a)') 'zone T=',trim(prefix)

      do i = 1, this%n
        do k = this%ia(i), this%ia(i+1)-1
          j = this%ja(k)
          write(fp,*) j, -i
        end do
      end do

      close(fp)

end subroutine Output_depen_tecplot

end module Common_Type_Dependency_mod
