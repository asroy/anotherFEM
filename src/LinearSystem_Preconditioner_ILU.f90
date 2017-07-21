module LinearSystem_Preconditioner_ILU_mod
      use KindDefinition_mod, only : DP
      use LinearSystem_Generic_mod
      use IO_mod
      use GlobalVariableForDebug_mod
      use Common_Type_CRS_mod
      use Common_Type_Dependency_mod
      use GMRES_mod
      implicit none
      save

      type, extends(preconditioner_t) :: ILU_t
        integer :: fillLevel
        type(CRS_t) :: ALU
        integer, allocatable :: iw(:)
        real(DP), allocatable :: tmat(:,:)
      contains
        procedure, pass(this) :: Reallocate_preconditioner => Reallocate_ILU
        procedure, pass(this) :: Calculate_preconditioner => Calculate_ILU
        procedure, pass(this) :: Apply_preconditioner => Apply_ILU
      end type ILU_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reallocate_ILU ( this, LS )
      implicit none
      class(ILU_t),            intent(out) :: this
      class(linear_system_t),  intent(in)  :: LS

      integer :: nnz_ilu
      integer, dimension(:), allocatable :: ia_ilu, iau_ilu, ja_ilu

      call this%ALU%Deallocate_CRS
      if( allocated( this%iw ) )    deallocate ( this%iw )
      if( allocated( this%tmat ) )  deallocate ( this%tmat )

      this%ndim = LS%ndim
      this%n    = LS%n

      this % fillLevel = 1

      call getIAJAFill(this%n, LS%A%nnz, LS%A%ia, LS%A%iau, LS%A%ja, nnz_ilu, ia_ilu, iau_ilu, ja_ilu, this%fillLevel)


      call this%ALU%Reallocate_CRS ( this%ndim, this%n, nnz_ilu, ia_ilu, ja_ilu )

      deallocate( ia_ilu, iau_ilu, ja_ilu )

     !call LS   % A   % output_CRS_tecplot( fpIO, 'output/A_depen.plt' )
     !call this % ALU % output_CRS_tecplot( fpIO, 'output/ALU_depen.plt' )

end subroutine Reallocate_ILU
!=====================================================================
!
!=====================================================================
subroutine Calculate_ILU ( this, LS )
      implicit none
      class(ILU_t), intent(inout) :: this
      class(linear_system_t),  intent(in)    :: LS

      integer :: n, ndim
      integer :: i, j, k, ii, jj, jstart, jend, kstart, kend, jcol, jcolnew, kcol

      integer :: icode

      n = this%n
      ndim = this%ndim

      ! copy the A into the ALU
      call this%ALU%Set_CRS_to_0

      do i = 1, n
        jstart = LS%A%ia(i)
        jend   = LS%A%ia(i+1) - 1
        kstart = this%ALU%ia(i)
        kend   = this%ALU%ia(i+1) - 1

        do j = jstart, jend
          jcol = LS%A%ja(j)
          jcolnew = jcol
          do k = kstart, kend ! New row
            kcol = this%ALU%ja(k) ! New column
            if ( kcol == jcolnew ) then  ! If new column matches the old column, put this value of A into Afill
              do ii = 1, ndim
              do jj = 1, ndim
                this%ALU%e(ii,jj,k) = LS%A%e(ii,jj,j)
              end do
              end do

              exit

            end if
          end do ! k loop
        end do ! jloop
      end do ! loop over rows to fill ALU


      ! calculate ALU
      call GBLKILU(this%n,this%ALU%nnz,this%ALU%ja,this%ALU%ia,this%ALU%e,this%ALU%iau,icode,this%ndim)

      if ( icode /= 0 ) then
        write(6,333)
333     format(1h ,'Error in BLKILU')
        stop
      end if

      ! for debug
     !call this%ALU%Set_CRS_to_unity

end subroutine Calculate_ILU
!=====================================================================
!
!=====================================================================
subroutine Apply_ILU ( this, x, y )
      implicit none
      class(ILU_t),                  intent(in)  :: this
      real(DP),      dimension(:,:), intent(out) :: x
      real(DP),      dimension(:),   intent(in)  :: y

      call GBLKSOL(this%n,this%ALU%nnz, y, x, this%ALU%e, this%ALU%ja, this%ALU%ia, this%ALU%iau, this%ndim)

end subroutine Apply_ILU

end module LinearSystem_Preconditioner_ILU_mod
