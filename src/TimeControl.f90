module TimeControl_mod
      use KindDefinition_mod, only : DP
      use Constant_mod
      use Common_Overload_mod
      save

! - time stage
      type :: time_control_t
        ! time residual
        logical :: Include_time
        integer :: ntimeStep_time_low
        integer :: ntimeStep_time_up
        real(DP) :: alph (ntimeStep_LOWEST:ntimeStep_UPPEST)    ! coeffcient for time discretization

        logical :: Use_dt_local
        real(DP) :: cfl                                         ! CFL number
        real(DP) :: dt_global                                   ! global dt for unsteady problem

      ! space residual
        logical :: Include_space
        integer :: itimeStep_space                              ! time step to be calculated for space residual/jacobian
        real(DP) :: ratio_space                                 ! ratio space residual to be added to residual
      contains
        generic, public :: time_discretization => time_discretization_real, time_discretization_dtype
        procedure, pass(this), private :: time_discretization_real
        procedure, pass(this), private :: time_discretization_dtype
      end type time_control_t

contains

!=====================================================================
!
!=====================================================================
function time_discretization_real ( this, q_ntime, dt, nq )  result ( dq_t )
      implicit none
      class(time_control_t),                                 intent(in) :: this
      integer,                                               intent(in) :: nq
      real(DP),              dimension(:,ntimeStep_LOWEST:), intent(in) :: q_ntime
      real(DP),                                              intent(in) :: dt
      real(DP)                                                          :: dq_t(nq)

      integer :: itime

      dq_t(1:nq) = 0.

      do itime = this % ntimeStep_time_low, this % ntimeStep_time_up
        dq_t(1:nq) = dq_t(1:nq) + this % alph(itime) * q_ntime(1:nq,itime)
      end do

      dq_t(1:nq) = dq_t(1:nq)/dt

end function time_discretization_real
!=====================================================================
!
!=====================================================================
function time_discretization_dtype ( this, q_ntime, dt, nq )  result ( dq_t )
      implicit none
      class(time_control_t),                                 intent(in) :: this
      integer,                                               intent(in) :: nq
      type(dType),           dimension(:,ntimeStep_LOWEST:), intent(in) :: q_ntime
      real(DP),                                              intent(in) :: dt
      type(dType)                                                       :: dq_t(nq)

      integer :: itime

      dq_t(1:nq) = 0.

      do itime = this % ntimeStep_time_low, this % ntimeStep_time_up
        dq_t(1:nq) = dq_t(1:nq) + this % alph(itime) * q_ntime(1:nq,itime)
      end do

      dq_t(1:nq) = dq_t(1:nq)/dt

end function time_discretization_dtype

end module TimeControl_mod
