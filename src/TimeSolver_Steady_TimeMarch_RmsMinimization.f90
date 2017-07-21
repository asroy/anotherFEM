module TimeSolver_Steady_TimeMarch_RmsMinimization_mod
      use KindDefinition_mod, only : DP
      use Common_Function_mod
      use IO_mod
      use TimeSolver_Generic_mod
      use GlobalVariableForDebug_mod
      save

      type, extends(time_solver_t) :: steady_time_march_rms_minimization_t
        logical :: Use_dt_local
        real(DP) :: dt_global

        real(DP) :: cfl_min
        real(DP) :: cfl_max

        integer :: iramp

        integer :: niteration
        real(DP) :: error
      contains
        procedure, pass(this) :: Reset_time_solver => Reset_steady_time_march_rms_minimization
        procedure, pass(this) :: Do_time_solver => Do_steady_time_march_rms_minimization
      end type steady_time_march_rms_minimization_t

contains

!=====================================================================
!
!=====================================================================
subroutine Reset_steady_time_march_rms_minimization ( this )
      implicit none
      class(steady_time_march_rms_minimization_t), intent(out) :: this

      this % Use_dt_local = .true.
      this % dt_global = 0.

      this % cfl_min = 1.
      this % cfl_max = 1.e6
      this % iramp = 100

      ! niteration, error
      this % niteration = 1000
      this % error = 1.e-14

end subroutine Reset_steady_time_march_rms_minimization
!=====================================================================
!
!=====================================================================
subroutine Do_steady_time_march_rms_minimization ( this, FEM, equation, linear_system, linear_solver, preconditioner )
      implicit none
      class(steady_time_march_rms_minimization_t),  intent(in)    :: this
      class(FEM_t),                                 intent(inout) :: FEM
      class(equation_t),                            intent(inout) :: equation
      class(linear_system_t),                       intent(inout) :: linear_system
      class(linear_solver_t),                       intent(inout) :: linear_solver
      class(preconditioner_t),                      intent(inout) :: preconditioner

      type(time_control_t) :: time_control
      real(DP), dimension(nequation_MAX) :: rhs_rms
      integer :: ndim, n
      integer :: iter
      real(DP), dimension(4) :: alf
      real(DP), dimension(0:3) :: omega, rms
      real(DP) :: omega_opt, rms_opt, omega_chosen, rms_chosen, omega_max
      real(DP) :: tmp, oa, ob, rmsa, rmsb
      real(DP) :: Amat(4,4), bvec(4)
      real(DP), parameter :: tol = 1.e-14

      ndim = linear_system % ndim
      n    = linear_system % n

      time_control % Include_time = .true.
      time_control % ntimeStep_time_low = 0
      time_control % ntimeStep_time_up  = 1
      time_control % alph(0:1) = [ -1., 1. ]
      time_control % Use_dt_local = .true.
      time_control % dt_global = 0.
      time_control % Include_space = .true.
      time_control % itimeStep_space = 1
      time_control % ratio_space = 1.


      do iter = 1, this % niteration
        if( iter == 1 ) then
          time_control % cfl = this % cfl_min
        end if

        call equation % Enforce_essential_boundary_condition ( FEM, 1 )

        call equation % Calculate_dt_local   ( FEM, time_control )
       !call equation % Calculate_dt_local_2 ( FEM, time_control )

        call FEM % Calculate_FEM_residual( linear_system % rhs, equation, time_control )
        linear_system % rhs(:,:) = - linear_system % rhs(:,:)

        rhs_rms = vector_rms ( linear_system % rhs, ndim, n )/n
        write(*,   200)  'Do_steady_time_march_rms_minimization:', iter, time_control % cfl, rhs_rms(1:ndim)
        write(fpd1,200)  'Do_steady_time_march_rms_minimization:', iter, time_control % cfl, rhs_rms(1:ndim)

        if( rhs_rms(1) < this % error ) exit

        call FEM % Calculate_FEM_jacobian( linear_system % A, equation, time_control )

        call linear_solver % Solve_linear_system( linear_system, preconditioner )


        ! line search RMS minimization
        ! maximn omega allowed
        omega_max = equation % max_omega_for_dq ( FEM, linear_system % x, 1 )

        ! omege0
        omega(0) = 0.
        rms(0) = sqrt( sum( rhs_rms(1:ndim)**2 ) )

        ! omega1
        omega(1) = 0.1
        call FEM % Update_FEM_solution( linear_system % x, omega(1), 1 )
        call FEM % Calculate_FEM_residual( linear_system % rhs, equation, time_control )

        rhs_rms = vector_rms ( linear_system % rhs, ndim, n )/n
        rms(1) = sqrt( sum( rhs_rms(1:ndim)**2 ) )

        call FEM % Update_FEM_solution( linear_system % x, -omega(1), 1 )

        if( isnan(rms(1)) ) then
          omega_chosen = 0.
          goto 1000
        end if

        ! omega2
        omega(2) = 0.5*( omega(1) + omega_max )
        call FEM % Update_FEM_solution( linear_system % x, omega(2), 1 )
        call FEM % Calculate_FEM_residual( linear_system % rhs, equation, time_control )

        rhs_rms = vector_rms ( linear_system % rhs, ndim, n )/n
        rms(2) = sqrt( sum( rhs_rms(1:ndim)**2 ) )

        call FEM % Update_FEM_solution( linear_system % x, -omega(2), 1 )

        if( isnan(rms(2)) ) then
          omega_chosen = 0.
          goto 1000
        end if

        ! omega3
        omega(3) = omega_max
        call FEM % Update_FEM_solution( linear_system % x, omega(3), 1 )
        call FEM % Calculate_FEM_residual( linear_system % rhs, equation, time_control )

        rhs_rms = vector_rms ( linear_system % rhs, ndim, n )/n
        rms(3) = sqrt( sum( rhs_rms(1:ndim)**2 ) )

        call FEM % Update_FEM_solution( linear_system % x, -omega(3), 1 )

        if( isnan(rms(3)) ) then
          omega_chosen = 0.
          goto 1000
        end if

        ! line search for best omega
        Amat(1,1:4) = [ 1., omega(0), omega(0)**2, omega(0)**3 ]
        Amat(2,1:4) = [ 1., omega(1), omega(1)**2, omega(1)**3 ]
        Amat(3,1:4) = [ 1., omega(2), omega(2)**2, omega(2)**3 ]
        Amat(4,1:4) = [ 1., omega(3), omega(3)**2, omega(3)**3 ]
        
        bvec(1:4) = [ rms(0), rms(1), rms(2), rms(3) ]

        alf(1:4) =  GaussJordan ( Amat, bvec, 4 )

        ! linear
        if( abs( alf(4) ) < tol .and. abs( alf(3) ) < tol ) then
          if( alf(2) >= 0. ) then
            omega_opt = omega(0)
          else
            omega_opt = omega(3)
          end if

        ! quadratic
        else if( abs( alf(4) ) < tol ) then
          omega_opt = - 0.5*alf(2)/alf(3)

        ! cubic
        else
          tmp = 4.*alf(3)**2 - 12.*alf(4)*alf(2)

          ! no real root, monotonic
          if( tmp <= 0. ) then
            if( rms(0) < rms(3) ) then
              omega_opt = omega(0)
            else
              omega_opt = omega(3)
            end if
          else
            oa = ( - 2.*alf(3) + sqrt(tmp) )/( 6.*alf(4)) 
            ob = ( - 2.*alf(3) - sqrt(tmp) )/( 6.*alf(4)) 

            rmsa = alf(1) + alf(2)*oa + alf(3)*oa**2 + alf(4)*oa**3
            rmsb = alf(1) + alf(2)*ob + alf(3)*ob**2 + alf(4)*ob**3

            if( rmsa < rmsb ) then
              omega_opt = oa
            else
              omega_opt = ob
            end if
          end if
        end if

        ! choose a omega
        omega_chosen = omega(minloc(rms,1))
        rms_chosen = minval(rms)

        if( omega_opt > omega(0) .and. omega_opt < omega(3) ) then
          rms_opt = alf(1) + alf(2)*omega_opt + alf(3)*omega_opt**2 + alf(4)*omega_opt**3
          
          if( rms_opt < rms_chosen ) then
            omega_chosen = omega_opt
            rms_chosen = rms_opt
          end if
        end if

1000    continue

        ! CFL
        if( omega_chosen < 0.1 ) then
          time_control % cfl =  time_control % cfl/10.
       !else if( omega_chosen > 1. ) then
        else if( omega_chosen > 0.6 ) then
          time_control % cfl = 2.*time_control % cfl
        else
          time_control % cfl = time_control % cfl
        end if

        if( time_control % cfl > this % cfl_max )  time_control % cfl = this % cfl_max

        write(*,100) 'omega', omega(0:3), omega_opt, omega_chosen
        write(*,100) 'rms  ', rms  (0:3),   rms_opt,   rms_chosen

        ! >>> debug
       !write(*,*) 'debug! Do_steady_time_march_rms_minimization'
       !if (iter < this % iramp ) then
       !  time_control % cfl = this % cfl_min + ( iter - 1 )*( this % cfl_max - this % cfl_min)/( this % iramp - 1 )
       !else
       !  time_control % cfl = this % cfl_max
       !end if
       !omega_chosen = 1.
        ! <<< debug


        call FEM % Update_FEM_solution( linear_system % x, omega_chosen, 1 )

        call FEM % Output_FEM_solution_tecplot( fpIO,'output/Do_steady_time_march_rms_mininization-q.plt' )

        call equation % Check_solution ( FEM, 1, fpIO, 'output/check_solution.plt' )

        call FEM % Update_FEM_solution_history ( 0, 1 )
      end do

100   format(a,20(es16.5))
200   format(a,i7,20(es16.5))

end subroutine Do_steady_time_march_rms_minimization

end module TimeSolver_Steady_TimeMarch_RmsMinimization_mod
