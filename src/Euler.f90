module Euler_mod
      use KindDefinition_mod, only : DP
      use Common_Overload_mod
      use Common_Function_mod
      use GlobalVariableForDebug_mod
      implicit none

      private

      public :: gama, gm1, gm2, gm3, gi, gm1i, ggm1i, gigm1i
      public :: conservative_variable_cq
      public :: cq_derivative_wrt_q
      public :: convective_flux_f, convective_flux_g, convective_flux_face
      public :: convective_jacobian_A, convective_jacobian_B, convective_jacobian_face
      public :: get_tau_inverse
      public :: roe_flux_face
      public :: shock_capture_flux_f_g
      public :: shock_capture_flux_fx_gy
      public :: shock_capture_flux_f_g_modified
      public :: shock_capture_flux_fx_gy_modified

      save

      real(DP), parameter :: gama = 1.4
      real(DP), parameter :: gm1 = gama - 1.
      real(DP), parameter :: gm2 = gama - 2.
      real(DP), parameter :: gm3 = gama - 3.
      real(DP), parameter :: gi = 1./gama
      real(DP), parameter :: gm1i = 1./( gama - 1. )
      real(DP), parameter :: ggm1i = gama/( gama - 1. )
      real(DP), parameter :: gigm1i = 1./gama/( gama - 1. )

      interface conservative_variable_cq
        module procedure conservative_variable_cq_real
        module procedure conservative_variable_cq_dtype
      end interface

      interface cq_derivative_wrt_q
        module procedure cq_derivative_wrt_q_real
        module procedure cq_derivative_wrt_q_dtype
      end interface

      interface convective_flux_face
        module procedure convective_flux_face_real
        module procedure convective_flux_face_dtype
      end interface

      interface convective_flux_f
        module procedure convective_flux_f_real
        module procedure convective_flux_f_dtype
      end interface

      interface convective_flux_g
        module procedure convective_flux_g_real
        module procedure convective_flux_g_dtype
      end interface

      interface convective_jacobian_face
        module procedure convective_jacobian_face_real
        module procedure convective_jacobian_face_dtype
      end interface

      interface convective_jacobian_A
        module procedure convective_jacobian_A_real
        module procedure convective_jacobian_A_dtype
      end interface

      interface convective_jacobian_B
        module procedure convective_jacobian_B_real
        module procedure convective_jacobian_B_dtype
      end interface

      interface get_tau_inverse
        module procedure get_tau_inverse_real
        module procedure get_tau_inverse_dtype
      end interface

      interface roe_flux_face
        module procedure roe_flux_face_real
        module procedure roe_flux_face_dtype
      end interface

      interface shock_capture_flux_f_g
        module procedure shock_capture_flux_f_g_real
        module procedure shock_capture_flux_f_g_dtype
      end interface

      interface shock_capture_flux_fx_gy
        module procedure shock_capture_flux_fx_gy_real
        module procedure shock_capture_flux_fx_gy_dtype
      end interface

      interface shock_capture_flux_f_g_modified
        module procedure shock_capture_flux_f_g_modified_real
        module procedure shock_capture_flux_f_g_modified_dtype
      end interface

      interface shock_capture_flux_fx_gy_modified
        module procedure shock_capture_flux_fx_gy_modified_real
        module procedure shock_capture_flux_fx_gy_modified_dtype
      end interface

contains
!=============================================================================
!
!=============================================================================
function conservative_variable_cq_real ( q )  result(cq)
      implicit none
      real(DP), dimension(4), intent(in)  :: q
      real(DP), dimension(4)              :: cq

      cq(1) = q(1)
      cq(2) = q(1)*q(2)
      cq(3) = q(1)*q(3)
      cq(4) = q(1)*( gigm1i*q(4) + 0.5*( q(2)**2 + q(3)**2 ) )

end function conservative_variable_cq_real
!=============================================================================
!
!=============================================================================
function conservative_variable_cq_dtype ( q )  result(cq)
      implicit none
      type(dType), dimension(4), intent(in)  :: q
      type(dType), dimension(4)              :: cq

      cq(1) = q(1)
      cq(2) = q(1)*q(2)
      cq(3) = q(1)*q(3)
      cq(4) = q(1)*( gigm1i*q(4) + 0.5*( q(2)**2 + q(3)**2 ) )

end function conservative_variable_cq_dtype
!=============================================================================
!
!=============================================================================
function cq_derivative_wrt_q_real ( q )  result(M)
      implicit none
      real(DP), dimension(4), intent(in)  :: q
      real(DP), dimension(4,4)            :: M

      M(1,1) = 1.                                      ;   M(1,2) = 0.        ;   M(1,3) = 0.        ;   M(1,4) = 0.
      M(2,1) = q(2)                                    ;   M(2,2) = q(1)      ;   M(2,3) = 0.        ;   M(2,4) = 0.
      M(3,1) = q(3)                                    ;   M(3,2) = 0.        ;   M(3,3) = q(1)      ;   M(3,4) = 0.
      M(4,1) = gigm1i*q(4) + 0.5*( q(2)**2 + q(3)**2 ) ;   M(4,2) = q(1)*q(2) ;   M(4,3) = q(1)*q(3) ;   M(4,4) = gigm1i*q(1)

end function cq_derivative_wrt_q_real
!=============================================================================
!
!=============================================================================
function cq_derivative_wrt_q_dtype ( q )  result(M)
      implicit none
      type(dType), dimension(4), intent(in)  :: q
      type(dType), dimension(4,4)            :: M

      M(1,1) = 1.                                      ;   M(1,2) = 0.        ;   M(1,3) = 0.        ;   M(1,4) = 0.
      M(2,1) = q(2)                                    ;   M(2,2) = q(1)      ;   M(2,3) = 0.        ;   M(2,4) = 0.
      M(3,1) = q(3)                                    ;   M(3,2) = 0.        ;   M(3,3) = q(1)      ;   M(3,4) = 0.
      M(4,1) = gigm1i*q(4) + 0.5*( q(2)**2 + q(3)**2 ) ;   M(4,2) = q(1)*q(2) ;   M(4,3) = q(1)*q(3) ;   M(4,4) = gigm1i*q(1)

end function cq_derivative_wrt_q_dtype
!=============================================================================
!
!=============================================================================
function convective_flux_face_real ( q, Sx, Sy )  result(flux)
      implicit none
      real(DP), dimension(4), intent(in)  :: q
      real(DP),               intent(in)  :: Sx, Sy
      real(DP), dimension(4)              :: flux

      real(DP) :: p, h

      p = gi*q(1)*q(4)
      h = gm1i*q(4) + 0.5*( q(2)**2 + q(3)**2 )

      flux(1) = q(2)*Sx + q(3)*Sy
      flux(2) = flux(1)*q(2) + p*Sx
      flux(3) = flux(1)*q(3) + p*Sy
      flux(4) = flux(1)*h

end function convective_flux_face_real
!=============================================================================
!
!=============================================================================
function convective_flux_face_dtype ( q, Sx, Sy )  result(flux)
      implicit none
      type(dType), dimension(4), intent(in)  :: q
      real(DP),                  intent(in)  :: Sx, Sy
      type(dType), dimension(4)              :: flux

      type(dType) :: p, h

      p = gi*q(1)*q(4)
      h = gm1i*q(4) + 0.5*( q(2)**2 + q(3)**2 )

      flux(1) = q(2)*Sx + q(3)*Sy
      flux(2) = flux(1)*q(2) + p*Sx
      flux(3) = flux(1)*q(3) + p*Sy
      flux(4) = flux(1)*h

end function convective_flux_face_dtype
!=============================================================================
!
!=============================================================================
function convective_flux_f_real ( q )  result(f)
      implicit none
      real(DP), dimension(4), intent(in)  :: q
      real(DP), dimension(4)              :: f

      real(DP) :: p, h

      p = gi*q(1)*q(4)
      h = gm1i*q(4) + 0.5*( q(2)**2 + q(3)**2 )

      f(1) = q(1)*q(2)
      f(2) = f(1)*q(2) + p
      f(3) = f(1)*q(3)
      f(4) = f(1)*h

end function convective_flux_f_real
!=============================================================================
!
!=============================================================================
function convective_flux_f_dtype ( q )  result(f)
      implicit none
      type(dType), dimension(4), intent(in)  :: q
      type(dType), dimension(4)              :: f

      type(dType) :: p, h

      p = gi*q(1)*q(4)
      h = gm1i*q(4) + 0.5*( q(2)**2 + q(3)**2 )

      f(1) = q(1)*q(2)
      f(2) = f(1)*q(2) + p
      f(3) = f(1)*q(3)
      f(4) = f(1)*h

end function convective_flux_f_dtype
!=============================================================================
!
!=============================================================================
function convective_flux_g_real ( q )  result(g)
      implicit none
      real(DP), dimension(4), intent(in) :: q
      real(DP), dimension(4)             :: g

      real(DP) :: p, h

      p = gi*q(1)*q(4)
      h = gm1i*q(4) + 0.5*( q(2)**2 + q(3)**2 )

      g(1) = q(1)*q(3)
      g(2) = g(1)*q(2)
      g(3) = g(1)*q(3) + p
      g(4) = g(1)*h

end function convective_flux_g_real
!=============================================================================
!
!=============================================================================
function convective_flux_g_dtype ( q )  result(g)
      implicit none
      type(dType), dimension(4), intent(in) :: q
      type(dType), dimension(4)             :: g

      type(dType) :: p, h

      p = gi*q(1)*q(4)
      h = gm1i*q(4) + 0.5*( q(2)**2 + q(3)**2 )

      g(1) = q(1)*q(3)
      g(2) = g(1)*q(2)
      g(3) = g(1)*q(3) + p
      g(4) = g(1)*h

end function convective_flux_g_dtype
!=============================================================================
!
! 2D convective flux Jacobian w.r.t. conservative variables, according to
!   J. Blazek's CFD book
!
!=============================================================================
function convective_jacobian_face_real ( q, Sx, Sy )  result(A)
      implicit none
      real(DP), dimension(4), intent(in)  :: q
      real(DP),               intent(in)  :: Sx, Sy
      real(DP), dimension(4,4)            :: A

      real(DP) :: ubar_S, ek, h, phi

      ubar_S = q(2)*Sx + q(3)*Sy
      ek   = 0.5*( q(2)**2 + q(3)**2 )
      h    = gm1i*q(4) + ek
      phi  = gm1*ek

      A(1,1) = 0.                   ;   A(1,2) = Sx                      ;   A(1,3) = Sy                      ;   A(1,4) = 0.
      A(2,1) = Sx*phi - q(2)*ubar_S ;   A(2,2) = ubar_S  - (gm2*Sx)*q(2) ;   A(2,3) = Sy*q(2) - (gm1*Sx)*q(3) ;   A(2,4) = gm1*Sx
      A(3,1) = Sy*phi - q(3)*ubar_S ;   A(3,2) = Sx*q(3) - (gm1*Sy)*q(2) ;   A(3,3) = ubar_S  - (gm2*Sy)*q(3) ;   A(3,4) = gm1*Sy
      A(4,1) = ubar_S*( phi - h )   ;   A(4,2) = Sx*h - gm1*q(2)*ubar_S  ;   A(4,3) = Sy*h - gm1*q(3)*ubar_S  ;   A(4,4) = gama*ubar_S

end function convective_jacobian_face_real
!=============================================================================
!
! 2D convective flux Jacobian w.r.t. conservative variables, according to
!   J. Blazek's CFD book
!
!=============================================================================
function convective_jacobian_face_dtype ( q, Sx, Sy )  result(A)
      implicit none
      type(dType), dimension(4), intent(in)  :: q
      real(DP),                  intent(in)  :: Sx, Sy
      type(dType), dimension(4,4)            :: A

      type(dType) :: ubar_S, ek, h, phi

      ubar_S = q(2)*Sx + q(3)*Sy
      ek   = 0.5*( q(2)**2 + q(3)**2 )
      h    = gm1i*q(4) + ek
      phi  = gm1*ek

      A(1,1) = 0.                   ;   A(1,2) = Sx                      ;   A(1,3) = Sy                      ;   A(1,4) = 0.
      A(2,1) = Sx*phi - q(2)*ubar_S ;   A(2,2) = ubar_S  - (gm2*Sx)*q(2) ;   A(2,3) = Sy*q(2) - (gm1*Sx)*q(3) ;   A(2,4) = gm1*Sx
      A(3,1) = Sy*phi - q(3)*ubar_S ;   A(3,2) = Sx*q(3) - (gm1*Sy)*q(2) ;   A(3,3) = ubar_S  - (gm2*Sy)*q(3) ;   A(3,4) = gm1*Sy
      A(4,1) = ubar_S*( phi - h )   ;   A(4,2) = Sx*h - gm1*q(2)*ubar_S  ;   A(4,3) = Sy*h - gm1*q(3)*ubar_S  ;   A(4,4) = gama*ubar_S

end function convective_jacobian_face_dtype
!=============================================================================
!
! 2D convective flux Jacobian w.r.t. conservative variables, according to
!   J. Blazek's CFD book
!
!=============================================================================
function convective_jacobian_A_real ( q )  result(A)
      implicit none
      real(DP), dimension(4), intent(in)  :: q
      real(DP), dimension(4,4)            :: A

      real(DP) :: ek, h, phi, q22, q23, q33

      q22 = q(2)**2
      q33 = q(3)**2
      q23 = q(2)*q(3)

      ek   = 0.5*( q22 + q33 )
      h    = gm1i*q(4) + ek
      phi  = gm1*ek

      A(1,1) = 0.               ;   A(1,2) =   1.          ;   A(1,3) =   0.       ;   A(1,4) = 0.
      A(2,1) = phi - q22        ;   A(2,2) = - gm3*q(2)    ;   A(2,3) = - gm1*q(3) ;   A(2,4) = gm1
      A(3,1) =     - q23        ;   A(3,2) =   q(3)        ;   A(3,3) =   q(2)     ;   A(3,4) = 0.
      A(4,1) = q(2)*( phi - h ) ;   A(4,2) =   h - gm1*q22 ;   A(4,3) = - gm1*q23  ;   A(4,4) = gama*q(2)

end function convective_jacobian_A_real
!=============================================================================
!
! 2D convective flux Jacobian w.r.t. conservative variables, according to
!   J. Blazek's CFD book
!
!=============================================================================
function convective_jacobian_A_dtype ( q )  result(A)
      implicit none
      type(dType), dimension(4), intent(in)  :: q
      type(dType), dimension(4,4)            :: A

      type(dType) :: ek, h, phi, q22, q23, q33

      q22 = q(2)**2
      q33 = q(3)**2
      q23 = q(2)*q(3)

      ek   = 0.5*( q22 + q33 )
      h    = gm1i*q(4) + ek
      phi  = gm1*ek

      A(1,1) = 0.               ;   A(1,2) =   1.          ;   A(1,3) =   0.       ;   A(1,4) = 0.
      A(2,1) = phi - q22        ;   A(2,2) = - gm3*q(2)    ;   A(2,3) = - gm1*q(3) ;   A(2,4) = gm1
      A(3,1) =     - q23        ;   A(3,2) =   q(3)        ;   A(3,3) =   q(2)     ;   A(3,4) = 0.
      A(4,1) = q(2)*( phi - h ) ;   A(4,2) =   h - gm1*q22 ;   A(4,3) = - gm1*q23  ;   A(4,4) = gama*q(2)

end function convective_jacobian_A_dtype
!=============================================================================
!
! 2D convective flux Jacobian w.r.t. conservative variables, according to
!   J. Blazek's CFD book
!
!=============================================================================
function convective_jacobian_B_real ( q )  result(B)
      implicit none
      real(DP), dimension(4), intent(in)  :: q
      real(DP), dimension(4,4)            :: B

      real(DP) :: ek, h, phi, q22, q23, q33

      q22 = q(2)**2
      q33 = q(3)**2
      q23 = q(2)*q(3)

      ek   = 0.5*( q22 + q33 )
      h    = gm1i*q(4) + ek
      phi  = gm1*ek

      B(1,1) = 0.               ;   B(1,2) =   0.       ;   B(1,3) =   1.           ;   B(1,4) = 0.
      B(2,1) =     - q23        ;   B(2,2) =   q(3)     ;   B(2,3) =   q(2)         ;   B(2,4) = 0.
      B(3,1) = phi - q33        ;   B(3,2) = - gm1*q(2) ;   B(3,3) = - gm3*q(3)     ;   B(3,4) = gm1
      B(4,1) = q(3)*( phi - h ) ;   B(4,2) = - gm1*q23  ;   B(4,3) =   h - gm1*q33  ;   B(4,4) = gama*q(3)

end function convective_jacobian_B_real
!=============================================================================
!
! 2D convective flux Jacobian w.r.t. conservative variables, according to
!   J. Blazek's CFD book
!
!=============================================================================
function convective_jacobian_B_dtype ( q )  result(B)
      implicit none
      type(dType), dimension(4), intent(in)  :: q
      type(dType), dimension(4,4)            :: B

      type(dType) :: ek, h, phi, q22, q23, q33

      q22 = q(2)**2
      q33 = q(3)**2
      q23 = q(2)*q(3)

      ek   = 0.5*( q22 + q33 )
      h    = gm1i*q(4) + ek
      phi  = gm1*ek

      B(1,1) = 0.               ;   B(1,2) =   0.       ;   B(1,3) =   1.           ;   B(1,4) = 0.
      B(2,1) =     - q23        ;   B(2,2) =   q(3)     ;   B(2,3) =   q(2)         ;   B(2,4) = 0.
      B(3,1) = phi - q33        ;   B(3,2) = - gm1*q(2) ;   B(3,3) = - gm3*q(3)     ;   B(3,4) = gm1
      B(4,1) = q(3)*( phi - h ) ;   B(4,2) = - gm1*q23  ;   B(4,3) =   h - gm1*q33  ;   B(4,4) = gama*q(3)

end function convective_jacobian_B_dtype
!=============================================================================
!
! Get the inverse of stablization matrix
!
! The transformation matrices (M, Mi), and eigenvector matrices (R, Ri) are
!   constructed following Antony Jameson's paper '(3D) Eigenvalues and Eigenvectors
!   for the Gas Dynamics Equations', with (R, Ri) being modified, due to the
!   lack of the 2D counterpart of 3D (R, Ri) given in the paper
!
!=============================================================================
function get_tau_inverse_real ( q, Sx, Sy, nfb )  result( taui )
      implicit none
      integer,                  intent(in) :: nfb
      real(DP), dimension(4),   intent(in) :: q
      real(DP), dimension(nfb), intent(in) :: Sx, Sy
      real(DP), dimension(4,4)             :: taui

      real(DP), dimension(4,4) :: M, Mi, Ahat, tmat
      real(DP), dimension(4,4) :: R, Ri
      real(DP), dimension(4)   :: eig

      real(DP), parameter :: epsilon = 0.1
      real(DP) :: c, c2, ci, c2i, ek, c_S, ubar_S, eig_epsilon
      real(DP) :: S, nx, ny

      integer :: i

      ! for debug
      integer :: ii


     !write(*,100) 'q', q(1:4)


      c2 = q(4)
      c = sqrt(c2)
      ci = 1./c
      c2i = ci**2
      ek = 0.5*( q(2)**2 + q(3)**2 )

      M(1,1) = 1.   ;   M(1,2) = 0.     ;   M(1,3) = 0.     ;   M(1,4) = 0.
      M(2,1) = q(2) ;   M(2,2) = c      ;   M(2,3) = 0.     ;   M(2,4) = 0.
      M(3,1) = q(3) ;   M(3,2) = 0.     ;   M(3,3) = c      ;   M(3,4) = 0.
      M(4,1) = ek   ;   M(4,2) = q(2)*c ;   M(4,3) = q(3)*c ;   M(4,4) = gm1i*c2

      Mi(1,1) =  1.         ;   Mi(1,2) =  0.           ;   Mi(1,3) =  0.           ;   Mi(1,4) = 0.
      Mi(2,1) = -q(2)*ci    ;   Mi(2,2) =  ci           ;   Mi(2,3) =  0.           ;   Mi(2,4) = 0.
      Mi(3,1) = -q(3)*ci    ;   Mi(3,2) =  0.           ;   Mi(3,3) =  ci           ;   Mi(3,4) = 0.
      Mi(4,1) =  gm1*ek*c2i ;   Mi(4,2) = -gm1*q(2)*c2i ;   Mi(4,3) = -gm1*q(3)*c2i ;   Mi(4,4) = gm1*c2i

      ! Ahat = Ri*eig*R
      Ahat(:,:) = 0.

      do i = 1, nfb
     !do i = 3, 3                         ! for debug
        S = sqrt( Sx(i)**2 + Sy(i)**2 )
        nx = Sx(i)/S
        ny = Sy(i)/S
        c_S = c*S
        ubar_S = q(2)*Sx(i) + q(3)*Sy(i)

        eig(1:2) = ubar_S
        eig(3)   = ubar_S + c_S
        eig(4)   = ubar_S - c_S

        where ( eig < 0. )  eig = -eig

        ! eignvalue smoothing
        eig_epsilon = epsilon*( abs(ubar_S) + c_S )
        where ( eig < eig_epsilon ) eig = 0.5*( eig**2/eig_epsilon + eig_epsilon )


        ! >>> debug
       !eig(1:2) = ubar_S
       !eig(3)   = ubar_S + c_S
       !eig(4)   = ubar_S - c_S

       !write(*,100) 'Sx, Sy', Sx(i), Sy(i)
       !write(*,100) 'eig', eig(1:4)
        ! <<< debug

        R(1,1) = 1. ;   R(1,2) =  0. ;   R(1,3) = 1. ;   R(1,4) =  1.
        R(2,1) = 0. ;   R(2,2) =  ny ;   R(2,3) = nx ;   R(2,4) = -nx
        R(3,1) = 0. ;   R(3,2) = -nx ;   R(3,3) = ny ;   R(3,4) = -ny
        R(4,1) = 0. ;   R(4,2) =  0. ;   R(4,3) = 1. ;   R(4,4) =  1.

        Ri(1,1) = 1. ;   Ri(1,2) =  0.     ;   Ri(1,3) =  0.     ;   Ri(1,4) = -1.
        Ri(2,1) = 0. ;   Ri(2,2) =  ny     ;   Ri(2,3) = -nx     ;   Ri(2,4) =  0.
        Ri(3,1) = 0. ;   Ri(3,2) =  0.5*nx ;   Ri(3,3) =  0.5*ny ;   Ri(3,4) =  0.5
        Ri(4,1) = 0. ;   Ri(4,2) = -0.5*nx ;   Ri(4,3) = -0.5*ny ;   Ri(4,4) =  0.5

        Ahat(1,1) = Ahat(1,1) + eig(1)*( R(1,1)*Ri(1,1) ) + eig(2)*( R(1,2)*Ri(2,1) ) + eig(3)*( R(1,3)*Ri(3,1) ) + eig(4)*( R(1,4)*Ri(4,1) )
        Ahat(1,2) = Ahat(1,2) + eig(1)*( R(1,1)*Ri(1,2) ) + eig(2)*( R(1,2)*Ri(2,2) ) + eig(3)*( R(1,3)*Ri(3,2) ) + eig(4)*( R(1,4)*Ri(4,2) )
        Ahat(1,3) = Ahat(1,3) + eig(1)*( R(1,1)*Ri(1,3) ) + eig(2)*( R(1,2)*Ri(2,3) ) + eig(3)*( R(1,3)*Ri(3,3) ) + eig(4)*( R(1,4)*Ri(4,3) )
        Ahat(1,4) = Ahat(1,4) + eig(1)*( R(1,1)*Ri(1,4) ) + eig(2)*( R(1,2)*Ri(2,4) ) + eig(3)*( R(1,3)*Ri(3,4) ) + eig(4)*( R(1,4)*Ri(4,4) )

        Ahat(2,1) = Ahat(2,1) + eig(1)*( R(2,1)*Ri(1,1) ) + eig(2)*( R(2,2)*Ri(2,1) ) + eig(3)*( R(2,3)*Ri(3,1) ) + eig(4)*( R(2,4)*Ri(4,1) )
        Ahat(2,2) = Ahat(2,2) + eig(1)*( R(2,1)*Ri(1,2) ) + eig(2)*( R(2,2)*Ri(2,2) ) + eig(3)*( R(2,3)*Ri(3,2) ) + eig(4)*( R(2,4)*Ri(4,2) )
        Ahat(2,3) = Ahat(2,3) + eig(1)*( R(2,1)*Ri(1,3) ) + eig(2)*( R(2,2)*Ri(2,3) ) + eig(3)*( R(2,3)*Ri(3,3) ) + eig(4)*( R(2,4)*Ri(4,3) )
        Ahat(2,4) = Ahat(2,4) + eig(1)*( R(2,1)*Ri(1,4) ) + eig(2)*( R(2,2)*Ri(2,4) ) + eig(3)*( R(2,3)*Ri(3,4) ) + eig(4)*( R(2,4)*Ri(4,4) )

        Ahat(3,1) = Ahat(3,1) + eig(1)*( R(3,1)*Ri(1,1) ) + eig(2)*( R(3,2)*Ri(2,1) ) + eig(3)*( R(3,3)*Ri(3,1) ) + eig(4)*( R(3,4)*Ri(4,1) )
        Ahat(3,2) = Ahat(3,2) + eig(1)*( R(3,1)*Ri(1,2) ) + eig(2)*( R(3,2)*Ri(2,2) ) + eig(3)*( R(3,3)*Ri(3,2) ) + eig(4)*( R(3,4)*Ri(4,2) )
        Ahat(3,3) = Ahat(3,3) + eig(1)*( R(3,1)*Ri(1,3) ) + eig(2)*( R(3,2)*Ri(2,3) ) + eig(3)*( R(3,3)*Ri(3,3) ) + eig(4)*( R(3,4)*Ri(4,3) )
        Ahat(3,4) = Ahat(3,4) + eig(1)*( R(3,1)*Ri(1,4) ) + eig(2)*( R(3,2)*Ri(2,4) ) + eig(3)*( R(3,3)*Ri(3,4) ) + eig(4)*( R(3,4)*Ri(4,4) )

        Ahat(4,1) = Ahat(4,1) + eig(1)*( R(4,1)*Ri(1,1) ) + eig(2)*( R(4,2)*Ri(2,1) ) + eig(3)*( R(4,3)*Ri(3,1) ) + eig(4)*( R(4,4)*Ri(4,1) )
        Ahat(4,2) = Ahat(4,2) + eig(1)*( R(4,1)*Ri(1,2) ) + eig(2)*( R(4,2)*Ri(2,2) ) + eig(3)*( R(4,3)*Ri(3,2) ) + eig(4)*( R(4,4)*Ri(4,2) )
        Ahat(4,3) = Ahat(4,3) + eig(1)*( R(4,1)*Ri(1,3) ) + eig(2)*( R(4,2)*Ri(2,3) ) + eig(3)*( R(4,3)*Ri(3,3) ) + eig(4)*( R(4,4)*Ri(4,3) )
        Ahat(4,4) = Ahat(4,4) + eig(1)*( R(4,1)*Ri(1,4) ) + eig(2)*( R(4,2)*Ri(2,4) ) + eig(3)*( R(4,3)*Ri(3,4) ) + eig(4)*( R(4,4)*Ri(4,4) )
      end do

      ! tmat = Ahat*Mi
      tmat(:,:) = 0.

      tmat(1,1) = Ahat(1,1)*Mi(1,1) + Ahat(1,2)*Mi(2,1) + Ahat(1,3)*Mi(3,1) + Ahat(1,4)*Mi(4,1)
      tmat(1,2) = Ahat(1,1)*Mi(1,2) + Ahat(1,2)*Mi(2,2) + Ahat(1,3)*Mi(3,2) + Ahat(1,4)*Mi(4,2)
      tmat(1,3) = Ahat(1,1)*Mi(1,3) + Ahat(1,2)*Mi(2,3) + Ahat(1,3)*Mi(3,3) + Ahat(1,4)*Mi(4,3)
      tmat(1,4) = Ahat(1,1)*Mi(1,4) + Ahat(1,2)*Mi(2,4) + Ahat(1,3)*Mi(3,4) + Ahat(1,4)*Mi(4,4)

      tmat(2,1) = Ahat(2,1)*Mi(1,1) + Ahat(2,2)*Mi(2,1) + Ahat(2,3)*Mi(3,1) + Ahat(2,4)*Mi(4,1)
      tmat(2,2) = Ahat(2,1)*Mi(1,2) + Ahat(2,2)*Mi(2,2) + Ahat(2,3)*Mi(3,2) + Ahat(2,4)*Mi(4,2)
      tmat(2,3) = Ahat(2,1)*Mi(1,3) + Ahat(2,2)*Mi(2,3) + Ahat(2,3)*Mi(3,3) + Ahat(2,4)*Mi(4,3)
      tmat(2,4) = Ahat(2,1)*Mi(1,4) + Ahat(2,2)*Mi(2,4) + Ahat(2,3)*Mi(3,4) + Ahat(2,4)*Mi(4,4)

      tmat(3,1) = Ahat(3,1)*Mi(1,1) + Ahat(3,2)*Mi(2,1) + Ahat(3,3)*Mi(3,1) + Ahat(3,4)*Mi(4,1)
      tmat(3,2) = Ahat(3,1)*Mi(1,2) + Ahat(3,2)*Mi(2,2) + Ahat(3,3)*Mi(3,2) + Ahat(3,4)*Mi(4,2)
      tmat(3,3) = Ahat(3,1)*Mi(1,3) + Ahat(3,2)*Mi(2,3) + Ahat(3,3)*Mi(3,3) + Ahat(3,4)*Mi(4,3)
      tmat(3,4) = Ahat(3,1)*Mi(1,4) + Ahat(3,2)*Mi(2,4) + Ahat(3,3)*Mi(3,4) + Ahat(3,4)*Mi(4,4)

      tmat(4,1) = Ahat(4,1)*Mi(1,1) + Ahat(4,2)*Mi(2,1) + Ahat(4,3)*Mi(3,1) + Ahat(4,4)*Mi(4,1)
      tmat(4,2) = Ahat(4,1)*Mi(1,2) + Ahat(4,2)*Mi(2,2) + Ahat(4,3)*Mi(3,2) + Ahat(4,4)*Mi(4,2)
      tmat(4,3) = Ahat(4,1)*Mi(1,3) + Ahat(4,2)*Mi(2,3) + Ahat(4,3)*Mi(3,3) + Ahat(4,4)*Mi(4,3)
      tmat(4,4) = Ahat(4,1)*Mi(1,4) + Ahat(4,2)*Mi(2,4) + Ahat(4,3)*Mi(3,4) + Ahat(4,4)*Mi(4,4)

      ! taui = M*tmat
      taui(1,1) = M(1,1)*tmat(1,1) + M(1,2)*tmat(2,1) + M(1,3)*tmat(3,1) + M(1,4)*tmat(4,1)
      taui(1,2) = M(1,1)*tmat(1,2) + M(1,2)*tmat(2,2) + M(1,3)*tmat(3,2) + M(1,4)*tmat(4,2)
      taui(1,3) = M(1,1)*tmat(1,3) + M(1,2)*tmat(2,3) + M(1,3)*tmat(3,3) + M(1,4)*tmat(4,3)
      taui(1,4) = M(1,1)*tmat(1,4) + M(1,2)*tmat(2,4) + M(1,3)*tmat(3,4) + M(1,4)*tmat(4,4)

      taui(2,1) = M(2,1)*tmat(1,1) + M(2,2)*tmat(2,1) + M(2,3)*tmat(3,1) + M(2,4)*tmat(4,1)
      taui(2,2) = M(2,1)*tmat(1,2) + M(2,2)*tmat(2,2) + M(2,3)*tmat(3,2) + M(2,4)*tmat(4,2)
      taui(2,3) = M(2,1)*tmat(1,3) + M(2,2)*tmat(2,3) + M(2,3)*tmat(3,3) + M(2,4)*tmat(4,3)
      taui(2,4) = M(2,1)*tmat(1,4) + M(2,2)*tmat(2,4) + M(2,3)*tmat(3,4) + M(2,4)*tmat(4,4)

      taui(3,1) = M(3,1)*tmat(1,1) + M(3,2)*tmat(2,1) + M(3,3)*tmat(3,1) + M(3,4)*tmat(4,1)
      taui(3,2) = M(3,1)*tmat(1,2) + M(3,2)*tmat(2,2) + M(3,3)*tmat(3,2) + M(3,4)*tmat(4,2)
      taui(3,3) = M(3,1)*tmat(1,3) + M(3,2)*tmat(2,3) + M(3,3)*tmat(3,3) + M(3,4)*tmat(4,3)
      taui(3,4) = M(3,1)*tmat(1,4) + M(3,2)*tmat(2,4) + M(3,3)*tmat(3,4) + M(3,4)*tmat(4,4)

      taui(4,1) = M(4,1)*tmat(1,1) + M(4,2)*tmat(2,1) + M(4,3)*tmat(3,1) + M(4,4)*tmat(4,1)
      taui(4,2) = M(4,1)*tmat(1,2) + M(4,2)*tmat(2,2) + M(4,3)*tmat(3,2) + M(4,4)*tmat(4,2)
      taui(4,3) = M(4,1)*tmat(1,3) + M(4,2)*tmat(2,3) + M(4,3)*tmat(3,3) + M(4,4)*tmat(4,3)
      taui(4,4) = M(4,1)*tmat(1,4) + M(4,2)*tmat(2,4) + M(4,3)*tmat(3,4) + M(4,4)*tmat(4,4)


      if (Do_print) then
        write(*,100) 'q', q(1:4)
        do ii = 1, 4
          write(*,100) 'Taui', Taui(ii,1:4)
        end do
      end if

100   format(a,20(es16.5))
end function get_tau_inverse_real
!=============================================================================
!
! Get the inverse of stablization matrix
!
! The transformation matrices (M, Mi), and eigenvector matrices (R, Ri) are
!   constructed following Antony Jameson's paper '(3D) Eigenvalues and Eigenvectors
!   for the Gas Dynamics Equations', with (R, Ri) being modified, due to the
!   lack of the 2D counterpart of 3D (R, Ri) given in the paper
!
!=============================================================================
function get_tau_inverse_dtype ( q, Sx, Sy, nfb )  result( taui )
      implicit none
      integer,                     intent(in) :: nfb
      type(dType), dimension(4),   intent(in) :: q
      real(DP),    dimension(nfb), intent(in) :: Sx, Sy
      type(dType), dimension(4,4)             :: taui

      type(dType), dimension(4,4) :: M, Mi, Ahat, tmat
      real(DP),    dimension(4,4) :: R, Ri
      type(dType), dimension(4)   :: eig

      real(DP), parameter :: epsilon = 0.1
      type(dType) :: c, c2, ci, c2i, ek, c_S, ubar_S, eig_epsilon
      real(DP)    :: S, nx, ny

      integer :: i

      c2 = q(4)
      c = sqrt(c2)
      ci = 1./c
      c2i = ci**2
      ek = 0.5*( q(2)**2 + q(3)**2 )

      M(1,1) = 1.   ;   M(1,2) = 0.     ;   M(1,3) = 0.     ;   M(1,4) = 0.
      M(2,1) = q(2) ;   M(2,2) = c      ;   M(2,3) = 0.     ;   M(2,4) = 0.
      M(3,1) = q(3) ;   M(3,2) = 0.     ;   M(3,3) = c      ;   M(3,4) = 0.
      M(4,1) = ek   ;   M(4,2) = q(2)*c ;   M(4,3) = q(3)*c ;   M(4,4) = gm1i*c2

      Mi(1,1) =  1.         ;   Mi(1,2) =  0.           ;   Mi(1,3) =  0.           ;   Mi(1,4) = 0.
      Mi(2,1) = -q(2)*ci    ;   Mi(2,2) =  ci           ;   Mi(2,3) =  0.           ;   Mi(2,4) = 0.
      Mi(3,1) = -q(3)*ci    ;   Mi(3,2) =  0.           ;   Mi(3,3) =  ci           ;   Mi(3,4) = 0.
      Mi(4,1) =  gm1*ek*c2i ;   Mi(4,2) = -gm1*q(2)*c2i ;   Mi(4,3) = -gm1*q(3)*c2i ;   Mi(4,4) = gm1*c2i

      ! Ahat = Ri*eig*R
      Ahat(:,:) = 0.

      do i = 1, nfb
        S = sqrt( Sx(i)**2 + Sy(i)**2 )
        nx = Sx(i)/S
        ny = Sy(i)/S
        c_S = c*S
        ubar_S = q(2)*Sx(i) + q(3)*Sy(i)

        eig(1:2) = ubar_S
        eig(3)   = ubar_S + c_S
        eig(4)   = ubar_S - c_S

        where ( eig < 0. )  eig = -eig

        ! eignvalue smoothing
        eig_epsilon = epsilon*( abs(ubar_S) + c_S )
        where ( eig < eig_epsilon ) eig = 0.5*( eig**2/eig_epsilon + eig_epsilon )

        R(1,1) = 1. ;   R(1,2) =  0. ;   R(1,3) = 1. ;   R(1,4) =  1.
        R(2,1) = 0. ;   R(2,2) =  ny ;   R(2,3) = nx ;   R(2,4) = -nx
        R(3,1) = 0. ;   R(3,2) = -nx ;   R(3,3) = ny ;   R(3,4) = -ny
        R(4,1) = 0. ;   R(4,2) =  0. ;   R(4,3) = 1. ;   R(4,4) =  1.

        Ri(1,1) = 1. ;   Ri(1,2) =  0.     ;   Ri(1,3) =  0.     ;   Ri(1,4) = -1.
        Ri(2,1) = 0. ;   Ri(2,2) =  ny     ;   Ri(2,3) = -nx     ;   Ri(2,4) =  0.
        Ri(3,1) = 0. ;   Ri(3,2) =  0.5*nx ;   Ri(3,3) =  0.5*ny ;   Ri(3,4) =  0.5
        Ri(4,1) = 0. ;   Ri(4,2) = -0.5*nx ;   Ri(4,3) = -0.5*ny ;   Ri(4,4) =  0.5

        Ahat(1,1) = Ahat(1,1) + eig(1)*( R(1,1)*Ri(1,1) ) + eig(2)*( R(1,2)*Ri(2,1) ) + eig(3)*( R(1,3)*Ri(3,1) ) + eig(4)*( R(1,4)*Ri(4,1) )
        Ahat(1,2) = Ahat(1,2) + eig(1)*( R(1,1)*Ri(1,2) ) + eig(2)*( R(1,2)*Ri(2,2) ) + eig(3)*( R(1,3)*Ri(3,2) ) + eig(4)*( R(1,4)*Ri(4,2) )
        Ahat(1,3) = Ahat(1,3) + eig(1)*( R(1,1)*Ri(1,3) ) + eig(2)*( R(1,2)*Ri(2,3) ) + eig(3)*( R(1,3)*Ri(3,3) ) + eig(4)*( R(1,4)*Ri(4,3) )
        Ahat(1,4) = Ahat(1,4) + eig(1)*( R(1,1)*Ri(1,4) ) + eig(2)*( R(1,2)*Ri(2,4) ) + eig(3)*( R(1,3)*Ri(3,4) ) + eig(4)*( R(1,4)*Ri(4,4) )

        Ahat(2,1) = Ahat(2,1) + eig(1)*( R(2,1)*Ri(1,1) ) + eig(2)*( R(2,2)*Ri(2,1) ) + eig(3)*( R(2,3)*Ri(3,1) ) + eig(4)*( R(2,4)*Ri(4,1) )
        Ahat(2,2) = Ahat(2,2) + eig(1)*( R(2,1)*Ri(1,2) ) + eig(2)*( R(2,2)*Ri(2,2) ) + eig(3)*( R(2,3)*Ri(3,2) ) + eig(4)*( R(2,4)*Ri(4,2) )
        Ahat(2,3) = Ahat(2,3) + eig(1)*( R(2,1)*Ri(1,3) ) + eig(2)*( R(2,2)*Ri(2,3) ) + eig(3)*( R(2,3)*Ri(3,3) ) + eig(4)*( R(2,4)*Ri(4,3) )
        Ahat(2,4) = Ahat(2,4) + eig(1)*( R(2,1)*Ri(1,4) ) + eig(2)*( R(2,2)*Ri(2,4) ) + eig(3)*( R(2,3)*Ri(3,4) ) + eig(4)*( R(2,4)*Ri(4,4) )

        Ahat(3,1) = Ahat(3,1) + eig(1)*( R(3,1)*Ri(1,1) ) + eig(2)*( R(3,2)*Ri(2,1) ) + eig(3)*( R(3,3)*Ri(3,1) ) + eig(4)*( R(3,4)*Ri(4,1) )
        Ahat(3,2) = Ahat(3,2) + eig(1)*( R(3,1)*Ri(1,2) ) + eig(2)*( R(3,2)*Ri(2,2) ) + eig(3)*( R(3,3)*Ri(3,2) ) + eig(4)*( R(3,4)*Ri(4,2) )
        Ahat(3,3) = Ahat(3,3) + eig(1)*( R(3,1)*Ri(1,3) ) + eig(2)*( R(3,2)*Ri(2,3) ) + eig(3)*( R(3,3)*Ri(3,3) ) + eig(4)*( R(3,4)*Ri(4,3) )
        Ahat(3,4) = Ahat(3,4) + eig(1)*( R(3,1)*Ri(1,4) ) + eig(2)*( R(3,2)*Ri(2,4) ) + eig(3)*( R(3,3)*Ri(3,4) ) + eig(4)*( R(3,4)*Ri(4,4) )

        Ahat(4,1) = Ahat(4,1) + eig(1)*( R(4,1)*Ri(1,1) ) + eig(2)*( R(4,2)*Ri(2,1) ) + eig(3)*( R(4,3)*Ri(3,1) ) + eig(4)*( R(4,4)*Ri(4,1) )
        Ahat(4,2) = Ahat(4,2) + eig(1)*( R(4,1)*Ri(1,2) ) + eig(2)*( R(4,2)*Ri(2,2) ) + eig(3)*( R(4,3)*Ri(3,2) ) + eig(4)*( R(4,4)*Ri(4,2) )
        Ahat(4,3) = Ahat(4,3) + eig(1)*( R(4,1)*Ri(1,3) ) + eig(2)*( R(4,2)*Ri(2,3) ) + eig(3)*( R(4,3)*Ri(3,3) ) + eig(4)*( R(4,4)*Ri(4,3) )
        Ahat(4,4) = Ahat(4,4) + eig(1)*( R(4,1)*Ri(1,4) ) + eig(2)*( R(4,2)*Ri(2,4) ) + eig(3)*( R(4,3)*Ri(3,4) ) + eig(4)*( R(4,4)*Ri(4,4) )
      end do

      ! tmat = Ahat*Mi
      tmat(:,:) = 0.

      tmat(1,1) = Ahat(1,1)*Mi(1,1) + Ahat(1,2)*Mi(2,1) + Ahat(1,3)*Mi(3,1) + Ahat(1,4)*Mi(4,1)
      tmat(1,2) = Ahat(1,1)*Mi(1,2) + Ahat(1,2)*Mi(2,2) + Ahat(1,3)*Mi(3,2) + Ahat(1,4)*Mi(4,2)
      tmat(1,3) = Ahat(1,1)*Mi(1,3) + Ahat(1,2)*Mi(2,3) + Ahat(1,3)*Mi(3,3) + Ahat(1,4)*Mi(4,3)
      tmat(1,4) = Ahat(1,1)*Mi(1,4) + Ahat(1,2)*Mi(2,4) + Ahat(1,3)*Mi(3,4) + Ahat(1,4)*Mi(4,4)

      tmat(2,1) = Ahat(2,1)*Mi(1,1) + Ahat(2,2)*Mi(2,1) + Ahat(2,3)*Mi(3,1) + Ahat(2,4)*Mi(4,1)
      tmat(2,2) = Ahat(2,1)*Mi(1,2) + Ahat(2,2)*Mi(2,2) + Ahat(2,3)*Mi(3,2) + Ahat(2,4)*Mi(4,2)
      tmat(2,3) = Ahat(2,1)*Mi(1,3) + Ahat(2,2)*Mi(2,3) + Ahat(2,3)*Mi(3,3) + Ahat(2,4)*Mi(4,3)
      tmat(2,4) = Ahat(2,1)*Mi(1,4) + Ahat(2,2)*Mi(2,4) + Ahat(2,3)*Mi(3,4) + Ahat(2,4)*Mi(4,4)

      tmat(3,1) = Ahat(3,1)*Mi(1,1) + Ahat(3,2)*Mi(2,1) + Ahat(3,3)*Mi(3,1) + Ahat(3,4)*Mi(4,1)
      tmat(3,2) = Ahat(3,1)*Mi(1,2) + Ahat(3,2)*Mi(2,2) + Ahat(3,3)*Mi(3,2) + Ahat(3,4)*Mi(4,2)
      tmat(3,3) = Ahat(3,1)*Mi(1,3) + Ahat(3,2)*Mi(2,3) + Ahat(3,3)*Mi(3,3) + Ahat(3,4)*Mi(4,3)
      tmat(3,4) = Ahat(3,1)*Mi(1,4) + Ahat(3,2)*Mi(2,4) + Ahat(3,3)*Mi(3,4) + Ahat(3,4)*Mi(4,4)

      tmat(4,1) = Ahat(4,1)*Mi(1,1) + Ahat(4,2)*Mi(2,1) + Ahat(4,3)*Mi(3,1) + Ahat(4,4)*Mi(4,1)
      tmat(4,2) = Ahat(4,1)*Mi(1,2) + Ahat(4,2)*Mi(2,2) + Ahat(4,3)*Mi(3,2) + Ahat(4,4)*Mi(4,2)
      tmat(4,3) = Ahat(4,1)*Mi(1,3) + Ahat(4,2)*Mi(2,3) + Ahat(4,3)*Mi(3,3) + Ahat(4,4)*Mi(4,3)
      tmat(4,4) = Ahat(4,1)*Mi(1,4) + Ahat(4,2)*Mi(2,4) + Ahat(4,3)*Mi(3,4) + Ahat(4,4)*Mi(4,4)

      ! taui = M*tmat
      taui(1,1) = M(1,1)*tmat(1,1) + M(1,2)*tmat(2,1) + M(1,3)*tmat(3,1) + M(1,4)*tmat(4,1)
      taui(1,2) = M(1,1)*tmat(1,2) + M(1,2)*tmat(2,2) + M(1,3)*tmat(3,2) + M(1,4)*tmat(4,2)
      taui(1,3) = M(1,1)*tmat(1,3) + M(1,2)*tmat(2,3) + M(1,3)*tmat(3,3) + M(1,4)*tmat(4,3)
      taui(1,4) = M(1,1)*tmat(1,4) + M(1,2)*tmat(2,4) + M(1,3)*tmat(3,4) + M(1,4)*tmat(4,4)

      taui(2,1) = M(2,1)*tmat(1,1) + M(2,2)*tmat(2,1) + M(2,3)*tmat(3,1) + M(2,4)*tmat(4,1)
      taui(2,2) = M(2,1)*tmat(1,2) + M(2,2)*tmat(2,2) + M(2,3)*tmat(3,2) + M(2,4)*tmat(4,2)
      taui(2,3) = M(2,1)*tmat(1,3) + M(2,2)*tmat(2,3) + M(2,3)*tmat(3,3) + M(2,4)*tmat(4,3)
      taui(2,4) = M(2,1)*tmat(1,4) + M(2,2)*tmat(2,4) + M(2,3)*tmat(3,4) + M(2,4)*tmat(4,4)

      taui(3,1) = M(3,1)*tmat(1,1) + M(3,2)*tmat(2,1) + M(3,3)*tmat(3,1) + M(3,4)*tmat(4,1)
      taui(3,2) = M(3,1)*tmat(1,2) + M(3,2)*tmat(2,2) + M(3,3)*tmat(3,2) + M(3,4)*tmat(4,2)
      taui(3,3) = M(3,1)*tmat(1,3) + M(3,2)*tmat(2,3) + M(3,3)*tmat(3,3) + M(3,4)*tmat(4,3)
      taui(3,4) = M(3,1)*tmat(1,4) + M(3,2)*tmat(2,4) + M(3,3)*tmat(3,4) + M(3,4)*tmat(4,4)

      taui(4,1) = M(4,1)*tmat(1,1) + M(4,2)*tmat(2,1) + M(4,3)*tmat(3,1) + M(4,4)*tmat(4,1)
      taui(4,2) = M(4,1)*tmat(1,2) + M(4,2)*tmat(2,2) + M(4,3)*tmat(3,2) + M(4,4)*tmat(4,2)
      taui(4,3) = M(4,1)*tmat(1,3) + M(4,2)*tmat(2,3) + M(4,3)*tmat(3,3) + M(4,4)*tmat(4,3)
      taui(4,4) = M(4,1)*tmat(1,4) + M(4,2)*tmat(2,4) + M(4,3)*tmat(3,4) + M(4,4)*tmat(4,4)

end function get_tau_inverse_dtype
!=============================================================================
!
!=============================================================================
function roe_flux_face_real ( q1, q2, Sx, Sy )  result(flux)
      implicit none
      real(DP), dimension(4),   intent(in) :: q1, q2
      real(DP),                 intent(in) :: Sx, Sy
      real(DP), dimension(4)               :: flux

      real(DP), dimension(4,4) :: M, Mi, Ahat, tmat
      real(DP), dimension(4,4) :: R, Ri
      real(DP), dimension(4)   :: eig

      real(DP), parameter :: epsilon = 0.1
      real(DP) :: h, c, c2, ci, c2i, ek, c_S, ubar_S, eig_epsilon
      real(DP) :: S, nx, ny
      real(DP), dimension(4) :: q, dcq, v1, v2, v3, v4, v5
      real(DP) :: h1, h2, sq_rho1, sq_rho2
      real(DP) :: coe1, coe2


      ! roe average
      h1 = gm1i*q1(4) + 0.5*( q1(2)**2 + q1(3)**2 )
      h2 = gm1i*q2(4) + 0.5*( q2(2)**2 + q2(3)**2 )

      sq_rho1 = sqrt(q1(1))
      sq_rho2 = sqrt(q2(1))
      coe1 = sq_rho1/(sq_rho1 + sq_rho2)
      coe2 = 1. - coe1

      q(1) = sq_rho1*sq_rho2
      q(2) = coe1*q1(1) + coe2*q2(1)
      q(3) = coe1*q1(3) + coe2*q2(3)
      h    = coe1*h1    + coe2*h2

      ek   = 0.5*( q(2)**2 + q(3)**2 )
      q(4) = gm1*( h - ek )

      c2 = q(4)
      c = sqrt(c2)
      ci = 1./c
      c2i = ci**2

      ! M, Mi
      M(1,1) = 1.   ;   M(1,2) = 0.     ;   M(1,3) = 0.     ;   M(1,4) = 0.
      M(2,1) = q(2) ;   M(2,2) = c      ;   M(2,3) = 0.     ;   M(2,4) = 0.
      M(3,1) = q(3) ;   M(3,2) = 0.     ;   M(3,3) = c      ;   M(3,4) = 0.
      M(4,1) = ek   ;   M(4,2) = q(2)*c ;   M(4,3) = q(3)*c ;   M(4,4) = gm1i*c2

      Mi(1,1) =  1.         ;   Mi(1,2) =  0.           ;   Mi(1,3) =  0.           ;   Mi(1,4) = 0.
      Mi(2,1) = -q(2)*ci    ;   Mi(2,2) =  ci           ;   Mi(2,3) =  0.           ;   Mi(2,4) = 0.
      Mi(3,1) = -q(3)*ci    ;   Mi(3,2) =  0.           ;   Mi(3,3) =  ci           ;   Mi(3,4) = 0.
      Mi(4,1) =  gm1*ek*c2i ;   Mi(4,2) = -gm1*q(2)*c2i ;   Mi(4,3) = -gm1*q(3)*c2i ;   Mi(4,4) = gm1*c2i


      ! eigenvalue
      S = sqrt( Sx**2 + Sy**2 )
      nx = Sx/S
      ny = Sy/S
      c_S = c*S
      ubar_S = q(2)*Sx+ q(3)*Sy

      eig(1:2) = ubar_S
      eig(3)   = ubar_S + c_S
      eig(4)   = ubar_S - c_S

      where ( eig < 0. )  eig = -eig

      ! eignvalue smoothing
      eig_epsilon = epsilon*( abs(ubar_S) + c_S )
      where ( eig < eig_epsilon ) eig = 0.5*( eig**2/eig_epsilon + eig_epsilon )


      ! R, Ri
      R(1,1) = 1. ;   R(1,2) =  0. ;   R(1,3) = 1. ;   R(1,4) =  1.
      R(2,1) = 0. ;   R(2,2) =  ny ;   R(2,3) = nx ;   R(2,4) = -nx
      R(3,1) = 0. ;   R(3,2) = -nx ;   R(3,3) = ny ;   R(3,4) = -ny
      R(4,1) = 0. ;   R(4,2) =  0. ;   R(4,3) = 1. ;   R(4,4) =  1.

      Ri(1,1) = 1. ;   Ri(1,2) =  0.     ;   Ri(1,3) =  0.     ;   Ri(1,4) = -1.
      Ri(2,1) = 0. ;   Ri(2,2) =  ny     ;   Ri(2,3) = -nx     ;   Ri(2,4) =  0.
      Ri(3,1) = 0. ;   Ri(3,2) =  0.5*nx ;   Ri(3,3) =  0.5*ny ;   Ri(3,4) =  0.5
      Ri(4,1) = 0. ;   Ri(4,2) = -0.5*nx ;   Ri(4,3) = -0.5*ny ;   Ri(4,4) =  0.5


      ! d_cq
      dcq(:) = conservative_variable_cq ( q2 ) - conservative_variable_cq ( q1 )

      v1(:) = matrix_vector_multiply( Mi, dcq, 4, 4 )
      v2(:) = matrix_vector_multiply( Ri, v1, 4, 4 )
      v3(:) = eig(:)*v2(:)
      v4(:) = matrix_vector_multiply( R, v3, 4, 4 )
      v5(:) = matrix_vector_multiply( M, v4, 4, 4 )

      flux(:) = 0.5*( convective_flux_face ( q1, Sx, Sy ) + convective_flux_face( q2, Sx, Sy ) - v5(:) )

end function roe_flux_face_real
!=============================================================================
!
!=============================================================================
function roe_flux_face_dtype ( q1, q2, Sx, Sy )  result(flux)
      implicit none
      type(dType), dimension(4),   intent(in) :: q1, q2
      real(DP),                    intent(in) :: Sx, Sy
      type(dType), dimension(4)               :: flux

      type(dType), dimension(4,4) :: M, Mi
      real(DP),    dimension(4,4) :: R, Ri
      type(dType), dimension(4)   :: eig

      real(DP), parameter :: epsilon = 0.1
      type(dType) :: h, c, c2, ci, c2i, ek, c_S, ubar_S, eig_epsilon
      type(dType) :: S, nx, ny
      type(dType), dimension(4) :: q, dcq, v1, v2, v3, v4, v5
      type(dType) :: h1, h2, sq_rho1, sq_rho2
      real(DP) :: coe1, coe2


      ! roe average
      h1 = gm1i*q1(4) + 0.5*( q1(2)**2 + q1(3)**2 )
      h2 = gm1i*q2(4) + 0.5*( q2(2)**2 + q2(3)**2 )

      sq_rho1 = sqrt(q1(1))
      sq_rho2 = sqrt(q2(1))
      coe1 = sq_rho1/(sq_rho1 + sq_rho2)
      coe2 = 1. - coe1

      q(1) = sq_rho1*sq_rho2
      q(2) = coe1*q1(1) + coe2*q2(1)
      q(3) = coe1*q1(3) + coe2*q2(3)
      h    = coe1*h1    + coe2*h2

      ek   = 0.5*( q(2)**2 + q(3)**2 )
      q(4) = gm1*( h - ek )

      c2 = q(4)
      c = sqrt(c2)
      ci = 1./c
      c2i = ci**2

      ! M, Mi
      M(1,1) = 1.   ;   M(1,2) = 0.     ;   M(1,3) = 0.     ;   M(1,4) = 0.
      M(2,1) = q(2) ;   M(2,2) = c      ;   M(2,3) = 0.     ;   M(2,4) = 0.
      M(3,1) = q(3) ;   M(3,2) = 0.     ;   M(3,3) = c      ;   M(3,4) = 0.
      M(4,1) = ek   ;   M(4,2) = q(2)*c ;   M(4,3) = q(3)*c ;   M(4,4) = gm1i*c2

      Mi(1,1) =  1.         ;   Mi(1,2) =  0.           ;   Mi(1,3) =  0.           ;   Mi(1,4) = 0.
      Mi(2,1) = -q(2)*ci    ;   Mi(2,2) =  ci           ;   Mi(2,3) =  0.           ;   Mi(2,4) = 0.
      Mi(3,1) = -q(3)*ci    ;   Mi(3,2) =  0.           ;   Mi(3,3) =  ci           ;   Mi(3,4) = 0.
      Mi(4,1) =  gm1*ek*c2i ;   Mi(4,2) = -gm1*q(2)*c2i ;   Mi(4,3) = -gm1*q(3)*c2i ;   Mi(4,4) = gm1*c2i


      ! eigenvalue
      S = sqrt( Sx**2 + Sy**2 )
      nx = Sx/S
      ny = Sy/S
      c_S = c*S
      ubar_S = q(2)*Sx+ q(3)*Sy

      eig(1:2) = ubar_S
      eig(3)   = ubar_S + c_S
      eig(4)   = ubar_S - c_S

      where ( eig < 0. )  eig = -eig

      ! eignvalue smoothing
      eig_epsilon = epsilon*( abs(ubar_S) + c_S )
      where ( eig < eig_epsilon ) eig = 0.5*( eig**2/eig_epsilon + eig_epsilon )


      ! R, Ri
      R(1,1) = 1. ;   R(1,2) =  0. ;   R(1,3) = 1. ;   R(1,4) =  1.
      R(2,1) = 0. ;   R(2,2) =  ny ;   R(2,3) = nx ;   R(2,4) = -nx
      R(3,1) = 0. ;   R(3,2) = -nx ;   R(3,3) = ny ;   R(3,4) = -ny
      R(4,1) = 0. ;   R(4,2) =  0. ;   R(4,3) = 1. ;   R(4,4) =  1.

      Ri(1,1) = 1. ;   Ri(1,2) =  0.     ;   Ri(1,3) =  0.     ;   Ri(1,4) = -1.
      Ri(2,1) = 0. ;   Ri(2,2) =  ny     ;   Ri(2,3) = -nx     ;   Ri(2,4) =  0.
      Ri(3,1) = 0. ;   Ri(3,2) =  0.5*nx ;   Ri(3,3) =  0.5*ny ;   Ri(3,4) =  0.5
      Ri(4,1) = 0. ;   Ri(4,2) = -0.5*nx ;   Ri(4,3) = -0.5*ny ;   Ri(4,4) =  0.5


      ! d_cq
      dcq(:) = conservative_variable_cq ( q2 ) - conservative_variable_cq ( q1 )

      v1(:) = matrix_vector_multiply( Mi, dcq, 4, 4 )
      v2(:) = matrix_vector_multiply( Ri, v1, 4, 4 )
      v3(:) = eig(:)*v2(:)
      v4(:) = matrix_vector_multiply( R, v3, 4, 4 )
      v5(:) = matrix_vector_multiply( M, v4, 4, 4 )

      flux(:) = 0.5*( convective_flux_face ( q1, Sx, Sy ) + convective_flux_face( q2, Sx, Sy ) - v5(:) )

end function roe_flux_face_dtype
!=============================================================================
!
! added dissipative flux for shock capturing, 
!   following Ryan Glasby and Erwin Taylor's paper on COFFE
!
!=============================================================================
subroutine shock_capture_flux_f_g_real ( f_shock, g_shock, q, qx, qy, h_length, k_shock )
      implicit none
      real(DP), dimension(4), intent(out) :: f_shock, g_shock
      real(DP), dimension(4), intent(in)  :: q, qx, qy
      real(DP),               intent(in)  :: h_length
      real(DP),               intent(in)  :: k_shock

      real(DP) :: p, px, py
      real(DP) :: ek, ekx, eky
      real(DP) :: h, hx, hy
      real(DP) :: e_shock
      real(DP) :: lamda
      real(DP) :: o1
      real(DP) :: o2
      real(DP) :: o3
      real(DP) :: o4
      real(DP) :: o5


      ! p
      p  = gi*q(1)*q(4)
      px = gi*( qx(1)*q(4) + q(1)*qx(4) )
      py = gi*( qy(1)*q(4) + q(1)*qy(4) )

      ! ek
      ek  = 0.5*( q(2)**2 + q(3)**2 )
      ekx = q(2)*qx(2) + q(3)*qx(3)
      eky = q(2)*qy(2) + q(3)*qy(3)

      ! h
      h   = gm1i*q(4) + ek
      hx  = gm1i*qx(4) + ekx
      hy  = gm1i*qy(4) + eky

      ! o1 = px**2 + py**2
      o1  = px**2 + py**2

      ! o2
      o2  = h_length*sqrt(o1)

      ! o3
      o3  = o2  + k_shock*p

      ! e_shock
      e_shock = o2/o3

      ! o4
      o4 = q(2)**2 + q(3)**2

      ! lamda
      lamda = sqrt( o4 )  + sqrt( q(4) )

      ! o5
      o5  = h_length*lamda*e_shock

      ! f_shock
      f_shock(1) = o5*  qx(1)
      f_shock(2) = o5*( qx(1)*q(2) + q(1)*qx(2) )
      f_shock(3) = o5*( qx(1)*q(3) + q(1)*qx(3) )
      f_shock(4) = o5*( qx(1)*h    + q(1)*hx    )

      ! g_shock
      g_shock(1) = o5*  qy(1)
      g_shock(2) = o5*( qy(1)*q(2) + q(1)*qy(2) )
      g_shock(3) = o5*( qy(1)*q(3) + q(1)*qy(3) )
      g_shock(4) = o5*( qy(1)*h    + q(1)*hy    )

end subroutine shock_capture_flux_f_g_real
!=============================================================================
!
! added dissipative flux for shock capturing, 
!   following Ryan Glasby and Erwin Taylor's paper on COFFE
!
!=============================================================================
subroutine shock_capture_flux_f_g_dtype ( f_shock, g_shock, q, qx, qy, h_length, k_shock )
      implicit none
      type(dType),  dimension(4), intent(out) :: f_shock, g_shock
      type(dType),  dimension(4), intent(in)  :: q, qx, qy
      real(DP),                   intent(in)  :: h_length
      real(DP),                   intent(in)  :: k_shock

      type(dType) :: p, px, py
      type(dType) :: ek, ekx, eky
      type(dType) :: h, hx, hy
      type(dType) :: e_shock
      type(dType) :: lamda
      type(dType) :: o1
      type(dType) :: o2
      type(dType) :: o3
      type(dType) :: o4
      type(dType) :: o5

      ! p
      p  = gi*q(1)*q(4)
      px = gi*( qx(1)*q(4) + q(1)*qx(4) )
      py = gi*( qy(1)*q(4) + q(1)*qy(4) )

      ! ek
      ek  = 0.5*( q(2)**2 + q(3)**2 )
      ekx = q(2)*qx(2) + q(3)*qx(3)
      eky = q(2)*qy(2) + q(3)*qy(3)

      ! h
      h   = gm1i*q(4) + ek
      hx  = gm1i*qx(4) + ekx
      hy  = gm1i*qy(4) + eky

      ! o1 = px**2 + py**2
      o1  = px**2 + py**2

      ! o2
      o2  = h_length*sqrt(o1)

      ! o3
      o3  = o2  + k_shock*p

      ! e_shock
      e_shock = o2/o3

      ! o4
      o4 = q(2)**2 + q(3)**2

      ! lamda
      lamda = sqrt( o4 )  + sqrt( q(4) )

      ! o5
      o5  = h_length*lamda*e_shock

      ! f_shock
      f_shock(1) = o5*  qx(1)
      f_shock(2) = o5*( qx(1)*q(2) + q(1)*qx(2) )
      f_shock(3) = o5*( qx(1)*q(3) + q(1)*qx(3) )
      f_shock(4) = o5*( qx(1)*h    + q(1)*hx    )

      ! g_shock
      g_shock(1) = o5*  qy(1)
      g_shock(2) = o5*( qy(1)*q(2) + q(1)*qy(2) )
      g_shock(3) = o5*( qy(1)*q(3) + q(1)*qy(3) )
      g_shock(4) = o5*( qy(1)*h    + q(1)*hy    )


end subroutine shock_capture_flux_f_g_dtype
!=============================================================================
!
! added dissipative flux for shock capturing, 
!   following Ryan Glasby and Erwin Taylor's paper on COFFE
!
!=============================================================================
subroutine shock_capture_flux_fx_gy_real ( fx_shock, gy_shock, q, qx, qy, qxx, qxy, qyy, h_length, k_shock )
      implicit none
      real(DP), dimension(4), intent(out) :: fx_shock, gy_shock
      real(DP), dimension(4), intent(in)  :: q, qx, qy, qxx, qxy, qyy
      real(DP),               intent(in)  :: h_length
      real(DP),               intent(in)  :: k_shock

      real(DP) :: p, px, py, pxx, pxy, pyy
      real(DP) :: ek, ekx, eky, ekxx, ekyy
      real(DP) :: h, hx, hy, hxx, hyy
      real(DP) :: e_shock, e_shock_x, e_shock_y
      real(DP) :: lamda, lamda_x, lamda_y
      real(DP) :: o1, o1x, o1y
      real(DP) :: o2, o2x, o2y
      real(DP) :: o3, o3x, o3y
      real(DP) :: o4, o4x, o4y
      real(DP) :: o5, o5x, o5y


      ! p
      p   = gi*q(1)*q(4)
      px  = gi*( qx(1)*q(4) + q(1)*qx(4) )
      py  = gi*( qy(1)*q(4) + q(1)*qy(4) )
      pxx = gi*( qxx(1)*q(4) + qx(1)*qx(4) + qx(1)*qx(4) + q(1)*qxx(4) )
      pxy = gi*( qxy(1)*q(4) + qx(1)*qy(4) + qy(1)*qx(4) + q(1)*qxy(4) )
      pyy = gi*( qyy(1)*q(4) + qy(1)*qy(4) + qy(1)*qy(4) + q(1)*qyy(4) )

      ! ek
      ek   = 0.5*( q(2)**2 + q(3)**2 )
      ekx  = q(2)*qx(2) + q(3)*qx(3)
      eky  = q(2)*qy(2) + q(3)*qy(3)
      ekxx = qx(2)*qx(2) + q(2)*qxx(2) + qx(3)*qx(3) + q(3)*qxx(3)
      ekyy = qy(2)*qy(2) + q(2)*qyy(2) + qy(3)*qy(3) + q(3)*qyy(3)

      ! h
      h   = gm1i*q  (4) + ek
      hx  = gm1i*qx (4) + ekx
      hy  = gm1i*qy (4) + eky
      hxx = gm1i*qxx(4) + ekxx
      hyy = gm1i*qyy(4) + ekyy

      ! o1 = px**2 + py**2
      o1  = px**2 + py**2
      o1x = 2.*( px*pxx + py*pxy )
      o1y = 2.*( px*pxy + py*pyy )

      ! o2
      o2  = h_length*sqrt(o1)
      o2x = h_length*0.5/sqrt(o1)*o1x
      o2y = h_length*0.5/sqrt(o1)*o1y

      ! o3
      o3  = o2  + k_shock*p
      o3x = o2x + k_shock*px
      o3y = o2y + k_shock*py

      ! e_shock
      e_shock = o2/o3
      e_shock_x = ( o2x*o3 - o2*o3x )/o3**2
      e_shock_y = ( o2y*o3 - o2*o3y )/o3**2

      ! o4
      o4 = q(2)**2 + q(3)**2
      o4x = 2.*( q(2)*qx(2) + q(3)*qx(3) )
      o4y = 2.*( q(2)*qy(2) + q(3)*qy(3) )

      ! lamda
      lamda = sqrt( o4 )  + sqrt( q(4) )
      lamda_x = 0.5/sqrt( o4 )*o4x + 0.5/sqrt( q(4) )*qx(4)
      lamda_y = 0.5/sqrt( o4 )*o4y + 0.5/sqrt( q(4) )*qy(4)

      ! o4
      o5  = h_length*lamda*e_shock
      o5x = h_length*( lamda_x*e_shock + lamda*e_shock_x )
      o5y = h_length*( lamda_y*e_shock + lamda*e_shock_y )

      ! fx_shock
      fx_shock(1) = o5x*  qx(1)                     + o5*  qxx(1)
      fx_shock(2) = o5x*( qx(1)*q(2) + q(1)*qx(2) ) + o5*( qxx(1)*q(2) + qx(1)*qx(2) + qx(1)*qx(2) + q(1)*qxx(2) )
      fx_shock(3) = o5x*( qx(1)*q(3) + q(1)*qx(3) ) + o5*( qxx(1)*q(3) + qx(1)*qx(3) + qx(1)*qx(3) + q(1)*qxx(3) )
      fx_shock(4) = o5x*( qx(1)*h    + q(1)*hx    ) + o5*( qxx(1)*h    + qx(1)*hx    + qx(1)*hx    + q(1)*hxx    )

      ! gy_shock
      gy_shock(1) = o5y*  qy(1)                     + o5*  qyy(1)
      gy_shock(2) = o5y*( qy(1)*q(2) + q(1)*qy(2) ) + o5*( qyy(1)*q(2) + qy(1)*qy(2) + qy(1)*qy(2) + q(1)*qyy(2) )
      gy_shock(3) = o5y*( qy(1)*q(3) + q(1)*qy(3) ) + o5*( qyy(1)*q(3) + qy(1)*qy(3) + qy(1)*qy(3) + q(1)*qyy(3) )
      gy_shock(4) = o5y*( qy(1)*h    + q(1)*hy    ) + o5*( qyy(1)*h    + qy(1)*hy    + qy(1)*hy    + q(1)*hyy    )

end subroutine shock_capture_flux_fx_gy_real
!=============================================================================
!
! added dissipative flux for shock capturing, 
!   following Ryan Glasby and Erwin Taylor's paper on COFFE
!
!=============================================================================
subroutine shock_capture_flux_fx_gy_dtype ( fx_shock, gy_shock, q, qx, qy, qxx, qxy, qyy, h_length, k_shock )
      implicit none
      type(dType), dimension(4), intent(out) :: fx_shock, gy_shock
      type(dType), dimension(4), intent(in)  :: q, qx, qy, qxx, qxy, qyy
      real(DP),                  intent(in)  :: h_length
      real(DP),                  intent(in)  :: k_shock

      type(dType) :: p, px, py, pxx, pxy, pyy
      type(dType) :: ek, ekx, eky, ekxx, ekyy
      type(dType) :: h, hx, hy, hxx, hyy
      type(dType) :: e_shock, e_shock_x, e_shock_y
      type(dType) :: lamda, lamda_x, lamda_y
      type(dType) :: o1, o1x, o1y
      type(dType) :: o2, o2x, o2y
      type(dType) :: o3, o3x, o3y
      type(dType) :: o4, o4x, o4y
      type(dType) :: o5, o5x, o5y


      ! p
      p   = gi*q(1)*q(4)
      px  = gi*( qx(1)*q(4) + q(1)*qx(4) )
      py  = gi*( qy(1)*q(4) + q(1)*qy(4) )
      pxx = gi*( qxx(1)*q(4) + qx(1)*qx(4) + qx(1)*qx(4) + q(1)*qxx(4) )
      pxy = gi*( qxy(1)*q(4) + qx(1)*qy(4) + qy(1)*qx(4) + q(1)*qxy(4) )
      pyy = gi*( qyy(1)*q(4) + qy(1)*qy(4) + qy(1)*qy(4) + q(1)*qyy(4) )

      ! ek
      ek   = 0.5*( q(2)**2 + q(3)**2 )
      ekx  = q(2)*qx(2) + q(3)*qx(3)
      eky  = q(2)*qy(2) + q(3)*qy(3)
      ekxx = qx(2)*qx(2) + q(2)*qxx(2) + qx(3)*qx(3) + q(3)*qxx(3)
      ekyy = qy(2)*qy(2) + q(2)*qyy(2) + qy(3)*qy(3) + q(3)*qyy(3)

      ! h
      h   = gm1i*q  (4) + ek
      hx  = gm1i*qx (4) + ekx
      hy  = gm1i*qy (4) + eky
      hxx = gm1i*qxx(4) + ekxx
      hyy = gm1i*qyy(4) + ekyy

      ! o1 = px**2 + py**2
      o1  = px**2 + py**2
      o1x = 2.*( px*pxx + py*pxy )
      o1y = 2.*( px*pxy + py*pyy )

      ! o2
      o2  = h_length*sqrt(o1)
      o2x = h_length*0.5/sqrt(o1)*o1x
      o2y = h_length*0.5/sqrt(o1)*o1y

      ! o3
      o3  = o2  + k_shock*p
      o3x = o2x + k_shock*px
      o3y = o2y + k_shock*py

      ! e_shock
      e_shock = o2/o3
      e_shock_x = ( o2x*o3 - o2*o3x )/o3**2
      e_shock_y = ( o2y*o3 - o2*o3y )/o3**2

      ! o4
      o4 = q(2)**2 + q(3)**2
      o4x = 2.*( q(2)*qx(2) + q(3)*qx(3) )
      o4y = 2.*( q(2)*qy(2) + q(3)*qy(3) )

      ! lamda
      lamda = sqrt( o4 )  + sqrt( q(4) )
      lamda_x = 0.5/sqrt( o4 )*o4x + 0.5/sqrt( q(4) )*qx(4)
      lamda_y = 0.5/sqrt( o4 )*o4y + 0.5/sqrt( q(4) )*qy(4)

      ! o4
      o5  = h_length*lamda*e_shock
      o5x = h_length*( lamda_x*e_shock + lamda*e_shock_x )
      o5y = h_length*( lamda_y*e_shock + lamda*e_shock_y )

      ! fx_shock
      fx_shock(1) = o5x*  qx(1)                     + o5*  qxx(1)
      fx_shock(2) = o5x*( qx(1)*q(2) + q(1)*qx(2) ) + o5*( qxx(1)*q(2) + qx(1)*qx(2) + qx(1)*qx(2) + q(1)*qxx(2) )
      fx_shock(3) = o5x*( qx(1)*q(3) + q(1)*qx(3) ) + o5*( qxx(1)*q(3) + qx(1)*qx(3) + qx(1)*qx(3) + q(1)*qxx(3) )
      fx_shock(4) = o5x*( qx(1)*h    + q(1)*hx    ) + o5*( qxx(1)*h    + qx(1)*hx    + qx(1)*hx    + q(1)*hxx    )

      ! gy_shock
      gy_shock(1) = o5y*  qy(1)                     + o5*  qyy(1)
      gy_shock(2) = o5y*( qy(1)*q(2) + q(1)*qy(2) ) + o5*( qyy(1)*q(2) + qy(1)*qy(2) + qy(1)*qy(2) + q(1)*qyy(2) )
      gy_shock(3) = o5y*( qy(1)*q(3) + q(1)*qy(3) ) + o5*( qyy(1)*q(3) + qy(1)*qy(3) + qy(1)*qy(3) + q(1)*qyy(3) )
      gy_shock(4) = o5y*( qy(1)*h    + q(1)*hy    ) + o5*( qyy(1)*h    + qy(1)*hy    + qy(1)*hy    + q(1)*hyy    )

end subroutine shock_capture_flux_fx_gy_dtype
!=============================================================================
!
! added dissipative flux for shock capturing, 
!   following Ryan Glasby and Erwin Taylor's paper on COFFE,
!   but modified artificial viscousity to avoid singularity when pressure is uniform
!
!=============================================================================
subroutine shock_capture_flux_f_g_modified_real ( f_shock, g_shock, q, qx, qy, h_length, k_shock, epsilon )
      implicit none
      real(DP), dimension(4), intent(out) :: f_shock, g_shock
      real(DP), dimension(4), intent(in)  :: q, qx, qy
      real(DP),               intent(in)  :: h_length
      real(DP),               intent(in)  :: k_shock
      real(DP),               intent(in)  :: epsilon

      real(DP) :: p, px, py
      real(DP) :: ek, ekx, eky
      real(DP) :: h, hx, hy
      real(DP) :: px_y_epsilon
      real(DP) :: e_shock
      real(DP) :: lamda
      real(DP) :: o1
      real(DP) :: o2
      real(DP) :: o3
      real(DP) :: o4
      real(DP) :: o5
      real(DP) :: o6


      ! p
      p  = gi*q(1)*q(4)
      px = gi*( qx(1)*q(4) + q(1)*qx(4) )
      py = gi*( qy(1)*q(4) + q(1)*qy(4) )

      ! ek
      ek  = 0.5*( q(2)**2 + q(3)**2 )
      ekx = q(2)*qx(2) + q(3)*qx(3)
      eky = q(2)*qy(2) + q(3)*qy(3)

      ! h
      h   = gm1i*q(4) + ek
      hx  = gm1i*qx(4) + ekx
      hy  = gm1i*qy(4) + eky

      ! h_length*( abs(px) + abs(py) )
      o1 = px
      o2 = py

      if( o1 < 0. )  o1 = - o1
      if( o2 < 0. )  o2 = - o2

      ! smoothing
      px_y_epsilon = epsilon*p/h_length
      if( o1 < px_y_epsilon )  o1 = ( 2.*o1**2 - o1**3/px_y_epsilon )/px_y_epsilon
      if( o2 < px_y_epsilon )  o2 = ( 2.*o2**2 - o2**3/px_y_epsilon )/px_y_epsilon

      o3  = h_length*( o1 + o2 )

      ! o4
      o4  = o3  + k_shock*p

      ! e_shock
      e_shock = o3/o4

      ! o4
      o5 = q(2)**2 + q(3)**2

      ! lamda
      lamda = sqrt( o5 )  + sqrt( q(4) )

      ! o5
      o6  = h_length*lamda*e_shock

      ! f_shock
      f_shock(1) = o6*  qx(1)
      f_shock(2) = o6*( qx(1)*q(2) + q(1)*qx(2) )
      f_shock(3) = o6*( qx(1)*q(3) + q(1)*qx(3) )
      f_shock(4) = o6*( qx(1)*h    + q(1)*hx    )

      ! g_shock
      g_shock(1) = o6*  qy(1)
      g_shock(2) = o6*( qy(1)*q(2) + q(1)*qy(2) )
      g_shock(3) = o6*( qy(1)*q(3) + q(1)*qy(3) )
      g_shock(4) = o6*( qy(1)*h    + q(1)*hy    )

end subroutine shock_capture_flux_f_g_modified_real
!=============================================================================
!
! added dissipative flux for shock capturing, 
!   following Ryan Glasby and Erwin Taylor's paper on COFFE,
!   but modified artificial viscousity to avoid singularity when pressure is uniform
!
!=============================================================================
subroutine shock_capture_flux_f_g_modified_dtype ( f_shock, g_shock, q, qx, qy, h_length, k_shock, epsilon )
      implicit none
      type(dType), dimension(4), intent(out) :: f_shock, g_shock
      type(dType), dimension(4), intent(in)  :: q, qx, qy
      real(DP),                  intent(in)  :: h_length
      real(DP),                  intent(in)  :: k_shock
      real(DP),                  intent(in)  :: epsilon

      type(dType) :: p, px, py
      type(dType) :: ek, ekx, eky
      type(dType) :: h, hx, hy
      real(DP)    :: px_y_epsilon
      type(dType) :: e_shock
      type(dType) :: lamda
      type(dType) :: o1
      type(dType) :: o2
      type(dType) :: o3
      type(dType) :: o4
      type(dType) :: o5
      type(dType) :: o6

      ! p
      p  = gi*q(1)*q(4)
      px = gi*( qx(1)*q(4) + q(1)*qx(4) )
      py = gi*( qy(1)*q(4) + q(1)*qy(4) )

      ! ek
      ek  = 0.5*( q(2)**2 + q(3)**2 )
      ekx = q(2)*qx(2) + q(3)*qx(3)
      eky = q(2)*qy(2) + q(3)*qy(3)

      ! h
      h   = gm1i*q(4) + ek
      hx  = gm1i*qx(4) + ekx
      hy  = gm1i*qy(4) + eky

      ! h_length*( abs(px) + abs(py) )
      o1 = px
      o2 = py

      if( o1 < 0. )  o1 = - o1
      if( o2 < 0. )  o2 = - o2

      ! smoothing
      px_y_epsilon = epsilon*p/h_length
      if( o1 < px_y_epsilon )  o1 = ( 2.*o1**2 - o1**3/px_y_epsilon )/px_y_epsilon
      if( o2 < px_y_epsilon )  o2 = ( 2.*o2**2 - o2**3/px_y_epsilon )/px_y_epsilon

      o3  = h_length*( o1 + o2 )

      ! o4
      o4  = o3  + k_shock*p

      ! e_shock
      e_shock = o3/o4

      ! o4
      o5 = q(2)**2 + q(3)**2

      ! lamda
      lamda = sqrt( o5 )  + sqrt( q(4) )

      ! o5
      o6  = h_length*lamda*e_shock

      ! f_shock
      f_shock(1) = o6*  qx(1)
      f_shock(2) = o6*( qx(1)*q(2) + q(1)*qx(2) )
      f_shock(3) = o6*( qx(1)*q(3) + q(1)*qx(3) )
      f_shock(4) = o6*( qx(1)*h    + q(1)*hx    )

      ! g_shock
      g_shock(1) = o6*  qy(1)
      g_shock(2) = o6*( qy(1)*q(2) + q(1)*qy(2) )
      g_shock(3) = o6*( qy(1)*q(3) + q(1)*qy(3) )
      g_shock(4) = o6*( qy(1)*h    + q(1)*hy    )


end subroutine shock_capture_flux_f_g_modified_dtype
!=============================================================================
!
! added dissipative flux for shock capturing, 
!   following Ryan Glasby and Erwin Taylor's paper on COFFE,
!   but modified artificial viscousity to avoid singularity when pressure is uniform
!
!=============================================================================
subroutine shock_capture_flux_fx_gy_modified_real ( fx_shock, gy_shock, q, qx, qy, qxx, qxy, qyy, h_length, k_shock, epsilon )
      implicit none
      real(DP), dimension(4), intent(out) :: fx_shock, gy_shock
      real(DP), dimension(4), intent(in)  :: q, qx, qy, qxx, qxy, qyy
      real(DP),               intent(in)  :: h_length
      real(DP),               intent(in)  :: k_shock
      real(DP),               intent(in)  :: epsilon

      real(DP) :: p, px, py, pxx, pxy, pyy
      real(DP) :: ek, ekx, eky, ekxx, ekyy
      real(DP) :: h, hx, hy, hxx, hyy
      real(DP) :: px_y_epsilon
      real(DP) :: e_shock, e_shock_x, e_shock_y
      real(DP) :: lamda, lamda_x, lamda_y
      real(DP) :: o1, o1x, o1y
      real(DP) :: o2, o2x, o2y
      real(DP) :: o3, o3x, o3y
      real(DP) :: o4, o4x, o4y
      real(DP) :: o5, o5x, o5y
      real(DP) :: o6, o6x, o6y
      real(DP) :: osave


      ! p
      p   = gi*q(1)*q(4)
      px  = gi*( qx(1)*q(4) + q(1)*qx(4) )
      py  = gi*( qy(1)*q(4) + q(1)*qy(4) )
      pxx = gi*( qxx(1)*q(4) + qx(1)*qx(4) + qx(1)*qx(4) + q(1)*qxx(4) )
      pxy = gi*( qxy(1)*q(4) + qx(1)*qy(4) + qy(1)*qx(4) + q(1)*qxy(4) )
      pyy = gi*( qyy(1)*q(4) + qy(1)*qy(4) + qy(1)*qy(4) + q(1)*qyy(4) )

      ! ek
      ek   = 0.5*( q(2)**2 + q(3)**2 )
      ekx  = q(2)*qx(2) + q(3)*qx(3)
      eky  = q(2)*qy(2) + q(3)*qy(3)
      ekxx = qx(2)*qx(2) + q(2)*qxx(2) + qx(3)*qx(3) + q(3)*qxx(3)
      ekyy = qy(2)*qy(2) + q(2)*qyy(2) + qy(3)*qy(3) + q(3)*qyy(3)

      ! h
      h   = gm1i*q  (4) + ek
      hx  = gm1i*qx (4) + ekx
      hy  = gm1i*qy (4) + eky
      hxx = gm1i*qxx(4) + ekxx
      hyy = gm1i*qyy(4) + ekyy

      ! h_length*( abs(px) + abs(py) )
      o1  = px
      o1x = pxx
      o1y = pxy

      o2  = py
      o2x = pxy
      o2y = pyy

      if( o1 < 0. ) then
        o1  = - o1
        o1x = - o1x
        o1y = - o1y
      end if

      if( o2 < 0. ) then
        o2  = - o2
        o2x = - o2x
        o2y = - o2y
      end if

      ! smoothing
      px_y_epsilon = epsilon*p/h_length

      if( o1 < px_y_epsilon ) then
        osave = o1
        o1  =     ( 2.*osave**2 -    osave**3/px_y_epsilon )/px_y_epsilon
        o1x = o1x*( 4.*osave    - 3.*osave**2/px_y_epsilon )/px_y_epsilon
        o1y = o1y*( 4.*osave    - 3.*osave**2/px_y_epsilon )/px_y_epsilon
      end if

      if( o2 < px_y_epsilon ) then
        osave = o2
        o2  =     ( 2.*osave**2 -    osave**3/px_y_epsilon )/px_y_epsilon
        o2x = o2x*( 4.*osave    - 3.*osave**2/px_y_epsilon )/px_y_epsilon
        o2y = o2y*( 4.*osave    - 3.*osave**2/px_y_epsilon )/px_y_epsilon
      end if

      o3  = h_length*( o1  + o2  )
      o3x = h_length*( o1x + o2x )
      o3y = h_length*( o1y + o2y )

      ! o4
      o4  = o3  + k_shock*p
      o4x = o3x + k_shock*px
      o4y = o3y + k_shock*py

      ! e_shock
      e_shock = o3/o4
      e_shock_x = ( o3x*o4 - o3*o4x )/o4**2
      e_shock_y = ( o3y*o4 - o3*o4y )/o4**2

      ! o4
      o5 = q(2)**2 + q(3)**2
      o5x = 2.*( q(2)*qx(2) + q(3)*qx(3) )
      o5y = 2.*( q(2)*qy(2) + q(3)*qy(3) )

      ! lamda
      lamda = sqrt( o5 )  + sqrt( q(4) )
      lamda_x = 0.5/sqrt( o5 )*o5x + 0.5/sqrt( q(4) )*qx(4)
      lamda_y = 0.5/sqrt( o5 )*o5y + 0.5/sqrt( q(4) )*qy(4)

      ! o4
      o6  = h_length*lamda*e_shock
      o6x = h_length*( lamda_x*e_shock + lamda*e_shock_x )
      o6y = h_length*( lamda_y*e_shock + lamda*e_shock_y )

      ! fx_shock
      fx_shock(1) = o6x*  qx(1)                     + o6*  qxx(1)
      fx_shock(2) = o6x*( qx(1)*q(2) + q(1)*qx(2) ) + o6*( qxx(1)*q(2) + qx(1)*qx(2) + qx(1)*qx(2) + q(1)*qxx(2) )
      fx_shock(3) = o6x*( qx(1)*q(3) + q(1)*qx(3) ) + o6*( qxx(1)*q(3) + qx(1)*qx(3) + qx(1)*qx(3) + q(1)*qxx(3) )
      fx_shock(4) = o6x*( qx(1)*h    + q(1)*hx    ) + o6*( qxx(1)*h    + qx(1)*hx    + qx(1)*hx    + q(1)*hxx    )

      ! gy_shock
      gy_shock(1) = o6y*  qy(1)                     + o6*  qyy(1)
      gy_shock(2) = o6y*( qy(1)*q(2) + q(1)*qy(2) ) + o6*( qyy(1)*q(2) + qy(1)*qy(2) + qy(1)*qy(2) + q(1)*qyy(2) )
      gy_shock(3) = o6y*( qy(1)*q(3) + q(1)*qy(3) ) + o6*( qyy(1)*q(3) + qy(1)*qy(3) + qy(1)*qy(3) + q(1)*qyy(3) )
      gy_shock(4) = o6y*( qy(1)*h    + q(1)*hy    ) + o6*( qyy(1)*h    + qy(1)*hy    + qy(1)*hy    + q(1)*hyy    )

end subroutine shock_capture_flux_fx_gy_modified_real
!=============================================================================
!
! added dissipative flux for shock capturing, 
!   following Ryan Glasby and Erwin Taylor's paper on COFFE
!   but modified artificial viscousity to avoid singularity when pressure is uniform
!
!=============================================================================
subroutine shock_capture_flux_fx_gy_modified_dtype ( fx_shock, gy_shock, q, qx, qy, qxx, qxy, qyy, h_length, k_shock, epsilon )
      implicit none
      type(dType), dimension(4), intent(out) :: fx_shock, gy_shock
      type(dType), dimension(4), intent(in)  :: q, qx, qy, qxx, qxy, qyy
      real(DP),                  intent(in)  :: h_length
      real(DP),                  intent(in)  :: k_shock
      real(DP),                  intent(in)  :: epsilon

      type(dType) :: p, px, py, pxx, pxy, pyy
      type(dType) :: ek, ekx, eky, ekxx, ekyy
      type(dType) :: h, hx, hy, hxx, hyy
      real(DP)    :: px_y_epsilon
      type(dType) :: e_shock, e_shock_x, e_shock_y
      type(dType) :: lamda, lamda_x, lamda_y
      type(dType) :: o1, o1x, o1y
      type(dType) :: o2, o2x, o2y
      type(dType) :: o3, o3x, o3y
      type(dType) :: o4, o4x, o4y
      type(dType) :: o5, o5x, o5y
      type(dType) :: o6, o6x, o6y
      type(dType) :: osave


      ! p
      p   = gi*q(1)*q(4)
      px  = gi*( qx(1)*q(4) + q(1)*qx(4) )
      py  = gi*( qy(1)*q(4) + q(1)*qy(4) )
      pxx = gi*( qxx(1)*q(4) + qx(1)*qx(4) + qx(1)*qx(4) + q(1)*qxx(4) )
      pxy = gi*( qxy(1)*q(4) + qx(1)*qy(4) + qy(1)*qx(4) + q(1)*qxy(4) )
      pyy = gi*( qyy(1)*q(4) + qy(1)*qy(4) + qy(1)*qy(4) + q(1)*qyy(4) )

      ! ek
      ek   = 0.5*( q(2)**2 + q(3)**2 )
      ekx  = q(2)*qx(2) + q(3)*qx(3)
      eky  = q(2)*qy(2) + q(3)*qy(3)
      ekxx = qx(2)*qx(2) + q(2)*qxx(2) + qx(3)*qx(3) + q(3)*qxx(3)
      ekyy = qy(2)*qy(2) + q(2)*qyy(2) + qy(3)*qy(3) + q(3)*qyy(3)

      ! h
      h   = gm1i*q  (4) + ek
      hx  = gm1i*qx (4) + ekx
      hy  = gm1i*qy (4) + eky
      hxx = gm1i*qxx(4) + ekxx
      hyy = gm1i*qyy(4) + ekyy

      ! h_length*( abs(px) + abs(py) )
      o1  = px
      o1x = pxx
      o1y = pxy

      o2  = py
      o2x = pxy
      o2y = pyy

      if( o1 < 0. ) then
        o1  = - o1
        o1x = - o1x
        o1y = - o1y
      end if

      if( o2 < 0. ) then
        o2  = - o2
        o2x = - o2x
        o2y = - o2y
      end if

      ! smoothing
      px_y_epsilon = epsilon*p/h_length

      if( o1 < px_y_epsilon ) then
        osave = o1
        o1  =     ( 2.*osave**2 -    osave**3/px_y_epsilon )/px_y_epsilon
        o1x = o1x*( 4.*osave    - 3.*osave**2/px_y_epsilon )/px_y_epsilon
        o1y = o1y*( 4.*osave    - 3.*osave**2/px_y_epsilon )/px_y_epsilon
      end if

      if( o2 < px_y_epsilon ) then
        osave = o2
        o2  =     ( 2.*osave**2 -    osave**3/px_y_epsilon )/px_y_epsilon
        o2x = o2x*( 4.*osave    - 3.*osave**2/px_y_epsilon )/px_y_epsilon
        o2y = o2y*( 4.*osave    - 3.*osave**2/px_y_epsilon )/px_y_epsilon
      end if

      o3  = h_length*( o1  + o2  )
      o3x = h_length*( o1x + o2x )
      o3y = h_length*( o1y + o2y )

      ! o4
      o4  = o3  + k_shock*p
      o4x = o3x + k_shock*px
      o4y = o3y + k_shock*py

      ! e_shock
      e_shock = o3/o4
      e_shock_x = ( o3x*o4 - o3*o4x )/o4**2
      e_shock_y = ( o3y*o4 - o3*o4y )/o4**2

      ! o4
      o5 = q(2)**2 + q(3)**2
      o5x = 2.*( q(2)*qx(2) + q(3)*qx(3) )
      o5y = 2.*( q(2)*qy(2) + q(3)*qy(3) )

      ! lamda
      lamda = sqrt( o5 )  + sqrt( q(4) )
      lamda_x = 0.5/sqrt( o5 )*o5x + 0.5/sqrt( q(4) )*qx(4)
      lamda_y = 0.5/sqrt( o5 )*o5y + 0.5/sqrt( q(4) )*qy(4)

      ! o4
      o6  = h_length*lamda*e_shock
      o6x = h_length*( lamda_x*e_shock + lamda*e_shock_x )
      o6y = h_length*( lamda_y*e_shock + lamda*e_shock_y )

      ! fx_shock
      fx_shock(1) = o6x*  qx(1)                     + o6*  qxx(1)
      fx_shock(2) = o6x*( qx(1)*q(2) + q(1)*qx(2) ) + o6*( qxx(1)*q(2) + qx(1)*qx(2) + qx(1)*qx(2) + q(1)*qxx(2) )
      fx_shock(3) = o6x*( qx(1)*q(3) + q(1)*qx(3) ) + o6*( qxx(1)*q(3) + qx(1)*qx(3) + qx(1)*qx(3) + q(1)*qxx(3) )
      fx_shock(4) = o6x*( qx(1)*h    + q(1)*hx    ) + o6*( qxx(1)*h    + qx(1)*hx    + qx(1)*hx    + q(1)*hxx    )

      ! gy_shock
      gy_shock(1) = o6y*  qy(1)                     + o6*  qyy(1)
      gy_shock(2) = o6y*( qy(1)*q(2) + q(1)*qy(2) ) + o6*( qyy(1)*q(2) + qy(1)*qy(2) + qy(1)*qy(2) + q(1)*qyy(2) )
      gy_shock(3) = o6y*( qy(1)*q(3) + q(1)*qy(3) ) + o6*( qyy(1)*q(3) + qy(1)*qy(3) + qy(1)*qy(3) + q(1)*qyy(3) )
      gy_shock(4) = o6y*( qy(1)*h    + q(1)*hy    ) + o6*( qyy(1)*h    + qy(1)*hy    + qy(1)*hy    + q(1)*hyy    )

end subroutine shock_capture_flux_fx_gy_modified_dtype

end module Euler_mod
