module FEMEntry_Quadrature_Gauss_mod
      use FEMEntry_Generic_mod
      use KindDefinition_mod, only : DP
      implicit none
      save

! - Lagranage mapping of element and edge
      type, extends (quadrature_2d_t) :: gauss_quadrature_2d_t
      contains
        procedure, pass(this) :: Reset_quadrature_2d => Reset_gauss_quadrature_2d
      end type gauss_quadrature_2d_t

      type, extends (quadrature_1d_t) :: gauss_quadrature_1d_t
      contains
        procedure, pass(this) :: Reset_quadrature_1d => Reset_gauss_quadrature_1d
      end type gauss_quadrature_1d_t

contains

!=====================================================================
!
! a way to determine (r,s,w) of each quadrature according to Pavel Solin's finite-element book:
!   p: polynomial order, (p+1)*(p+2)/2 is the number of terms of the complete polynomials of order p
!   N: number of quadratures, each quadrature has 3 unknow (r,s,w)
! number of unknown should be no less than number of equations, which is also the number of terms of the complete polynomials
!   3*N >= (p+1)*(p+2)/2
! using the criterion:
!   p          :   1;   2;   3;   4;   5;   6;   7;   8;   9;
!   N (minimum):   1;   2;   4;   5;   7;  10;  12;  15;  17;
!   N (actual) :   1;   3;   4;   7;   7;  12;  12;  16;  16;
!
!=====================================================================
subroutine Reset_gauss_quadrature_2d( this, polynomial_order )
      implicit none
      class(gauss_quadrature_2d_t), intent(inout) :: this
      integer,                      intent(in)    :: polynomial_order

      this%polynomial_order = polynomial_order

      select case( polynomial_order )
      case(0:1)
        this%nqd = 1
        this%rqd(1) = 1./3.;  this%sqd(1) = 1./3.;  this%wqd(1) = 1.
      case(2)
        this%nqd = 3
        this%rqd(1) = 1./6.;  this%sqd(1) = 1./6.;  this%wqd(1) = 1./3.
        this%rqd(2) = 1./6.;  this%sqd(2) = 4./6.;  this%wqd(2) = 1./3.
        this%rqd(3) = 4./6.;  this%sqd(3) = 1./6.;  this%wqd(3) = 1./3.
      case(3)
        this%nqd = 4
        this%rqd(1) = 1./3.;  this%sqd(1) = 1./3.;  this%wqd(1) = -0.5625
        this%rqd(2) = 0.6  ;  this%sqd(2) = 0.2  ;  this%wqd(2) =  0.52083333333333333
        this%rqd(3) = 0.2  ;  this%sqd(3) = 0.6  ;  this%wqd(3) =  0.52083333333333333
        this%rqd(4) = 0.2  ;  this%sqd(4) = 0.2  ;  this%wqd(4) =  0.52083333333333333
      case(4:5)
        this%nqd = 7
        this%rqd(1) = 1./3.            ;  this%sqd(1) = 1./3.            ;  this%wqd(1) = 0.225
        this%rqd(2) = 0.797426985353087;  this%sqd(2) = 0.101286507323456;  this%wqd(2) = 0.125939180544827
        this%rqd(3) = 0.101286507323456;  this%sqd(3) = 0.797426985353087;  this%wqd(3) = 0.125939180544827
        this%rqd(4) = 0.101286507323456;  this%sqd(4) = 0.101286507323456;  this%wqd(4) = 0.125939180544827
        this%rqd(5) = 0.059715871789770;  this%sqd(5) = 0.470142064105115;  this%wqd(5) = 0.132394152788506
        this%rqd(6) = 0.470142064105115;  this%sqd(6) = 0.059715871789770;  this%wqd(6) = 0.132394152788506
        this%rqd(7) = 0.470142064105115;  this%sqd(7) = 0.470142064105115;  this%wqd(7) = 0.132394152788506
      case(6:7)
        this%nqd = 12
        this%rqd( 1) = 0.2492867451709100;  this%sqd( 1) = 0.2492867451709100;  this%wqd( 1) = 0.1167862757263788
        this%rqd( 2) = 0.2492867451709100;  this%sqd( 2) = 0.5014265096581790;  this%wqd( 2) = 0.1167862757263788
        this%rqd( 3) = 0.5014265096581790;  this%sqd( 3) = 0.2492867451709100;  this%wqd( 3) = 0.1167862757263788
        this%rqd( 4) = 0.0630890144915020;  this%sqd( 4) = 0.0630890144915020;  this%wqd( 4) = 0.0508449063702069
        this%rqd( 5) = 0.0630890144915020;  this%sqd( 5) = 0.8738219710169960;  this%wqd( 5) = 0.0508449063702069
        this%rqd( 6) = 0.8738219710169960;  this%sqd( 6) = 0.0630890144915020;  this%wqd( 6) = 0.0508449063702069
        this%rqd( 7) = 0.3103524510337840;  this%sqd( 7) = 0.6365024991213990;  this%wqd( 7) = 0.0828510756183738
        this%rqd( 8) = 0.6365024991213990;  this%sqd( 8) = 0.0531450498448170;  this%wqd( 8) = 0.0828510756183738
        this%rqd( 9) = 0.0531450498448170;  this%sqd( 9) = 0.3103524510337840;  this%wqd( 9) = 0.0828510756183738
        this%rqd(10) = 0.3103524510337840;  this%sqd(10) = 0.0531450498448170;  this%wqd(10) = 0.0828510756183738
        this%rqd(11) = 0.6365024991213990;  this%sqd(11) = 0.3103524510337840;  this%wqd(11) = 0.0828510756183738
        this%rqd(12) = 0.0531450498448170;  this%sqd(12) = 0.6365024991213990;  this%wqd(12) = 0.0828510756183738
     !case(8:9)
      case(8: )
        this%nqd = 16
        this%rqd( 1) = 0.3333333333333335;  this%sqd( 1) = 0.3333333333333335;  this%wqd( 1) = 0.1443156076777871
        this%rqd( 2) = 0.4592925882927230;  this%sqd( 2) = 0.4592925882927230;  this%wqd( 2) = 0.0950916342672850
        this%rqd( 3) = 0.4592925882927230;  this%sqd( 3) = 0.0814148234145540;  this%wqd( 3) = 0.0950916342672850
        this%rqd( 4) = 0.0814148234145540;  this%sqd( 4) = 0.4592925882927230;  this%wqd( 4) = 0.0950916342672850
        this%rqd( 5) = 0.1705693077517600;  this%sqd( 5) = 0.1705693077517600;  this%wqd( 5) = 0.1032173705347180
        this%rqd( 6) = 0.1705693077517600;  this%sqd( 6) = 0.6588613844964800;  this%wqd( 6) = 0.1032173705347180
        this%rqd( 7) = 0.6588613844964800;  this%sqd( 7) = 0.1705693077517600;  this%wqd( 7) = 0.1032173705347180
        this%rqd( 8) = 0.0505472283170310;  this%sqd( 8) = 0.0505472283170310;  this%wqd( 8) = 0.0324584976231980
        this%rqd( 9) = 0.0505472283170310;  this%sqd( 9) = 0.8989055433659380;  this%wqd( 9) = 0.0324584976231975
        this%rqd(10) = 0.8989055433659380;  this%sqd(10) = 0.0505472283170310;  this%wqd(10) = 0.0324584976231980
        this%rqd(11) = 0.2631128296346380;  this%sqd(11) = 0.7284923929554040;  this%wqd(11) = 0.0272303141744350
        this%rqd(12) = 0.7284923929554040;  this%sqd(12) = 0.0083947774099580;  this%wqd(12) = 0.0272303141744350
        this%rqd(13) = 0.0083947774099580;  this%sqd(13) = 0.2631128296346380;  this%wqd(13) = 0.0272303141744350
        this%rqd(14) = 0.2631128296346380;  this%sqd(14) = 0.0083947774099580;  this%wqd(14) = 0.0272303141744350
        this%rqd(15) = 0.7284923929554040;  this%sqd(15) = 0.2631128296346380;  this%wqd(15) = 0.0272303141744350
        this%rqd(16) = 0.0083947774099580;  this%sqd(16) = 0.7284923929554040;  this%wqd(16) = 0.0272303141744350
      case default
        write(*,*) 'wrong! gauss_quadrature_2d_t%Reset_gauss_qudrature_2d'
        stop
      end select

end subroutine Reset_gauss_quadrature_2d
!=====================================================================
!
! a way to determine (r,w) of each quadrature according to Pavel Solin's finite-element book:
!   p: polynomial order, p+1 is the number of terms of the complete polynomials of order p
!   N: number of quadratures, each quadrature has 2 unknow (r,w)
! number of unknown should be no less than number of equations, which is also the number of terms of the complete polynomials
!   2*N >= p+1
! using the criterion:
!   p:   1;   2;   3;   4;   5;   6;   7;   8;   9;
!   N:   1;   2;   2;   3;   3;   4;   4;   5;   5;
!
!=====================================================================
subroutine Reset_gauss_quadrature_1d( this, polynomial_order )
      implicit none
      class(gauss_quadrature_1d_t), intent(inout) :: this
      integer,                      intent(in)    :: polynomial_order

      this%polynomial_order = polynomial_order

      select case( polynomial_order )
      case(0:1)
        this%nqd = 1
        this%rqd(1) = 0.5;  this%wqd(1) = 1.0
      case(2:3)
        this%nqd = 2
        this%rqd(1) = 0.5*(1. - sqrt(1./3.));  this%wqd(1) = 0.5
        this%rqd(2) = 0.5*(1. + sqrt(1./3.));  this%wqd(2) = 0.5
      case(4:5)
        this%nqd = 3
        this%rqd(1) = 0.5*(1. - sqrt(3./5.));  this%wqd(1) = 5./18.
        this%rqd(2) = 0.5                   ;  this%wqd(2) = 8./18.
        this%rqd(3) = 0.5*(1. + sqrt(3./5.));  this%wqd(3) = 5./18.
      case(6: )
     !case(6:7)
        this%nqd = 4
        this%rqd(1) = 0.0694318442029738;  this%wqd(1) = 0.1739274225687269
        this%rqd(2) = 0.3300094782075718;  this%wqd(2) = 0.3260725774312730
        this%rqd(3) = 0.6699905217924282;  this%wqd(3) = 0.3260725774312730
        this%rqd(4) = 0.9305681557970262;  this%wqd(4) = 0.1739274225687269
     !case(8: )
     !  this%nqd = 5
     !  this%rqd(1) = 0.5*( 1. - 1./3.*sqrt( 5. + 2.*sqrt(10./7.) ) );  this%wqd(1) = ( 322. - 13.*sqrt(70.) )/1800.
     !  this%rqd(2) = 0.5*( 1. - 1./3.*sqrt( 5. - 2.*sqrt(10./7.) ) );  this%wqd(2) = ( 322. + 13.*sqrt(70.) )/1800.
     !  this%rqd(3) = 0.5                                            ;  this%wqd(3) = 128./450.
     !  this%rqd(4) = 0.5*( 1. + 1./3.*sqrt( 5. - 2.*sqrt(10./7.) ) );  this%wqd(4) = ( 322. + 13.*sqrt(70.) )/1800.
     !  this%rqd(5) = 0.5*( 1. + 1./3.*sqrt( 5. + 2.*sqrt(10./7.) ) );  this%wqd(5) = ( 322. - 13.*sqrt(70.) )/1800.
      case default
        write(*,*) 'wrong! gauss_quadrature_1d_t%Reset_gauss_qudrature_1d', polynomial_order
        stop
      end select

end subroutine Reset_gauss_quadrature_1d

end module FEMEntry_Quadrature_Gauss_mod
