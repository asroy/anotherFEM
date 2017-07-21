!
! Functions for overloading operators
! Note that we designate a type of dType that has a value and derivatives
! For the variables below, we let "real" be a variable of type dType
! An integer is just an integer but a "const" is a real-valued constant
!
module Common_Overload_mod
      use KindDefinition_mod, only : dp
      implicit none
      save

      integer, parameter:: maxOverloadLength = 120
      integer :: overloadLength = -66

      type :: dType
        real(dp) :: value
        real(dp), dimension(maxOverloadLength) :: deriv
      end type dType

      type :: testType
        real(dp) :: value
        real(dp), dimension(:), allocatable :: deriv
      end type testType

      interface assignment(=)
        module procedure dTypeEqReal
        module procedure RealEqdType
        module procedure dTypeEqInteger
      end interface

      interface operator(*)
        module procedure real_real_multiply
        module procedure integer_real_multiply
        module procedure real_integer_multiply
        module procedure const_real_multiply
        module procedure real_const_multiply
      end interface

      interface operator(/)
        module procedure real_real_divide
        module procedure integer_real_divide
        module procedure real_integer_divide
        module procedure const_real_divide
        module procedure real_const_divide
      end interface

      interface operator(+)
        module procedure real_real_add
        module procedure integer_real_add
        module procedure real_integer_add
        module procedure const_real_add
        module procedure real_const_add
        module procedure keep_sign
      end interface

      interface operator(-)
        module procedure real_real_subtract
        module procedure integer_real_subtract
        module procedure real_integer_subtract
        module procedure const_real_subtract
        module procedure real_const_subtract
        module procedure change_sign
      end interface

      interface operator(**)
        module procedure real_real_exponent
        module procedure integer_real_exponent
        module procedure real_integer_exponent
        module procedure const_real_exponent
        module procedure real_const_exponent
      end interface

      interface operator( < )
        module procedure real_real_lt
        module procedure integer_real_lt
        module procedure real_integer_lt
        module procedure const_real_lt
        module procedure real_const_lt
      end interface

      interface operator( <= )
        module procedure real_real_le
        module procedure integer_real_le
        module procedure real_integer_le
        module procedure const_real_le
        module procedure real_const_le
      end interface

      interface operator( > )
        module procedure real_real_gt
        module procedure integer_real_gt
        module procedure real_integer_gt
        module procedure const_real_gt
        module procedure real_const_gt
      end interface

      interface operator( >= )
        module procedure real_real_ge
        module procedure integer_real_ge
        module procedure real_integer_ge
        module procedure const_real_ge
        module procedure real_const_ge
      end interface

      interface operator( == )
        module procedure real_real_eq
        module procedure integer_real_eq
        module procedure real_integer_eq
        module procedure const_real_eq
        module procedure real_const_eq
      end interface

      interface operator( /= )
        module procedure real_real_ne
        module procedure integer_real_ne
        module procedure real_integer_ne
        module procedure const_real_ne
        module procedure real_const_ne
      end interface

      interface sqrt
        module procedure real_sqrt
      end interface

      interface abs
        module procedure real_abs
      end interface

      interface exp
        module procedure real_exp
      end interface

      interface log
        module procedure real_log
      end interface

      interface max
        module procedure real_real_max
        module procedure real_const_max
        module procedure const_real_max
      end interface

      interface min
        module procedure real_real_min
        module procedure real_const_min
        module procedure const_real_min
      end interface

      interface sin
        module procedure real_sin
      end interface

      interface cos
        module procedure real_cos
      end interface

contains

!======================== Initialize ===================================
!
! Initializes variables of type dType
! If index is zero, all derivatives are zero
! If index is not zero, this implies that the variable is an independent
! variable so that derivative is set to unity
!
! We should only need to call this for the independent variables
!
!=======================================================================
      subroutine initialize(v,index,value)
        type (dType) :: v
        integer :: index
        real(dp) :: value
!
! Check overloadLength
!
      if(overloadLength.eq.-66)then
        write(6,'("overloadLength needs to be set")')
        stop
      end if
      if(overloadLength.gt.maxOverloadLength)then
        write(6,'("overloadLength exceeds the maximum")')
!        write(6,overloadLength)
        print * , overloadLength
        stop
      end if
!
! Set the value and the derivatives
!
        v%value = value
        v%deriv = 0.
        if(index.ne.0)then
          v%deriv(index) = 1.0
        end if
      end subroutine initialize

!============================ dTYPEEQREAL ==============================
!
! Sets a variable of type (dType) to a real number
!
!=======================================================================
        elemental subroutine dTypeEqReal(v1,r1)
          type (dType), intent(out) :: v1
          real(dp), intent(in) :: r1
          v1%value = r1
          v1%deriv(1:overloadLength) = 0.
        end subroutine dTypeEqReal

!============================ dTYPEEQREAL ==============================
!
! Sets a variable of type (dType) to a real number
!
!=======================================================================
        elemental subroutine RealEqdType(r1,v1)
          real(dp), intent(out) :: r1
          type (dType), intent(in) :: v1
          r1 = v1%value
        end subroutine RealEqdType

!============================ dTYPEEQINTEGER ===========================
!
! Sets a variable of type (dType) to an integer
!
!=======================================================================
        elemental subroutine dTypeEqInteger(v1,i2)
          type (dType), intent(out) :: v1
          integer, intent(in) :: i2
          v1%value = real(i2)
          v1%deriv(1:overloadLength) = 0.
        end subroutine dTypeEqInteger
!
! MULTIPLY
!
!======================== Real_Real_Multiply ===========================
!
! Multiplies two real numbers v1*v2
! d(v1*v2) = v1*d(v2) + v2*d(v1)
!
!=======================================================================
      elemental function real_real_multiply(v1,v2)
        type (dtype) :: real_real_multiply
        type (dType), intent(in) :: v1,v2
        real_real_multiply%value = v1%value * v2%value
        real_real_multiply%deriv(1:overloadLength) = &
          (v1%value * v2%deriv(1:overloadLength)) &
          + (v2%value * v1%deriv(1:overloadLength))
      end function real_real_multiply

!======================== Integer_Real_Multiply ========================
!
! Multiplies integer*real i1*v2
! d(i1*v2) = i1*d(v2)
!
!=======================================================================
      elemental function integer_real_multiply(i1,v2)
        type (dtype) :: integer_real_multiply
        integer, intent(in) :: i1
        type (dType), intent(in) :: v2
        integer_real_multiply%value = i1 * v2%value
        integer_real_multiply%deriv(1:overloadLength) = (i1 * v2%deriv(1:overloadLength))
      end function integer_real_multiply

!======================== Real_Integer_Multiply ========================
!
! Multiplies real*integer v1*i2
! d(v1*i2) = i2*d(v1)
!
!=======================================================================
      elemental function real_integer_multiply(v1,i2)
        type (dtype) :: real_integer_multiply
        type (dType), intent(in) :: v1
        integer, intent(in) :: i2
        real_integer_multiply%value = v1%value * i2
        real_integer_multiply%deriv(1:overloadLength) = (i2 * v1%deriv(1:overloadLength))
      end function real_integer_multiply

!======================== Const_Real_Multiply ==========================
!
! Multiplies real_constant*real r1*v2
! d(r1*v2) = r1*d(v2)
!
!=======================================================================
      elemental function const_real_multiply(r1,v2)
        type (dtype) :: const_real_multiply
        real(dp), intent(in) :: r1
        type (dType), intent(in) :: v2
        const_real_multiply%value = r1 * v2%value
        const_real_multiply%deriv(1:overloadLength) = (r1 * v2%deriv(1:overloadLength))
      end function const_real_multiply

!======================== Real_Const_Multiply ========================
!
! Multiplies real*real_constant v1*r2
! d(v1*r2) = r2*d(v1)
!
!=======================================================================
      elemental function real_const_multiply(v1,r2)
        type (dtype) :: real_const_multiply
        type (dType), intent(in) :: v1
        real(dp), intent(in) :: r2
        real_const_multiply%value = v1%value * r2
        real_const_multiply%deriv(1:overloadLength) = (r2 * v1%deriv(1:overloadLength))
      end function real_const_multiply
!
! DIVIDE
!
!======================== Real_Real_Divide =============================
!
! Divides two real numbers v1/v2
! d(v1/v2) = (v2*d(v1) - v1*d(v2))/(v2*v2)
!
!=======================================================================
      elemental function real_real_divide(v1,v2)
        type (dtype) :: real_real_divide
        type (dType), intent(in) :: v1,v2
        real_real_divide%value = v1%value / v2%value
        real_real_divide%deriv(1:overloadLength) = &
         ((v2%value * v1%deriv(1:overloadLength)) &
       - (v1%value * v2%deriv(1:overloadLength)))/(v2%value*v2%value)
      end function real_real_divide

!======================== Integer_Real_Divide ==========================
!
! Divides two real numbers i1/v2
! d(i1/v2) =  -i1*d(v2))/(v2*v2)
!
!=======================================================================
      elemental function integer_real_divide(i1,v2)
        type (dtype) :: integer_real_divide
        integer, intent(in) :: i1
        type (dType), intent(in) :: v2
        integer_real_divide%value = real(i1)/ v2%value
        integer_real_divide%deriv(1:overloadLength) =  -(i1 * v2%deriv(1:overloadLength))/(v2%value*v2%value)
      end function integer_real_divide

!======================== Real_Integer_Divide ==========================
!
! Divides two real numbers v1/i2
! d(v1/i2) = (i2*d(v1))/(v2*v2)
!
!=======================================================================
      elemental function real_integer_divide(v1,i2)
        type (dtype) :: real_integer_divide
        type (dType), intent(in) :: v1
        integer, intent(in) :: i2
        real_integer_divide%value = v1%value / real(i2)
        real_integer_divide%deriv(1:overloadLength) = (i2 * v1%deriv(1:overloadLength))/(real(i2)*real(i2))
      end function real_integer_divide

!======================== Const_Real_Divide ============================
!
! Divides real_const/real numbers r1/v2
! d(r1/v2) =  -r1*d(v2))/(v2*v2)
!
!=======================================================================
      elemental function const_real_divide(r1,v2)
        type (dtype) :: const_real_divide
        real(dp), intent(in) :: r1
        type (dType), intent(in) :: v2
        const_real_divide%value = r1/ v2%value
        const_real_divide%deriv(1:overloadLength) =  -(r1 * v2%deriv(1:overloadLength))/(v2%value*v2%value)
      end function const_real_divide

!======================== Real_Const_Divide ============================
!
! Divides real/real_const numbers v1/r2
! d(v1/r2) = (r2*d(v1))/(v2*v2)
!
!=======================================================================
      elemental function real_const_divide(v1,r2)
        type (dtype) :: real_const_divide
        type (dType), intent(in) :: v1
        real(dp), intent(in) :: r2
        real_const_divide%value = v1%value / r2
        real_const_divide%deriv(1:overloadLength) = (r2 * v1%deriv(1:overloadLength))/(r2*r2)
      end function real_const_divide
!
! ADD
!
!======================== Real_Real_Add ================================
!
! Add two real numbers v1 + v2
! d(v1+v2) = d(v1) + d(v2)
!
!=======================================================================
      elemental function real_real_add(v1,v2)
        type (dtype) :: real_real_add
        type (dType), intent(in) :: v1,v2
        real_real_add%value = v1%value + v2%value
        real_real_add%deriv(1:overloadLength) = v1%deriv(1:overloadLength) + v2%deriv(1:overloadLength)
      end function real_real_add

!======================== Integer_Real_Add =============================
!
! Add integer + real i1 + v2
! d(i1+v2) = d(v2)
!
!=======================================================================
      elemental function integer_real_add(i1,v2)
        type (dtype) :: integer_real_add
        integer, intent(in) :: i1
        type (dType), intent(in) :: v2
        integer_real_add%value = real(i1) + v2%value
        integer_real_add%deriv(1:overloadLength) = v2%deriv(1:overloadLength)
      end function integer_real_add

!======================== Real_Integer_Add =============================
!
! Add real + integer v1 + i2
! d(v1 + i2) = d(v1)
!
!=======================================================================
      elemental function real_integer_add(v1,i2)
        type (dtype) :: real_integer_add
        type (dType), intent(in) :: v1
        integer, intent(in) :: i2
        real_integer_add%value = v1%value + real(i2)
        real_integer_add%deriv(1:overloadLength) = v1%deriv(1:overloadLength)
      end function real_integer_add

!======================== Const_Real_Add ===============================
!
! Add real_constant + real r1 + v2
! d(r1+v2) = d(v2)
!
!=======================================================================
      elemental function const_real_add(r1,v2)
        type (dtype) :: const_real_add
        real(dp), intent(in) :: r1
        type (dType), intent(in) :: v2
        const_real_add%value = r1 + v2%value
        const_real_add%deriv(1:overloadLength) = v2%deriv(1:overloadLength)
      end function const_real_add

!======================== Real_Const_Add ===============================
!
! Add real + real_constant v1+v2
! d(v1 + r2) = d(v1)
!
!=======================================================================
      elemental function real_const_add(v1,r2)
        type (dtype) :: real_const_add
        type (dType), intent(in) :: v1
        real(dp), intent(in) :: r2
        real_const_add%value = v1%value + r2
        real_const_add%deriv(1:overloadLength) = v1%deriv(1:overloadLength)
      end function real_const_add

!============================ Keep_Sign ================================
!
! Handles the case where v1 = +v
!
!=======================================================================
      elemental function keep_sign(v1)
        type (dtype) :: keep_sign
        type (dType), intent(in) :: v1
        keep_sign%value = v1%value
        keep_sign%deriv(1:overloadLength) = v1%deriv(1:overloadLength)
      end function keep_sign
!
! SUBTRACT
!
!======================== Real_Real_Subtract ===========================
!
! Subtract two real numbers v1 - v2
! d(v1 - v2) = d(v1) - d(v2)
!
!=======================================================================
      elemental function real_real_subtract(v1,v2)
        type (dtype) :: real_real_subtract
        type (dType), intent(in) :: v1,v2
        real_real_subtract%value = v1%value - v2%value
        real_real_subtract%deriv(1:overloadLength) = v1%deriv(1:overloadLength) - v2%deriv(1:overloadLength)
      end function real_real_subtract

!======================== Integer_Real_Subtract ========================
!
! Subtract integer - real i1 - v2
! d(i1-v2) = -d(v2)
!
!=======================================================================
      elemental function integer_real_subtract(i1,v2)
        type (dtype) :: integer_real_subtract
        integer, intent(in) :: i1
        type (dType), intent(in) :: v2
        integer_real_subtract%value = real(i1) - v2%value
        integer_real_subtract%deriv(1:overloadLength) = -v2%deriv(1:overloadLength)
      end function integer_real_subtract

!======================== Real_Integer_Subtract ========================
!
! Subtract real - integer v1 - i2
! d(v1-i2) = d(v1)
!
!=======================================================================
      elemental function real_integer_subtract(v1,i2)
        type (dtype) :: real_integer_subtract
        type (dType), intent(in) :: v1
        integer, intent(in) :: i2
        real_integer_subtract%value = v1%value - real(i2)
        real_integer_subtract%deriv(1:overloadLength) = v1%deriv(1:overloadLength)
      end function real_integer_subtract

!======================== Const_Real_Subtract ==========================
!
! Subtract real_constant - real r1 - v2
! d(r1-v2) = -d(v2)
!
!=======================================================================
      elemental function const_real_subtract(r1,v2)
        type (dtype) :: const_real_subtract
        real(dp), intent(in) :: r1
        type (dType), intent(in) :: v2
        const_real_subtract%value = r1 - v2%value
        const_real_subtract%deriv(1:overloadLength) = -v2%deriv(1:overloadLength)
      end function const_real_subtract

!======================== Real_Const_Subtract ==========================
!
! Subtract real - real_constant  v1 - r2
! d(v1 - r2) = d(v1)
!
!=======================================================================
      elemental function real_const_subtract(v1,r2)
        type (dtype) :: real_const_subtract
        type (dType), intent(in) :: v1
        real(dp), intent(in) :: r2
        real_const_subtract%value = v1%value - r2
        real_const_subtract%deriv(1:overloadLength) = v1%deriv(1:overloadLength)
      end function real_const_subtract

!============================ Change_Sign ==============================
!
! Subtract real - real_constant  v1 - r2
! d(v1 - r2) = d(v1)
!
!=======================================================================
      elemental function change_sign(v1)
        type (dtype) :: change_sign
        type (dType), intent(in) :: v1
        change_sign%value = -v1%value
        change_sign%deriv(1:overloadLength) = -v1%deriv(1:overloadLength)
      end function change_sign
!
! EXPONENTIAL
!
!======================== Real_Integer_Exponent ========================
!
! Exponent real**integer v1**i2
! d(v1**i2) = i2*v1**(i2-1)*d(v1)
!
!=======================================================================
      elemental function real_integer_exponent(v1,i2)
        type (dtype) :: real_integer_exponent
        type (dType), intent(in) :: v1
        integer, intent(in) :: i2
        real_integer_exponent%value = v1%value**i2
        real_integer_exponent%deriv(1:overloadLength) = i2*v1%value**(i2-1)*v1%deriv(1:overloadLength)
      end function real_integer_exponent

!======================== Integer_Real_Exponent ================================
!
! Exponent integer**real  phi = i1**v2
! Note that if i1 is negative, this is not a defined operation
! This is derived using natural logarithms
! phi = i1**v2
! log(phi) = v2*log(i1)
! 1./phi*d(phi)/d(chi) = log(i1)*d(v2)/d(chi)
! d(phi)/d(chi) = phi*log(i1)*d(v2)/d(chi)
!
!=======================================================================
      elemental function integer_real_exponent(i1,v2)
        type (dtype) :: integer_real_exponent
        type (dType), intent(in) :: v2
        integer, intent(in) :: i1
        integer_real_exponent%value = i1**v2%value
        if(i1.gt.0.)then
          integer_real_exponent%deriv(1:overloadLength) = integer_real_exponent%value*log(real(i1))*v2%deriv(1:overloadLength)
        end if
      end function integer_real_exponent

!======================== Real_Real_Exponent ================================
!
! Exponent real**real  phi = x1**x2
! Note that if x1 is negative, this is not a defined operation
! This is derived using natural logarithms
! phi = x1**x2
! log(phi) = x2*log(x1)
! 1./phi*d(phi)/d(chi) = x2*(1/x1*d(x1)/d(chi)) + log(x1)*d(x2)/d(chi)
! d(phi)/d(chi) = phi*(x2/x1*d(x1)/d(chi)) + log(x1)*d(x2)/d(chi))
!
!=======================================================================
      elemental function real_real_exponent(v1,v2)
        type (dtype) :: real_real_exponent
        type (dType), intent(in) :: v1,v2
        real_real_exponent%value = v1%value**v2%value
!       real_real_exponent%deriv(1:overloadLength) = v2%value*v1%value**(v2%value-1.)*v1%deriv(1:overloadLength)
        if(v1%value.gt.0.)then
          real_real_exponent%deriv(1:overloadLength) = &
          real_real_exponent%value*(v2%value/v1%value*v1%deriv(1:overloadLength)&
           + log(v1%value)*v2%deriv(1:overloadLength))
        end if
      end function real_real_exponent

!======================== Real_Const_Exponent ==========================
!
! Exponent real**real_constant  v1**r2
! d(v1**r2) = r2*x1**(r2-1)*d(v1)
!
!=======================================================================
      elemental function real_const_exponent(v1,r2)
        type (dtype) :: real_const_exponent
        type (dType), intent(in) :: v1
        real(dp), intent(in) :: r2
        real_const_exponent%value = v1%value**r2
        real_const_exponent%deriv(1:overloadLength) = r2*v1%value**(r2-1.)*v1%deriv(1:overloadLength)
      end function real_const_exponent

!======================== Const_Real_Exponent ================================
!
! Exponent real**real  phi = r1**v2
! Note that if r1 is negative or zero, this is not a defined operation
! This is derived using natural logarithms
! phi = r1**v2
! log(phi) = x2*log(r1)
! 1./phi*d(phi)/d(chi) = log(r1)*d(v2)/d(chi)
! d(phi)/d(chi) = phi*log(r1)*d(v2)/d(chi)
!
!=======================================================================
      elemental function const_real_exponent(r1,v2)
        type (dtype) :: const_real_exponent
        type (dType), intent(in) :: v2
        real(dp), intent(in) :: r1
        const_real_exponent%value = r1**v2%value
        if(r1.gt.0.)then
          const_real_exponent%deriv(1:overloadLength) = const_real_exponent%value*log(r1)*v2%deriv(1:overloadLength)
        end if
      end function const_real_exponent
!
! Comparisons
!
! .LT.
!======================== Real_Real_LT =================================
!
! Compare v1 < v2
!
!=======================================================================
      elemental function real_real_lt(v1,v2)
        logical :: real_real_lt
        type (dType), intent(in) :: v1,v2
        real_real_lt = .false.
        if(v1%value.lt.v2%value)real_real_lt = .true.
      end function real_real_lt

!======================== Real_Integer_LT ==============================
!
! Compare v1 < i2
!
!=======================================================================
      elemental function real_integer_lt(v1,i2)
        logical :: real_integer_lt
        type (dType), intent(in) :: v1
        integer, intent(in) :: i2
        real_integer_lt = .false.
        if(v1%value.lt.real(i2))real_integer_lt = .true.
      end function real_integer_lt

!======================== Integer_Real_LT ==============================
!
! Compare i1 < v2
!
!=======================================================================
      elemental function integer_real_lt(i1,v2)
        logical :: integer_real_lt
        integer, intent(in) :: i1
        type (dType), intent(in) :: v2
        integer_real_lt = .false.
        if(real(i1).lt.v2%value)integer_real_lt = .true.
      end function integer_real_lt

!======================== Real_Const_LT =================================
!
! Compare v1 < r2
!
!=======================================================================
      elemental function real_const_lt(v1,r2)
        logical :: real_const_lt
        type (dType), intent(in) :: v1
        real(dp), intent(in) :: r2
        real_const_lt = .false.
        if(v1%value.lt.r2)real_const_lt = .true.
      end function real_const_lt

!======================== Const_Real_LT ================================
!
! Compare r1 < v2
!
!=======================================================================
      elemental function const_real_lt(r1,v2)
        logical :: const_real_lt
        real(dp), intent(in) :: r1
        type (dType), intent(in) :: v2
        const_real_lt = .false.
        if(r1.lt.v2%value)const_real_lt = .true.
      end function const_real_lt

! .LE.
!======================== Real_Real_LE =================================
!
! Compare v1 <= v2
!
!=======================================================================
      elemental function real_real_le(v1,v2)
        logical :: real_real_le
        type (dType), intent(in) :: v1,v2
        real_real_le = .false.
        if(v1%value.le.v2%value)real_real_le = .true.
      end function real_real_le

!======================== Real_Integer_LE ==============================
!
! Compare v1 <= i2
!
!=======================================================================
      elemental function real_integer_le(v1,i2)
        logical :: real_integer_le
        type (dType), intent(in) :: v1
        integer, intent(in) :: i2
        real_integer_le = .false.
        if(v1%value.le.real(i2))real_integer_le = .true.
      end function real_integer_le

!======================== Integer_Real_LE ==============================
!
! Compare i1 <= v2
!
!=======================================================================
      elemental function integer_real_le(i1,v2)
        logical :: integer_real_le
        integer, intent(in) :: i1
        type (dType), intent(in) :: v2
        integer_real_le = .false.
        if(real(i1).le.v2%value)integer_real_le = .true.
      end function integer_real_le

!======================== Real_Const_LE =================================
!
! Compare v1 <= r2
!
!=======================================================================
      elemental function real_const_le(v1,r2)
        logical :: real_const_le
        type (dType), intent(in) :: v1
        real(dp), intent(in) :: r2
        real_const_le = .false.
        if(v1%value.le.r2)real_const_le = .true.
      end function real_const_le

!======================== Const_Real_LE ================================
!
! Compare r1 <= v2
!
!=======================================================================
      elemental function const_real_le(r1,v2)
        logical :: const_real_le
        real(dp), intent(in) :: r1
        type (dType), intent(in) :: v2
        const_real_le = .false.
        if(r1.le.v2%value)const_real_le = .true.
      end function const_real_le
!
! .GT.
!======================== Real_Real_GT =================================
!
! Compare v1 > v2
!
!=======================================================================
      elemental function real_real_gt(v1,v2)
        logical :: real_real_gt
        type (dType), intent(in) :: v1,v2
        real_real_gt = .false.
        if(v1%value.gt.v2%value)real_real_gt = .true.
      end function real_real_gt

!======================== Real_Integer_GT ==============================
!
! Compare v1 > i2
!
!=======================================================================
      elemental function real_integer_gt(v1,i2)
        logical :: real_integer_gt
        type (dType), intent(in) :: v1
        integer, intent(in) :: i2
        real_integer_gt = .false.
        if(v1%value.gt.real(i2))real_integer_gt = .true.
      end function real_integer_gt

!======================== Integer_Real_GT ==============================
!
! Compare i1 > v2
!
!=======================================================================
      elemental function integer_real_gt(i1,v2)
        logical :: integer_real_gt
        integer, intent(in) :: i1
        type (dType), intent(in) :: v2
        integer_real_gt = .false.
        if(real(i1).gt.v2%value)integer_real_gt = .true.
      end function integer_real_gt

!======================== Real_Const_GT =================================
!
! Compare v1 > r2
!
!=======================================================================
      elemental function real_const_gt(v1,r2)
        logical :: real_const_gt
        type (dType), intent(in) :: v1
        real(dp), intent(in) :: r2
        real_const_gt = .false.
        if(v1%value.gt.r2)real_const_gt = .true.
      end function real_const_gt

!======================== Const_Real_GT ================================
!
! Compare r1 > v2
!
!=======================================================================
      elemental function const_real_gt(r1,v2)
        logical :: const_real_gt
        real(dp), intent(in) :: r1
        type (dType), intent(in) :: v2
        const_real_gt = .false.
        if(r1.gt.v2%value)const_real_gt = .true.
      end function const_real_gt
!
! .GE.
!======================== Real_Real_GE =================================
!
! Compare v1 >= v2
!
!=======================================================================
      elemental function real_real_ge(v1,v2)
        logical :: real_real_ge
        type (dType), intent(in) :: v1,v2
        real_real_ge = .false.
        if(v1%value.ge.v2%value)real_real_ge = .true.
      end function real_real_ge

!======================== Real_Integer_GE ==============================
!
! Compare v1 >= i2
!
!=======================================================================
      elemental function real_integer_ge(v1,i2)
        logical :: real_integer_ge
        type (dType), intent(in) :: v1
        integer, intent(in) :: i2
        real_integer_ge = .false.
        if(v1%value.ge.real(i2))real_integer_ge = .true.
      end function real_integer_ge

!======================== Integer_Real_GE ==============================
!
! Compare i1 >= v2
!
!=======================================================================
      elemental function integer_real_ge(i1,v2)
        logical :: integer_real_ge
        integer, intent(in) :: i1
        type (dType), intent(in) :: v2
        integer_real_ge = .false.
        if(real(i1).ge.v2%value)integer_real_ge = .true.
      end function integer_real_ge

!======================== Real_Const_GE =================================
!
! Compare v1 >= r2
!
!=======================================================================
      elemental function real_const_ge(v1,r2)
        logical :: real_const_ge
        type (dType), intent(in) :: v1
        real(dp), intent(in) :: r2
        real_const_ge = .false.
        if(v1%value.ge.r2)real_const_ge = .true.
      end function real_const_ge

!======================== Const_Real_GE ================================
!
! Compare r1 >= v2
!
!=======================================================================
      elemental function const_real_ge(r1,v2)
        logical :: const_real_ge
        real(dp), intent(in) :: r1
        type (dType), intent(in) :: v2
        const_real_ge = .false.
        if(r1.ge.v2%value)const_real_ge = .true.
      end function const_real_ge

! .EQ.
!======================== Real_Real_EQ =================================
!
! Compare v1 == v2
!
!=======================================================================
      elemental function real_real_eq(v1,v2)
        logical :: real_real_eq
        type (dType), intent(in) :: v1,v2
        real_real_eq = .false.
        if(v1%value.eq.v2%value)real_real_eq = .true.
      end function real_real_eq

!======================== Real_Integer_EQ ==============================
!
! Compare v1 == i2
!
!=======================================================================
      elemental function real_integer_eq(v1,i2)
        logical :: real_integer_eq
        type (dType), intent(in) :: v1
        integer, intent(in) :: i2
        real_integer_eq = .false.
        if(v1%value.eq.real(i2))real_integer_eq = .true.
      end function real_integer_eq

!======================== Integer_Real_EQ ==============================
!
! Compare i1 == v2
!
!=======================================================================
      elemental function integer_real_eq(i1,v2)
        logical :: integer_real_eq
        integer, intent(in) :: i1
        type (dType), intent(in) :: v2
        integer_real_eq = .false.
        if(real(i1).eq.v2%value)integer_real_eq = .true.
      end function integer_real_eq

!======================== Real_Const_EQ =================================
!
! Compare v1 == r2
!
!=======================================================================
      elemental function real_const_eq(v1,r2)
        logical :: real_const_eq
        type (dType), intent(in) :: v1
        real(dp), intent(in) :: r2
        real_const_eq = .false.
        if(v1%value.eq.r2)real_const_eq = .true.
      end function real_const_eq

!======================== Const_Real_EQ ================================
!
! Compare r1 == v2
!
!=======================================================================
      elemental function const_real_eq(r1,v2)
        logical :: const_real_eq
        real(dp), intent(in) :: r1
        type (dType), intent(in) :: v2
        const_real_eq = .false.
        if(r1.eq.v2%value)const_real_eq = .true.
      end function const_real_eq
!
! .NE.
!======================== Real_Real_NE =================================
!
! Compare v1 /= v2
!
!=======================================================================
      elemental function real_real_ne(v1,v2)
        logical :: real_real_ne
        type (dType), intent(in) :: v1,v2
        real_real_ne = .false.
        if(v1%value.ne.v2%value)real_real_ne = .true.
      end function real_real_ne

!======================== Real_Integer_NE ==============================
!
! Compare v1 /= i2
!
!=======================================================================
      elemental function real_integer_ne(v1,i2)
        logical :: real_integer_ne
        type (dType), intent(in) :: v1
        integer, intent(in) :: i2
        real_integer_ne = .false.
        if(v1%value.ne.real(i2))real_integer_ne = .true.
      end function real_integer_ne

!======================== Integer_Real_NE ==============================
!
! Compare i1 /= v2
!
!=======================================================================
      elemental function integer_real_ne(i1,v2)
        logical :: integer_real_ne
        integer, intent(in) :: i1
        type (dType), intent(in) :: v2
        integer_real_ne = .false.
        if(real(i1).ne.v2%value)integer_real_ne = .true.
      end function integer_real_ne

!======================== Real_Const_NE =================================
!
! Compare v1 /= r2
!
!=======================================================================
      elemental function real_const_ne(v1,r2)
        logical :: real_const_ne
        type (dType), intent(in) :: v1
        real(dp), intent(in) :: r2
        real_const_ne = .false.
        if(v1%value.ne.r2)real_const_ne = .true.
      end function real_const_ne

!======================== Const_Real_NE ================================
!
! Compare r1 /= v2
!
!=======================================================================
      elemental function const_real_ne(r1,v2)
        logical :: const_real_ne
        real(dp), intent(in) :: r1
        type (dType), intent(in) :: v2
        const_real_ne = .false.
        if(r1.ne.v2%value)const_real_ne = .true.
      end function const_real_ne
!
! SQRT(real)
!
!============================ Real_SQRT ================================
!
! sqrt(real)
!
!=======================================================================
      elemental function real_sqrt(v1)
        type (dtype) :: real_sqrt
        type (dType), intent(in) :: v1
        real_sqrt%value = sqrt(v1%value)
        if(v1%value.eq.0)then
          real_sqrt%deriv(1:overloadLength) = 0.
        else
          real_sqrt%deriv(1:overloadLength) = 0.5*v1%deriv(1:overloadLength)/sqrt(v1%value)
        end if
      end function real_sqrt
!
! ABS(real)
!
!============================ Real_ABS =================================
!
! abs(real)
!
!=======================================================================
      elemental function real_abs(v1)
        type (dtype) :: real_abs
        type (dType), intent(in) :: v1
        real_abs%value = abs(v1%value)
        if(v1%value.ge.0)then
          real_abs%value = v1%value
          real_abs%deriv(1:overloadLength) = v1%deriv(1:overloadLength)
        else
          real_abs%value = -v1%value
          real_abs%deriv(1:overloadLength) = -v1%deriv(1:overloadLength)
        end if
      end function real_abs
!
! EXP(real)
!
!============================ Real_EXP =================================
!
! exp(real)
!
!=======================================================================
      elemental function real_exp(v1)
        type (dtype) :: real_exp
        type (dType), intent(in) :: v1
        real_exp%value = exp(v1%value)
        real_exp%deriv(1:overloadLength) = real_exp%value*v1%deriv(1:overloadLength)
      end function real_exp
!
! LOG(real)
!
!============================ Real_LOG =================================
!
! log(real)
!
!=======================================================================
      elemental function real_log(v1)
        type (dtype) :: real_log
        type (dType), intent(in) :: v1
        real_log%value = log(v1%value)
        real_log%deriv(1:overloadLength) = v1%deriv(1:overloadLength)/v1%value
      end function real_log
!
! Minimum
!
!============================ Real_Real_Min =============================
!
! min(real,real)
!
!=======================================================================
     elemental function real_real_min(v1,v2)
       type (dtype) :: real_real_min
       type (dType), intent(in) :: v1,v2
       if(v1%value > v2%value) then
         real_real_min%value = v2%value
         real_real_min%deriv(1:overloadLength) = v2%deriv(1:overloadLength)
       else
         real_real_min%value = v1%value
         real_real_min%deriv(1:overloadLength) = v1%deriv(1:overloadLength)
       endif
     end function real_real_min

!============================ Real_Const_Min ===========================
!
! min(real,const)
!
!=======================================================================
     elemental function real_const_min(v1,r2)
       type (dtype) :: real_const_min
       type (dType), intent(in) :: v1
       real(dp), intent(in) :: r2
       if(v1%value > r2) then
         real_const_min%value = r2
         real_const_min%deriv(1:overloadLength) = 0.
       else
         real_const_min%value = v1%value
         real_const_min%deriv(1:overloadLength) = v1%deriv(1:overloadLength)
       endif
     end function real_const_min

!========================== Const_Real_Min =============================
!
! min(const,real)
!
!=======================================================================
     elemental function const_real_min(r1,v2)
       type (dtype) :: const_real_min
       type (dType), intent(in) :: v2
       real(dp), intent(in) :: r1
       if(r1 > v2%value) then
         const_real_min%value = v2%value
         const_real_min%deriv(1:overloadLength) = v2%deriv(1:overloadLength)
       else
         const_real_min%value = r1
         const_real_min%deriv(1:overloadLength) = 0.
       endif
     end function const_real_min
!
! Maximum
!
!============================ Real_Real_Max =============================
!
! max(real,real)
!
!=======================================================================
     elemental function real_real_max(v1,v2)
       type (dtype) :: real_real_max
       type (dType), intent(in) :: v1,v2
       if(v1%value < v2%value) then
         real_real_max%value = v2%value
         real_real_max%deriv(1:overloadLength) = v2%deriv(1:overloadLength)
       else
         real_real_max%value = v1%value
         real_real_max%deriv(1:overloadLength) = v1%deriv(1:overloadLength)
       endif
     end function real_real_max

!============================ Real_Const_Max ===========================
!
! max(real,const)
!
!=======================================================================
     elemental function real_const_max(v1,r2)
       type (dtype) :: real_const_max
       type (dType), intent(in) :: v1
       real(dp), intent(in) :: r2
       if(v1%value < r2) then
         real_const_max%value = r2
         real_const_max%deriv(1:overloadLength) = 0.
       else
         real_const_max%value = v1%value
         real_const_max%deriv(1:overloadLength) = v1%deriv(1:overloadLength)
       endif
     end function real_const_max

!========================== Const_Real_Min =============================
!
! max(const,real)
!
!=======================================================================
     elemental function const_real_max(r1,v2)
       type (dtype) :: const_real_max
       type (dType), intent(in) :: v2
       real(dp), intent(in) :: r1
       if(r1 < v2%value) then
         const_real_max%value = v2%value
         const_real_max%deriv(1:overloadLength) = v2%deriv(1:overloadLength)
       else
         const_real_max%value = r1
         const_real_max%deriv(1:overloadLength) = 0.
       endif
     end function const_real_max
!
! SIN(real)
!
!============================ Real_Sin =================================
!
! sin(real)
!
!=======================================================================
      elemental function real_sin(v1)
        type (dtype) :: real_sin
        type (dType), intent(in) :: v1
        real_sin%value = sin(v1%value)
        real_sin%deriv(1:overloadLength) = cos(v1%value)*v1%deriv(1:overloadLength)
      end function real_sin
!
! COS(real)
!
!============================ Real_Cos =================================
!
! cos(real)
!
!=======================================================================
      elemental function real_cos(v1)
        type (dtype) :: real_cos
        type (dType), intent(in) :: v1
        real_cos%value = cos(v1%value)
        real_cos%deriv(1:overloadLength) = -sin(v1%value)*v1%deriv(1:overloadLength)
      end function real_cos

end module Common_Overload_mod

