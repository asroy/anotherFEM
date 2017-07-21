module Constant_mod
      use KindDefinition_mod, only : DP
      implicit none
      save

      integer, parameter :: ntimeStep_LOWEST = -1, ntimeStep_UPPEST = 1
      integer, parameter :: nequation_MAX = 4
      integer, parameter :: nfb_MAX  = 10
      integer, parameter :: nfbnz_MAX  = 4
      integer, parameter :: nnd_MAX = 10
      integer, parameter :: nqd_2d_MAX = 16
      integer, parameter :: nqd_1d_MAX = 5
      integer, parameter :: nbasis_order_MAX = 3

end module Constant_mod
