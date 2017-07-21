module KindDefinition_mod
      implicit none
      private

      public :: DP

! double precision (IEEE 754)
      integer, parameter :: SYSTEM_R8 = selected_real_kind(15, 307)
      integer, parameter :: DP = SYSTEM_R8

end module KindDefinition_mod
