module GlobalVariableForDebug_mod
      use KindDefinition_mod, only : DP
      save

      integer, parameter :: fpd1 = 201, fpd2 = 202, fpd3 = 203, fpd4 = 204

      integer, parameter :: fpIO = 222

      logical :: Do_print = .false.

end module GlobalVariableForDebug_mod
