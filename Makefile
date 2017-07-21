include ../make.local

CODE_DIR = ../src

FC = $(F_COMPILER)

include $(CODE_DIR)/make.rules

define build_obj
$(FC) -c $< $(FFLAGS)
endef

define build_exec
$(FC) -o $@ $+ $(LDFLAGS)
endef


clean_obj::
	rm -f *.o *.mod

clean:: clean_obj clean_exec

