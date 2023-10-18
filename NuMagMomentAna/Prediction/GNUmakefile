LIB_TYPE    := shared
LIB         := lib$(PACKAGE)Prediction
LIBCXXFILES := $(wildcard *.cxx)

include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_boost.mk
include SoftRelTools/arch_spec_eigen.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_ifdhc.mk
include SoftRelTools/arch_spec_stan.mk
include SoftRelTools/arch_spec_srproxy.mk
include SoftRelTools/arch_spec_osclib.mk
include SoftRelTools/arch_spec_cafanacore.mk

override LIBLIBS += $(LOADLIBES) -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR)