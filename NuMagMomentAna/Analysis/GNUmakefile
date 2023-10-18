LIB_TYPE    := shared
LIB         := lib$(PACKAGE)Analysis
LIBCXXFILES := $(wildcard *.cxx)

include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_stan.mk
include SoftRelTools/arch_spec_srproxy.mk
include SoftRelTools/arch_spec_cafanacore.mk

# Don't include the arch_specs, only need the headers, don't want the libs
override CPPFLAGS +=-I$(GENIE_INC) -I$(GENIE_INC)/GENIE -I$(NUGEN_INC)/ -I$(NOVARWGT_INC)

override CPPFLAGS += -fno-var-tracking-assignments

override LIBLIBS += $(LOADLIBES) -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR)
