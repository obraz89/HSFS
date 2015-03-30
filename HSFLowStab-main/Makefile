#
# HSFlowStab makefile

# Project : HSFLowStab-main

#
# Author: Andrey V. Novikov
# Modified : A.Obraz
#

OUTLIB := HSFLowStab-main

#-------------------------

topDIR = ../..
plugDIR = ../../components
libDIR = $(topDIR)/lib

include $(topDIR)/build_linux/Make-conf

INTDIR := $(topDIR)/$(OUTDIR)/obj/$(OUTLIB)
OUTDIR := $(topDIR)/$(OUTDIR)

CFLAGS += -I. -I$(topDIR)/shared -I$(plugDIR)/SmallMat -I$(plugDIR)/PhysCommon -I$(plugDIR)/Profile
CFLAGS += $(WXFLAGS) -vec_report3
CFLAGS += -I$(libDIR)/mpi/include
LFLAGS += -L$(OUTDIR) -Wl,-no-undefined -lshared.so -lSmallMat.so -lPhysCommon.so -lProfile.so
LFLAGS += -L$(libDIR)/mpi/lib -lmpi.a

SRCDIRS := . 
OBJS := $(wildcard $(addsuffix /*.cpp, $(SRCDIRS) ) )
OBJS := $(notdir $(OBJS))
OBJS := $(addprefix  $(INTDIR)/, $(OBJS:.cpp=.o))

VPATH := $(SRCDIRS)
#-------------------------

$(OUTDIR)/$(OUTLIB).exe : $(OBJS) $(OUTDIR)
	$(CC) -o $@ $(LFLAGS) $(OBJS)

$(INTDIR)/%.o : %.cpp $(INTDIR)
	$(CC) -o $@ -fPIC $(CFLAGS) -c $<


# Create output and intermed dirs
$(OUTDIR):
	mkdir -p $(OUTDIR)
$(INTDIR):
	mkdir -p $(INTDIR)

