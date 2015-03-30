#
# HSFlowStab makefile

# Project : MF.CGNS2D

#
# Author: Andrey V. Novikov
# Modified : A.Obraz
#

OUTLIB := MF.CGNS2D

#-------------------------

topDIR = ../..
plugDIR = ../../components
libDIR = $(topDIR)/lib

include $(topDIR)/build_linux/Make-conf

INTDIR := $(topDIR)/$(OUTDIR)/obj/$(OUTLIB)
OUTDIR := $(topDIR)/$(OUTDIR)

CFLAGS += -I. -I$(topDIR)/shared -I$(plugDIR)/SmallMat -I$(plugDIR)/PhysCommon -I$(plugDIR)/Profile -I$(plugDIR)/MF.CGNS-shared
CFLAGS += -I$(libDIR)/cgns
CFLAGS += $(WXFLAGS) -vec_report3
LFLAGS += -L$(OUTDIR) -Wl,-no-undefined -lshared.so -lSmallMat.so -lPhysCommon.so -lProfile.so -lMF.CGNS-shared.a
LFLAGS += -L$(libDIR)/cgns -llibhdf5.a -llibcgns.a
LFLAGS += -L$(libDIR)/mpi -llibmpi.a

SRCDIRS := . 
OBJS := $(wildcard $(addsuffix /*.cpp, $(SRCDIRS) ) )
OBJS := $(notdir $(OBJS))
OBJS := $(addprefix  $(INTDIR)/, $(OBJS:.cpp=.o))

VPATH := $(SRCDIRS)
#-------------------------

$(OUTDIR)/$(OUTLIB).so : $(OBJS) $(OUTDIR)
	$(CC) -o $@ -shared $(LFLAGS) $(OBJS)

$(INTDIR)/%.o : %.cpp $(INTDIR)
	$(CC) -o $@ -fPIC $(CFLAGS) -c $<


# Create output and intermed dirs
$(OUTDIR):
	mkdir -p $(OUTDIR)
$(INTDIR):
	mkdir -p $(INTDIR)

