#
# HSFlowStab makefile

# Project : MF.CGNS-shared

#
# Author: Andrey V. Novikov
# Modified : A.Obraz
#

OUTLIB := MF.CGNS-shared

#-------------------------

topDIR = ../..
libDIR = $(topDIR)/lib
plugDIR = ../../components

include $(topDIR)/build_linux/Make-conf

INTDIR := $(topDIR)/$(OUTDIR)/obj/$(OUTLIB)
OUTDIR := $(topDIR)/$(OUTDIR)

CFLAGS += -I. -I$(topDIR)/shared -I$(plugDIR)/SmallMat -I$(plugDIR)/PhysCommon 
CFLAGS += -I$(libDIR)/cgns/$(CGNS_ARCH)/include
CFLAGS += $(WXFLAGS) -vec_report3
#LFLAGS += -L$(OUTDIR) -Wl,-no-undefined -lshared.so -lSmallMat.so -lPhysCommon.so

SRCDIRS := . 
OBJS := $(wildcard $(addsuffix /*.cpp, $(SRCDIRS) ) )
OBJS := $(notdir $(OBJS))
OBJS := $(addprefix  $(INTDIR)/, $(OBJS:.cpp=.o))

VPATH := $(SRCDIRS)
#-------------------------

$(OUTDIR)/$(OUTLIB).a : $(OBJS) $(OUTDIR)
	ar rcs $@ $(OBJS)

$(INTDIR)/%.o : %.cpp $(INTDIR)
	$(CC) -o $@ -fPIC $(CFLAGS) -c $<

# Create output and intermed dirs
$(OUTDIR):
	mkdir -p $(OUTDIR)
$(INTDIR):
	mkdir -p $(INTDIR)

