#
# HSFlowStab makefile

# Project : PhysCommon

#
# Author: Andrey V. Novikov
# Modified : A.Obraz
#

OUTLIB := PhysCommon

#-------------------------

topDIR = ../..
plugDIR = ../../components

include $(topDIR)/build_linux/Make-conf

INTDIR := $(topDIR)/$(OUTDIR)/obj/$(OUTLIB)
OUTDIR := $(topDIR)/$(OUTDIR)

CFLAGS += -I. -DPHYSCOMMON_EXPORT -I$(topDIR)/shared -I$(plugDIR)/SmallMat $(WXFLAGS) -vec_report3
LFLAGS += -L$(OUTDIR) -Wl,-no-undefined -lshared.so -lSmallMat.so

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

