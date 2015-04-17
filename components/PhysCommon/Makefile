#
# HSFlowStab makefile

# Project : PhysCommon

#
# Author: Andrey V. Novikov
# Modified : A.Obraz
#

OUTLIB := libhsfs_physcommon

#-------------------------

topDIR = $(realpath ../..)
plugDIR = $(topDIR)/components

include $(topDIR)/build_linux/Make-conf

INTDIR := $(topDIR)/$(OUTDIR)/obj/$(OUTLIB)
OUTDIR := $(topDIR)/$(OUTDIR)

CFLAGS += -I. -DPHYSCOMMON_EXPORT -I$(topDIR)/shared -I$(plugDIR)/SmallMat $(WXFLAGS) -vec_report3
#MYLIBS := shared.so SmallMat.so
#MYLIBS := $(addprefix $(OUTDIR)/, $(MYLIBS))
LFLAGS += -L$(OUTDIR) -Wl,-no-undefined -lhsfs_shared -lhsfs_smallmat $(WXLIBS)

SRCDIRS := . 
OBJS := $(wildcard $(addsuffix /*.cpp, $(SRCDIRS) ) )
OBJS := $(notdir $(OBJS))
OBJS := $(addprefix  $(INTDIR)/, $(OBJS:.cpp=.o))

VPATH := $(SRCDIRS)
#-------------------------

$(OUTDIR)/$(OUTLIB).so : $(OBJS) $(OUTDIR)
	$(CC) -o $@ -shared $(OBJS) $(LFLAGS)

$(INTDIR)/%.o : %.cpp $(INTDIR)
	$(CC) -o $@ -fPIC $(CFLAGS) -c $<


# Create output and intermed dirs
$(OUTDIR):
	mkdir -p $(OUTDIR)
$(INTDIR):
	mkdir -p $(INTDIR)

