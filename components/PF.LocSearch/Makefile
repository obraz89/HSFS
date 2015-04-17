#
# HSFlowStab makefile

# Project : PF.GlobSearch_Lapack

#
# Author: Andrey V. Novikov
# Modified : A.Obraz
#

OUTLIB := PF.LocSearch

#-------------------------

topDIR = $(realpath ../..)
plugDIR = $(topDIR)/components
libDIR = $(topDIR)/lib

include $(topDIR)/build_linux/Make-conf

INTDIR := $(topDIR)/$(OUTDIR)/obj/$(OUTLIB)
OUTDIR := $(topDIR)/$(OUTDIR)

CFLAGS += -I. -I$(topDIR)/shared -I$(plugDIR)/SmallMat -I$(plugDIR)/PhysCommon -I$(plugDIR)/Profile
CFLAGS += $(WXFLAGS) -vec_report3
#MY_LIBS := shared.so SmallMat.so PhysCommon.so Profile.so
#MY_LIBS := $(addprefix $(OUTDIR)/, $(MY_LIBS))
LFLAGS += -L$(OUTDIR) -lhsfs_shared -lhsfs_smallmat -lhsfs_physcommon -lhsfs_profile $(WXLIBS)

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

