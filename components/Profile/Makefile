#
# HSFlowStab makefile

# Project : Profile

#
# Author: Andrey V. Novikov
# Modified : A.Obraz
#

OUTLIB := libhsfs_profile

#-------------------------

topDIR = $(realpath ../..)
plugDIR = $(topDIR)/components

include $(topDIR)/build_linux/Make-conf

INTDIR := $(topDIR)/$(OUTDIR)/obj/$(OUTLIB)
OUTDIR := $(topDIR)/$(OUTDIR)

CFLAGS += -I$(topDIR)/shared -I$(plugDIR)/SmallMat -I$(plugDIR)/PhysCommon -I. $(WXFLAGS) -vec_report3
CFLAGS += -I$(MPI_INC_DIR)

#MYLIBS := shared.so SmallMat.so PhysCommon.so
#MYLIBS := $(addprefix $(OUTDIR)/, $(MYLIBS))

LFLAGS += -L$(OUTDIR) -Wl,-no-undefined -lhsfs_shared -lhsfs_smallmat -lhsfs_physcommon $(WXLIBS)
LFLAGS += -L$(MPI_LIB_DIR) -lmpi
SRCDIRS := . 
OBJS := $(wildcard $(addsuffix /*.cpp, $(SRCDIRS) ) )
OBJS := $(notdir $(OBJS))
OBJS := $(addprefix  $(INTDIR)/, $(OBJS:.cpp=.o))

VPATH := $(SRCDIRS)
#-------------------------

$(OUTDIR)/$(OUTLIB).so : $(OBJS) $(OUTDIR)
	$(CC) -o $@ -shared $(OBJS) $(LFLAGS)

$(INTDIR)/%.o : %.cpp $(INTDIR)
	$(CC) -o $@ -fPIC $(CFLAGS) -c -H $<


# Create output and intermed dirs
$(OUTDIR):
	mkdir -p $(OUTDIR)
$(INTDIR):
	mkdir -p $(INTDIR)

