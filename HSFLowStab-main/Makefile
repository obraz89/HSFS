#
# HSFlowStab makefile

# Project : HSFlowStab-main

#
# Author: Andrey V. Novikov
# Modified : A.Obraz
#

OUTLIB := HSFlowStab-main

#-------------------------

topDIR = $(realpath ..)
plugDIR = $(topDIR)/components
libDIR = $(topDIR)/lib

include $(topDIR)/build_linux/Make-conf

INTDIR := $(topDIR)/$(OUTDIR)/obj/$(OUTLIB)
OUTDIR := $(topDIR)/$(OUTDIR)

CFLAGS += -I. -I$(topDIR)/shared -I$(plugDIR)/SmallMat -I$(plugDIR)/PhysCommon -I$(plugDIR)/Profile -I$(CGNS_DIR)/include_hdf5
CFLAGS += $(MPI_INC_DIR)
CFLAGS += $(WXFLAGS) -vec_report3
#MY_LIBS := shared.so SmallMat.so PhysCommon.so Profile.so
#MY_LIBS := $(addprefix $(OUTDIR)/, $(MY_LIBS))
LFLAGS += -L$(OUTDIR) -lhsfs_shared -lhsfs_smallmat -lhsfs_physcommon -lhsfs_profile $(WXLIBS)
LFLAGS += $(MPI_LFLAGS) $(HDF5_LFLAGS)

SRCDIRS := . 
OBJS := $(wildcard $(addsuffix /*.cpp, $(SRCDIRS) ) )
OBJS := $(notdir $(OBJS))
OBJS := $(addprefix  $(INTDIR)/, $(OBJS:.cpp=.o))

VPATH := $(SRCDIRS)
#-------------------------

$(OUTDIR)/$(OUTLIB).exe : $(OBJS) $(OUTDIR)
	$(CC) -o $@ $(OBJS) $(LFLAGS)

$(INTDIR)/%.o : %.cpp $(INTDIR)
	$(CC) -o $@ -fPIC $(CFLAGS) -c $<


# Create output and intermed dirs
$(OUTDIR):
	mkdir -p $(OUTDIR)
$(INTDIR):
	mkdir -p $(INTDIR)

