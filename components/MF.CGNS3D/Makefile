#
# HSFlowStab makefile

# Project : MF.CGNS3D

#
# Author: Andrey V. Novikov
# Modified : A.Obraz
#

OUTLIB := MF.CGNS3D

#-------------------------

topDIR = $(realpath ../..)
plugDIR = $(topDIR)/components
libDIR = $(topDIR)/lib

include $(topDIR)/build_linux/Make-conf

INTDIR := $(topDIR)/$(OUTDIR)/obj/$(OUTLIB)
OUTDIR := $(topDIR)/$(OUTDIR)

CFLAGS += -I. -I$(topDIR)/shared -I$(plugDIR)/SmallMat -I$(plugDIR)/PhysCommon -I$(plugDIR)/Profile -I$(plugDIR)/MF.CGNS-shared
CFLAGS += -I$(CGNS_INC_DIR)
CFLAGS += $(WXFLAGS) -vec_report3
#MY_LIBS_D := shared.so SmallMat.so PhysCommon.so Profile.so
#MY_LIBS_D := $(addprefix $(OUTDIR)/, $(MY_LIBS_D))
#MY_LIBS_S := MF.CGNS-shared.a
#MY_LIBS_S := $(addprefix $(OUTDIR)/, $(MY_LIBS_S))
#LFLAGS += -Wl,-Bstatic $(libDIR)/cgns/$(CGNS_ARCH)/lib/libhdf5.a $(libDIR)/cgns/$(CGNS_ARCH)/lib/libcgns.a
LFLAGS +=  -L$(OUTDIR) -lhsfs_shared -lhsfs_smallmat -lhsfs_physcommon -lhsfs_profile $(CGNS_LFLAGS)  $(WXLIBS)
SRCDIRS := ../MF.CGNS-shared .
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

