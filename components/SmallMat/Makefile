#
# HSFlowStab makefile
#
# Author: Andrey V. Novikov
# Modified : A.Obraz
#

OUTLIB := libhsfs_smallmat

#-------------------------

topDIR = $(realpath ../..)

include $(topDIR)/build_linux/Make-conf

INTDIR := $(topDIR)/$(OUTDIR)/obj/$(OUTLIB)
OUTDIR := $(topDIR)/$(OUTDIR)

CFLAGS += -I. -I$(topDIR)/shared $(WXFLAGS) -I$(MKL_INC_DIR) -DSMALLMAT_EXPORT -vec_report3
#MY_LIBS := shared.so
#MY_LIBS := $(addprefix $(OUTDIR)/, $(MY_LIBS))
LFLAGS += -L$(OUTDIR) -Wl,-no-undefined -lhsfs_shared $(WXLIBS) $(MKL_LFLAGS)

SRCDIRS := . 
OBJS := $(wildcard $(addsuffix /*.cpp, $(SRCDIRS) ) )
OBJS := $(notdir $(OBJS))
#OBJS := test.cpp
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

