#
# HSFlow makefile
#
# Author: Andrey V. Novikov
#

OUTLIB := libhsfs_shared

#-------------------------

topDIR = $(realpath ../)

include $(topDIR)/build_linux/Make-conf

INTDIR := $(topDIR)/$(OUTDIR)/obj/$(OUTLIB)
OUTDIR := $(topDIR)/$(OUTDIR)

CFLAGS += -I. $(WXFLAGS) -vec_report3 -DSHARED_EXPORT
#LFLAGS += -Wl,-no-undefined

SRCDIRS := . 
#OBJS := $(wildcard $(addsuffix /*.cpp, $(SRCDIRS) ) )
#OBJS := $(notdir $(OBJS))
OBJS := common_data.cpp gen_exception.cpp io_helpers.cpp log.cpp
OBJS := $(addprefix  $(INTDIR)/, $(OBJS:.cpp=.o))

VPATH := $(SRCDIRS)
#-------------------------

$(OUTDIR)/$(OUTLIB).so : $(OBJS) $(OUTDIR)
	$(CC) -o $@ -shared $(OBJS) $(WXLIBS)

$(INTDIR)/%.o : %.cpp $(INTDIR)
	$(CC) -o $@ -fPIC $(CFLAGS) -c $<


# Create output and intermed dirs
$(OUTDIR):
	mkdir -p $(OUTDIR)
$(INTDIR):
	mkdir -p $(INTDIR)

