#
# HSFlowStab makefile
#
# Author: Andrey V. Novikov
# Modified : A.Obraz
#

OUTLIB := SmallMat

#-------------------------

topDIR = ../..

include $(topDIR)/build_linux/Make-conf

INTDIR := $(topDIR)/$(OUTDIR)/obj/$(OUTLIB)
OUTDIR := $(topDIR)/$(OUTDIR)/plugins

CFLAGS += -I$(topDIR)/shared -I. $(WXFLAGS) -vec_report3
LFLAGS += -L(OUTDIR) -Wl,-no-undefined -lshared.so

SRCDIRS := . 
OBJS := $(wildcard $(addsuffix /*.cpp, $(SRCDIRS) ) )
OBJS := $(notdir $(OBJS))
OBJS := $(addprefix  $(INTDIR)/, $(OBJS:.cpp=.o))

VPATH := $(SRCDIRS)
#-------------------------

$(OUTDIR)/$(OUTLIB).so : $(OBJS) $(OUTDIR) $(OUTDIR)/shared.so
	$(CC) -o $@ -shared $(LFLAGS) $(OBJS)

$(INTDIR)/%.o : %.cpp $(INTDIR)
	$(CC) -o $@ -fPIC $(CFLAGS) -c $<


# Create output and intermed dirs
$(OUTDIR):
	mkdir -p $(OUTDIR)
$(INTDIR):
	mkdir -p $(INTDIR)

