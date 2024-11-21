F77           = f90

FFLAGS        = +O3 +U77

LINKER	      = f90

LIBFLAGS_A    =  -L/home/local/lib  -lkempo1 -lcgraph -loption +U77

LIBFLAGS_B    = -L/home/local/lib -lkempo1 -lcX11 -loption -lX11

LIBFLAGS_C    = -L/home/staff/omura/lib -lkempo1 -lcgraphp  -loption -ltcs


OBJS	      = bfield.o \
		charge.o \
		chargv.o \
		chkprm.o \
		curntv.o \
		currnt.o \
		ecrrct.o \
		efield.o \
		energy.o \
		renrgy.o \
		fldplt.o \
		inital.o \
		initam.o \
		initan.o \
		input.o \
		kspplt.o \
		main.o \
		phsplt.o \
		pltprm.o \
		positn.o \
		reader.o \
		renorm.o \
		spectr.o \
		vdsplt.o \
		velcty.o \
		rvlcty.o \
		bvlcty.o \
		ringds.o \
		rexch.o \
		vinit.o \
		qsort.o \
		writer.o 

PROGRAM_A     = a

PROGRAM_B     = b

PROGRAM_C     = c

.f.o:
		$(F77) $(FFLAGS) -c $*.f

all:		$(PROGRAM_A)

$(PROGRAM_A):   $(OBJS)
		$(LINKER) $(OBJS) -o $@ $(LIBFLAGS_A)

$(PROGRAM_B):   $(OBJS)
		$(LINKER) $(OBJS) -o $@ $(LIBFLAGS_B)

$(PROGRAM_C):   $(OBJS)
		$(LINKER) $(OBJS) -o $@ $(LIBFLAGS_C)

clean:;		rm -f $(OBJS) $(PROGRAM_A) $(PROGRAM_B) $(PROGRAM_C) core cgraph.dat cgraph.aux
