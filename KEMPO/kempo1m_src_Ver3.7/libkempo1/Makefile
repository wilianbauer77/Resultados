DESTDIR = /home/wilian/.local/lib
F77           =  f95
FFLAGS        =  +O3
CFLAGS  = 
OBJS    = \
etrans.o	prmplt.o	skfft.o		strnd1.o	xaxis1.o \
four1.o		qlook.o		strnda.o	strndm.o	xaxis2.o \
gnumbr.o	qlook2.o	yaxis1.o 	unrnda.o	rfftcv.o \
ascale.o	lxaxis.o	realft.o	unrndm.o	yaxis2.o \
dplot.o		lyaxis.o	rkfft.o		unrnd1.o	wkfft.o \
enumbr.o	maxmin.o	rndoff.o	wkplot.o	
LDFLAGS = 
PROGRAM = libkempo1.a

.f.o:
		$(F77) $(FFLAGS) -c $*.f

all:$(PROGRAM)

$(PROGRAM): $(OBJS) $(CSRCS)
	ar rv $(PROGRAM) $(OBJS)

install:	$(PROGRAM)
	mkdir -p /home/wilian/.local/lib
	mv $(PROGRAM) $(DESTDIR)
clean:;
	rm -f *.o $(PROGRAM)
###
