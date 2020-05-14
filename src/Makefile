CC ?= gcc
CCPP ?= g++

LNAME = satskymap
L = lib$(LNAME).a
LIBSO = lib$(LNAME).so

# HEALPix C lib
HPXDIR = ./HealpixC
HPXLIB = -L$(HPXDIR) -lchealpix
HPXINC = -I$(HPXDIR)

CFLAGS = -Wall -O3 -fPIC $(HPXINC) # -DDEBUG

# Using shared lib
# LDLIBS = -L. -l$(LNAME) -lm

# Using static lib
LDLIBS = $(HPXLIB) $(L) -lm

DESTLIB = /usr/local/lib
BINDIR = ${HOME}/bin

VER = v02b_`/bin/date +"%Y%m%d"`
dist = $(LNAME)_$(VER)


LIBOBJS = sgp.o sgp4.o sgp8.o sdp4.o sdp8.o deep.o basics.o get_el.o common.o observe.o date2mjd.o mjd2date.o skysep_h.o lmst_hr.o dechalat2alt.o lpsun_radec.o earthtilt.o zenith_site_ra.o
LIBFILES = $(LIBOBJS:.o=.c)

EXE = sat_skymap
MAIN = $(EXE:=.c)

EXTRAFILES = Makefile norad.h norad_in.h observe.h const_def.h sat_skymap_def.h stations.txt default.tle tle_retrive.sh $(HPXDIR)

all: lib libso libhpx $(EXE)

exe: $(EXE)
lib: $(L)
libso: $(LIBSO)

$(L): $(LIBOBJS)
	ar rv $@ $?
	ranlib $@

$(LIBSO): $(LIBOBJS)
	$(CC) $(CFLAGS) -shared -o $@ $(LIBOBJS) $(LIBS)

libhpx:
	@ cd $(HPXDIR) ;  $(MAKE) WITHOUT_CFITSIO=1

dist:
	tar zcvf $(dist).tar.gz $(LIBFILES) $(MAIN) $(EXTRAFILES)

tar:
	tar cvf $(dist).tar $(LIBFILES) $(MAIN) $(EXTRAFILES)

install: $(L) $(LIBSO)
	cp -a $(L) $(LIBSO) $(DESTLIB)
	mv $(EXE) $(BINDIR)

clean:
	rm -f *.o $(L) $(LIBSO) $(EXE)
	@ cd $(HPXDIR) ;  $(MAKE) clean
