#IGNORE:
#CMPLR = f77
CMPLR  = gfortran
CPLUS  = g++
FFLAGS = -lstdc++
LFLAGS = const_bse.h zdata.h 
#ls
MODEL  = GENEVA
#K2STR  = yes

EXETAIL = .bonn
ifeq ($(MODEL),GENEVA)
	CFLAGS  = -DGENEVAMODEL
	EXETAIL = .geneva
endif
ifeq ($(MODEL),MARIGO)
	CFLAGS  = -DMARIGOMODEL
	EXETAIL = .marigo
endif
ifeq ($(MODEL),KINUGAWA)
	CFLAGS  = -DKINUGAWACE
	EXETAIL = .bonn
endif
ifeq ($(K2STR),yes)
	CFLAGS += -DK2_YOSHIDA
	K2TAIL  = .k2str
endif
SSEEXE = sse$(EXETAIL)$(K2TAIL)
BSEEXE = bse$(EXETAIL)$(K2TAIL)
POPEXE = popbin2$(EXETAIL)$(K2TAIL)

.f.o:
	$(CMPLR) -c $(FFLAGS) $<

SRCE1 = \
sse.f deltat.f evolv1.f hrdiag.f kick.f mlwind.f mrenv.f \
ran3.f star.f zcnsts.f zfuncs.f
SRCE0 = empfuncs.cpp phasefuncs.cpp

OBJT1 = $(SRCE1:.f=.o)
OBJT0 = $(SRCE0:.cpp=.o)

sse:    $(OBJT1) $(OBJT0) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT1) $(OBJT0) -o sse 
	mv sse $(SSEEXE)

empfuncs.o: empfuncs.cpp
	$(CPLUS) -c empfuncs.cpp $(CFLAGS)

phasefuncs.o: phasefuncs.cpp
	$(CPLUS) -c phasefuncs.cpp

SRCE2 = \
bse.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f merger.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f
 
OBJT2 = $(SRCE2:.f=.o)

bse:    $(OBJT2) $(OBJT0) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT2) $(OBJT0) -o bse 
	mv bse $(BSEEXE)

SRCE3 = \
popbin.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f merger.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f
 
OBJT3 = $(SRCE3:.f=.o)

popbin: $(OBJT3) $(OBJT0) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT3) $(OBJT0) -o popbin

SRCE4 = \
popbin2.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f merger.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f
 
OBJT4 = $(SRCE4:.f=.o)

popbin2: $(OBJT4) $(OBJT0) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT4) $(OBJT0) -o popbin2
	mv popbin2 $(POPEXE)

work:
	mkdir temp
	ln -s ../ffbonn temp/
	ln -s ../popbin2 temp/
	ln -s ../popbin2.sh temp/

clean:
	rm -f *~ *.o fort.99 errorwrite.txt binary.dat search.out binaries.out

distclean:
	make clean
	rm -f popbin
	rm -f sse.bonn* bse.bonn* popbin2.bonn*
	rm -f sse.geneva* bse.geneva* popbin2.geneva*

