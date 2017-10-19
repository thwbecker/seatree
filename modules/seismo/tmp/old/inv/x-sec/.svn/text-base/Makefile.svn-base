include ../../makefile.system

all: ../$(ARCH)/cross_sec

clean:
	rm ../$(ARCH)/cross_sec.o ../$(ARCH)/param.o

../$(ARCH)/cross_sec:  ../$(ARCH)/cross_sec.o ../$(ARCH)/param.o
	$(F77) $(FFLAGS) $(ADDFLAGS) -o ../$(ARCH)/cross_sec ../$(ARCH)/cross_sec.o ../$(ARCH)/param.o

../$(ARCH)/cross_sec.o: cross_sec.f
	$(F77) $(FFLAGS) $(ADDFLAGS) -c -o ../$(ARCH)/cross_sec.o cross_sec.f

#param.o: param.f
#	ifc $(FLAGSSUSE) -c -o param.o param.f

../$(ARCH)/param.o: ../../mat/param.f
	$(F77) $(FFLAGS) $(ADDFLAGS)  -c -o ../$(ARCH)/param.o ../../mat/param.f
