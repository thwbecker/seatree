#
# run the inversion
#
PARA_FILE = ../global_para.dat 

include ../makefile.include

all: tools invert 


tools: bin/invray bin/mod2rms bin/dumbcurv bin/dumbder

invert: lcurve.dat

clean:
	rm maps/* hitcount.dat hitcount.xyz lcurve.txt  fort.1 *~ invray bin/*


#
# binaries
#
bin/invray: invray.f
	$(F77) $(FFLAGS) $(F_EXT_SOURCE_FLAG) invray.f -o bin/invray

bin/mod2rms: mod2rms.f
	$(F77) $(FFLAGS) $(F_EXT_SOURCE_FLAG) mod2rms.f -o bin/mod2rms

bin/dumbcurv: curves/dumbcurv.f
	$(F77) $(FFLAGS) $(F_EXT_SOURCE_FLAG) curves/dumbcurv.f -o bin/dumbcurv

bin/dumbder: curves/dumbder.f
	$(F77) $(FFLAGS) $(F_EXT_SOURCE_FLAG) curves/dumbder.f -o bin/dumbder
#
#
#
lcurve.dat:	$(PARA_FILE) invray.f 
	./drive_inversion $(PARA_FILE)
