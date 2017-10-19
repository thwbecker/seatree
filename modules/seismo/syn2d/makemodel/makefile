#
# make a model file
#
PARA_FILE = ../global_para.dat


include ../makefile.include 

all: tools model

tools: bin/chkbd

model: chkbd.xyz


clean_all: clean clean_bin

clean:
	rm chkbd.px chkbd.xyz chkbd.xyz.grd chkbd.xyz.ps *~ bin/*

#
# binary to create checkboard
#
bin/chkbd: chkbd.f
	$(F77) $(FFLAGS) $(F_EXT_SOURCE_FLAG) chkbd.f -o bin/chkbd


#
# generate output model file
#
# xmax and dx are arguments to chkbd script
chkbd.xyz:	chkbd.f $(PARA_FILE) ../pm1.cpt 
	./drive_makemodel $(PARA_FILE) ../pm1.cpt

image.xyz:	chkbd.f $(PARA_FILE) ../pm1.cpt image.pgm
	./drive_makemodel $(PARA_FILE) ../pm1.cpt



