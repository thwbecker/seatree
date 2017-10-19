#
# makefile for the data generator
#

# global parameter file
PARA_FILE = ../global_para.dat

include ../makefile.include

all: tools data

tools: bin/noisydatamaker bin/shootray bin/shootray_sr bin/make_sr

data:  rays.ind


#
# automatic data generater
rays.ind: $(PARA_FILE) bin/noisydatamaker bin/shootray  bin/shootray_sr bin/make_sr
	./drive_makedata $(PARA_FILE)


clean_all: clean clean_bin

clean:
	rm rays.* receivers.txt sources.txt  *~ *.o rayc.ps  bin/*

OBJS = bin/

#
# binaries
#
bin/noisydatamaker: $(OBJS)/noisydatamaker.o $(OBJS)/ran.o
	$(F77) $(OBJS)/ran.o $(OBJS)/noisydatamaker.o -o bin/noisydatamaker

bin/shootray: $(OBJS)/shootray.o $(OBJS)/ran.o
	$(F77) $(OBJS)/ran.o $(OBJS)/shootray.o -o bin/shootray

bin/make_sr: $(OBJS)/make_sr.o $(OBJS)/ran_c.o
	$(CC) $(OBJS)/ran_c.o $(OBJS)/make_sr.o -o bin/make_sr -lm 

bin/shootray_sr: $(OBJS)/shootray_sr.o $(OBJS)/ran.o
	$(F77) $(OBJS)/ran.o $(OBJS)/shootray_sr.o -o bin/shootray_sr

$(OBJS)/ran.o: ran.f
	$(F77) ran.f -c $(FFLAGS) -o $(OBJS)/ran.o

$(OBJS)/shootray.o: shootray.f
	$(F77) shootray.f -c $(FFLAGS) $(F_EXT_SOURCE_FLAG) -o $(OBJS)/shootray.o

#$(OBJS)/make_sr.o: make_sr.f
#	$(F77) make_sr.f -c $(FFLAGS) $(F_EXT_SOURCE_FLAG) -o $(OBJS)/make_sr.o

$(OBJS)/make_sr.o: make_sr.c
	$(CC) make_sr.c -c $(CFLAGS) -o $(OBJS)/make_sr.o

$(OBJS)/ran_c.o: ran_c.c
	$(CC) ran_c.c -c $(CFLAGS) -o $(OBJS)/ran_c.o

$(OBJS)/shootray_sr.o: shootray_sr.f
	$(F77) shootray_sr.f -c $(FFLAGS) $(F_EXT_SOURCE_FLAG) -o $(OBJS)/shootray_sr.o

$(OBJS)/noisydatamaker.o: noisydatamaker.f
	$(F77)  noisydatamaker.f -c $(FFLAGS) $(F_EXT_SOURCE_FLAG) -o $(OBJS)/noisydatamaker.o
