FLAGSSUSE= -ffixed-line-length-none -O3
FC=gfortran
SRCDIR=.

all: make.mk mapview_3d
	rm *.o

#EXECUTABLE FILE
mapview_3d: mapview_3d.f param.o
	$(FC) $(FLAGSSUSE) -o mapview_3d mapview_3d.f param.o

#OBJECT

param.o: param.f
	$(FC) $(FLAGSSUSE) -c -o param.o param.f

clean:
	rm *.o
