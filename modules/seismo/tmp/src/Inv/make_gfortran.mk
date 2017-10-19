FLAGSSUSE= -ffixed-line-length-none -O3
FC=gfortran
SRCDIR=.

all: make.mk joint_lsqr_vx_isodamp vox2roughness
	rm *.o

#EXECUTABLE FILES
joint_lsqr_vx_isodamp: main.o joint_lsqr_vx_isodamp.o gradamp.o readmatrix.o
	$(FC) $(FLAGSSUSE) -o joint_lsqr_vx_isodamp main.o \
	joint_lsqr_vx_isodamp.o gradamp.o readmatrix.o

vox2roughness: vox2roughness.f
	$(FC) $(FLAGSSUSE) -o vox2roughness vox2roughness.f

#OBJECTS

main.o: main.f90
	$(FC) $(FLAGSSUSE) -c -o main.o main.f90

joint_lsqr_vx_isodamp.o: joint_lsqr_vx_isodamp.f
	$(FC) $(FLAGSSUSE) -c -o joint_lsqr_vx_isodamp.o \
	joint_lsqr_vx_isodamp.f

gradamp.o: gradamp.f
	$(FC) $(FLAGSSUSE) -c -o gradamp.o \
	gradamp.f

readmatrix.o: readmatrix.f
	$(FC) $(FLAGSSUSE) -c -o readmatrix.o \
	readmatrix.f

clean:
	rm *.o
