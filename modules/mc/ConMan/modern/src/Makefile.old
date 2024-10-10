#
# makefile for the program conman on 
# 

COMMON=ConMan.o adjust_velocity_rhs_for_bc.o elminp.o error.o \
	form_temp_matrix.o form_velocity_rhs.o \
	form_temp_rhs_explicit.o form_velocity_stiffness_matrix.o \
	form_temp_rhs_impicit.o velo_stiff_matrix_assemble.o \
	form_temp_mass_matrix.o \
        time_driver.o dyna_P.o genien.o genienIP.o \
        genmat.o genshg.o genshp.o \
        input.o pickdt.o printd.o \
        rdflux.o ddiff.o vadd.o batchelor.o deos.o \
        intemp.o stress.o fluxke.o geoid.o \
        output.o masflx.o \
	projct2.o prj4q.o grdpar.o ppmprt.o \
	print_reg.o prtstr_reg.o print_compbm_data.o print.o \
        rheol.o output_rheol.o

LIBS= ../libs/libsolver.a ../libs/libgrid.a ../libs/libutils.a 

PROGRAM=../conman
#

FC = gfortran
FFLAGS = -O3 -I. -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 
LDFLAGS = -I. -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8

conman: $(COMMON)  
	(cd solver.src ; make)
	(cd grid.src ; make)
	(cd utils.src ; make)
	$(FC) $(LDFLAGS) -o $(PROGRAM) $(COMMON) $(LIBS)

clean:
	rm -rf core $(COMMON) $(PROGRAM)

veryclean:
	rm -rf core $(COMMON) ../libs/lib* *.mod $(PROGRAM)
	(cd solver.src ; make clean)
	(cd grid.src ; make clean)
	(cd utils.src ; make clean)
