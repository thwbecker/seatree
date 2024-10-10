#
# makefile for the program conman on 
# 

COMMON=ConMan.o adjust_velocity_rhs_for_bc.o elminp.o error.o \
	form_temp_matrix.o form_velocity_rhs.o \
	form_temp_rhs_explicit.o form_velocity_stiffness_matrix.o \
	form_temp_rhs_implicit.o velo_stiff_matrix_assemble.o \
	form_temp_mass_matrix.o \
        time_driver.o dyna_P.o genien.o genienIP.o \
        genmat.o genshg.o genshp.o input.o pickdt.o printd.o \
        rdflux.o ddiff.o vadd.o deos.o \
        intemp.o stress.o fluxke.o geoid.o \
        output.o masflx.o projct2.o prj4q.o grdpar.o ppmprt.o \
	print_compbm_data.o print.o rheol.o output_rheol.o 

SOLVER=./solver.src/back.o ./solver.src/coldot.o ./solver.src/factor.o \
	./solver.src/unback.o ./solver.src/unfact.o ./solver.src/zfactor.o

UTILS=./utils.src/fftsubs.o ./utils.src/mydate.o \
	./utils.src/mytime.o ./utils.src/timer.o

GRID=./grid.src/clear.o ./grid.src/genBC.o ./grid.src/genfl.o \
	./grid.src/genfl12.o ./grid.src/geniBC.o ./grid.src/genNBC.o \
	./grid.src/iclear_mov.o ./grid.src/move.o ./grid.src/smove.o

SUBDUCT=./subduct.src/fault.o ./subduct.src/batchelor.o

PROGRAM=../conman
#
#
# GNU Fortran (GCC) 9.3.1 20200408 (Red Hat 9.3.1-2)
# WARNING: -Ofast made large runs core dump
#

FC = gfortran
FFLAGS = -O3 -I. -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 
LDFLAGS = -I. -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8

conman: $(COMMON) $(SOLVER) $(UTILS) $(GRID) $(SUBDUCT)
	$(FC) $(LDFLAGS) -o $(PROGRAM) $(COMMON) $(SOLVER) $(UTILS) $(GRID) $(SUBDUCT)

clean:
	rm -rf core $(COMMON) $(SOLVER) $(UTILS) $(GRID) $(PROGRAM) $(SUBDUCT)

