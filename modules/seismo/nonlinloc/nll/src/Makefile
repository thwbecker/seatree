# Makefile for NonLinLoc software package
#
# Invocation:
#     Solaris: make all
#     Linux:   make -R all


#BINDIR=${MYBIN}
BINDIR=../bin/

INCLUDE_DIR=./
CCFLAGS = $(CFLAGS)
#CCFLAGS=-O3 -Wall
#CCFLAGS=-O3 -pg -Wall
#CCFLAGS=-g  -Wall

# Linux
TIME3D_CCFLAGS=

# Custom builds
CUSTOM=
#CUSTOM=ETH

GRID_LIB_OBJS=GridLib.o util.o geo.o octtree.o nrutil.o nrmatrix.o ran1.o map_project.o

all : ${BINDIR}/NLLoc ${BINDIR}/Vel2Grid ${BINDIR}/Grid2Time ${BINDIR}/Grid2GMT ${BINDIR}/LocSum \
	${BINDIR}/Time2EQ ${BINDIR}/PhsAssoc ${BINDIR}/hypoe2hyp ${BINDIR}/fpfit2hyp ${BINDIR}/oct2grid \
	${BINDIR}/read_angles

distrib : NLLoc Vel2Grid Grid2Time Grid2GMT LocSum Time2EQ PhsAssoc hypoe2hyp fpfit2hyp oct2grid


# --------------------------------------------------------------------------
# NLLoc
#
OBJS0=NLLocLib.o ${GRID_LIB_OBJS} calc_crust_corr.o velmod.o GridMemLib.o phaselist.o loclist.o
OPTIONS0=
# Custom builds
ifeq ($(CUSTOM), ETH)
	OBJS1=$(OBJS0) eth_functions.o get_region_name_nr.o magnitude.o new_sedlib.o complex.o
	OPTIONS=$(OPTIONS0) -D CUSTOM_ETH
else
	OBJS1=$(OBJS0)
	OPTIONS=$(OPTIONS0)
endif
PVER=1
#NLLoc : ${BINDIR}/NLLoc
${BINDIR}/NLLoc : NLLoc_main.o NLLoc${PVER}.o ${OBJS1}
	$(CC) NLLoc_main.o NLLoc${PVER}.o ${OBJS1} ${CCFLAGS} -o ${BINDIR}/NLLoc -lm
NLLoc_main.o : NLLoc_main.c NLLocLib.h GridLib.h GridMemLib.h ${INCLUDE_DIR}nrutil.h
	$(CC) -c ${CCFLAGS}  NLLoc_main.c  $(OPTIONS)
NLLoc${PVER}.o : NLLoc${PVER}.c NLLocLib.h GridLib.h GridMemLib.h ${INCLUDE_DIR}nrutil.h
	$(CC) -c ${CCFLAGS}  NLLoc${PVER}.c  $(OPTIONS)
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# NLDiffLoc
#
OBJS100=NLLocLib.o ${GRID_LIB_OBJS} calc_crust_corr.o velmod.o GridMemLib.o
#NLDiffLoc : ${BINDIR}/NLDiffLoc
${BINDIR}/NLDiffLoc : NLDiffLoc.o ${OBJS100}
	$(CC) NLDiffLoc.o ${OBJS100} ${CCFLAGS} -o ${BINDIR}/NLDiffLoc -lm
NLDiffLoc.o : NLDiffLoc.c NLDiffLoc.h NLLocLib.h GridLib.h GridMemLib.h ${INCLUDE_DIR}nrutil.h
	$(CC) -c ${CCFLAGS}  NLDiffLoc.c  $(OPTIONS)
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# Vel2Grid
#
OBJS2=${GRID_LIB_OBJS} velmod.o
PVER=1
#Vel2Grid : ${BINDIR}/Vel2Grid
${BINDIR}/Vel2Grid : Vel2Grid${PVER}.o ${OBJS2}
	$(CC) Vel2Grid${PVER}.o  ${OBJS2} ${CCFLAGS} -o ${BINDIR}/Vel2Grid -lm
Vel2Grid${PVER}.o : Vel2Grid${PVER}.c GridLib.h
	$(CC) -c Vel2Grid${PVER}.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# Grid2Time
#
OBJS3=${GRID_LIB_OBJS} Time_3d.o
PVER=1
#Grid2Time : ${BINDIR}/Grid2Time
${BINDIR}/Grid2Time : Grid2Time${PVER}.o ${OBJS3}
	$(CC) Grid2Time${PVER}.o ${OBJS3} ${CCFLAGS}  \
		-o ${BINDIR}/Grid2Time -lm
Grid2Time${PVER}.o : Grid2Time${PVER}.c GridLib.h
	$(CC) -c ${CCFLAGS} Grid2Time${PVER}.c
Time_3d.o : Time_3d.c
	$(CC) ${TIME3D_CCFLAGS} -c -DNO_IEEE_PROTOCOL Time_3d.c
#	$(CC) -c -DNO_IEEE_PROTOCOL Time_3d.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# Grid2GMT
#
OBJS4=Grid2GMT.o ${GRID_LIB_OBJS} vector.o GridGraphLib.o
#Grid2GMT : ${BINDIR}/Grid2GMT
${BINDIR}/Grid2GMT : ${OBJS4}
	$(CC) ${OBJS4} ${CCFLAGS} -o ${BINDIR}/Grid2GMT -lm
Grid2GMT.o : Grid2GMT.c GridLib.h GridGraphLib.h
	$(CC) -c Grid2GMT.c
GridGraphLib.o : GridGraphLib.c GridLib.h GridGraphLib.h
	$(CC) -c GridGraphLib.c
#
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# LocSum
#
OBJS5=LocSum.o ${GRID_LIB_OBJS}
#LocSum : ${BINDIR}/LocSum
${BINDIR}/LocSum : ${OBJS5}
	$(CC) ${OBJS5} ${CCFLAGS} -o ${BINDIR}/LocSum -lm
LocSum.o : LocSum.c GridLib.h
	$(CC) -c LocSum.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# Time2EQ
#
OBJS6=${GRID_LIB_OBJS}
PVER=1
#Time2EQ : ${BINDIR}/Time2EQ
${BINDIR}/Time2EQ : Time2EQ${PVER}.o ${OBJS6}
	$(CC) Time2EQ${PVER}.o ${OBJS6} ${CCFLAGS} -o ${BINDIR}/Time2EQ -lm

Time2EQ${PVER}.o : Time2EQ${PVER}.c GridLib.h ${IN!=CLUDE_DIR}ran1.h
	$(CC) -c ${CCFLAGS}  Time2EQ${PVER}.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# hypoe2hyp
#
OBJS7=hypoe2hyp.o ${GRID_LIB_OBJS}
#hypoe2hyp : ${BINDIR}/hypoe2hyp
${BINDIR}/hypoe2hyp : ${OBJS7}
	$(CC) ${OBJS7} ${CCFLAGS} -o ${BINDIR}/hypoe2hyp -lm
hypoe2hyp.o : hypoe2hyp.c GridLib.h
	$(CC) -c hypoe2hyp.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# fpfit2hyp
#
OBJS8=fpfit2hyp.o ${GRID_LIB_OBJS}
#fpfit2hyp : ${BINDIR}/fpfit2hyp
${BINDIR}/fpfit2hyp : ${OBJS8}
	$(CC) ${OBJS8} ${CCFLAGS} -o ${BINDIR}/fpfit2hyp -lm
fpfit2hyp.o : fpfit2hyp.c GridLib.h
	$(CC) -c fpfit2hyp.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# PhsAssoc
#
OBJS9=PhsAssoc.o ${GRID_LIB_OBJS} calc_crust_corr.o
#PhsAssoc : ${BINDIR}/PhsAssoc
${BINDIR}/PhsAssoc : ${OBJS9}
	$(CC) ${OBJS9} ${CCFLAGS} -o ${BINDIR}/PhsAssoc -lm
PhsAssoc.o : PhsAssoc.c GridLib.h
	$(CC) -c PhsAssoc.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# oct2grid
#
OBJS10=oct2grid.o ${GRID_LIB_OBJS}
#oct2grid : ${BINDIR}/oct2grid
${BINDIR}/oct2grid : ${OBJS10}
	$(CC) ${OBJS10} ${CCFLAGS} -o ${BINDIR}/oct2grid -lm
oct2grid.o : oct2grid.c GridLib.h
	$(CC) -c ${CCFLAGS} oct2grid.c
# --------------------------------------------------------------------------


# read_agnles
${BINDIR}/read_angles : read_angles.o
	$(CC) read_angles.o ${CCFLAGS} -o ${BINDIR}/read_angles -lm




# --------------------------------------------------------------------------
# Librarires
#
NLLocLib.o : NLLocLib.c NLLocLib.h GridLib.h
	$(CC) -c ${CCFLAGS}  NLLocLib.c

GridLib.o : GridLib.c GridLib.h  geometry.h nrutil.h util.h geo.h ${INCLUDE_DIR}octtree.h
	$(CC) -c ${CCFLAGS}  GridLib.c

GridMemLib.o : GridMemLib.c GridMemLib.h  GridLib.h
	$(CC) -c ${CCFLAGS}  GridMemLib.c

phaselist.o : phaselist.c phaseloclist.h  GridLib.h
	$(CC) -c ${CCFLAGS}  phaselist.c

loclist.o : loclist.c phaseloclist.h  GridLib.h
	$(CC) -c ${CCFLAGS}  loclist.c

util.o : ${INCLUDE_DIR}util.c ${INCLUDE_DIR}util.h
	$(CC) -c ${CCFLAGS}  ${INCLUDE_DIR}util.c

octtree.o : ${INCLUDE_DIR}octtree.c ${INCLUDE_DIR}octtree.h
	$(CC) -c ${CCFLAGS}  ${INCLUDE_DIR}octtree.c

nrutil.o : ${INCLUDE_DIR}nrutil.c ${INCLUDE_DIR}nrutil.h
	$(CC) -c ${CCFLAGS}  ${INCLUDE_DIR}nrutil.c

nrmatrix.o : ${INCLUDE_DIR}nrmatrix.c
	$(CC) -c ${CCFLAGS}  ${INCLUDE_DIR}nrmatrix.c

ran1.o : ${INCLUDE_DIR}ran1.c ${INCLUDE_DIR}ran1.h
	$(CC) -c ${CCFLAGS}  ${INCLUDE_DIR}ran1.c

vector.o : ${INCLUDE_DIR}vector.c ${INCLUDE_DIR}vector.h
	$(CC) -c ${CCFLAGS}  ${INCLUDE_DIR}vector.c

geo.o : ${INCLUDE_DIR}geo.c ${INCLUDE_DIR}geo.h
	$(CC) -c ${CCFLAGS}  ${INCLUDE_DIR}geo.c

velmod.o : ${INCLUDE_DIR}velmod.c
	$(CC) -c ${CCFLAGS}  ${INCLUDE_DIR}velmod.c

calc_crust_corr.o :   GridLib.h ${INCLUDE_DIR}calc_crust_corr.c ${INCLUDE_DIR}crust_corr_model.h  ${INCLUDE_DIR}crust_type.h  ${INCLUDE_DIR}crust_type_key.h
	$(CC) -c ${CCFLAGS}  ${INCLUDE_DIR}calc_crust_corr.c

GMT_INCLUDE=./
map_project.o : ${GMT_INCLUDE}map_project.c
	$(CC) -c ${CCFLAGS}  ${GMT_INCLUDE}map_project.c
#
# --------------------------------------------------------------------------

clean :
	rm *.o



# --------------------------------------------------------------------------
# Custom Librarires
#

# ETH ----------------------------------------------------------------------
ETH_DIR=custom_eth/
eth_functions.o : $(ETH_DIR)eth_functions.c $(ETH_DIR)eth_functions.h $(ETH_DIR)new_sedlib.h $(ETH_DIR)complex.h
	$(CC) -c ${CCFLAGS}  $(ETH_DIR)eth_functions.c
get_region_name_nr.o : $(ETH_DIR)get_region_name_nr.c $(ETH_DIR)get_region_name_nr.h $(ETH_DIR)new_sedlib.h
	$(CC) -c ${CCFLAGS}  $(ETH_DIR)get_region_name_nr.c
magnitude.o : $(ETH_DIR)magnitude.c $(ETH_DIR)new_sedlib.h $(ETH_DIR)complex.h
	$(CC) -c ${CCFLAGS}  $(ETH_DIR)magnitude.c
new_sedlib.o : $(ETH_DIR)new_sedlib.c $(ETH_DIR)new_sedlib.h
	$(CC) -c ${CCFLAGS}  $(ETH_DIR)new_sedlib.c  -D UNIX
complex.o : $(ETH_DIR)complex.c $(ETH_DIR)complex.h
	$(CC) -c ${CCFLAGS}  $(ETH_DIR)complex.c  -D UNIX
#gselib.o : $(ETH_DIR)gselib.c $(ETH_DIR)gselib.h
#	$(CC) -c ${CCFLAGS}  $(ETH_DIR)gselib.c  -D UNIX

clean_eth :
	rm eth_functions.o get_region_name_nr.o new_sedlib.o complex.o magnitude.o


#
# --------------------------------------------------------------------------

