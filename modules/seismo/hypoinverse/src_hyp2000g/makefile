BDIR = ../bin

all: hyp install

#  HYP	    HYPOINVERSE EARTHQUAKE LOCATION PROGRAM.
hyp  :  hyp.o hybda.o hycmd.o hysta.o hydel.o hyate.o hyfmc.o hycal.o \
hyxmc.o hycrh.o hycrt.o hystl.o hyopen.o hyinit.o hyphs.o hysou.o \
hycin.o hycout.o hytrl.o hyloc.o hysol.o hysvd.o hytra.o hytrh.o hytrt.o \
hymag.o hyrep.o hylst.o hysum.o hyinp.o hypro.o medwt.o utmcal.o \
init_event.o hydelt.o hyedit.o hytime.o hybeg.o hypref.o hymagp.o \
hyfile.o allsubs.o
	$(F77) $(FFLAGS) -g hyp.o hybda.o hycmd.o hysta.o hydel.o hyate.o hyfmc.o hypref.o \
	hyxmc.o hycrh.o hycrt.o hystl.o hyopen.o hyinit.o hyphs.o hycal.o \
	hycin.o hytrl.o hyloc.o hysol.o hysvd.o hytra.o hytrh.o hytrt.o \
	hymag.o hyrep.o hylst.o hysum.o hyinp.o hypro.o medwt.o utmcal.o \
	init_event.o hydelt.o hyedit.o hytime.o hybeg.o hycout.o \
	hymagp.o hysou.o hyfile.o \
	-Bstatic allsubs.o -o hyp

# The -Bstatic option loads all fortran system subroutines into the executable
# so hyp can be copied to a similar system without fortran to run without
# the fortran system libraries.

# line for loading all "common" subroutines in a library using f77:
#	-Bstatic -L/home/calnet/klein/subs -lfs -o hyp
# this is no longer used because it is not portable across systems & compilers

#  allsubs - collection of general subroutines into one file
#  This is the $(F77) $(FFLAGS) version of readq subroutine (no q format)

allsubs.o  :  allsubs.f 
	$(F77) $(FFLAGS) -c allsubs.f

#  HYP	   * MAIN HYPOINVERSE PROGRAM.
hyp.o  :  hyp.for common.inc
	$(F77) $(FFLAGS) -c hyp.for

#  HYBDA   * BLOCK DATA INITIALIZATION OF COMMON.
hybda.o  :  hybda.for common.inc
	$(F77) $(FFLAGS) -c hybda.for

#  HYBEG   * INITIALIZATION OF OTHER VARIABLES.
hybeg.o  :  hybeg.f common.inc
	$(F77) $(FFLAGS) -c hybeg.f

#  HYCMD   * GETS AND PROCESSES COMMANDS.
hycmd.o  :  hycmd.for common.inc
	$(F77) $(FFLAGS) -c hycmd.for

#  HYSTA   * READS IN STATIONS.
hysta.o  :  hysta.for common.inc
	$(F77) $(FFLAGS) -c hysta.for

#  HYDEL   * READS IN STATION DELAYS (FOR MULTIPLE MODELS).
hydel.o  :  hydel.for common.inc
	$(F77) $(FFLAGS) -c hydel.for

#  HYATE   * READS IN STATION ATTENUATION HISTORY.
hyate.o  :  hyate.for common.inc
	$(F77) $(FFLAGS) -c hyate.for

#  HYCAL   * READS IN STATION CAL FACTOR HISTORY.
hycal.o  :  hycal.for common.inc
	$(F77) $(FFLAGS) -c hycal.for

#  HYFMC   * READS IN STATION FMAG CORRECTIONS.
hyfmc.o  :  hyfmc.for common.inc
	$(F77) $(FFLAGS) -c hyfmc.for

#  HYFILE   * DETERMINES INPUT FILE TYPE.
hyfile.o  :  hyfile.for
	$(F77) $(FFLAGS) -c hyfile.for

#  HYXMC   * READS IN STATION XMAG CORRECTIONS.
hyxmc.o  :  hyxmc.for common.inc
	$(F77) $(FFLAGS) -c hyxmc.for

#  HYCRH   * READS IN HOMOGENOUS LAYER CRUSTAL MODELS.
hycrh.o  :  hycrh.for common.inc
	$(F77) $(FFLAGS) -c hycrh.for

#  HYCRT   * READS IN TRAVEL-TIME-TABLE CRUSTAL MODELS.
hycrt.o  :  hycrt.for common.inc
	$(F77) $(FFLAGS) -c hycrt.for

#  HYSTL   * OUTPUTS STATIONS, CRUST & PARAMETERS TO PRINT FILE.
hystl.o  :  hystl.for common.inc
	$(F77) $(FFLAGS) -c hystl.for

#  HYOPEN  * OPENS FILES FOR LOCATION RUN.
hyopen.o  :  hyopen.for common.inc
	$(F77) $(FFLAGS) -c hyopen.for

#  HYINIT  * INITIALIZES SOME VARIABLES FOR LOCATION RUN.
hyinit.o  :  hyinit.for common.inc
	$(F77) $(FFLAGS) -c hyinit.for

#  HYPHS   * READS IN PHASE DATA FOR ONE EVENT.
hyphs.o  :  hyphs.for common.inc
	$(F77) $(FFLAGS) -c hyphs.for

#  HYCIN   * INPUTS PHASE DATA FROM CUSP MEM FILES (ALTERNATE TO HYPHS)
hycin.o  :  hycin.f
	$(F77) $(FFLAGS) -c hycin.f

#  HYCOUT   * OUTPUTS PHASE DATA TO CUSP MEM FILES (ALTERNATE TO HYARC)
hycout.o  :  hycout.f
	$(F77) $(FFLAGS) -c hycout.f

#  HYTRL   * SETS TRIAL HYPOCENTER.
hytrl.o  :  hytrl.for common.inc
	$(F77) $(FFLAGS) -c hytrl.for

#  HYLOC   * LOCATES ONE EVENT.
hyloc.o  :  hyloc.for common.inc
	$(F77) $(FFLAGS) -c hyloc.for

#  HYSOL   * PERFORMS INVERSION FOR ONE ITERATION.
hysol.o  :  hysol.for common.inc
	$(F77) $(FFLAGS) -c hysol.for

#  HYSVD     CANNED SINGLE-VALUE-DECOMPOSITION ROUTINE.
hysvd.o  :  hysvd.f
	$(F77) $(FFLAGS) -c hysvd.f

#  HYTRA   * MANAGE CRUST MODEL CHOICE & AVERAGING.
hytra.o  :  hytra.for common.inc
	$(F77) $(FFLAGS) -c hytra.for

#  HYTRH   * CALC TRAVEL TIMES AND DERIVS FOR HOMO LAYER MODEL.
hytrh.o  :  hytrh.for common.inc
	$(F77) $(FFLAGS) -c hytrh.for

#  HYTRT   * CALC TRAVEL TIMES AND DERIVS FROM TRAV-TIME TABLE.
hytrt.o  :  hytrt.for common.inc
	$(F77) $(FFLAGS) -c hytrt.for

#  HYMAG   * COMPUTES DUR & MAX AMP MAGNITUDE.
hymag.o  :  hymag.for common.inc
	$(F77) $(FFLAGS) -c hymag.for

#  HYMAGP   * COMPUTES P AMP MAGNITUDE.
hymagp.o  :  hymagp.for common.inc
	$(F77) $(FFLAGS) -c hymagp.for

#  HYPREF  * SELECTS PREFERRED MAGNITUDE.
hypref.o  :  hypref.for common.inc
	$(F77) $(FFLAGS) -c hypref.for

#  HYREP   * REPORTS A LOCATION ON THE TERMINAL.
hyrep.o  :  hyrep.for common.inc
	$(F77) $(FFLAGS) -c hyrep.for

#  HYSOU   * TABULATES MOST COMMON DATA SOURCES.
hysou.o  :  hysou.for common.inc
	$(F77) $(FFLAGS) -c hysou.for

#  HYLST   * OUTPUTS DATA BY STATION TO PRINT & ARCHIVE FILES.
hylst.o  :  hylst.for common.inc
	$(F77) $(FFLAGS) -c hylst.for

#  HYSUM   * OUTPUTS SUMMARY RECORD (FINAL LOCATION).
hysum.o  :  hysum.for common.inc
	$(F77) $(FFLAGS) -c hysum.for

#  HYINP     FOR INTERACTIVE ENTRY OF PHASE DATA.
hyinp.o  :  hyinp.for
	$(F77) $(FFLAGS) -c hyinp.for

#  HYPRO   * INTERACTIVELY PROCESSES A SERIES OF EVENTS.
hypro.o  :  hypro.for common.inc
	$(F77) $(FFLAGS) -c hypro.for

#  MEDWT     COMPUTES THE WEIGHTED MEDIAN OF A SERIES OF MAGNITUDES.
medwt.o  :  medwt.for
	$(F77) $(FFLAGS) -c medwt.for

#  UTMCAL    CALCULATE STATION DISTANCES ON A UTM GRID
utmcal.o :  utmcal.for
	$(F77) $(FFLAGS) -c utmcal.for

#  INIT_EVENT  INITIALIZES CUSP FOR READING OF MEM DATA FILES (VAX ONLY)
init_event.o  :  init_event.f
	$(F77) $(FFLAGS) -c init_event.f

#  HYDELETE  DELETES FILES IN INTERACTIVE PROCESSING
hydelt.o  :  hydelt.f
	$(F77) $(FFLAGS) -c hydelt.f

#  HYEDTIOR  RUNS AN EDTIOR WITHIN A PROCESS
hyedit.o  :  hyedit.f
	$(F77) $(FFLAGS) -c hyedit.f

#  HYTIME    GETS CURRENT SYSTEM TIME FOR LABELING PRINT FILE
hytime.o  :  hytime.f
	$(F77) $(FFLAGS) -c hytime.f

clean:	
	rm *.o hyp

install:
	cp hyp $(BDIR)



#  ALSO USES THESE GENERAL PURPOSE SUBROUTINES FROM THE SUBROUTINE LIBRARY
#  KLAS      ASSIGNS A NAME AND NUMBER TO AN EVENT BASED ON LOCATION.
#  KLASS (NET 1), BOX2 (NET2) - USED BY KLAS.
#  ISITIN    DETERMINES IF EVENT IS IN A PARTICULAR REGION.
#  UPSTR     CONVERTS A STRING TO UPPER CASE.
#  JASK	     INTERACTIVE PROMPT & ENTRY OF AN INTEGER.
#  ASKC      INTERACTIVE PROMPT AND ENTRY OF A STRING.
#  ASKR      INTERACTIVE PROMPT AND ENTRY OF A REAL VALUE.
#  LASK      INTERACTIVE PROMPT AND ENTRY OF A LOGICAL VALUE.
#  LENG      DETERMINES THE NON-BLANK LENGTH OF A STRING.

