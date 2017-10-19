      SUBROUTINE HYBEG
C--BEGIN HYPOINVERSE BY INITIALIZING CERTAIN STRINGS.
C--THIS IS EQUIVALENT TO INITIALIZATION IN A DATA STATEMENT.
      INCLUDE 'common.inc'

C--SUN/UNIX
       TERMIN='/dev/tty'
       INFILE(1)='hypinst'
       INFILE(0)='/home/calnet/klein/hypfiles/cal2000.hyp'
       STAFIL='/home/calnet/klein/hypfiles/all2.sta'
C       ATNFIL='/home/calnet/klein/hypfiles/all2.atn'
C       FMCFIL='/home/calnet/klein/hypfiles/all2.fmc'
C       XMCFIL='/home/calnet/klein/hypfiles/all2.xmc'
       BSTAFL='/home/calnet/klein/hypfiles/allsta2.bin'
       BCRUFL='/home/calnet/klein/hypfiles/multmod2.bin'
       ATNFIL=' '
       FMCFIL=' '
       XMCFIL=' '

C--VAX
C      TERMIN='TT:'
C      INFILE(1)='HYPINST.'
C      INFILE(0)='HOME:[KLEIN.HYPFILES]CAL2000.HYP'
CC      INFILE(0)='HYPOINV$MODELS'            !CUSP LOGICAL NAME

C--STATION DATA
      DO I=1,MAXSTA
        JCEXP(I)=0
        JLMOD(I)=.FALSE.
      END DO

C--PHASE DATA
      DO I=1,MAXPHS
        KRMK6(I)=' '
      END DO
      RETURN
      END
