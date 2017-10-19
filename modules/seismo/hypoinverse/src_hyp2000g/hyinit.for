      SUBROUTINE HYINIT
C--CALLED BY HYPOINV TO INITIALIZE SOME VALUES BEFORE EACH LOCATION RUN
      INCLUDE 'common.inc'

C--INITIALIZE THE SHADOW RECORDS
      LENSHA=0
      SHADO=' '
      DO I=1,MSHA
        LSHA1(I)=0
        SHAD1(I)=' '
      END DO
      DO K=1,MAXPHS
        KLSHA(K)=0
          KPAMP(K)=0
        KSHAD(K)=' '
      END DO

C--INITIALIZE SOME VALUES
      REMK=' '
      INUM=0
      LTBIG=.FALSE.
C--SET FLAG USED TO INDICATE END OF PHASE CARD FILE
      KEND=0

C--SET FLAGS FOR EACH STATION FOR MAGNITUDES TO BE COMPUTED
C--FMAG
      DO J=1,JSTA
        IF (NCPF1.LT.0) THEN
          JFM1(J)=.TRUE.
        ELSE
          JFM1(J)=.FALSE.
          DO I=1,NCPF1
            IF (COMPF1(I)(1:NCOMP) .EQ. JCOMP3(J)(1:NCOMP)) THEN
              JFM1(J)=.TRUE.
              GOTO 11
            END IF
          END DO
        END IF

C--FMAG2
11      IF (NCPF2.LT.0) THEN
          JFM2(J)=.TRUE.
        ELSE
          JFM2(J)=.FALSE.
          DO I=1,NCPF2
            IF (COMPF2(I)(1:NCOMP) .EQ. JCOMP3(J)(1:NCOMP)) THEN
              JFM2(J)=.TRUE.
              GOTO 12
            END IF
          END DO
        END IF

C--XMAG
12      IF (LXCH) THEN
C--CHOOSE XMAG BY COMPONENT
          IF (NCPX1.LT.0) THEN
            JXM1(J)=.TRUE.
          ELSE
            JXM1(J)=.FALSE.
            DO I=1,NCPX1
              IF (COMPX1(I)(1:NCOMP) .EQ. JCOMP3(J)(1:NCOMP)) THEN
                JXM1(J)=.TRUE.
                GOTO 13
              END IF
            END DO
          END IF

        ELSE
C--CHOOSE XMAG BY INSTRUMENT TYPE
          IF (NXTYP1.LT.0) THEN
            JXM1(J)=.TRUE.
          ELSE
            JXM1(J)=.FALSE.
            DO I=1,NXTYP1
              IF (IXTYP1(I) .EQ. JTYPE(J)) THEN
                JXM1(J)=.TRUE.
                GOTO 13
              END IF
            END DO
          END IF
        END IF

C--XMAG2
13      IF (LXCH) THEN
C--CHOOSE XMAG2 BY COMPONENT
          IF (NCPX2.LT.0) THEN
            JXM2(J)=.TRUE.
          ELSE
            JXM2(J)=.FALSE.
            DO I=1,NCPX2
              IF (COMPX2(I)(1:NCOMP) .EQ. JCOMP3(J)(1:NCOMP)) THEN
                JXM2(J)=.TRUE.
                GOTO 14
              END IF
            END DO
          END IF

        ELSE
C--CHOOSE XMAG2 BY INSTRUMENT TYPE
          IF (NXTYP2.LT.0) THEN
            JXM2(J)=.TRUE.
          ELSE
            JXM2(J)=.FALSE.
            DO I=1,NXTYP2
              IF (IXTYP2(I) .EQ. JTYPE(J)) THEN
                JXM2(J)=.TRUE.
                GOTO 14
              END IF
            END DO
          END IF
        END IF

C--PAMAG
14      IF (NCPP1.LT.0) THEN
          JPM1(J)=.TRUE.
        ELSE
          JPM1(J)=.FALSE.
          DO I=1,NCPP1
            IF (COMPP1(I)(1:NCOMP) .EQ. JCOMP3(J)(1:NCOMP)) THEN
              JPM1(J)=.TRUE.
              GOTO 15
            END IF
          END DO
        END IF

C--PAMAG2
15      IF (NCPP2.LT.0) THEN
          JPM2(J)=.TRUE.
        ELSE
          JPM2(J)=.FALSE.
          DO I=1,NCPP2
            IF (COMPP2(I)(1:NCOMP) .EQ. JCOMP3(J)(1:NCOMP)) THEN
              JPM2(J)=.TRUE.
              GOTO 16
            END IF
          END DO
        END IF

      END DO
16    RETURN
      END
