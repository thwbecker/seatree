      SUBROUTINE HYSOU
C--TABULATES THE MOST COMMON DATA SOURCE CODES FOR HYPOINVERSE PRIOR TO OUTPUT
C--THIS CODE USED TO BE IN HYLST
      INCLUDE 'common.inc'
      DIMENSION NSOU(5),NSOUF(5),NSOUX(5),NSOUF2(5),NSOUX2(5)
      DIMENSION NSOUP(5),NSOUP2(5)
      CHARACTER CH(5)*1,CHF(5)*1,CHX(5)*1,CHF2(5)*1,CHX2(5)*1
      CHARACTER CHP(5)*1,CHP2(5)*1

C--GET THE PRIMARY DATA SOURCES FOR WEIGHTED PHASES, FMAGS, XMAGS
      DO I=1,5
        NSOU(I)=0
        NSOUF(I)=0
        NSOUX(I)=0
        NSOUF2(I)=0
        NSOUX2(I)=0
        NSOUP(I)=0
        NSOUP2(I)=0
      END DO
      XMSOU=' '
      FMSOU=' '
      XMSOU2=' '
      FMSOU2=' '
      PSOUR=' '
      PSOUR2=' '
      SOUCOD=' '

C--CHECK ALL PHASES FOR PRIMARY DATA SOURCE
C--KEEP COUNTS OF 7 DIFFERENT SOURCE CODES
C--P OR S ARRIVAL TIME
      DO 11 IM=1,M
        IF (W(IM).LT..1) GOTO 11
C--DETERMINE STATION INDEX & REMOVE S FLAG
        K=IND(IM)
        KPS=K/10000
        K=K-10000*KPS
        DO I=1,5
          IF (NSOU(I).EQ.0) THEN
            NSOU(I)=1
            CH(I)=KSOU(K)
            GOTO 11
          END IF
          IF (KSOU(K).EQ.CH(I)) THEN
            NSOU(I)=NSOU(I)+1
            GOTO 11
          END IF
        END DO
11    CONTINUE

C--COUNT FMAG AND XMAG SOURCE CODES FROM WEIGHTED READINGS FOR ALL STATIONS
      DO 15 K=1,KSTA
        J=KINDX(K)
C--COUNT FIRST FMAG SOURCE CODES
        IF (KFWT(K).GT.3 .OR. KFMAG(K).EQ.0 .OR. .NOT.JFM1(J)) GOTO 12
        DO I=1,5
          IF (NSOUF(I).EQ.0) THEN
            NSOUF(I)=1
            CHF(I)=KSOU(K)
            GOTO 12
          END IF
          IF (KSOU(K).EQ.CHF(I)) THEN
            NSOUF(I)=NSOUF(I)+1
            GOTO 12
          END IF
        END DO

C--COUNT FIRST XMAG SOURCE CODES 
12      IF (KXWT(K).GT.3 .OR. KXMAG(K).EQ.0 .OR. .NOT.JXM1(J)) GOTO 13
        DO I=1,5
          IF (NSOUX(I).EQ.0) THEN
            NSOUX(I)=1
            CHX(I)=KSOU(K)
            GOTO 13
          END IF
          IF (KSOU(K).EQ.CHX(I)) THEN
            NSOUX(I)=NSOUX(I)+1
            GOTO 13
          END IF
        END DO

C--COUNT SECOND FMAG SOURCE CODES
13      IF (KFWT(K).GT.3 .OR. KFMAG(K).EQ.0 .OR. .NOT.JFM2(J)) GOTO 14
        DO I=1,5
          IF (NSOUF2(I).EQ.0) THEN
            NSOUF2(I)=1
            CHF2(I)=KSOU(K)
            GOTO 14
          END IF
          IF (KSOU(K).EQ.CHF2(I)) THEN
            NSOUF2(I)=NSOUF2(I)+1
            GOTO 14
          END IF
        END DO

C--COUNT SECOND XMAG SOURCE CODES 
14      IF (KXWT(K).GT.3 .OR. KXMAG(K).EQ.0 .OR. .NOT.JXM2(J)) GOTO 15
        DO I=1,5
          IF (NSOUX2(I).EQ.0) THEN
            NSOUX2(I)=1
            CHX2(I)=KSOU(K)
            GOTO 15
          END IF
          IF (KSOU(K).EQ.CHX2(I)) THEN
            NSOUX2(I)=NSOUX2(I)+1
            GOTO 15
          END IF
        END DO
15    CONTINUE

C--COUNT PAMAG SOURCE CODES FROM WEIGHTED READINGS FOR ALL STATIONS
      IF (LPMAG) THEN
      DO 17 K=1,KSTA
        J=KINDX(K)
C--COUNT PRIMARY P AMP MAG SOURCE CODES 
        IF (PAWT(K).LT..1 .OR. KPAMP(K).EQ.0 .OR. .NOT.JPM1(J)) GOTO 16
        DO I=1,5
          IF (NSOUP(I).EQ.0) THEN
            NSOUP(I)=1
            CHP(I)=KSOU(K)
            GOTO 16
          END IF
          IF (KSOU(K).EQ.CHP(I)) THEN
            NSOUP(I)=NSOUP(I)+1
            GOTO 16
          END IF
        END DO

C--COUNT SECONDARY P AMP MAG SOURCE CODES 
16      IF (PAWT(K).LT..1 .OR. KPAMP(K).EQ.0 .OR. .NOT.JPM2(J)) GOTO 17
        DO I=1,5
          IF (NSOUP2(I).EQ.0) THEN
            NSOUP2(I)=1
            CHP2(I)=KSOU(K)
            GOTO 17
          END IF
          IF (KSOU(K).EQ.CHP2(I)) THEN
            NSOUP2(I)=NSOUP2(I)+1
            GOTO 17
          END IF
        END DO

17    CONTINUE
      END IF

C--REPORT THE MOST FREQUENT DATA SOURCE CODES
      K1=0
      K2=0
      KF2=0
      KX2=0
      KP1=0
      KP2=0
      IM=0
      DO I=1,5
        IF (NSOU(I).GT.IM) THEN
          IM=NSOU(I)
          SOUCOD=CH(I)
        END IF

        IF (NSOUF(I).GT.K1) THEN
          K1=NSOUF(I)
          FMSOU=CHF(I)
        END IF

        IF (NSOUX(I).GT.K2) THEN
          K2=NSOUX(I)
          XMSOU=CHX(I)
        END IF

        IF (NSOUF2(I).GT.KF2) THEN
          KF2=NSOUF2(I)
          FMSOU2=CHF2(I)
        END IF

        IF (NSOUX2(I).GT.KX2) THEN
          KX2=NSOUX2(I)
          XMSOU2=CHX2(I)
        END IF

        IF (NSOUP(I).GT.KP1) THEN
          KP1=NSOUP(I)
          PSOUR=CHP(I)
        END IF

        IF (NSOUP2(I).GT.KP2) THEN
          KP2=NSOUP2(I)
          PSOUR2=CHP2(I)
        END IF
      END DO

      RETURN
      END
