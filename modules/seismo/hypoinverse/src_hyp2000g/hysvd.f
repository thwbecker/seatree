      SUBROUTINE  HYSVD  (A, S, V, MMAX, NMAX, M, N, P, WITHU, WITHV,T)

C--NOTE THAT FOR USE IN HYPOINVERSE:
C-(1) THE TEMPORARY T IS DIMENSIONED T(MMAX) IN THE CALLING ROUTINE.
C-(2) THE NUMBER OF COLUMNS OF A IS NPMAX = NMAX + PMAX. SINCE PMAX
C-IS 1 AND NMAX IS 4, NPMAX = 5.
      PARAMETER (NPMAX=5)

      INTEGER    MMAX, NMAX, M, N, P
      REAL       R, W, CS, SN, TOL, F, X, EPS, G, T, Y
      REAL       ETA, H, Q, Z
      INTEGER    I, J, K, L, L1, N1, NP
      LOGICAL    WITHU, WITHV
      DIMENSION       S(NMAX), V(NMAX,NMAX), T(MMAX)
      DIMENSION A(MMAX,NPMAX)

C     ------------------------------------------------------------------

C     THIS IS A TRANSLATION OF A CDC 6600 FORTRAN PROGRAM TO IBM 360
C     FORTRAN IV.  THIS SUBROUTINE USES SHORT PRECISION ARITHMETIC.
C     A LONG PRECISION VERSION IS AVAILABLE UNDER THE NAME 'DSVD'.

C     THIS SUBROUTINE REPLACES EARLIER SUBROUTINES WITH THE SAME NAME,
C    689   6       &   &0&  S OF A COMPLEX AR&THMETIC PROGRAM, PUBLISHED
C     AS ALGORITHM 358.  THIS CURRENT PROGRAM IS FASTER, MORE ACCURATE
C     AND LESS OBSCURE IN DESCRIBING ITS CAPABILITIES.

C     ORIGINAL PROGRAMMER=  R. C. SINGLETON
C     360 VERSION BY=       J. G. LEWIS
C     LAST REVISION OF THIS SUBROUTINE=  4 DECEMBER 1973

C     ------------------------------------------------------------------

C     ADDITIONAL SUBROUTINE NEEDED=  ROTATE
C--THE HYROT SUBROUTINE WAS ELIMINATED TO PERMIT USE OF VIRTUAL ARRAYS

C     ------------------------------------------------------------------


C     THIS SUBROUTINE COMPUTES THE SINGULAR VALUE DECOMPOSITION
C     OF A REAL M*N MATRIX A, I.E. IT COMPUTES MATRICES U, S, AND V
C     SUCH THAT
C
C                  A = U * S * VT ,
C     WHERE
C              U IS AN M*N MATRIX AND UT*U = I, (UT=TRANSPOSE
C                                                    OF U),
C              V IS AN N*N MATRIX AND VT*V = I, (VT=TRANSPOSE
C                                                    OF V),
C        AND   S IS AN N*N DIAGONAL MATRIX.
C
C     DESCRIPTION OF PARAMETERS=
C
C     A = REAL ARRAY. A CONTAINS THE MATRIX TO BE DECOMPOSED.
C         THE ORIGINAL DATA ARE LOST.  IF WITHV=.TRUE., THEN
C         THE MATRIX U IS COMPUTED AND STORED IN THE ARRAY A.
C
C     MMAX = INTEGER VARIABLE.  THE NUMBER OF ROWS IN THE
C            ARRAY A.
C
C     NMAX = INTEGER VARIABLE.  THE NUMBER OF ROWS IN THE
C            ARRAY V.
C
C     M,N = INTEGER VARIABLES.  THE NUMBER OF ROWS AND COLUMNS
C           IN THE MATRIX STORED IN A.  (N<=M<=100.  IF IT IS
C           NECESSARY TO SOLVE A LARGER PROBLEM, THEN THE
C           AMOUNT OF STORAGE ALLOCATED TO THE ARRAY T MUST
C           BE INCREASED ACCORDINGLY.)  IF M<N , THEN EITHER
C           TRANSPOSE THE MATRIX A OR ADD ROWS OF ZEROS TO
C           INCREASE M TO N.
C
C     P = INTEGER VARIABLE.  IF P>0, THEN COLUMNS N+1, . . . ,
C         N+P OF A ARE ASSUMED TO CONTAIN THE COLUMNS OF AN M*P
C         MATRIX B.  THIS MATRIX IS MULTIPLIED BY UT, AND UPON
C         EXIT, A CONTAINS IN THESE SAME COLUMNS THE N*P MATRIX
C         UT*B. (P>=0)
C
C     WITHU, WITHV = LOGICAL VARIABLES.  IF WITHU=.TRUE., THEN
C         THE MATRIX U IS COMPUTED AND STORED IN THE ARRAY A.
C         IF WITHV=.TRUE., THEN THE MATRIX V IS COMPUTED AND
C         STORED IN THE ARRAY V.
C
C     S = REAL ARRAY.  S(1), . . . , S(N) CONTAIN THE DIAGONAL
C         ELEMENTS OF THE MATRIX S ORDERED SO THAN S(I)>=S(I+1),
C         I=1, . . . , N-1.
C
C     V = REAL ARRAY.  V CONTAINS THE MATRIX V.  IF WITHU
C         AND WITHV ARE NOT BOTH =.TRUE., THEN THE ACTUAL
C         PARAMETER CORRESPONDING TO A AND V MAY BE THE SAME.
C
C     THIS SUBROUTINE IS A REAL VERSION OF A FORTRAN SUBROUTINE
C     BY BUSINGER AND GOLUB, ALGORITHM 358=  SINGULAR VALUE
C     DECOMPOSITION OF A COMPLEX MATRIX, COMM. ACM, V. 12,
C     NO. 10, PP. 564-565 (OCT. 1969).
C     WITH REVISIONS BY RC SINGLETON, MAY 1972.
C     ------------------------------------------------------------------


C    ECLIPSE CONSTANTS
      SAVE ETA,TOL
      DATA ETA,TOL/1.9E-7,1.1E-30/

C     DATA ETA,TOL /1.5E-15,1.E-250/


C     ETA (16**-6) AND TOL (16**-59) ARE MACHINE DEPENDENT CONSTANTS
C     FOR IBM 360/370 COMPUTERS (SHORT FORM ARITHMETIC).
C     ETA IS THE MACHINE EPSILON (RELATIVE ACCURACY)!
C     TOL IS THE SMALLEST REPRESENTABLE REAL DIVIDED BY ETA.

      NP = N + P
      N1 = N + 1

C     HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM
      G = 0.0
      EPS = 0.0
      L = 1
   10 T(L) = G
      K = L
      L = L + 1

C     ELIMINATION OF A(I,K), I=K+1, . . . , M
      S(K) = 0.0
      Z = 0.0
      DO 20 I = K,M
   20   Z = Z + A(I,K)**2
      IF (Z.LT.TOL) GOTO 50
      G = SQRT(Z)
      F = A(K,K)
      IF (F.GE.0.0) G = - G
      S(K) = G
      H = G * (F - G)
      A(K,K) = F - G
      IF (K.EQ.NP) GOTO 50
      DO 40 J = L,NP
         F = 0
         DO 30 I = K,M
   30       F = F + A(I,K)*A(I,J)
         F = F/H
         DO 40 I = K,M
   40       A(I,J) = A(I,J) + F*A(I,K)

C     ELIMINATION OF A(K,J), J=K+2, . . . , N
   50 EPS = AMAX1(EPS,ABS(S(K)) + ABS(T(K)))
      IF (K.EQ.N) GOTO 100
      G = 0.0
      Z = 0.0
      DO 60 J = L,N
   60    Z = Z + A(K,J)**2
      IF (Z.LT.TOL) GOTO 10
      G = SQRT(Z)
      F = A(K,L)
      IF (F.GE.0.0) G = - G
      H = G * (F - G)
      A(K,L) = F - G
      DO 70 J = L,N
   70    T(J) = A(K,J)/H
      DO 90 I = L,M
         F = 0
         DO 80 J = L,N
   80       F = F + A(K,J)*A(I,J)
         DO 90 J = L,N
   90       A(I,J) = A(I,J) + F*T(J)

      GOTO 10

C     TOLERANCE FOR NEGLIGIBLE ELEMENTS
  100 EPS = EPS*ETA

C     ACCUMULATION OF TRANSFORMATIONS
      IF (.NOT.WITHV) GOTO 160
      K = N
      GOTO 140
  110    IF (T(L).EQ.0.0) GOTO 140
         H = A(K,L)*T(L)
         DO 130 J = L,N
            Q = 0
            DO 120 I = L,N
  120          Q = Q + A(K,I)*V(I,J)
            Q = Q/H
            DO 130 I = L,N
  130          V(I,J) = V(I,J) + Q*A(K,I)
  140    DO 150 J = 1,N
  150       V(K,J) = 0
         V(K,K) = 1.0
         L = K
         K = K - 1
         IF (K.NE.0) GOTO 110

  160  K = N
       IF (.NOT.WITHU) GOTO 230
       G = S(N)
       IF (G.NE.0.0) G = 1.0/G
       GO TO 210
  170    DO 180 J = L,N
  180       A(K,J) = 0
         G = S(K)
         IF (G.EQ.0.0) GOTO 210
         H = A(K,K)*G
         DO 200 J = L,N
            Q = 0
            DO 190 I = L,M
  190          Q = Q + A(I,K)*A(I,J)
            Q = Q/H
            DO 200 I = K,M
  200          A(I,J) = A(I,J) + Q*A(I,K)
         G = 1.0/G
  210    DO 220 J = K,M
  220       A(J,K) = A(J,K)*G
         A(K,K) = A(K,K) + 1.0
         L = K
         K = K - 1
         IF (K.NE.0) GOTO 170

C     QR DIAGONALIZATION
      K = N

C     TEST FOR SPLIT
  230    L = K
  240       IF (ABS(T(L)).LE.EPS) GOTO 290
            L = L - 1
            IF (ABS(S(L)).GT.EPS) GOTO 240

C     CANCELLATION
         CS = 0.0
         SN = 1.0
         L1 = L
         L = L + 1
         DO 280 I = L,K
            F = SN*T(I)
            T(I) = CS*T(I)
            IF (ABS(F).LE.EPS) GOTO 290
            H = S(I)
            W = SQRT(ABS(F*F + H*H))
            S(I) = W
            CS = H/W
            SN = - F/W
C           IF (WITHU) CALL HYROT(A(1,L1), A(1,I), CS, SN, M)
            IF (WITHU) THEN
              DO 260 JJ=1,M
                XX=A(JJ,L1)
                YY=A(JJ,I)
                A(JJ,L1)=XX*CS + YY*SN
260             A(JJ,I)= YY*CS - XX*SN
            END IF
            IF (NP.EQ.N) GOTO 280
            DO 270 J = N1,NP
               Q = A(L1,J)
               R = A(I,J)
               A(L1,J) = Q*CS + R*SN
 270           A(I,J) = R*CS - Q*SN
  280     CONTINUE

C     TEST FOR CONVERGENCE
  290    W = S(K)
         IF (L.EQ.K) GOTO 360

C     ORIGIN SHIFT
         X = S(L)
         Y = S(K-1)
         G = T(K-1)
         H = T(K)
         F = ((Y - W)*(Y + W) + (G - H)*(G + H))/(2.0*H*Y)
         G = SQRT(F*F + 1.0)
         IF (F.LT.0.0) G = - G
         F = ((X - W)*(X + W) + (Y/(F + G) - H)*H)/X

C     QR STEP
         CS = 1.0
         SN = 1.0
         L1 = L + 1
         DO 350 I = L1,K
            G = T(I)
            Y = S(I)
            H = SN*G
            G = CS*G
            W = SQRT(H*H + F*F)
            T(I-1) = W
            CS = DIV (F,W,TOL)
            SN = DIV (H,W,TOL)
            F = X*CS + G*SN
            G = G*CS - X*SN
            H = Y*SN
            Y = Y*CS
C           IF (WITHV) CALL HYROT(V(1,I-1), V(1,I), CS, SN, N)
          IF (WITHV) THEN
            DO 320 JJ=1,N
            XX=V(JJ,I-1)
            YY=V(JJ,I)
            V(JJ,I-1)=XX*CS + YY*SN
320            V(JJ,I)=  YY*CS - XX*SN
          END IF
            W = SQRT(ABS(H*H + F*F))
            S(I-1) = W
            CS = DIV (F,W,TOL)
            SN = DIV (H,W,TOL)
            F = CS*G + SN*Y
            X = CS*Y - SN*G
C           IF (WITHU) CALL HYROT(A(1,I-1), A(1,I), CS, SN, M)
          IF (WITHU) THEN
            DO 330 JJ=1,M
            XX=A(JJ,I-1)
            YY=A(JJ,I)
            A(JJ,I-1)=XX*CS + YY*SN
330            A(JJ,I)=  YY*CS - XX*SN
          END IF
            IF (N.EQ.NP) GOTO 350
            DO 340 J = N1,NP
               Q = A(I-1,J)
               R = A(I,J)
               A(I-1,J) = Q*CS + R*SN
 340  A(I,J) = R*CS - Q*SN
  350       CONTINUE

         T(L) = 0.0
         T(K) = F
         S(K) = X
         GOTO 230

C     CONVERGENCE
  360    IF (W.GE.0.0) GOTO 380
         S(K) = - W
         IF (.NOT.WITHV) GOTO 380
         DO 370 J = 1,N
  370       V(J,K) = - V(J,K)
  380    K = K - 1
         IF (K.NE.0) GO TO 230

C     SORT SINGULAR VALUES
      DO 450 K = 1,N
         G = -1.0
         DO 390 I = K,N
            IF (S(I).LT.G) GOTO 390
            G = S(I)
            J = I
  390       CONTINUE
      IF (J .EQ. K) GO TO 450
      S(J) = S(K)
         S(K) = G
         IF (.NOT.WITHV) GOTO 410
         DO 400 I = 1,N
            Q = V(I,J)
            V(I,J) = V(I,K)
  400       V(I,K) = Q
  410    IF (.NOT.WITHU) GOTO 430
         DO 420 I = 1,M
            Q = A(I,J)
            A(I,J) = A(I,K)
  420       A(I,K) = Q
  430    IF (N.EQ.NP) GOTO 450
         DO 440 I = N1,NP
            Q = A(J,I)
            A(J,I) = A(K,I)
  440       A(K,I) = Q
  450    CONTINUE

      RETURN
      END
C      SUBROUTINE   HYROT  (X, Y, CS, SN, N)
C      INTEGER N
C      REAL    X(N), Y(N), CS, SN
C      REAL    XX
C      INTEGER J
C      DO 10 J = 1, N
C         XX = X(J)
C         X(J) = XX*CS + Y(J)*SN
C   10    Y(J) = Y(J)*CS - XX*SN
C      RETURN
C      END

      FUNCTION DIV (A,B,TOL)
C--AVOID DIVION BY 0, AND UNDEFINED 0/0 DIVIOSN
C--DIV RETURNS 0 FOR 0/0 INPUT
      IF (A.EQ.0.) THEN
        DIV=0.
        RETURN
      END IF
      IF (B.EQ.0.) THEN
        DIV=SIGN(TOL,A)
        RETURN
      END IF
      DIV=A/B
      RETURN
      END
