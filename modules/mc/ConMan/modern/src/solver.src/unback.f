      SUBROUTINE UNBACK(A,C,B,IDIAG,NEQ)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),C(*),B(*),IDIAG(*)
C
C     FORWARD REDUCTION
C
      JR=0
      DO 100 J=1,NEQ
      JD=IDIAG(J)
      JH=JD-JR
      IF(JH .LE. 1) GO TO 50
      IS=J+1-JH
      B(J)=B(J)-COLDOT(C(JR+1),B(IS),JH-1)
   50 JR=JD
  100 CONTINUE
C
C     BACK SUBSTITUTION
C
      J=NEQ
      JD=IDIAG(J)
  200 CONTINUE
      B(J)=B(J)/A(JD)
      D=B(J)
      J=J-1
      IF(J.EQ.0) RETURN
      JR=IDIAG(J)
      JH=JD-JR
      IF(JH .LE. 1) GO TO 400
      IS=J-JH+2
      K=JR-IS+1
      DO 300 I=IS,J
      B(I)=B(I)-A(I+K)*D
  300 CONTINUE
  400 CONTINUE
      JD=JR
      GO TO 200
      END
