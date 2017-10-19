	subroutine rsoinc(A,N,IDX)
c--sorts array A in order of increasing entry value
c--integer array idx remembers how entries have been swapped
	implicit real*8 (a-h,o-z)
      DIMENSION A(1),IDX(1) 
      IF (N.EQ.1) GO TO 65
      IF (N.LE.0) GO TO 60
      DO 1 I = 1,N
      IDX(I) = I
    1 CONTINUE
      N2 = N/2
      N21 = N2 + 2
      ICT=1 
      I=2 
   11 N1=N21-I
      NN=N
      IK=N1 
   15 C=A(IK) 
      IC=IDX(IK)
  100 JK=2*IK 
      IF (JK.GT.NN) GO TO 140 
      IF (JK.EQ.NN) GO TO 120 
       IF (A(JK+1).LE.A(JK)) GO TO 120
      JK=JK+1 
  120 IF (A(JK).LE. C) GO TO 140
      A(IK)=A(JK) 
      IDX(IK)=IDX(JK) 
      IK=JK 
      GO TO 100 
  140 A(IK)=C 
      IDX(IK)=IC
      GO TO (3,45) ,ICT 
    3 IF (I.GE.N2) GO TO 35 
      I=I+1 
      GO TO 11
   35 ICT=2 
      NP2=N+2 
      I=2 
   37 N1=NP2-I
      NN=N1 
      IK=1
      GO TO 15
  45  CONTINUE
      T = A(1)
      A(1) = A(N1)
      A(N1) = T 
      IT = IDX(1) 
      IDX(1) = IDX(N1)
      IDX(N1) = IT
      IF (I.GE.N) GO TO 55
      I=I+1 
      GO TO 37
   55 RETURN
c60	continue
   60 WRITE(16,500)
  500 FORMAT('ERROR RETURN FROM SORTD1 - N LESS THAN OR EQUAL TO 1')
      STOP
   65 IDX(1)=1
      RETURN
      END

