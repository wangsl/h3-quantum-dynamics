
C     $Id$

      SubRoutine ForwardLegendreTransform(
     $     CPsi, LPsi, NR, NTheta, NLeg, WLegP)
      Implicit None
      Integer NR, NTheta, NLeg
      Complex*16 CPsi(NR, NTheta), LPsi(NR, NLeg)
      Real*8 WLegP(NTheta, NLeg)

      Integer I, K, L

      Call DCZeros(NR*NLeg, LPsi)
      
C$OMP Parallel Do
C$OMP& Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(L, K, I)
      Do L = 1, NLeg
         Do K = 1, NTheta
            Do I = 1, NR
               LPsi(I,L) = LPsi(I,L) + CPsi(I,K)*WLegP(K,L)
            EndDo
         EndDo
      EndDo
C$OMP End Parallel Do
      
      Return
      End

      SubRoutine BackwardLegendreTransform(
     $     CPsi, LPsi, NR, NTheta, NLeg, LegP)
      Implicit None
      Integer NR, NTheta, NLeg
      Complex*16 CPsi(NR, NTheta), LPsi(NR, NLeg)
      Real*8 LegP(NTheta, NLeg)
      
      Integer I, K, L
      
      Call DCZeros(NR*NTheta, CPsi)
      
C$OMP Parallel Do
C$OMP& Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(L, K, I)
      Do K = 1, NTheta
         Do L = 1, NLeg
            Do I = 1, NR
               CPsi(I,K) = CPsi(I,K) + LPsi(I,L)*LegP(K,L)
            EndDo
         EndDo
      EndDo
C$OMP End Parallel Do

      Return
      End
      
#if 0
C     C = A*B
      SubRoutine MyTestC_AB(A, B, C, N1, N2, N3)
      Implicit None
      Integer N1, N2, N3
      Real*8 A(N1,N2), B(N2,N3), C(N1,N3)

      Integer I, J, K

      Call DZeros(C, N1*N3)

C$OMP Parallel Do
C$OMP& Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(I, J, K)
      Do J = 1, N3
         Do K = 1, N2
            Do I = 1, N1
               C(I,J) = C(I,J) + A(I,K)*B(K,J)
            EndDo
         EndDo
      EndDo
C$OMP End Parallel Do

      End
#endif
