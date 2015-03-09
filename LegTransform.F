

C     $Id$

      SubRoutine ForwardLegendreTransform(
     $     CPsi, LPsi, NR, NTheta, NLeg, W, LegP)
      Implicit None
      Integer NR, NTheta, NLeg
      Complex*16 CPsi(NR, NTheta), LPsi(NR, NLeg)
      Real*8 W(NTheta), LegP(NTheta, NLeg)
      
      Integer I, K, L
      Real*8 WP, F
      
      Call AClear(2*NR*NLeg, LPsi)

C$OMP Parallel Do
C$OMP& Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(L, K, I, WP)
      Do L = 1, NLeg
         Do K = 1, NTheta
            WP = W(K)*LegP(K,L)
            Do I = 1, NR
               LPsi(I,L) = LPsi(I,L) + WP*CPsi(I,K)
            EndDo
         EndDo
      EndDo
C$OMP End Parallel Do

C$OMP Parallel Do
C$OMP& Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(L, I, F)
      Do L = 1, NLeg
         F = L - 0.5D0
         Do I = 1, NR
            LPsi(I,L) = F*LPsi(I,L)
         EndDo
      EndDo
C$OMP End Parallel Do
      
      Return
      End

      
      SubRoutine BackwardLegendreTransform(
     $     CPsi, LPsi, NR, NTheta, NLeg, W, LegP)
      Implicit None
      Integer NR, NTheta, NLeg
      Complex*16 CPsi(NR, NTheta), LPsi(NR, NLeg)
      Real*8 W(NTheta), LegP(NTheta, NLeg)
      
      Integer I, K, L
      
      Call AClear(2*NR*NTheta, CPsi)
      
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

      
