

C     $Id$

      SubRoutine ForwardLegendreTransform(
     $     CPsi, LPsi, NR, NTheta, NLeg, W, LegP)
      Implicit None
      Integer NR, NTheta, NLeg
      Complex*16 CPsi(NR, NTheta), LPsi(NR, NLeg)
      Real*8 W(NTheta), LegP(NTheta, NLeg)
      
      Real*8 WLegP(NTheta, NLeg)
      Integer I, K, L
      Real*8 WP, F

C$OMP Parallel Do
C$OMP& Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(L, K, F)
      Do L = 1, NLeg
         F = L - 0.5D0
         Do K = 1, NTheta
            WLegP(K,L) = W(K)*LegP(K,L)*F
         EndDo
      EndDo
C$OMP End Parallel Do

      Call DZeros(2*NR*NLeg, LPsi)

C$OMP Parallel Do
C$OMP& Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(L, K, I, WP)
      Do L = 1, NLeg
         Do K = 1, NTheta
            WP = WLegP(K,L)
            Do I = 1, NR
               LPsi(I,L) = LPsi(I,L) + CPsi(I,K)*WP
            EndDo
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
      
      Call DZeros(2*NR*NTheta, CPsi)
      
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

      
