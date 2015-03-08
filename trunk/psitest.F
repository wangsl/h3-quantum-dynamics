
C     $Id$


#if 0
      SubRoutine PsiTest(psi, n1, n2, nTheta, dr1, dr2, 
     $     x, w)
      Implicit None
      Integer n1, n2, nTheta
      Complex*16 psi(n1, n2, nTheta)
      Real*8 dr1, dr2, x(*), w(*)

      Real*8 S, SK
      Character*120 Line
      Integer I, J, K

      S = 0.0D0
C$OMP Parallel Do 
C$OMP& Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(I, J, K, SK)
C$OMP& Reduction(+:S)
      Do K = 1, NTheta
         SK = 0.0D0
         Do J = 1, N2
            Do I = 1, N1
               SK = SK + Conjg(Psi(I,J,K))*Psi(I,J,K)
            EndDo
         EndDo
         S = S + w(k)*SK
      EndDo
C$OMP End Parallel Do

      S = S*dr1*dr2
      
      Write(Line, '(F20.15)') S
      Call mexPrintf(Line//achar(10))
      
      End

#else

      SubRoutine PsiTest(Psi, n1, n2, nTheta, dr1, dr2, 
     $     x, w)
      Implicit None
      Integer n1, n2, nTheta
      Complex*16 psi(n1*n2, nTheta)
      Real*8 dr1, dr2, x(*), w(*)
      
      Real*8 S, SK
      Character*120 Line
      Integer J, K, N1N2
      
      N1N2 = N1*N2
      
      S = 0.0D0
C$OMP Parallel Do 
C$OMP& Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(J, K, SK)
C$OMP& Reduction(+:S)
      Do K = 1, NTheta
         SK = 0.0D0
         Do J = 1, N1N2
            SK = SK + Conjg(Psi(J,K))*Psi(J,K)
         EndDo
         S = S + w(k)*SK
      EndDo
C$OMP End Parallel Do
      
      S = S*dr1*dr2
      
      Write(Line, '(F20.15)') S
      Call mexPrintf(Line//achar(10))
      
      End
#endif

