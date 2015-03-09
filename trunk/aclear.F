

C     $Id$

      Subroutine AClear(N, A)
      Implicit Real*8(A-H,O-Z)
      
      Real*8 A(*)
      Save Zero
      Data Zero/0.0D0/
      
C$OMP Parallel Do 
C$OMP& Schedule(Static, 1)
C$OMP& Default(Shared)
C$OMP& Private(I)
      Do I = 1, N
         A(I) = Zero
      EndDo
C$OMP End Parallel Do

      Return
      End

