! PROGRAM 01 - DIFFUSION PROCESS
! Partial Differential Equation: 
! dU/dt - a * d^2(U)/dx^2 = 0
! Solution using Finite Element Method

Program Pr
      Implicit none
      
      INTEGER, parameter:: IO = 12 ! input-output unit
      INTEGER NX,NT,I,J,ID,m,ISCHEME
      REAL,ALLOCATABLE :: U(:),UN(:),X(:)
      REAL L,h,VNM,dt,t,Time
      REAL C0, C1, a
            
      WRITE(*,*) 'Read input file' 
      OPEN(IO,FILE='Input.txt')
      READ(IO,*) L
      READ(IO,*) NX
      READ(IO,*) Time
      READ(IO,*) VNM
      READ(IO,*) a
      READ(IO,*) C0, C1
      READ(IO,*) ISCHEME
      CLOSE(IO)
      
      ALLOCATE(U(1:NX),UN(1:NX),X(1:NX))      
      Time = (L*L/16)/a 
      h= L/NX
      dt= VNM*h*h/a
      NT= Time/dt
      
      WRITE(*,*) 'L=',L, 'h=', h, 'NX=', NX
      WRITE(*,*) 'VNM=', VNM, 'dt=', dt, 'Time=', Time, 'NT=', NT
      !pause 1

      X(1)=0.0
      DO I=2, NX-1
            X(I)=X(I-1)+h
      END DO
      X(NX)=L
      
      U(:)=0.0
      UN(:)=0.0
          
      CALL InitValue(NX, U)
      CALL BoundValue(NX, U, C0, C1)
      
      OPEN(IO,FILE='Res.dat')
     
!-------------------------  Solve equation ------------------
	DO I=1, NT
		UN(1) = 0
		UN(NX) = 1
		DO J=2, NX - 1
			UN(J) = U(J) + VNM*(U(J-1)-2*U(J)+U(J+1))
		ENDDO
		U=UN
     	        UN(:)=0.0
	ENDDO


!---------------------------Output results--------------------                       

      call Output(NX, U, X, IO)
      close(IO)

!-------------------------Plot results-------------------------
      
      open(unit = 10, file = "Res.dat", status = "unknown")      
      call system("gnuplot -p plot_res.plt")

End program

!----------------------- Set Initial Value -----------------      
      SUBROUTINE InitValue(NX, U)
            IMPLICIT NONE
            INTEGER NX
            REAL U(NX)
            U(:)=0.0

      END SUBROUTINE

!----------------------- Set Boundary Condition ------------            
      SUBROUTINE BoundValue(NX, U, C0, C1)
            IMPLICIT NONE
            INTEGER NX
            REAL U(NX), C0, C1
            U(1) = C0
            U(NX) = C1
      END SUBROUTINE

!----------------------- 
      SUBROUTINE Output(NX, U, X, IO)
            IMPLICIT NONE
	      INTEGER I, NX, IO
	      REAL U(NX), X(NX)
	      DO I=1, NX
		      WRITE(IO,*) U(I), 'X =', X(I)
	      END DO
      END SUBROUTINE

