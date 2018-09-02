!-----------------------------------------------------------------
! This module computes the solution of a sparse 
! symmetric/unsymmetric linear system : Ax = b
! A is a sparse matrix 
!----------------------------------------------------------------- 
! This module comes from Quantsolver project 
!
! Julien Versaci 
!-----------------------------------------------------------------

module linearSystemSolver
implicit none 

!..     Internal solver memory pointer 
!        INTEGER*8 pt(64)

!..     All other variables 
!        INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
!        INTEGER iparm(64)
!        REAL*8  dparm(64) 

!        INTEGER i, j, idum, solver
!        REAL*8  waltime1, waltime2, ddum, normb, normr

!        DATA  nrhs /1/, maxfct /1/, mnum /1/

! External Subroutines ..
      EXTERNAL         DGESV

contains 

!.. ----------------------------------------------------------------        
!.. unsymmetric linear system
!.. ----------------------------------------------------------------

!        subroutine unsymmetric_linear_solver(a,ia,ja,b,x,Nt) 
!        REAL*8 :: a(:)
!        REAL*8 :: b(:)
!        REAL*8 :: x(:)
!        INTEGER :: ia(:)
!        INTEGER :: ja(:)
!        INTEGER :: Nt
             

!  .. Setup Pardiso control parameters und initialize the solvers     
!     internal adress pointers. This is only necessary for the FIRST   
!     call of the PARDISO solver.                                     
!     
!      mtype     = 11      ! unsymmetric matrix
!      solver    = 10     ! use sparse direct method

!      call pardisoinit(pt, mtype, solver, iparm, dparm, error)

!  .. Numbers of Processors ( value of OMP_NUM_THREADS )
!      iparm(3) = 1


!..   Reordering and Symbolic Factorization, This step also allocates
!     all memory that is necessary for the factorization 
!
!      phase     = 11      ! only reordering and symbolic factorization
!      msglvl    = 1       ! with statistical information

      
!      CALL pardiso (pt, maxfct, mnum, mtype, phase, Nt, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error,dparm)
     
!      WRITE(*,*) 'Reordering completed ... '

!      IF (error .NE. 0) THEN
!        WRITE(*,*) 'The following ERROR was detected: ', error
!        STOP
!      END IF

!      WRITE(*,*) 'Number of nonzeros in factors   = ',iparm(18)
!      WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

!.. Factorization.
!      phase     = 22  ! only factorization
!      CALL pardiso (pt, maxfct, mnum, mtype, phase, Nt, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error,dparm) 

!      WRITE(*,*) 'Factorization completed ...  '
!      IF (error .NE. 0) THEN
!         WRITE(*,*) 'The following ERROR was detected: ', error
!        STOP
!      ENDIF

!.. Back substitution and iterative refinement
!      phase     = 33  ! only solve
!      iparm(8)  = 1   ! max numbers of iterative refinement steps

!      CALL pardiso (pt, maxfct, mnum, mtype, phase, Nt, a, ia, ja, idum, nrhs, iparm, msglvl, b, x, error, dparm) 
!      WRITE(*,*) 'Solve completed ... '
!     
!      WRITE(*,*) 'The solution of the system is '
!      DO i = 1, Nt
!        WRITE(*,*) ' x(',i,') = ', x(i)
!      END DO

!.. Back substitution and solution with A^Tx=b  and iterative refinement
!     phase     = 33  ! only solve
!      iparm(8)  = 1   ! max numbers of iterative refinement steps
!      iparm(12) = 1   ! Solving with transpose matrix
!      do i = 1, Nt
!        b(i) = 1.d0
!      end do
!
!      CALL pardiso (pt, maxfct, mnum, mtype, phase, Nt, a, ia, ja, idum, nrhs, iparm, msglvl, b, x, error, dparm)
!      WRITE(*,*) 'Solve completed ... '
    
!      WRITE(*,*) 'The solution of the system is '
!      DO i = 1, Nt
!        WRITE(*,*) ' x(',i,') = ', x(i)
!      END DO


!.. Termination and release of memory
!      phase     = -1           ! release internal memory
!      CALL pardiso (pt, maxfct, mnum, mtype, phase, Nt, ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm)

! end subroutine


 subroutine lapack_solver(matrix,cl_tab,Nt,NRHS) 
      ! Parameters ..
      INTEGER, INTENT(IN) :: NRHS, Nt
      REAL*8 :: matrix(:,:)
      REAL*8 :: cl_tab(:)
      INTEGER :: LDA , LDB 
      ! local scalar
      INTEGER ::  INFO 
      ! local array 
      INTEGER :: IPIV(Nt)

     LDA = Nt
     LDB = Nt

      !!  Executable Statements ..
      WRITE(*,*)'DGESV Example Program Results'
!
!     Solve the equations A*X = B.
!
      CALL DGESV( Nt, NRHS, matrix, LDA, IPIV, cl_tab, LDB, INFO )

  
    end subroutine 


end module 

