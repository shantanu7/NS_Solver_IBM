! Source code for a 2D Poisson Equation solver.
! Written by Shantanu Bailoor at The George Washington University.
!
! Capabilities:
!--------------
! 1. Code can handle Dirichlet and Neumann Boundary conditions (homogenous/non-homogenous)
! 2. Iterative schemes implemented: Jacobi, Successive Over-relaxation, Conjugate Gradient
! 3. Output file is cell-centered, Tecplot-compatible, ASCII data file (out.plt).
 
 SUBROUTINE Poisson_Solver
 
 USE grid_params
 USE global
 USE flo_vars
 USE Poisson_cal_params, ONLY: converged, U
 
 IMPLICIT NONE 
 
 CALL get_coeffs 
 CALL initialize_Poisson
 
 U = Pr
 DO WHILE (.NOT. converged)
	 CALL impose_BC
	 CALL iterate_Poisson
 END DO 
 
 IF (iProblem .EQ. TaylorGreen) CALL compute_global_error
 
 END SUBROUTINE Poisson_Solver 
 !-----------------------------------------------------------------------! 
 
 
 SUBROUTINE allocate_Poisson_memory
 
 USE grid_params
 USE global
 USE Poisson_cal_params
 IMPLICIT NONE
 
 INTEGER :: i, j, k, seed
 
 WRITE(*,'(A)', ADVANCE = 'NO') ' ALLOCATING MEMORY FOR POISSON SOLVER ...'
 
 ALLOCATE( U(nci,ncj), Uold(nci,ncj), Uex(nci-2,ncj-2), f(nci,ncj) )
 
 !Initialize solution variables
 !-----------------------------!
 U = 0.; Uold = 0.; Uex = 0; f = 0.
 
 !~CALL CPU_TIME(tstart)
 
 iter = 1
 
 !~WRITE(*,*) ' DONE!'
 !~PRINT*,''

 !~WRITE(*,'(A)', ADVANCE = 'NO') ' ALLOCATING MEMORY FOR CALCULATION VARIABLES ...'	
 
 ALLOCATE(aP(nci,ncj), aE(nci,ncj), aW(nci,ncj), aS(nci,ncj), aN(nci,ncj))
 
 IF (Sol_Method .EQ. ConGrad) THEN
	ALLOCATE(A_CG(nci-1,ncj-1), b_CG(nci-1,ncj-1))
	ALLOCATE(phi_CG(nci-1,ncj-1), phiO_CG(nci-1,ncj-1), resid_CG(nci-1,ncj-1))
	ALLOCATE(Ar_CG(nci-1,ncj-1))
	
	A_CG = 0.; b_CG = 0.; resid_CG = 0.; Ar_CG = 0.; 
	phi_CG = 0.; phiO_CG = 0.
 END IF
 
 WRITE(*,*) ' DONE!'
 PRINT*,''
 
 END SUBROUTINE allocate_Poisson_memory
 !-----------------------------------------------------------------------!
 
 
 SUBROUTINE get_coeffs
 
 USE global
 USE grid_params
 USE Poisson_cal_params
 IMPLICIT NONE
 
 INTEGER :: i, j, ki, kj, ki1, kj1
 
 !Matrix form
 INTEGER :: iEp, iWp, iNp, iSp
 
 REAL(KIND = 8) :: dxce, dxcw, dycn, dycs, dxew, dyns
 REAL(KIND = 8) :: angX, angY
  
 DO i = 2, nci-1
	DO j = 2, ncj-1
	
	! First set the RHS vector
	 !~angX = 2*pi*freq*xc(i,j); angY = 2*pi*freq*yc(i,j)
	 !~IF (Mode .EQ. 1) THEN
		!~f(i,j) = -2.*((2.*pi*freq)**2)*(SIN(ANGX)*SIN(ANGY))
	 !~ELSEIF (Mode .EQ. 2) THEN
		!~f(i,j) = -2.*((2.*pi*freq)**2)*(COS(ANGX)*COS(ANGY))
	 !~END IF
	
	! Set coefficients for A matrix elements
	 dxce = xc(i+1,j)-xc(i,j); dxcw = xc(i,j)-xc(i-1,j);
	 dycn = yc(i,j+1)-yc(i,j); dycs = yc(i,j)-yc(i,j-1);
	 !~dxew = 0.5*( (x(i+1,j)-x(i,j))+(x(i+1,j+1)-x(i,j+1)) )
	 !~dyns = 0.5*( (y(i,j+1)-y(i,j))+(y(i+1,j+1)-y(i+1,j)) )
	 dxew = ux(i+1,j)-ux(i,j)
	 dyns = vy(i,j+1)-vy(i,j)
	 
	 aE(i,j) = 1./(dxce*dxew)
	 aW(i,j) = 1./(dxcw*dxew)
	 aN(i,j) = 1./(dycn*dyns) 
	 aS(i,j) = 1./(dycs*dyns)
	 aP(i,j) = -1.*(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j))
	 
	 !~PRINT*,'i,j:',i,j
	 !~PRINT*,aP(i,j), aE(i,j), aW(i,j),aN(i,j), aS(i,j)
	 !~PRINT*,''
	 
	 IF (Sol_Method .EQ. ConGrad) THEN
		ki = i-1; kj = j-1
		b_CG(ki,kj) = f(i,j)
	 END IF

	END DO
 END DO	

 END SUBROUTINE get_coeffs
 !-----------------------------------------------------------------------!
 
 
 SUBROUTINE initialize_Poisson
 
 USE grid_params
 USE global
 USE flo_vars
 USE Poisson_cal_params
 IMPLICIT NONE
 
 U = Pr; Uold = Pr
 iter = 1; resid0 = 0.
 converged = .FALSE.
 
 END SUBROUTINE initialize_Poisson 
 !-----------------------------------------------------------------------!
 
 
 SUBROUTINE impose_BC
 
 USE grid_params
 USE global
 USE Poisson_cal_params
 IMPLICIT NONE
 
 INTEGER :: i, j, ki, kj, nbc, sumBC, ind
 
 REAL :: xy1, xy2, UBC, ANGX, ANGY
 
 sumBC = 0
 
 IF (BoundCond(1) .EQ. WALL) THEN
	sumBC = sumBC + 1
	U(1,:) = U(2,:)! + BCVal(1)*(xc(2,:)-xc(1,:))
 ELSEIF (BoundCond(1) .EQ. PERIODIC) THEN
	sumBC = sumBC + 1
	U(1,:) = U(nci-1,:)	
 ELSE !Dirichlet -- Fix
	xy1 = x(2,1)-xc(1,1); xy2 = xc(2,1)-x(2,1)
	U(1,:) = ((xy1+xy2)*BCVal(1)-xy1*U(2,:))/xy2
 END IF
 
 IF (BoundCond(2) .EQ. WALL) THEN
	sumBC = sumBC + 1
	U(nci,:) = U(nci-1,:)! + BCVal(2)*(xc(nci,:)-xc(nci-1,:))
 ELSEIF (BoundCond(2) .EQ. PERIODIC) THEN
	sumBC = sumBC + 1
	U(nci,:) = U(2,:)
 ELSE	
	xy2 = x(ngi-1,1)-xc(nci-1,1); xy1 = xc(nci,1)-x(ngi-1,1)
	U(nci,:) = ((xy1+xy2)*BCVal(2)-xy1*U(nci-1,:))/xy2
 END IF
 
 IF (BoundCond(3) .EQ. WALL) THEN
	sumBC = sumBC + 1
	U(:,1) = U(:,2)! + BCVal(3)*(yc(:,2)-yc(:,1))
 ELSEIF (BoundCond(3) .EQ. PERIODIC) THEN
	sumBC = sumBC + 1
	U(:,1) = U(:,ncj-1)
 ELSE
	xy1 = y(1,2)-yc(1,1); xy2 = yc(1,2)-y(1,2)
	U(:,1) = ((xy1+xy2)*BCVal(3)-xy1*U(:,2))/xy2
 END IF
 
 IF (BoundCond(4) .EQ. WALL) THEN
	sumBC = sumBC + 1
	U(:,ncj) = U(:,ncj-1)! + BCVal(4)*(yc(:,ncj)-yc(:,ncj-1))
 ELSEIF (BoundCond(4) .EQ. PERIODIC) THEN
	sumBC = sumBC + 1
	U(:,ncj) = U(:,2)
 ELSE
	xy2 = y(1,ngj-1)-yc(1,ncj-1); xy1 = yc(1,ncj)-y(1,ngj-1)
	U(:,ncj) = ((xy1+xy2)*BCVal(4)-xy1*U(:,ncj-1))/xy2
 END IF
 
 !Handle four corners
 !-------------------!
 !~U(1,1) = 0.5*(U(1,2)+U(2,1))
 !~U(nci,1) = 0.5*(U(nci,2)+U(nci-1,1))
 !~U(1,ncj) = 0.5*(U(2,ncj)+U(1,ncj-1))
 !~U(nci,ncj) = 0.5*(U(nci-1,ncj)+U(nci,ncj-1))
 
 !Special case to avoid singularity
 !---------------------------------!
 IF ( sumBC .EQ. 4 ) THEN
	
	ind = INT(nci/2)
	IF (iProblem .EQ. TaylorGreen) THEN
		U(ind,1) = -1./4.*exp(-4.*time)*(COS(2.*xc(ind,1))+COS(2.*yc(ind,1)))
	ELSEIF (iProblem .EQ. DrivenCavity) THEN
		xy1 = y(1,2)-yc(1,1); xy2 = yc(1,2)-y(1,2)
		U(ind,1) = -1.*U(ind,2)
	ELSE
		U(ind,1) = -1.*U(ind,2)
	END IF	
	
 END IF
 
 END SUBROUTINE impose_BC 
 !-----------------------------------------------------------------------!
 
 
 SUBROUTINE iterate_Poisson
 
 USE grid_params
 USE global
 USE Poisson_cal_params
 IMPLICIT NONE
 
 INTEGER :: i, j, ki, kj
 REAL(KIND = 8) :: summ
 
SELECT CASE (Sol_Method)

	CASE (JACOBI)
	DO i = 2, nci-1
		DO j = 2, ncj-1
			U(i,j) = (f(i,j)-aE(i,j)*Uold(i+1,j)-aW(i,j)*Uold(i-1,j)- & 
					aN(i,j)*Uold(i,j+1)-aS(i,j)*Uold(i,j-1))/aP(i,j)
		END DO
	END DO

	CASE (SOR)
	DO i = 2, nci-1
		DO j = 2, ncj-1
			U(i,j) = (omega*f(i,j)- omega*(aE(i,j)*Uold(i+1,j) + &
					aN(i,j)*Uold(i,j+1) + aW(i,j)*U(i-1,j) + &
					aS(i,j)*U(i,j-1)) -(omega-1)*aP(i,j)*Uold(i,j))/aP(i,j)
		END DO
	END DO	
		
	CASE (ConGrad)
	DO i = 2, nci-1
		DO j = 2, ncj-1
			ki = i-1; kj = j-1
			phiO_CG(ki,kj) = Uold(i,j)
		END DO
	END DO 

 END SELECT

 CALL compute_residual

 IF (Sol_Method .EQ. Congrad) THEN
	 DO i = 2, nci-1
		DO j = 2, ncj-1
			ki = i-1; kj = j-1
			phi_CG(ki,kj) = phiO_CG(ki,kj) + alpha*resid_CG(ki,kj)
			U(i,j) = phi_CG(ki,kj)
		END DO
	 END DO 
 END IF
 
 Uold = U
 
 END SUBROUTINE iterate_Poisson
 !-----------------------------------------------------------------------!


 SUBROUTINE compute_residual
 
 USE grid_params
 USE global
 USE Poisson_cal_params
 USE flo_vars
 IMPLICIT NONE
 
 INTEGER, PARAMETER :: fReisd = 102
 INTEGER :: i, j, ki, kj, fail = 0
 INTEGER :: resMaxLoc(2)
 
 REAL(KIND = 8), PARAMETER :: residLimit = 1E5
 REAL(KIND = 8) :: resid, residMax, temp
 
 converged = .FALSE.
 residG = 0.; residMax = 0.
 
 IF ( (Sol_Method .EQ. JACOBI) .OR. (Sol_Method .EQ. SOR) ) THEN
	 resid = 0.
	 DO i = 2, nci-1
		DO j = 2, ncj-1
			resid = f(i,j)-aE(i,j)*U(i+1,j)-aW(i,j)*U(i-1,j)- & 
					 aN(i,j)*U(i,j+1)-aS(i,j)*U(i,j-1)-aP(i,j)*U(i,j)
					 
			residG = residG + resid**2
			IF (ABS(resid) .GT. residMax) THEN
				residMax = resid
				resMaxLoc(:) = (/i,j/)
			END IF
		END DO
	 END DO
 ELSEIF (Sol_Method .EQ. ConGrad) THEN
	DO i = 2, nci-1
		DO j = 2, ncj-1
			ki = i-1; kj = j-1
			IF ( (iter .EQ. 1) .OR. (MOD(iter,1) .EQ. 0) ) THEN
				resid_CG(ki,kj) = b_CG(ki,kj) - ( aE(i,j)*U(i+1,j) + & 
							 aW(i,j)*U(i-1,j) + aN(i,j)*U(i,j+1) + &
							 aS(i,j)*U(i,j-1) + aP(i,j)*U(i,j))
			ELSE
				resid_CG(ki,kj) = resid_CG(ki,kj) - alpha*Ar_CG(ki,kj)
			END IF
			residG = residG + resid_CG(ki,kj)**2
		END DO
	END DO
	CALL compute_alpha
 END IF
 
 residG = SQRT(residG)/SQRT(REAL(nci*ncj))
 IF (iter .EQ. 1) resid0 = residG
 
 !~IF (iter .EQ. 1) THEN
	!~OPEN(fReisd, FILE = 'convergence.dat')
	!~WRITE(fReisd,'(A)') 'VARIABLES = iter, residG'
 !~END IF
 
 !WRITE(fReisd,'(2X,I5,2X,E14.6)') iter, residG
 
 !~IF (iter .GE. 20000) THEN
 IF (MOD(iter,1000) .EQ. 0) THEN
	WRITE(*,'(A,I7,5X,A,E14.7)')' >> ITER: ', iter, 'RESIDUAL NORM: ', residG
	WRITE(*,'(A,2I3,2X,E14.7)')'    RESIDMAX AT:', resMaxLoc, residMax
 END IF
 !~END IF 
 !~STOP

 IF (residG .LE. residMon .AND. iter .NE. 1) THEN
	converged = .TRUE.
	fail = 0
 ELSE IF (residG/resid0 .GE. residLimit .AND. iter .GE. 1000) THEN
	converged = .TRUE.
	fail = 1
 END IF
 
 IF (converged .EQV. .TRUE.) THEN
	pr = U
	!~CLOSE(fReisd)
	!~CALL CPU_TIME(tend)
	!~minute = INT((tend-tstart)/60); sec = (tend-tstart)-60.*REAL(minute)
	!~PRINT*,''
	!~PRINT*,'!--------------------------- SOLUTION REPORT ---------------------------!'
	IF (fail .EQ. 0) THEN
		IF (iterG .EQ. 1 .OR. MOD(iterG,1) .EQ. 0 .OR. iterG .EQ. itermax) THEN
			WRITE(*,'(A,I7)')' SOLUTION CONVERGED AT ITERATION: ', iter
		END IF
	ELSE
		!~IF (iterG .EQ. 1 .OR. MOD(iterG,1) .EQ. 0 .OR. iterG .EQ. itermax) THEN
			WRITE(*,'(A,I7)')' SOLUTION DIVERGED AT ITERATION: ', iter
			WRITE(*,'(A)')' ABORTING!'
			STOP
		!~END IF	
	END IF
	!~WRITE(*,'(A,I2.2,A,F7.4,A)') ' SOLUTION TIME: ',minute, &
	!~		' MINUTES, ', sec, ' SECONDS'
	!~WRITE(*,'(A,E14.7)')' GLOBAL RESIDUAL NORM: ', residG
 END IF
 
 iter = iter + 1
 
 END SUBROUTINE compute_residual
 !-----------------------------------------------------------------------!
 
 
 SUBROUTINE compute_alpha
 
 USE grid_params
 USE global
 USE Poisson_cal_params
 IMPLICIT NONE
 
 INTERFACE
 
	 REAL(KIND = 8) FUNCTION calc_inner_prod (V1, V2)
	 REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: V1, V2
	 END FUNCTION calc_inner_prod 
	 
	 SUBROUTINE mat_prod(M1, M2, M3)
	 REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: M1, M2
	 REAL(KIND = 8), DIMENSION(:,:), INTENT(OUT) :: M3
	 END SUBROUTINE mat_prod
 
 END INTERFACE
 
 REAL(KIND = 8) :: rE, rW, rN, rS, rP
 REAL(KIND = 8) :: num, den
 
 INTEGER :: i, j, ki, kj
  
  num = 0.; den = 0.
  
  DO i = 2, nci-1
	DO j = 2, ncj-1
		rE = 0.; rW = 0.; rN = 0.; rS = 0.
		ki = i-1; kj = j-1
		
		!Assign neighboring residuals:
		!-----------------------------
		
		!Not left boundary
		IF (ki .NE. 1) rW = resid_CG(ki-1,kj)
		!Not right boundary
		IF (ki .NE. ncx) rE = resid_CG(ki+1,kj)
		!Not bottom boundary
		IF (kj .NE. 1) rS = resid_CG(ki,kj-1)
		!Not top boundary
		IF (kj .NE. ncy) rN = resid_CG(ki,kj+1)
		
		rP = resid_CG(ki,kj)
		
		Ar_CG(ki,kj) = aE(i,j)*rE + aW(i,j)*rW + aN(i,j)*rN + &
					   aS(i,j)*rS + aP(i,j)*rP
					   
		!Compute numerator and denominator for alpha
		num = num + rP**2
		den = den + rP*Ar_CG(ki,kj)
	END DO
  END DO
  
  alpha = num/den
 
 END SUBROUTINE compute_alpha
 !-----------------------------------------------------------------------!
 
 
 SUBROUTINE mat_prod(M1, M2, M3)
 IMPLICIT NONE

 REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: M1, M2
 REAL(KIND = 8), DIMENSION(:,:), INTENT(OUT) :: M3
 REAL(KIND = 8) :: sum1

 INTEGER :: p1, q1, p2, q2, p3, q3
 INTEGER :: i, j, k

 p1 = SIZE(M1(:,1)); q1 = SIZE(M1(1,:));
 p2 = SIZE(M2(:,1)); q2 = SIZE(M2(1,:));
 p3 = p1; q3 = q2

 IF (p2 .ne. q1) THEN
     PRINT*,''
     PRINT*,'ERROR: INCOMPATIBLE MATRIX DIMENSIONS'
     PRINT*,''
     STOP
 END IF

 DO i = 1, p1
     DO j = 1, q2
         sum1 = 0.
         DO k = 1, p2
	     sum1 = sum1 + M1(i,k)*M2(k,j) 
         END DO
         M3(i,j) = sum1
     END DO    
 END DO

 END SUBROUTINE mat_prod
 !-----------------------------------------------------------------------!
 
 
 REAL(KIND = 8) FUNCTION calc_inner_prod (V1, V2)
 
 REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: V1, V2
 REAL(KIND = 8) :: temp
 INTEGER :: n1, n2, i
 
 n1 = SIZE(V1); n2 = SIZE(V2); temp = 0.
 
 IF (n1 .EQ. n2) THEN
	 DO i = 1, n1
		temp = temp + V1(i)*V2(i)
		calc_inner_prod = temp
	 END DO
 ELSE
	 calc_inner_prod = 0.
	 WRITE(*,'(A)') ' ERROR! VECTOR DIMENSIONS FOR INNER PRODUCT DO NOT MATCH!'
	 WRITE(*,'(A,I5)') ' n1:', n1
	 WRITE(*,'(A,I5)') ' n2:', n2
	 WRITE(*,'(A)') ' ABORTING!'
	 PRINT*,''
	 STOP
 END IF
 
 END FUNCTION calc_inner_prod
 !-----------------------------------------------------------------------!
 
 
 SUBROUTINE compute_global_error
 
 USE grid_params
 USE global
 USE Poisson_cal_params
 IMPLICIT NONE
 
 INTEGER :: i, j
 
 REAL (KIND = 8) :: ErrNorm, Error, MaxErr
 REAL (KIND = 8) :: AngX, AngY
 
 MaxErr = 0.; ErrNorm = 0.
 
 DO i = 2, nci-1
	DO j = 2, ncj-1
		Uex(i-1,j-1) = -1./4.*exp(-4.*time)*(COS(2.*xc(i,j))+COS(2.*yc(i,j)))
		Error = Uex(i-1,j-1) - U(i,j)
		ErrNorm = ErrNorm + Error**2
		IF (ABS(Error) .GT. ABS(MaxErr)) MaxErr = Error
	END DO
 END DO
 
 ErrNorm = SQRT(ErrNorm)/(REAL(nci*ncj))
 
 IF (MOD(iterG,1) .EQ. 0 .OR. iterG .EQ. itermax) THEN
	 WRITE(*,'(A,E14.7)')' GLOBAL ERROR NORM: ', ErrNorm
	 WRITE(*,'(A,E14.7)')' MAX ERROR: ', MaxErr
	 !PRINT*,''
 END IF	 
 
 END SUBROUTINE compute_global_error
 !-----------------------------------------------------------------------!
 
 
 SUBROUTINE free_Poisson_memory
 
 USE grid_params
 USE global
 USE Poisson_cal_params
 
 WRITE(*,'(A)', ADVANCE = 'NO') ' FREEING UP USED MEMORY ...'
 
 !DEALLOCATE (x, y, xc, yc)
 DEALLOCATE (U, Uold, Uex, f)
 DEALLOCATE (aP, aE, aW, aN, aS)
 
 IF (Sol_Method .EQ. ConGrad) THEN
	DEALLOCATE (A_CG, Ar_CG, b_CG, phi_CG, phiO_CG, resid_CG)
 END IF
 
 WRITE(*,'(A)') ' DONE!'
 PRINT*,''
 
 END SUBROUTINE free_Poisson_memory
 !-----------------------------------------------------------------------!
