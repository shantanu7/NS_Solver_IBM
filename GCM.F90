 SUBROUTINE create_GCM
 
 USE grid_params
 USE global
 USE GCM
 IMPLICIT NONE
 
 INTEGER, PARAMETER :: PRESS = 1, XVEL = 2, YVEL = 3
 INTEGER, ALLOCATABLE, DIMENSION(:,:) :: tempGC, tempGmap, nb14
 INTEGER :: ivars, is, ie, js, je
 INTEGER :: i, j, i1, j1, ki, kj, ind, found
 INTEGER :: nem, iloc, jloc, iGC, jGC, ileft ,iright, jbot, jtop
 INTEGER :: v1, v2 
 INTEGER :: nGC1
 
 REAL(KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: xtemp, ytemp, iBlankT
 REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: xBIt, yBIt, xIPt, yIPt
 REAL(KIND = 8) :: x1, x2, y1, y2
 REAL(KIND = 8) :: xx, yy, dArea, dyn, dys, dxe, dxw, dprojXY
 REAL(KIND = 8) :: vect(2), vectL(2), vectR(2), prodL, prodR, lenL, lenR
 REAL(KIND = 8) :: dist, distMin, inProd, lenProj, lenBI, diaglen, lenIP
 REAl(KIND = 8) :: ratio

 ALLOCATE(gmap(nci,ncj), gmapU(nci+1,ncj), gmapV(nci,ncj+1))
 ALLOCATE(NormElX(nElemMarker), NormElY(nElemMarker))
 ALLOCATE(xElemCent(nElemMarker), yElemCent(nElemMarker))
 ALLOCATE(iBlank(nci,ncj), iBlanku(nci+1,ncj), iBlankv(nci,ncj+1))
 ALLOCATE(AreaBodyMarker(nElemMarker), nbMarker(nElemMarker,2))
 ALLOCATE(LSFMarker(nElemMarker,4),pBodyMarker(nElemMarker)) 
 
 gmap = 0; gmapU = 0; gmapV = 0
 iBlank = 0; iBlanku = 0; iBlankv = 0
 
 WRITE(*,'(A)', ADVANCE = 'NO') ' CREATING GCM DATA STRUCTURES...'
 
 DO ivars = 1, nvarsGC
 
	IF (ivars .EQ. PRESS) THEN
		is = 1; ie = nci
		js = 1; je = ncj
		ALLOCATE(xtemp(ie-is+1, je-js+1), ytemp(ie-is+1, je-js+1))
		xtemp = xc; ytemp = yc
	ELSEIF (ivars .EQ. XVEL) THEN
		is = 1; ie = nci+1
		js = 1; je = ncj
		ALLOCATE(xtemp(ie-is+1, je-js+1), ytemp(ie-is+1, je-js+1))
		xtemp = ux; ytemp = uy
	ELSEIF (ivars .EQ. YVEL) THEN
		is = 1; ie = nci
		js = 1; je = ncj+1
		ALLOCATE(xtemp(ie-is+1, je-js+1), ytemp(ie-is+1, je-js+1))
		xtemp = vx; ytemp = vy
	END IF
	
	ALLOCATE(tempGmap(ie-is+1,je-js+1))
	ALLOCATE(tempGC((ie-is+1)*(je-js+1),2))
	tempGmap = 0; nGC1 = 0
 
	!Compute face normals
	DO nem = 1, nElemMarker
	
		v1 = ElemMarkerConnex(nem,1); v2 = ElemMarkerConnex(nem,2);
		x1 = xBodyMarker(v1); y1 = yBodyMarker(v1)
		x2 = xBodyMarker(v2); y2 = yBodyMarker(v2)
		dArea = SQRT((x2-x1)**2+(y2-y1)**2)
		AreaBodyMarker(nem) = dArea
		normElX(nem) = -1.*(y2-y1)/dArea
		normElY(nem) = (x2-x1)/dArea
		xElemCent(nem) = 0.5*(x1+x2); yElemCent(nem) = 0.5*(y1+y2)
		
		xx = xElemCent(nem); yy = yElemCent(nem)
	
		DO j = js+1, je-1
			dyn = ABS(ytemp(1,j+1)-ytemp(1,j)); 
			dys = ABS(ytemp(1,j)-ytemp(1,j-1))
			IF ( (yy .GT. ytemp(1,j)-0.5*dys) .AND. &
				 (yy .LE. ytemp(1,j)+0.5*dyn)) THEN
				jloc = j
				!Determine top and bottom cells
			END IF
		END DO
		
		DO i = is+1, ie-1
			dxe = ABS(xtemp(i+1,jloc)-xtemp(i,jloc)); 
			dxw = ABS(xtemp(i,jloc)-xtemp(i-1,jloc))
			IF ( (xx .GT. xtemp(i,jloc)-dxw/2.) .AND. &
				 (xx .LE. xtemp(i,jloc)+dxe/2.) ) THEN
				iloc = i
			END IF
		END DO
		
		IF (ABS(normElX(nem)) .GT. ABS(normElY(nem))) THEN
		
			!Determine porjected coordinate at ytemp(iloc,jloc)
			dprojXY = xx + (x2-x1)/(y2-y1)*(ytemp(iloc,jloc)-yy)
			
			!Determine left and right cells
			IF (dprojXY .LT. xtemp(iloc,jloc)) THEN
				ileft = iloc-1; iright = iloc
			ELSE
				ileft = iloc; iright = iloc+1
			END IF
			
			vectL(1) = xtemp(ileft,jloc)-xx; vectL(2) = ytemp(ileft,jloc)-yy
			lenL = SQRT((xtemp(ileft,jloc)-xx)**2+(ytemp(ileft,jloc)-yy)**2)
			vectL(1) = vectL(1)/lenL
			vectL(2) = vectL(2)/lenL
			
			vectR(1) = xtemp(iright,jloc)-xx; vectR(2) = ytemp(iright,jloc)-yy
			lenR = SQRT((xtemp(iright,jloc)-xx)**2+(ytemp(iright,jloc)-yy)**2)
			vectR(1) = vectR(1)/lenR
			vectR(2) = vectR(2)/lenR
				
			!Compute inner products to determine position of cells wrt solid
			prodL = normElX(nem)*vectL(1) + normElY(nem)*vectL(2)
			prodR = normElX(nem)*vectR(1) + normElY(nem)*vectR(2)

			IF (prodL .LT. prodR) THEN
				!Left cell is in the solid
				iGC = ileft
			ELSEIF (prodR .LT. prodL) THEN
				!Right cell is in the solid
				iGC = iright
			END IF
			jGC = jloc
			!~PRINT*,'INNER CELL AT: ',xtemp(iloc,jloc), ytemp(iloc,jloc)
			
		ELSE
		
			!Determine porjected coordinate at xtemp(iloc,jloc)
			dprojXY = yy + (y2-y1)/(x2-x1)*(xtemp(iloc,jloc)-xx)
			
			!Determine bottom and top cells
			IF (dprojXY .LT. ytemp(iloc,jloc)) THEN
				jbot = jloc-1; jtop = jloc
			ELSE
				jbot = jloc; jtop = jloc+1
			END IF
			
			vectL(1) = xtemp(iloc,jbot)-xx; vectL(2) = ytemp(iloc,jbot)-yy
			lenL = SQRT((xtemp(iloc,jbot)-xx)**2+(ytemp(iloc,jbot)-yy)**2)
			vectL(1) = vectL(1)/lenL
			vectL(2) = vectL(2)/lenL
			
			vectR(1) = xtemp(iloc,jtop)-xx; vectR(2) = ytemp(iloc,jtop)-yy
			lenR = SQRT((xtemp(iloc,jtop)-xx)**2+(ytemp(iloc,jtop)-yy)**2)
			vectR(1) = vectR(1)/lenR
			vectR(2) = vectR(2)/lenR		
			
			!Compute inner products to determine position of cells wrt solid
			prodL = normElX(nem)*vectL(1) + normElY(nem)*vectL(2)
			prodR = normElX(nem)*vectR(1) + normElY(nem)*vectR(2)

			IF (prodL .LT. prodR) THEN
				!Bottom cell is in the solid
				jGC = jbot
			ELSEIF (prodR .LT. prodL) THEN
				!Top cell is in the solid
				jGC = jtop
			END IF
			iGC = iloc
			!~PRINT*,'INNER CELL AT: ',xtemp(iloc,jloc), ytemp(iloc,jloc)
			
		END IF
		
		IF (iVerbose .EQ. 1) THEN
			PRINT*,'FOR elem:',nem
			PRINT*,'normElemX:', normElX(nem)
			PRINT*,'normElemY:', normElY(nem)
			PRINT*,'dprojXY:', dprojXY
			PRINT*,'prodL', prodL
			PRINT*,'prodR', prodR
			PRINT*,'XX', xx, 'YY', yy
			PRINT*,'Bounding cell',iloc-1, jloc-1
			PRINT*,'Ghost cell',iGC-1, jGC-1
		END IF
		
		!Get only unique cells
		IF (tempGmap(iGC,jGC) .NE. 1) THEN
			nGC1 = nGC1 + 1
			tempGC(nGC1,1) = iGC
			tempGC(nGC1,2) = jGC
			tempGmap(iGC,jGC) = 1
		END IF
		
	END DO
	
	IF (ivars .EQ. PRESS) THEN
		gmap = tempGmap
		nGC = nGC1
		ALLOCATE(GCp(nGC1,2))
		GCp(:,1) = tempGC(1:nGC1,1)
		GCp(:,2) = tempGC(1:nGC1,2)
	ELSEIF (ivars .EQ. XVEL) THEN	
		gmapU = tempGmap
		nGCu = nGC1
		ALLOCATE(GCu(nGC1,2))
		GCu(:,1) = tempGC(1:nGC1,1)
		GCu(:,2) = tempGC(1:nGC1,2)
	ELSEIF (ivars .EQ. YVEL) THEN	
		gmapV = tempGmap
		nGCv = nGC1
		ALLOCATE(GCv(nGC1,2))
		GCv(:,1) = tempGC(1:nGC1,1)
		GCv(:,2) = tempGC(1:nGC1,2)
	END IF
	
	!~DEALLOCATE(xtemp, ytemp)
	!~DEALLOCATE(tempGC, tempGmap)
	
	!Create iBlank
	 ALLOCATE( iBlankT(ie-is+1,je-js+1) )
	 DO i = is+1, ie-1
		DO j = js+1, je-1
			!IF (tempGmap(i,j) .EQ. 0) THEN
				distMin = 1.0E10
				DO nem = 1, nElemMarker
					!Calculate distance between cell center and element center				
					dist = SQRT( (xtemp(i,j)-xElemCent(nem))**2 + &
					(ytemp(i,j)-yElemCent(nem))**2 )
					IF (dist .LT. distMin) THEN
						distMin = dist
						iloc = nem
					END IF	
				END DO
				
				vect(1) = xElemCent(iloc)-xtemp(i,j); 
				vect(2) = yElemCent(iloc)-ytemp(i,j)
				inProd = normElX(iloc)*vect(1) + normElY(iloc)*vect(2)
				
				IF (inProd .GE. 0.) iBlankT(i,j) = 1
			!END IF	
		END DO
	 END DO
	 
	 
	 !Calculate body intercept, image point and interpolating stencil
	 IF (ivars .EQ. PRESS) THEN
		 iBlank = iBlankT
		 ALLOCATE(xBI(nGC), yBI(nGC), xIP(nGC), yIP(nGC), prIP(nGC))
		 ALLOCATE(nb4v1(nGC,2))
	 ELSEIF (ivars .EQ. XVEL) THEN	 
		 iBlanku = iBlankT
		 ALLOCATE(xBIu(nGCu), yBIu(nGCu), xIPu(nGCu), yIPu(nGCu), uIP(nGCu))
		 ALLOCATE(nb4v1u(nGCu,2))
	 ELSEIF (ivars .EQ. YVEL) THEN
		 iBlankv = iBlankT
		 ALLOCATE(xBIv(nGCv), yBIv(nGCv), xIPv(nGCv), yIPv(nGCv), vIP(nGCv))
		 ALLOCATE(nb4v1v(nGCv,2))
	 END IF	
	 
	 !Allocate temp vars
	 ALLOCATE(xBIt(nGC1), yBIt(nGC1), xIPt(nGC1), yIPt(nGC1), nb14(nGC1,2))
	 
	 DO ind = 1, nGC1
		i = tempGC(ind,1); j = tempGC(ind,2)
		!Find closest Lagrangian element center
		distMin = 1.0E10
		DO iloc = 1, nElemMarker		
			dist = SQRT( (xtemp(i,j)-xElemCent(iloc))**2 + (ytemp(i,j)-yElemCent(iloc))**2 )
			IF (dist .LT. distMin) THEN
				distMin = dist
				nem = iloc
			END IF
			xx = xElemCent(nem); yy = yElemCent(nem)
		END DO
		
		!Calculate length of BI
		lenBI = (xx-xtemp(i,j))*normElX(nem) + (yy-ytemp(i,j))*normElY(nem)
		!Calculate coordinates of BI
		xBIt(ind) = xtemp(i,j) + lenBI*normElX(nem)
		yBIt(ind) = ytemp(i,j) + lenBI*normElY(nem)
		
		!Calculate diagonal length of cell
		IF (ivars .EQ. PRESS) THEN
			diaglen = SQRT((x(i+1,j+1)-x(i,j))**2+(y(i+1,j+1)-y(i,j))**2)
		ELSEIF (ivars .EQ. XVEL) THEN
			diaglen = SQRT((vx(i+1,j)-ux(i,j))**2+(vy(i+1,j)-uy(i,j))**2)
		ELSEIF (ivars .EQ. YVEL) THEN
			diaglen = SQRT((ux(i,j+1)-vx(i,j))**2+(uy(i,j+1)-vy(i,j))**2)
		END IF
		
		!Project to image point
		ratio = lenBI/diaglen
		IF(ratio .LT. 0.20) THEN
			lenIP = 1.1*diaglen
		ELSEIF(ratio .LT. 0.45) THEN
			lenIP = 0.75*diaglen
		ELSE
			lenIP = lenBI
		END IF
		
		xIPt(ind) = xBIt(ind) + lenIP*normElX(nem)
		yIPt(ind) = yBIt(ind) + lenIP*normElY(nem)
		
		!Debug stuff
		!-----------
		!~PRINT*,''
		!~PRINT*,'iVARS:', ivars
		!~!~PRINT*,'FOR EL:', i-1, j-1
		!~PRINT*,'FOR EL:', i, j
		!~PRINT*,'RATIO:',lenBI/diaglen
		!~PRINT*,'NEM:', nem
		!~PRINT*,'xBI:', xBIt(ind)
		!~PRINT*,'yBI:', yBIt(ind)
		!~PRINT*,'xIP:', xIPt(ind)
		!~PRINT*,'yIP:', yIPt(ind)
		!~PRINT*,'SUM:', SQRT(xBIt(ind)**2+yBIt(ind)**2)
		
		!Find four nodes surrounding image point
		distMin = 1.0E10
		DO ki = i-3, i+2
			DO kj = j-3, j+2
				dist = SQRT((xtemp(ki,kj)-xIPt(ind))**2+ & 
							(ytemp(ki,kj)-yIPt(ind))**2)
				IF (dist .LT. distMin) THEN
					distMin = dist
					i1 = ki; j1 = kj
				END IF
			END DO
		END DO
		IF (xIPt(ind) .LT. xtemp(i1,j1)) THEN
			is = i1-1; ie = i1
		ELSE
			is = i1; ie = i1+1
		END IF
		
		IF (yIPt(ind) .LT. ytemp(i1,j1)) THEN
			js = j1-1; je = j1
		ELSE
			js = j1; je = j1+1
		END IF
		
		nb14(ind,1) = is; nb14(ind,2) = js
		
		!~PRINT*,'FOUR NODES SURROUNDING IP:'
		!~PRINT*,is, js
		!~PRINT*,is+1, js
		!~PRINT*,is+1, js+1
		!~PRINT*,is, js+1
		!~PRINT*,is-1, js-1
		!~PRINT*,is+1-1, js-1
		!~PRINT*,is+1-1, js+1-1
		!~PRINT*,is-1, js+1-1
		
		!Assign values to global variables
		IF (ivars .EQ. PRESS) THEN
			xBI = xBIt; yBI = yBIt
			xIP = xIPt; yIP = yIPt
			nb4v1 = nb14
		ELSEIF (ivars .EQ. XVEL) THEN
			xBIu = xBIt; yBIu = yBIt
			xIPu = xIPt; yIPu = yIPt
			nb4v1u = nb14
		ELSEIF (ivars .EQ. YVEL) THEN
			xBIv = xBIt; yBIv = yBIt
			xIPv = xIPt; yIPv = yIPt
			nb4v1v = nb14
		END IF
		
	 END DO
	 
	 CALL calc_coeffs_nb(ivars)
	 
	 !Deallocate memory
	 !-----------------!
	 DEALLOCATE(xtemp, ytemp)
	 DEALLOCATE(tempGC, tempGmap)
	 DEALLOCATE(iBlankT)
	 DEALLOCATE(xBIt, yBIt, xIPt, yIPt)
	 DEALLOCATE(nb14)

 END DO
 
 !Process neighbor cells for marker elements
 !------------------------------------------!
 DO i = 1, nElemMarker
	x1 = xElemCent(i); y1 = yElemCent(i)
	DO ki = 2, nci-1
		DO kj = 2, ncj-1
			IF (  (x1.GT.xc(ki,kj)) .AND. (x1.LE.xc(ki+1,kj)) & 
			.AND. (y1.GT.yc(ki,kj)) .AND. (y1.LE.yc(ki,kj+1))) THEN
			nbMarker(i,1) = ki; nbmarker(i,2) = kj
			END IF
		END DO
	END DO
	!~PRINT*,'FOUR NODES SURROUNDING MARKER:', i
	!~PRINT*,'Xecent, Yecent', x1, y1
	!~PRINT*,nbMarker(i,1), nbMarker(i,2)
	!~PRINT*,nbMarker(i,1)+1, nbMarker(i,2)
	!~PRINT*,nbMarker(i,1)+1, nbMarker(i,2)+1
	!~PRINT*,nbMarker(i,1), nbMarker(i,2)+1
 END DO
 
 !~CALL compute_LSF_flo
 CALL compute_LSF_marker
 
 WRITE(*,'(A)') ' DONE!'
 PRINT*,''
 
 END SUBROUTINE create_GCM
 !----------------------------------------------------------------------------!


 SUBROUTINE compute_LSF_marker
 
 USE grid_params
 USE global
 USE flo_vars
 USE GCM
 IMPLICIT NONE

 INTEGER :: ind
 
 REAL(KIND = 8) :: x1, x2, x3, x4, y1, y2, y3, y4
 REAL(KIND = 8) :: Lx, Ly
 REAL(KIND = 8) :: xi, eta
 
 DO ind = 1, nbodyMarker
	x1 = xc(nbMarker(ind,1), nbMarker(ind,2))
	y1 = yc(nbMarker(ind,1), nbMarker(ind,2))
	
	x2 = xc(nbMarker(ind,1)+1, nbMarker(ind,2))
	y2 = yc(nbMarker(ind,1)+1, nbMarker(ind,2))
	
	x3 = xc(nbMarker(ind,1)+1, nbMarker(ind,2)+1)
	y3 = yc(nbMarker(ind,1)+1, nbMarker(ind,2)+1)
	
	x4 = xc(nbMarker(ind,1), nbMarker(ind,2)+1)
	y4 = yc(nbMarker(ind,1), nbMarker(ind,2)+1)
	
	Lx = x2-x1; Ly = y4-y1
		
	xi = (xElemCent(ind)-x1)/Lx; eta = (yElemCent(ind)-y1)/Ly
	
	LSFMarker(ind,1) = (1.-xi)*(1.-eta)
	LSFMarker(ind,2) = (xi)*(1.-eta)
	LSFMarker(ind,3) = (xi)*(eta)
	LSFMarker(ind,4) = (1.-xi)*(eta)
 END DO
 
 END SUBROUTINE compute_LSF_marker
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE calc_coeffs_nb (ivars)
 
 USE grid_params
 USE global
 USE flo_vars
 USE GCM
 IMPLICIT NONE
 
 INTEGER, INTENT (IN) :: ivars
 INTEGER, PARAMETER :: PRESS = 1, XVEL = 2, YVEL = 3
 INTEGER, ALLOCATABLE, DIMENSION (:,:) :: nb14
 INTEGER :: ind, nGC1
 INTEGER :: is, ie, js, je
 
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: xtemp, ytemp
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:) :: xIPt, yIPt
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: phit
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:) :: tempvar
 REAL(KIND = 8) :: x1, x2, x3, x4, y1, y2, y3, y4, xpip, ypip
 REAL(KIND = 8) :: phi(4), phip, q1, q2 
 
 IF (ivars .EQ. PRESS) THEN
	nGC1 = nGC
	is = 1; ie = nci
	js = 1; je = ncj
	
	ALLOCATE(xtemp(ie-is+1,je-js+1), ytemp(ie-is+1,je-js+1))
	ALLOCATE(tempvar(nGC1))
	ALLOCATE(nb14(nGC1,2), xIPt(nGC1), yIPt(nGC1))
	ALLOCATE(phit(ie-is+1,je-js+1))
	xtemp = xc; ytemp = yc
	phit = pr
	nb14 = nb4v1; xIPt = xIP; yIPt = yIP
 ELSEIF (ivars .EQ. XVEL) THEN
	nGC1 = nGCu
	is = 1; ie = nci+1
	js = 1; je = ncj
	
	ALLOCATE(xtemp(ie-is+1,je-js+1), ytemp(ie-is+1,je-js+1))
	ALLOCATE(tempvar(nGC1))
	ALLOCATE(nb14(nGC1,2), xIPt(nGC1), yIPt(nGC1))
	ALLOCATE(phit(ie-is+1,je-js+1))
	xtemp = ux; ytemp = uy
	phit = velx
	nb14 = nb4v1u; xIPt = xIPu; yIPt = yIPu
 ELSEIF (ivars .EQ. YVEL) THEN
	nGC1 = nGCv
	is = 1; ie = nci
	js = 1; je = ncj+1
	
	ALLOCATE(xtemp(ie-is+1,je-js+1), ytemp(ie-is+1,je-js+1))
	ALLOCATE(tempvar(nGC1))
	ALLOCATE(nb14(nGC1,2), xIPt(nGC1), yIPt(nGC1))
	ALLOCATE(phit(ie-is+1,je-js+1))
	xtemp = vx; ytemp = vy
	phit = vely
	nb14 = nb4v1v; xIPt = xIPv; yIPt = yIPv
 END IF
 

 DO ind = 1, nGC1
	x1 = xtemp(nb14(ind,1),nb14(ind,2))
	y1 = ytemp(nb14(ind,1),nb14(ind,2))
	
	x2 = xtemp(nb14(ind,1)+1,nb14(ind,2))
	y2 = ytemp(nb14(ind,1)+1,nb14(ind,2))
	
	x3 = xtemp(nb14(ind,1)+1,nb14(ind,2)+1)
	y3 = ytemp(nb14(ind,1)+1,nb14(ind,2)+1)
	
	x4 = xtemp(nb14(ind,1),nb14(ind,2)+1)
	y4 = ytemp(nb14(ind,1),nb14(ind,2)+1)

	xpip = xIPt(ind); ypip = yIPt(ind)
	
	!Assign phi values
	phi(1) = phit(nb14(ind,1),nb14(ind,2))
	phi(2) = phit(nb14(ind,1)+1,nb14(ind,2))
	phi(3) = phit(nb14(ind,1)+1,nb14(ind,2)+1)
	phi(4) = phit(nb14(ind,1),nb14(ind,2)+1)
	
	!Compute phi at image point
	q1 = phi(1) + (phi(4)-phi(1))/(y4-y1)*(ypip-y1)
	q2 = phi(2) + (phi(3)-phi(2))/(y3-y2)*(ypip-y2)
	phip = q1 + (q2-q1)/(x2-x1)*(xpip-x1)
	tempvar(ind) = phip
 END DO 
 
 IF (ivars .EQ. PRESS) THEN
	prIP = tempvar
 ELSEIF (ivars .EQ. XVEL) THEN
	uIP = tempvar
 ELSEIF (ivars .EQ. YVEL) THEN
	vIP = tempvar
 END IF
 
 DEALLOCATE(xtemp, ytemp, xIPt, yIPt, nb14, tempvar)
 
 END SUBROUTINE calc_coeffs_nb
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE impose_IBM_bc
 
 USE grid_params
 USE global
 USE flo_vars
 USE GCM
 IMPLICIT NONE
 
 INTEGER, PARAMETER :: PRESS = 1, XVEL = 2, YVEL = 3
 INTEGER :: ind, i, j
 
 REAL (KIND = 8) :: xxGC, yyGC, xxBI, yyBI, xxIP, yyIP
 REAL (KIND = 8) :: d1, d2
 REAL (KIND = 8), PARAMETER :: uBIbc = 0., vBIbc = 0., prBIbc = 0.
 
 !Set velocity boundary conditions
 !X-velocity
 !CALL calc_coeffs_nb(XVEL)
 DO ind = 1, nGCu
	i = GCu(ind,1); j = GCu(ind,2)
	xxGC = ux(i,j); yyGC = uy(i,j)
	xxBI = xBIu(ind); yyBI = yBIu(ind)
	xxIP = xIPu(ind); yyIP = yIPu(ind)
	d1 = SQRT( (xxGC-xxBI)**2+(yyGC-yyBI)**2 )
	d2 = SQRT( (xxIP-xxBI)**2+(yyIP-yyBI)**2 )
	velxH(i,j) = ((d1+d2)*uBIbc - d1*uIP(ind))/d2
 END DO
 
 !Y-velocity
 !CALL calc_coeffs_nb(YVEL)
 DO ind = 1, nGCv
	i = GCv(ind,1); j = GCv(ind,2)
	xxGC = vx(i,j); yyGC = vy(i,j)
	xxBI = xBIv(ind); yyBI = yBIv(ind)
	xxIP = xIPv(ind); yyIP = yIPv(ind)
	d1 = SQRT( (xxGC-xxBI)**2+(yyGC-yyBI)**2 )
	d2 = SQRT( (xxIP-xxBI)**2+(yyIP-yyBI)**2 )
	velyH(i,j) = ((d1+d2)*vBIbc - d1*vIP(ind))/d2
 END DO
 
 !Set pressure boundary conditions
 !CALL calc_coeffs_nb(PRESS)
 DO ind = 1, nGC
	i = GCp(ind,1); j = GCp(ind,2)
	xxGC = xc(i,j); yyGC = yc(i,j)
	xxBI = xBI(ind); yyBI = yBI(ind)
	xxIP = xIP(ind); yyIP = yIP(ind)
	d1 = SQRT( (xxGC-xxBI)**2+(yyGC-yyBI)**2 )
	d2 = SQRT( (xxIP-xxBI)**2+(yyIP-yyBI)**2 )
	pr(i,j) = prIP(ind) - (d1+d2)*prBIbc	
 END DO
 
 END SUBROUTINE impose_IBM_bc
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE CDCL
 
 USE grid_params
 USE global
 USE flo_vars
 USE GCM
 IMPLICIT NONE
 
 INTEGER :: ind, ib, i, j
 INTEGER :: indnb(4,2)
 
 REAL(KIND = 8) :: sum_prx, sum_pry, CX, CY
 
 sum_prx = 0.; sum_pry = 0.
 
 DO ind = 1, nElemMarker
	!Compute pressure on ElemCent
	indnb(1,1) = nbMarker(ind,1); indnb(1,2) = nbMarker(ind,2)
	indnb(2,1) = nbMarker(ind,1)+1; indnb(2,2) = nbMarker(ind,2)
	indnb(3,1) = nbMarker(ind,1)+1; indnb(3,2) = nbMarker(ind,2)+1
	indnb(4,1) = nbMarker(ind,1); indnb(4,2) = nbMarker(ind,2)+1
	
	pBodyMarker(ind) = 0.
	!~PRINT*,''
	!~PRINT*,'nEM:', ind
	!~PRINT*,'press:'
	DO ib = 1, 4
		pBodyMarker(ind) = pBodyMarker(ind) + & 
		pr(indnb(ib,1),indnb(ib,2))*LSFMarker(ind,ib)
		!~PRINT*,pr(indnb(ib,1),indnb(ib,2)),LSFMarker(ind,ib)
	END DO
	!~PRINT*,'pr', pBodymarker(ind)
	sum_prx = sum_prx + pBodyMarker(ind)*-1.*normElX(ind)
	sum_pry = sum_pry + pBodyMarker(ind)*normElY(ind)
 END DO
 
 CX = sum_prx/(0.5*rho*(uinf**2+vinf**2))
 CY = sum_pry/(0.5*rho*(uinf**2+vinf**2))
 
 WRITE(CDCLbase,'(I6,3(2X,E14.7))') iterG, time, CX, CY
 
 END SUBROUTINE CDCL
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE free_GCM_vars
 
 USE GCM
 
 DEALLOCATE(xBodyMarker, yBodyMarker, xElemCent, yElemCent)
 DEALLOCATE(AreaBodyMarker)
 DEALLOCATE(pBodyMarker)
 DEALLOCATE(normElX, normElY, xBI, yBI, xIP, yIP)
 DEALLOCATE(ElemMarkerConnex)
 DEALLOCATE(gmap, gmapU, gmapV, gmapOut, GCp, GCu, GCv, iBlank, iBlankOut)
 DEALLOCATE(nb4v1, nb4v1u, nb4v1v, nbMarker)
 DEALLOCATE(LSFMarker)
 
 END SUBROUTINE free_GCM_vars
