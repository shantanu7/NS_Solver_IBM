 !Program to solve Navier-Stoke's equation
 !Written by Shantanu Bailoor at The George Washington University.
 !
 !Capabilities:
 !-------------
 ! 1. Code can simulate Taylor-Green Vortex, Lid Driven Cavity 
 !	  and Flow over stationary cylinder problems
 ! 2. Boundary conditions implemented: 
 !    Wall 				= 1
 !	  Inlet 			= 2 
 !	  Zero Grad Outlet	= 3
 !    Periodic 			= 4
 !	  Convective Outlet = 5
 !	  Symmetric			= 6
 ! 3. Ghost-cell based immersed boundary method implemented for FSI.
 ! 4. Coupled with SOR-based Poisson solver
 ! 5. Output file is cell-centered, Tecplot-compatible, ASCII data file (NS_out*.plt).
 
 PROGRAM NS_Solver

 USE grid_params
 USE global, ONLY: iterG, itermax, ntec, iRestart
 USE GCM, ONLY: iIntlBody
 IMPLICIT NONE

 CALL read_input

 IF (gen_opt .EQ. 0) THEN
	CALL make_grid 
 ELSEIF (gen_opt .EQ. 1) THEN
	CALL read_grid
 END IF	
 
 CALL allocate_flo_vars
 
 IF (iRestart .EQ. 0) THEN
	CALL init_flo_vars
 ELSE
	CALL read_restart_data
 END IF	
 
 CALL allocate_Poisson_memory
 
 CALL init_mon_files
 
 DO WHILE (iterG .LE. itermax)	 	 
	 CALL iterate
	 
	 IF (MOD(iterG,ntec) .EQ. 0) THEN
		 CALL out_data
		 CALL out_restart_data
		 IF (iIntlBody .EQ. 1) CALL out_marker_data
	 END IF
 END DO
 
 CALL close_mon_files
 CALL free_Poisson_memory
 CALL free_flo_vars
 CALL free_grid_params
 IF (iIntlBody .EQ. 1) CALL free_GCM_vars
 
 END PROGRAM NS_Solver
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE read_input
 
 USE grid_params
 USE global
 USE Poisson_cal_params
 USE flo_vars
 USE GCM
 IMPLICIT NONE
 
 INTEGER, PARAMETER :: finp = 101, fmon = 102, fmarker = 103
 INTEGER :: i
 
 WRITE(*,'(A)', ADVANCE = 'NO')' READING INPUT FILES...'
 
 OPEN(finp, FILE = 'NavSto.inp', ACTION = 'READ')
 READ(finp,*)
 READ(finp,*)
 READ(finp,*)
 READ(finp,*)Xmin, Xmax, Ymin, Ymax
 READ(finp,*)
 READ(finp,*)ngx, ngy, gen_opt
 
 READ(finp,*)
 READ(finp,*)
 READ(finp,*)Tmin, Tmax, dt, ntec
 
 READ(finp,*)
 READ(finp,*)
 READ(finp,*)iProblem, iRestart, nRead, iIntlBody, iCDCL
 READ(finp,*)
 READ(finp,*)Sol_method, time_scheme
 
 READ(finp,*)
 READ(finp,*)
 READ(finp,*) BoundCond(:)
 READ(finp,*)
 READ(finp,*) BCVal(:)
 
 READ(finp,*)
 READ(finp,*)
 READ(finp,*) residMon, Var_ex
 
 READ(finp,*)
 READ(finp,*)
 READ(finp,*) omega 

 READ(finp,*)
 READ(finp,*)
 READ(finp,*) iVerbose
 CLOSE(finp)
 
 WRITE(*,'(A)')' DONE!'
 PRINT*,''
 
 CALL disp_simulation_message
 
 WRITE(*,'(A)', ADVANCE = 'NO')' READING PROBE INPUT...'
 
 OPEN(fmon, FILE = 'mon.inp', ACTION = 'READ')
 
 READ(fmon,*)
 READ(fmon,*)nmon
 READ(fmon,*)
 
 ALLOCATE(xmon(nmon), ymon(nmon))
 ALLOCATE(imon(nmon), jmon(nmon))
 
 DO i = 1, nmon
	READ(fmon,*) xmon(i), ymon(i)
 END DO
 
 CLOSE(fmon)
 
 WRITE(*,'(A)')' DONE!'
 PRINT*,''
 
 IF (iIntlBody .EQ. 1) THEN
	WRITE(*,'(A)', ADVANCE = 'NO')' READING MARKER FILE...'
	OPEN(fmarker, FILE = 'canonical_body_in.dat', ACTION = 'READ')
	READ(fmarker,*)nBodyMarker
	ALLOCATE (xBodyMarker(nBodyMarker), yBodyMarker(nBodyMarker))
	DO i = 1, nBodyMarker
		READ(fmarker,*)xBodyMarker(i), yBodyMarker(i)
	END DO
	
	READ(fmarker,*)
	READ(fmarker,*) nElemMarker
	ALLOCATE (ElemMarkerConnex(nElemMarker,2))
	
	DO i = 1, nElemMarker
		READ(fmarker,*)ElemMarkerConnex(i,:)
	END DO
	CLOSE(fmarker)
	WRITE(*,'(A)')' DONE!'
	PRINT*,''
	
	WRITE(*,'(A,I5)') ' NUMBER OF MARKER POINTS READ: ', nBodyMarker
	PRINT*,''
 END IF
 
 END SUBROUTINE read_input
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE make_grid
 
 USE grid_params
 USE global
 IMPLICIT NONE
 
 INTEGER, PARAMETER :: fgrid = 101
 INTEGER :: i, j
 
 time = Tmin
 itermax = INT((Tmax-Tmin)/dt); iterG = 0
 
 !Define global grid parameters
 ngi = ngx + 2; ngj = ngy + 2
 ncx = ngx - 1; ncy = ngy - 1
 nci = ncx + 2; ncj = ncy + 2
 nodes = ngi*ngj; elem = nci*ncj
 
 IF (iProblem .EQ. TaylorGreen) THEN
	 Xmin = Xmin*2.*pi; Xmax = Xmax*2.*pi
	 Ymin = Ymin*2.*pi; Ymax = Ymax*2.*pi
 END IF
 
 dx = (Xmax-Xmin)/REAL(ncx); dy = (Ymax-Ymin)/REAL(ncy)
 XminG = Xmin - dx; XmaxG = Xmax + dx
 YminG = Ymin - dy; YmaxG = Ymax + dy
 
 WRITE(*,'(A)', ADVANCE = 'NO')' ALLOCATING MEMORY FOR GRID VARIABLES...'
 
 ALLOCATE ( x(ngi,ngj), y(ngi,ngj), xc(nci,ncj), yc(nci,ncj) )
 ALLOCATE ( xgm(ngi,ngj), ygm(ngi,ngj) )
 ALLOCATE ( px(nci,ncj), py(nci,ncj), ux(nci+1,ncj), uy(nci+1,ncj), &
			vx(nci,ncj+1), vy(nci,ncj+1) )

 ALLOCATE ( dxue(nci+1,ncj), dxuw(nci+1,ncj), dyun(nci+1,ncj), dyus(nci+1,ncj) )
 ALLOCATE ( dxve(nci,ncj+1), dxvw(nci,ncj+1), dyvn(nci,ncj+1), dyvs(nci,ncj+1) )
 
 WRITE(*,'(A)')' DONE!'
 PRINT*,''
 
 WRITE(*,'(A)', ADVANCE = 'NO')' CREATING STAGGERED GRID...'
 
 !Define position for grid points
 DO i = 1, ngi
	xgm(i,:) = Xmin + ((i-1)-1)*dx
	x(i,:) = Xmin + ((i-1)-1)*dx
 END DO
 DO j = 1, ngj
	ygm(:,j) = Ymin + ((j-1)-1)*dy
	y(:,j) = Ymin + ((j-1)-1)*dy
 END DO
 
 !Define position for cell centers
 DO i = 1, nci
	xc(i,:) = Xmin - dx/2. + (i-1)*dx
 END DO
 DO j = 1, ncj
	yc(:,j) = Ymin - dy/2. + (j-1)*dy
 END DO
 
 !Define position for pressure nodes
 DO i = 1, nci
	px(i,:) = Xmin - dx/2. + (i-1)*dx
 END DO
 DO j = 1, ncj
	py(:,j) = Ymin - dy/2. + (j-1)*dy
 END DO
 
 !Define position for x-velocity nodes
 DO i = 1, nci+1
	ux(i,:) = Xmin + (i-2)*dx
	IF (i .NE. 1) THEN
		dxue(i-1,:) = ux(i,:) - ux(i-1,:)
		dxuw(i,:) = ux(i,:) - ux(i-1,:)
	END IF
 END DO
 DO j = 1, ncj
	uy(:,j) = Ymin - dy/2. + (j-1)*dy
	IF (j .NE. 1) THEN
		dyun(:,j-1) = uy(:,j) - uy(:,j-1)
		dyus(:,j) = uy(:,j) - uy(:,j-1)
	END IF
 END DO
 
 !Define position for y-velocity nodes
 DO i = 1, nci
	vx(i,:) = Xmin -dx/2. + (i-1)*dx
	IF (i .NE. 1) THEN
		dxve(i-1,:) = vx(i,:) - vx(i-1,:)
		dxvw(i,:) = vx(i,:) - vx(i-1,:)
	END IF
 END DO
 DO j = 1, ncj+1
	vy(:,j) = Ymin + (j-2)*dy
	IF (j .NE. 1) THEN
		dyvn(:,j-1) = vy(:,j) - vy(:,j-1)
		dyvs(:,j) = vy(:,j) - vy(:,j-1)
	END IF
 END DO
 
 !Create Tecplot compatible grid file
 OPEN(fgrid, FILE = 'grid.plt')
	WRITE(fgrid,'(A)') ' TITLE = "STAGGERED GRID ARRANGEMENT"'
	WRITE(fgrid,'(A)') ' VARIABLES = "X", "Y"'
	WRITE(fgrid,71) ' ZONE T = "MESH", I = ',ngi,', J = ',ngj,', F = POINT'
	DO j = 1, ngj
		DO i = 1, ngi
			WRITE(fgrid,'(F12.9,5X,F12.9)') xgm(i,j), ygm(i,j)
		END DO
	END DO
	
	WRITE(fgrid,71) ' ZONE T = "P-NODES", I = ',nci,', J = ',ncj,', F = POINT'
	DO j = 1, ncj
		DO i = 1, nci
			WRITE(fgrid,'(F12.9,5X,F12.9)') px(i,j), py(i,j)
		END DO
	END DO
	
	WRITE(fgrid,71) ' ZONE T = "U-NODES", I = ',nci+1,', J = ',ncj,', F = POINT'
	DO j = 1, ncj
		DO i = 1, nci+1
			WRITE(fgrid,'(F12.9,5X,F12.9)') ux(i,j), uy(i,j)
		END DO
	END DO
	
	WRITE(fgrid,71) ' ZONE T = "V-NODES", I = ',nci,', J = ',ncj+1,', F = POINT'
	DO j = 1, ncj+1
		DO i = 1, nci
			WRITE(fgrid,'(F12.9,5X,F12.9)') vx(i,j), vy(i,j)
		END DO
	END DO
	
 CLOSE(fgrid) 
 
 WRITE(*,'(A)')' DONE!'
 PRINT*,''
 
 71 FORMAT (A,I3,A,I3,A)
 
 END SUBROUTINE make_grid
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE read_grid
 
 USE grid_params
 USE global
 IMPLICIT NONE
 
 INTEGER, PARAMETER :: fgrid = 101, Agrid = 102
 INTEGER :: i, j
 
 time = Tmin
 itermax = INT((Tmax-Tmin)/dt); iterG = 0
 
 OPEN (Agrid, FILE='Agrid.dat', ACTION = 'READ')

 READ(Agrid,*)
 READ(Agrid,*) ngx, ngy
 
 !Define global grid parameters
 ngi = ngx + 2; ngj = ngy + 2
 ncx = ngx - 1; ncy = ngy - 1
 nci = ncx + 2; ncj = ncy + 2
 nodes = ngi*ngj; elem = nci*ncj
 
 WRITE(*,'(A)', ADVANCE = 'NO')' ALLOCATING MEMORY FOR GRID VARIABLES...'
 
 ALLOCATE ( x(ngi,ngj), y(ngi,ngj), xc(nci,ncj), yc(nci,ncj) )
 ALLOCATE ( xgm(ngi,ngj), ygm(ngi,ngj) )
 ALLOCATE ( px(nci,ncj), py(nci,ncj), ux(nci+1,ncj), uy(nci+1,ncj), &
			vx(nci,ncj+1), vy(nci,ncj+1) )

 ALLOCATE ( dxue(nci+1,ncj), dxuw(nci+1,ncj), dyun(nci+1,ncj), dyus(nci+1,ncj) )
 ALLOCATE ( dxve(nci,ncj+1), dxvw(nci,ncj+1), dyvn(nci,ncj+1), dyvs(nci,ncj+1) )
 
 WRITE(*,'(A)')' DONE!'
 PRINT*,''
 
 WRITE(*,'(A)', ADVANCE = 'NO')' READING GRID COORDINATES...'
 
 x = 0.; y = 0.
 
 READ(Agrid,*)x(2:ngi-1,1)
 READ(Agrid,*)y(1,2:ngj-1)
 
 !Create ghost layers
 x(1,1) = x(2,1) - (x(3,1)-x(2,1))
 x(ngi,1) = x(ngi-1,1) + (x(ngi-1,1)-x(ngi-2,1))
 
 y(1,1) = y(1,2) - (y(1,3)-y(1,2))
 y(1,ngj) = y(1,ngj-1) + (y(1,ngj-1)-y(1,ngj-2))
 
 Xmin = x(2,1); XminG = x(1,1)
 Ymin = y(1,2); YminG = y(1,1)
 
 DO i = 2, ngi
	y(i,:) = y(1,:)
 END DO
 
 DO j = 2, ngj
	x(:,j) = x(:,1)
 END DO
 
 WRITE(*,'(A)')' DONE!'
 PRINT*,''
 
 CLOSE (Agrid)
 
  IF (iProblem .EQ. TaylorGreen) THEN
	 x = x*2.*pi; y = y*2.*pi; 
 END IF
 
 !Other needed processing for grids
 xgm = x; ygm = y
 
 !Cell centers
 DO i = 1, nci
	xc(i,:) = 0.5*(x(i,:)+x(i+1,:))
 END DO
 
 DO j = 1, ncj
	yc(:,j) = 0.5*(y(:,j)+y(:,j+1))
 END DO
 
 !Define position for pressure nodes
 px = xc; py = yc
 
 !Define position for x-velocity nodes
 ux = x
 uy = yc
 
 DO i = 1, nci+1
	IF (i .NE. 1) THEN
		dxue(i-1,:) = ux(i,:) - ux(i-1,:)
		dxuw(i,:) = ux(i,:) - ux(i-1,:)
	END IF
 END DO
 DO j = 1, ncj
	IF (j .NE. 1) THEN
		dyun(:,j-1) = uy(:,j) - uy(:,j-1)
		dyus(:,j) = uy(:,j) - uy(:,j-1)
	END IF
 END DO
 
 !Define position for y-velocity nodes
 vx = xc
 vy = y 
 
 DO i = 1, nci
	IF (i .NE. 1) THEN
		dxve(i-1,:) = vx(i,:) - vx(i-1,:)
		dxvw(i,:) = vx(i,:) - vx(i-1,:)
	END IF
 END DO
 DO j = 1, ncj+1
	IF (j .NE. 1) THEN
		dyvn(:,j-1) = vy(:,j) - vy(:,j-1)
		dyvs(:,j) = vy(:,j) - vy(:,j-1)
	END IF
 END DO 
 
 !Create Tecplot compatible grid file
 OPEN(fgrid, FILE = 'grid.plt')
	WRITE(fgrid,'(A)') ' TITLE = "STAGGERED GRID ARRANGEMENT"'
	WRITE(fgrid,'(A)') ' VARIABLES = "X", "Y"'
	WRITE(fgrid,71) ' ZONE T = "MESH", I = ',ngi,', J = ',ngj,', F = POINT'
	DO j = 1, ngj
		DO i = 1, ngi
			WRITE(fgrid,'(E13.6,5X,E13.6)') xgm(i,j), ygm(i,j)
		END DO
	END DO
	
	WRITE(fgrid,71) ' ZONE T = "P-NODES", I = ',nci,', J = ',ncj,', F = POINT'
	DO j = 1, ncj
		DO i = 1, nci
			WRITE(fgrid,'(E13.6,5X,E13.6)') px(i,j), py(i,j)
		END DO
	END DO
	
	WRITE(fgrid,71) ' ZONE T = "U-NODES", I = ',nci+1,', J = ',ncj,', F = POINT'
	DO j = 1, ncj
		DO i = 1, nci+1
			WRITE(fgrid,'(E13.6,5X,E13.6)') ux(i,j), uy(i,j)
		END DO
	END DO
	
	WRITE(fgrid,71) ' ZONE T = "V-NODES", I = ',nci,', J = ',ncj+1,', F = POINT'
	DO j = 1, ncj+1
		DO i = 1, nci
			WRITE(fgrid,'(E13.6,5X,E13.6)') vx(i,j), vy(i,j)
		END DO
	END DO
	
 CLOSE(fgrid) 
 
 WRITE(*,'(A)')' DONE!'
 PRINT*,''
 
 71 FORMAT (A,I3,A,I3,A)
 
 END SUBROUTINE read_grid
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE allocate_flo_vars
 
 USE flo_vars
 USE grid_params
 IMPLICIT NONE
 
 ALLOCATE (velx(nci+1,ncj), vely(nci,ncj+1), pr(nci,ncj))
 ALLOCATE (velxH(nci+1,ncj), velyH(nci,ncj+1))
 ALLOCATE (Fc_u(nci+1,ncj), Fc_v(nci,ncj+1))
 ALLOCATE (Fv_u(nci+1,ncj), Fv_v(nci,ncj+1))
 ALLOCATE (Fp_u(nci+1,ncj), Fp_v(nci,ncj+1))
 ALLOCATE (Gi_u(nci+1,ncj), Gi_v(nci,ncj+1)) 
 
 END SUBROUTINE allocate_flo_vars
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE init_flo_vars
 
 USE grid_params
 USE global
 USE flo_vars
 USE GCM
 IMPLICIT NONE
 
 REAL (KIND = 8) :: xi, yi
 
 INTEGER :: i, j, is, ie, js, je
 
 iToggle = 0
 
 !Initialize flow fields
 velx = 0.; vely = 0.; pr = 0.
 velxH = 0.; velyH = 0.
 IF (iProblem .EQ. TaylorGreen) THEN
	 !X-velocity
	 is = 1; ie = nci+1; js = 1; je = ncj
	 DO i = is, ie
		DO j = js, je
			xi = ux(i,j); yi = uy(i,j)
			velx(i,j) = -1.*exp(-2.*time)*COS(xi)*SIN(yi)
		END DO
	 END DO
	 !Y-velocity
	 is = 1; ie = nci; js = 1; je = ncj+1
	 DO i = is, ie
		DO j = js, je
			xi = vx(i,j); yi = vy(i,j)
			vely(i,j) = exp(-2.*time)*SIN(xi)*COS(yi)
		END DO
	 END DO
	 !Pressure
	 is = 1; ie = nci; js = 1; je = ncj
	 DO i = is, ie
		DO j = js, je
			xi = px(i,j); yi = py(i,j)
			pr(i,j) = -1./4.*exp(-4.*time)*(COS(2.*xi)+COS(2.*yi))
		END DO
	 END DO
 ELSEIF	 (iProblem .EQ. Cylinder) THEN
	velx = MAXVAL(ABS(BCVal)); vely = 0.
 END IF	 
 
 !Initialize flux variables
 Fc_u = 0.; Fc_v = 0.
 Fv_u = 0.; Fv_v = 0.
 Fp_u = 0.; Fp_v = 0.
 Gi_u = 0.; Gi_v = 0.
 
 uinf = 0.; vinf = 0.
 
 IF (iIntlBody .EQ. 1) CALL create_GCM
 
 CALL out_data
 
 END SUBROUTINE init_flo_vars 
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE init_mon_files
 
 USE grid_params
 USE global
 IMPLICIT NONE
 
 INTEGER :: i, j, ind
 
 REAL(KIND = 8) :: dist, distMin
 
 CHARACTER*50 :: fname, fname_ex, fCDCL
 
 !Find position of closest cell-center to each probe
 DO ind = 1, nmon
 
	distMin = 1.0E10
	DO i = 1, nci
		DO j = 1, ncj
			dist = SQRT((xc(i,j)-xmon(ind))**2 + (yc(i,j)-ymon(ind))**2)
			IF (dist .LT. distMin) THEN
				distMin = dist
				imon(ind) = i; jmon(ind) = j
			END IF
		END DO
	END DO
	WRITE(*,'(A,I1,A,2F7.3)')'MON ', ind,' FOUND AT: ', & 
	xc(imon(ind),jmon(ind)), yc(imon(ind),jmon(ind))
 
	!Create mon files
	WRITE(fname,'(A,I3.3,A)')'mon', ind, '.dat'
	IF (iProblem .EQ. TaylorGreen) WRITE(fname_ex,'(A,I3.3,A)')'mon_ex', ind, '.dat'
	
	IF (iRestart .EQ. 0) THEN
		OPEN(monBase+ind, FILE = TRIM(fname))
		WRITE(monBase+ind,'(A)')'VARIABLES = iter, time, u, v, p'
		IF (iProblem .EQ. TaylorGreen) THEN
			OPEN(monBaseEx+ind, FILE = TRIM(fname_ex))
			WRITE(monBaseEx+ind,'(A)')'VARIABLES = iter, time, u, v, p'
		END IF
	ELSEIF (iRestart .EQ. 1) THEN
		OPEN(monBase+ind, FILE = TRIM(fname), POSITION = 'APPEND')
		IF (iProblem .EQ. TaylorGreen) THEN
			OPEN(monBaseEx+ind, FILE = TRIM(fname_ex), POSITION = 'APPEND')
		END IF
	END IF	
	
 END DO 
 
 PRINT*,''
 
 IF (iProblem .EQ. Cylinder .AND. iCDCL .EQ. 1) THEN
	 WRITE(*,'(A)', ADVANCE = 'NO')' OPENING CDCL FILE...'
	 WRITE(fCDCL,'(A)')'CDCL.dat'
	 IF (iRestart .EQ. 0) THEN
		OPEN(CDCLbase, FILE = TRIM(fCDCL))
		WRITE(CDCLbase,*)'VARIABLES = "iter", "time", "CX", "CY"'
	 ELSE
		OPEN(CDCLbase, FILE = TRIM(fCDCL), POSITION = 'APPEND')
	 END IF	 
	 WRITE(*,'(A)')' DONE!'
	 PRINT*,''
 END IF
 
 END SUBROUTINE init_mon_files
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE set_velBC
 
 USE grid_params
 USE global
 USE flo_vars
 USE cal_params
 IMPLICIT NONE
 
 INTEGER :: i, j
 
 REAL (KIND = 8) :: xy1, xy2, UBC, VBC
 
 !Left Boundary:
 !-------------!
 IF (BoundCond(1) .EQ. WALL ) THEN
	UBC = 0.; VBC = BCVal(1)
	!X-velocity
	xy1 = ux(2,1)-ux(1,1); xy2 = ux(3,1)-ux(2,1)
	velx(1,:) = ((xy1+xy2)*UBC-xy1*velx(3,:))/xy2
	velx(2,:) = UBC
	!Y-velocity
	xy1 = x(2,1)-vx(1,1); xy2 = vx(2,1)-x(2,1)
	vely(1,:) = ((xy1+xy2)*VBC-xy1*vely(2,:))/xy2
 ELSEIF	(BoundCond(1) .EQ. INLET) THEN
	UBC = BCVal(1); VBC = 0.
	!X-velocity
	!~xy1 = ux(2,1)-ux(1,1); xy2 = ux(3,1)-ux(2,1)
	!~velx(1,:) = ((xy1+xy2)*UBC-xy1*velx(3,:))/xy2
	velx(2,:) = UBC
	!Y-velocity
	xy1 = x(2,1)-vx(1,1); xy2 = vx(2,1)-x(2,1)
	vely(1,:) = ((xy1+xy2)*VBC-xy1*vely(2,:))/xy2
	IF (iProblem .EQ. Cylinder) uinf = UBC
 ELSEIF (BoundCond(1) .EQ. OUTLET) THEN
	UBC = BCVal(1); VBC = 0.
	!X-velocity
	velx(1,:) = velx(2,:) + UBC*(ux(2,:)-ux(1,:))
	!Y-velocity
	xy1 = x(2,1)-vx(1,1); xy2 = vx(2,1)-x(2,1)
	vely(1,:) = ((xy1+xy2)*VBC-xy1*vely(2,:))/xy2
 ELSEIF (BoundCond(1) .EQ. PERIODIC) THEN
	!X-velocity
	velx(1,:) = velx(nci-1,:)
	!Y-velocity
	vely(1,:) = vely(nci-1,:)
	!Pressure
	pr(1,:) = pr(nci-1,:)
 ELSEIF (BoundCond(1) .EQ. CONVECT) THEN
	velx(1,:) = BCVal(1)
 ELSEIF (BoundCond(1) .EQ. SYMM) THEN
	velx(1,:) = -1.*velx(3,:)
	vely(1,:) = vely(2,:)
 END IF
 
 !Right Boundary:
 !--------------!
 IF (BoundCond(2) .EQ. WALL ) THEN
	UBC = 0.; VBC = BCVal(2)
	!X-velocity
	xy1 = ux(nci+1,1)-ux(nci,1); xy2 = ux(nci,1)-ux(nci-1,1)
	velx(nci+1,:) = ((xy1+xy2)*UBC-xy1*velx(nci-1,:))/xy2
	velx(nci,:) = UBC
	!Y-velocity
	xy1 = vx(nci,1)-x(ngi-1,1); xy2 = x(ngi-1,1)-vx(nci-1,1)
	vely(nci,:) = ((xy1+xy2)*VBC-xy1*vely(nci-1,:))/xy2
 ELSEIF	(BoundCond(2) .EQ. INLET) THEN
	UBC = BCVal(2); VBC = 0.
	!X-velocity
	xy1 = ux(nci+1,1)-ux(nci,1); xy2 = ux(nci,1)-ux(nci-1,1)
	velx(nci+1,:) = ((xy1+xy2)*UBC-xy1*velx(nci-1,:))/xy2
	!Y-velocity
	xy1 = vx(nci,1)-x(ngi-1,1); xy2 = x(ngi-1,1)-vx(nci-1,1)
	vely(nci,:) = ((xy1+xy2)*VBC-xy1*vely(nci-1,:))/xy2
	IF (iProblem .EQ. Cylinder) uinf = UBC
 ELSEIF (BoundCond(2) .EQ. OUTLET) THEN
	UBC = BCVal(2); VBC = 0.
	!X-velocity
	velx(nci+1,:) = velx(nci,:) + UBC*(ux(nci+1,:)-ux(nci,:))
	!Y-velocity
	xy1 = vx(nci,1)-x(ngi-1,1); xy2 = x(ngi-1,1)-vx(nci-1,1)
	vely(nci,:) = ((xy1+xy2)*VBC-xy1*vely(nci-1,:))/xy2
 ELSEIF (BoundCond(2) .EQ. PERIODIC) THEN
	!X-velocity
	velx(nci+1,:) = velx(3,:)
	!Y-velocity
	vely(nci,:) = vely(2,:)
	!Pressure
	pr(nci,:) = pr(2,:)
 ELSEIF (BoundCond(2) .EQ. CONVECT) THEN
	DO j =1, ncj
		velx(nci+1,j) = velx(nci+1,j) + dt2*BCVal(2)* &
		(velx(nci-1,j)-velx(nci,j))/(ux(nci,j)-ux(nci-1,j))
	END DO
 ELSEIF (BoundCond(2) .EQ. SYMM) THEN
	velx(nci+1,:) = -1.*velx(nci-1,:)	
	vely(nci,:) = vely(nci-1,:)
 END IF
 
 !Bottom Boundary:
 !---------------!
 IF (BoundCond(3) .EQ. WALL ) THEN
	UBC = BCVal(3); VBC = 0.
	!X-velocity
	xy1 = y(1,2)-uy(1,1); xy2 = uy(1,2)-x(1,2)
	velx(:,1) = ((xy1+xy2)*UBC-xy1*velx(:,2))/xy2
	!Y-velocity
	xy1 = vy(1,2)-vy(1,1); xy2 = vy(1,3)-vy(1,2)
	vely(:,1) = ((xy1+xy2)*VBC-xy1*vely(:,3))/xy2
	vely(:,2) = VBC
 ELSEIF	(BoundCond(3) .EQ. INLET) THEN
	UBC = 0.; VBC = BCVal(3)
	!X-velocity
	xy1 = y(1,2)-uy(1,1); xy2 = uy(1,2)-x(1,2)
	velx(:,1) = ((xy1+xy2)*UBC-xy1*velx(:,2))/xy2
	!Y-velocity
	xy1 = vy(1,2)-vy(1,1); xy2 = vy(1,3)-vy(1,2)
	vely(:,1) = ((xy1+xy2)*VBC-xy1*vely(:,3))/xy2
	IF (iProblem .EQ. Cylinder) vinf = VBC
 ELSEIF (BoundCond(3) .EQ. OUTLET) THEN
	UBC = 0.; VBC = BCVal(3)
	!X-velocity
	xy1 = y(1,2)-uy(1,1); xy2 = uy(1,2)-x(1,2)
	velx(:,1) = ((xy1+xy2)*UBC-xy1*velx(:,2))/xy2
	!Y-velocity
	vely(:,1) = vely(:,2) + VBC*(vy(:,2)-vy(:,1))
 ELSEIF (BoundCond(3) .EQ. PERIODIC) THEN
	!X-velocity
	velx(:,1) = velx(:,ncj-1)
	!Y-velocity
	vely(:,1) = vely(:,ncj-1)
	!Pressure
	pr(:,1) = pr(:,ncj-1)
 ELSEIF (BoundCond(3) .EQ. CONVECT) THEN
	velx(:,1) = BCVal(3)
 ELSEIF (BoundCond(3) .EQ. SYMM) THEN
	!vely(:,1) = -1.*vely(:,3)
	vely(:,1) = 0.; vely(:,2) = 0.
	velx(:,1) = velx(:,2)	
	!velx(:,1) = BCVal(1)
 END IF
 
 !Top Boundary:
 !------------!
 IF (BoundCond(4) .EQ. WALL ) THEN
	UBC = BCVal(4); VBC = 0.
	!X-velocity
	xy1 = uy(1,ncj)-y(1,ngj-1); xy2 = y(1,ngj-1)-uy(1,ncj-1)
	velx(:,ncj) = ((xy1+xy2)*UBC-xy1*velx(:,ncj-1))/xy2
	!Y-velocity
	xy1 = vy(1,ncj+1)-vy(1,ncj); xy2 = vy(1,ncj)-vy(1,ncj-1)
	vely(:,ncj+1) = ((xy1+xy2)*VBC-xy1*vely(:,ncj-1))/xy2
	vely(:,ncj) = VBC
 ELSEIF	(BoundCond(4) .EQ. INLET) THEN
	UBC = 0.; VBC = BCVal(4)
	!X-velocity
	xy1 = uy(1,ncj)-y(1,ngj-1); xy2 = y(1,ngj-1)-uy(1,ncj-1)
	velx(:,ncj) = ((xy1+xy2)*UBC-xy1*velx(:,ncj-1))/xy2
	!Y-velocity
	xy1 = vy(1,ncj+1)-vy(1,ncj); xy2 = vy(1,ncj)-vy(1,ncj-1)
	vely(:,ncj+1) = ((xy1+xy2)*VBC-xy1*vely(:,ncj-1))/xy2
	IF (iProblem .EQ. Cylinder) vinf = VBC
 ELSEIF (BoundCond(4) .EQ. OUTLET) THEN
	UBC = 0.; VBC = BCVal(4)
	!X-velocity
	xy1 = uy(1,ncj)-y(1,ngj-1); xy2 = y(1,ngj-1)-uy(1,ncj-1)
	velx(:,ncj) = ((xy1+xy2)*UBC-xy1*velx(:,ncj-1))/xy2
	!Y-velocity
	vely(:,ncj) = vely(:,ncj-1) + VBC*(vy(:,ncj)-vy(:,ncj-1))
 ELSEIF (BoundCond(4) .EQ. PERIODIC) THEN
	!X-velocity
	velx(:,ncj) = velx(:,2)
	!Y-velocity
	vely(:,ncj+1) = vely(:,3)
	!Pressure
	pr(:,ncj) = pr(:,2)
 ELSEIF (BoundCond(4) .EQ. CONVECT) THEN
	velx(:,ncj+1) = BCVal(4)	
 ELSEIF (BoundCond(4) .EQ. SYMM) THEN
	!vely(:,ncj+1) = -1.*vely(:,ncj-1)
	vely(:,ncj+1) = 0.; vely(:,ncj) = 0.
	velx(:,ncj) = velx(:,ncj-1)
	!velx(:,ncj) = BCVal(1)
 END IF
 
 !~IF (iVerbose .EQ. 1 .AND. iProblem .EQ. DrivenCavity) THEN
 !~PRINT*,'vel at top:', MAXVAL(0.5*(velx(:,ncj)+velx(:,ncj-1))), MAXVAL(vely(:,ncj))
 !~END IF
 
 END SUBROUTINE set_velBC
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE iterate
 
 USE grid_params
 USE global
 USE flo_vars
 USE cal_params
 USE Poisson_cal_params
 USE GCM
 IMPLICIT NONE
 
 INTEGER, PARAMETER :: X_MOM = 1, Y_MOM = 2
 INTEGER :: stage 
 INTEGER :: i, j, is, ie, js, je
 
 REAL (KIND = 8) :: divH, dpdx, dpdy, div, divex, duex, dvex, diffmax
 
 iterG = iterG + 1
 !time = time + dt
 
 IF (iterG .EQ. 1 .OR. MOD(iterG,1) .EQ. 0 .OR. iterG .EQ. itermax) THEN
	WRITE(*,'(A,I6,A,F8.5)')' >> ITER: ', iterG, ' TIME: ', time
 END IF
 
 SELECT CASE (time_scheme)
	CASE (RK3)
		DO stage = 1, 3
			dt2 = dtFactor(stage)*dt
			time = time + dt2
			
			CALL set_velBC
			
			!Compute Predicted x-velocities
			is = 2; ie = nci; js = 2; je = ncj-1
			CALL compute_flux(X_MOM)
			DO i = is, ie
				DO j = js, je
					Gi_u(i,j) = GiFactor(stage)*Gi_u(i,j) &
								+ Fc_u(i,j) + Fv_u(i,j)
					velxH(i,j) = velx(i,j) + dt2*Gi_u(i,j)
				END DO
			END DO
			
			!Compute Predicted y-velocities
			is = 2; ie = nci-1; js = 2; je = ncj
			CALL compute_flux(Y_MOM)
			DO i = is, ie
				DO j = js, je
					Gi_v(i,j) = GiFactor(stage)*Gi_v(i,j) &
								+ Fc_v(i,j) + Fv_v(i,j)
					velyH(i,j) = vely(i,j) + dt2*Gi_v(i,j)
				END DO
			END DO
			
			IF (iIntlBody .EQ. 1) CALL impose_IBM_bc
			
			!Compute divergence of Predicted velocity
			is = 2; ie = nci-1; js = 2; je = ncj-1
			diffmax = 0.
			DO i = is, ie
				DO j = js, je
					divH = (velxH(i+1,j)-velxH(i,j))/(ux(i+1,j)-ux(i,j)) &
						 + (velyH(i,j+1)-velyH(i,j))/(vy(i,j+1)-vy(i,j))
					f(i,j) = divH/dt2
				END DO
			END DO
			
			!Solve Poisson equation
			CALL Poisson_Solver
			IF (iterG .EQ. 1 .OR. MOD(iterG,1) .EQ. 0 .OR. iterG .EQ. itermax) THEN
				WRITE(*,'(A,I1,A)')' SOLVING POISSON EQUATION FOR STAGE ', stage, ' COMPLETE!'
			END IF
			
			is = 2; ie = nci; js = 2; je = ncj-1
			DO i = is, ie
				DO j = js, je
					dpdx = (pr(i,j)-pr(i-1,j))/(px(i,j)-px(i-1,j))
					velx(i,j) = velxH(i,j) - dt2*dpdx/rho
				END DO
			END DO
			is = 2; ie = nci-1; js = 2; je = ncj
			DO i = is, ie
				DO j = js, je
					dpdy = (pr(i,j)-pr(i,j-1))/(py(i,j)-py(i,j-1))
					vely(i,j) = velyH(i,j) - dt2*dpdy/rho
				END DO
			END DO
			IF (iterG .EQ. 1 .OR. MOD(iterG,1) .EQ. 0 .OR. iterG .EQ. itermax) THEN
				is = 2; ie = nci-1; js = 2; je = ncj-1
				div = 0.
				DO i = is, ie
					DO j = js, je
						div = div + (velx(i+1,j)-velx(i,j))/(ux(i+1,j)-ux(i,j)) &
							 + (vely(i,j+1)-vely(i,j))/(vy(i,j+1)-vy(i,j))
					END DO
				END DO
				PRINT*,'TOTAL DIV: ', div
				!IF (iVerbose .EQ. 1) PRINT*,'Umax:', MAXVAL(velx(2:nci,2:ncj-1))
				PRINT*,''
			END IF	
		END DO	
 END SELECT
 
 CALL out_mon_data
 IF (iCDCL .EQ. 1) CALL CDCL
 
 IF (iterG .EQ. 1 .OR. MOD(iterG,1) .EQ. 0 .OR. iterG .EQ. itermax) THEN
	PRINT*,'!-----------------------------------------------!'
	PRINT*,''
END IF

 END SUBROUTINE iterate 
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE compute_flux(eqn)
 
 USE grid_params
 USE global
 USE flo_vars
 USE cal_params
 IMPLICIT NONE
 
 INTEGER, INTENT(IN) :: eqn
 INTEGER, PARAMETER :: X_MOM = 1, Y_MOM = 2
 INTEGER :: i, j, is, ie, js, je
 
 REAL (KIND = 8) :: ue, uw, un, us, up
 REAL (KIND = 8) :: ve, vw, vn, vs, vp
 REAL (KIND = 8) :: dxe, dxw, dyn, dys, dxew, dyns
 
 SELECT CASE (eqn)
	CASE (X_MOM)
		is = 2; ie = nci; js = 2; je = ncj-1
		DO i = is, ie
			DO j = js, je
				!Assign required velocities at cell boundaries
				up = velx(i,j)
				!~ue = 0.5*(velx(i,j)+velx(i+1,j))
				!~uw = 0.5*(velx(i-1,j)+velx(i,j))
				!~un = 0.5*(velx(i,j)+velx(i,j+1))
				!~us = 0.5*(velx(i,j-1)+velx(i,j))
				!~vn = 0.5*(vely(i,j+1)+vely(i-1,j+1))
				!~vs = 0.5*(vely(i,j)+vely(i-1,j))
				
				ue = velx(i,j) + (xc(i,j)-ux(i,j))* &
				(velx(i+1,j)-velx(i,j))/(ux(i+1,j)-ux(i,j))
				uw = velx(i-1,j) + (xc(i-1,j)-ux(i-1,j))* &
				(velx(i,j)-velx(i-1,j))/(ux(i,j)-ux(i-1,j))
				un = velx(i,j) + (y(i,j+1)-uy(i,j))* &
				(velx(i,j+1)-velx(i,j))/(uy(i,j+1)-uy(i,j))
				us = velx(i,j-1) + (y(i,j)-uy(i,j-1))* &
				(velx(i,j)-velx(i,j-1))/(uy(i,j)-uy(i,j-1))
				vn = vely(i-1,j+1) + (x(i,j+1)-vx(i-1,j+1))* &
				(vely(i,j+1)-vely(i-1,j+1))/(vx(i,j+1)-vx(i-1,j+1))
				vs = vely(i-1,j) + (x(i,j)-vx(i-1,j))* &
				(vely(i,j)-vely(i-1,j))/(vx(i,j)-vx(i-1,j))
				
				!Assign required distances
				dxe = dxue(i,j); dxw = dxuw(i,j);
				dyn = dyun(i,j); dys = dyus(i,j);
				!~dxew = 0.5*(ux(i+1,j)+ux(i,j)) - 0.5*(ux(i,j)+ux(i-1,j))
				!~dyns = 0.5*(uy(i,j+1)+uy(i,j)) - 0.5*(uy(i,j)+uy(i,j-1))
				dxew = xc(i,j)-xc(i-1,j)
				dyns = y(i,j+1)-y(i,j)
				
				!Compute flux terms
				Fc_u(i,j) = -1.*( (ue**2/dxe - uw**2/dxw) &
								+ (un*vn/dyn - us*vs/dys))

				Fv_u(i,j) = nu*( ((velx(i+1,j)-velx(i,j))/dxe - &
								(velx(i,j)-velx(i-1,j))/dxw)/dxew + &
								((velx(i,j+1)-velx(i,j))/dyn - &
								(velx(i,j)-velx(i,j-1))/dys)/dyns )

				Fp_u(i,j) = -1./rho*(pr(i,j)-pr(i-1,j))/dxew
			END DO
		END DO
	CASE (Y_MOM)
		is = 2; ie = nci-1; js = 2; je = ncj
		DO i = is, ie
			DO j = js, je
				!Assign required velocities at cell boundaries
				vp = vely(i,j)
				!~ve = 0.5*(vely(i,j)+vely(i+1,j))
				!~vw = 0.5*(vely(i-1,j)+vely(i,j))
				!~vn = 0.5*(vely(i,j)+vely(i,j+1))
				!~vs = 0.5*(vely(i,j-1)+vely(i,j))
				!~ue = 0.5*(velx(i+1,j)+velx(i+1,j-1))
				!~uw = 0.5*(velx(i,j)+velx(i,j-1))
				
				ve = vely(i,j) + (x(i+1,j)-vx(i,j))* &
				(vely(i+1,j)-vely(i,j))/(vx(i+1,j)-vx(i,j))
				vw = vely(i-1,j) + (x(i,j)-vx(i-1,j))* &
				(vely(i,j)-vely(i-1,j))/(vx(i,j)-vx(i-1,j))
				vn = vely(i,j) + (yc(i,j)-vy(i,j))* &
				(vely(i,j+1)-vely(i,j))/(vy(i,j+1)-vy(i,j))
				vs = vely(i,j-1) + (yc(i,j-1)-vy(i,j-1))* &
				(vely(i,j)-vely(i,j-1))/(vy(i,j)-vy(i,j-1))
				ue = velx(i+1,j-1) + (y(i+1,j)-uy(i+1,j-1))* &
				(velx(i+1,j)-velx(i+1,j-1))/(uy(i+1,j)-uy(i+1,j-1))
				uw = velx(i,j-1) + (y(i,j)-uy(i,j-1))* &
				(velx(i,j)-velx(i,j-1))/(uy(i,j)-uy(i,j-1))
				
				!Assign required distances
				dxe = dxve(i,j); dxw = dxvw(i,j);
				dyn = dyvn(i,j); dys = dyvs(i,j);
				!~dxew = 0.5*(vx(i+1,j)+vx(i,j)) - 0.5*(vx(i,j)+vx(i-1,j))
				!~dyns = 0.5*(vy(i,j+1)+vy(i,j)) - 0.5*(vy(i,j)+vy(i,j-1))
				dxew = x(i+1,j)-x(i,j)
				dyns = yc(i,j)-yc(i,j-1)
				
				!Compute flux terms
				Fc_v(i,j) = -1.*( (vn**2/dyn - vs**2/dys) &
								+ (ue*ve/dxe - uw*vw/dxw))

				Fv_v(i,j) = nu*( ((vely(i+1,j)-vely(i,j))/dxe - &
								(vely(i,j)-vely(i-1,j))/dxw)/dxew + &
								((vely(i,j+1)-vely(i,j))/dyn - &
								(vely(i,j)-vely(i,j-1))/dys)/dyns )

				Fp_v(i,j) = -1./rho*(pr(i,j)-pr(i,j-1))/dyns
			END DO
		END DO
 END SELECT 
 
 END SUBROUTINE compute_flux
 !----------------------------------------------------------------------------!
  
 
 SUBROUTINE out_mon_data
 
 USE grid_params
 USE global
 USE flo_vars
 IMPLICIT NONE
 
 INTEGER :: i, j, ind, ki, kj
 
 REAL(KIND = 8) :: umon, vmon, pmon, uex, vex, pex
 
 DO ind = 1, nmon
 
	 !Interpolate vars on cell-centers
	 ki = imon(ind); kj = jmon(ind)
	 umon = 0.5*(velx(ki+1,kj) + velx(ki,kj))
	 vmon = 0.5*(vely(ki,kj+1) + vely(ki,kj))
	 pmon = pr(ki,kj)
	 !Calculate exact values of parameters
	 IF (iProblem .EQ. TaylorGreen) THEN
		 uex = -1.*exp(-2.*time)*(COS(xc(ki,kj))*SIN(yc(ki,kj)))
		 vex = exp(-2.*time)*(SIN(xc(ki,kj))*COS(yc(ki,kj)))
		 pex = -1./4.*exp(-4.*time)*(COS(2.*xc(ki,kj))+COS(2.*yc(ki,kj)))
	 END IF	 
	 WRITE(monBase+ind,71) iterG, time, umon, vmon, pmon
	 IF (iProblem .EQ. TaylorGreen) THEN
		WRITE(monBaseEx+ind,71) iterG, time, uex, vex, pex
	 END IF
 
 END DO
 
 71 FORMAT (I5,2X,F7.3,3(2X,E13.6))
 
 END SUBROUTINE out_mon_data
 !----------------------------------------------------------------------------!
  
 
 SUBROUTINE out_data
 
 USE grid_params
 USE global
 USE flo_vars
 USE GCM
 IMPLICIT NONE
 
 INTEGER, PARAMETER :: fout = 101
 INTEGER :: i, j, kx, ky
 INTEGER :: a, b, c, d
 
 CHARACTER (LEN = 50) :: fname_out
 
 WRITE(*,'(A)', ADVANCE = 'NO') ' WRITING TECPLOT-COMPATIBLE OUTPUT FILES ...'
 
 ALLOCATE(Xout(ngx,ngy), Yout(ngx,ngy)) 
 ALLOCATE(Uout(nci-2,ncj-2), Vout(nci-2,ncj-2), Pout(nci-2,ncj-2), &
 Vor(nci-2,ncj-2), Phi_ex(nci-2,ncj-2))
 IF (iIntlBody .EQ. 1) ALLOCATE(iBlankOut(nci-2,ncj-2))
 
 CALL compute_out_vars
 
 !Create output file name
 WRITE(fname_out,'("NS_out.",I6.6,".plt")') iterG
 
 OPEN(fout, FILE = TRIM(fname_out))
 WRITE(fout,'(A)')' TITLE = "NAVIER-STOKES OUTPUT"'
 IF (iProblem .EQ. TaylorGreen) THEN
	 IF (Var_ex .EQ. 1) THEN
		 WRITE(fout,'(A)')' VARIABLES = "X", "Y", "U", "Uex", "V", "P", "Vor"'
		 WRITE(fout,'(A,I4,A,I4,A)')' ZONE I = ',ngx,', J = ', ngy, &
									', DATAPACKING = BLOCK, VARLOCATION=(3=CELLCENTERED, &
										   4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED, &
										   7=CELLCENTERED)'

		 WRITE(fout,'(10(E15.7,2X))') Xout, Yout, Uout, Phi_ex, Vout, Pout, Vor	
	 ELSEIF (Var_ex .EQ. 2) THEN
		 WRITE(fout,'(A)')' VARIABLES = "X", "Y", "U", "V", "Vex", "P", "Vor"'
		 WRITE(fout,'(A,I4,A,I4,A)')' ZONE I = ',ngx,', J = ', ngy, &
									', DATAPACKING = BLOCK, VARLOCATION=(3=CELLCENTERED, &
										   4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED, &
										   7=CELLCENTERED)'

		 WRITE(fout,'(10(E15.7,2X))') Xout, Yout, Uout, Vout, Phi_ex, Pout, Vor	 
	 ELSEIF (Var_ex .EQ. 3) THEN
		 WRITE(fout,'(A)')' VARIABLES = "X", "Y", "U", "V", "P", "Pex", "Vor"'
		 WRITE(fout,'(A,I4,A,I4,A)')' ZONE I = ',ngx,', J = ', ngy, &
									', DATAPACKING = BLOCK, VARLOCATION=(3=CELLCENTERED, &
										   4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED, &
										   7=CELLCENTERED)'

		 WRITE(fout,'(10(E15.7,2X))') Xout, Yout, Uout, Vout, Pout, Phi_ex, Vor	
	 END IF
 ELSEIF (iProblem .EQ. DrivenCavity) THEN	 
	 WRITE(fout,'(A)')' VARIABLES = "X", "Y", "U", "V", "P", "Vor"'
	 WRITE(fout,'(A,I4,A,I4,A)')' ZONE I = ',ngx,', J = ', ngy, &
									', DATAPACKING = BLOCK, VARLOCATION=(3=CELLCENTERED, &
										   4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED)'
										   
	 WRITE(fout,'(10(E15.7,2X))') Xout, Yout, Uout, Vout, Pout, Vor
 ELSEIF (iProblem .EQ. Cylinder) THEN
	 WRITE(fout,'(A)')' VARIABLES = "X", "Y", "U", "V", "P", "Vor", "BL"'
	 WRITE(fout,'(A,I4,A,I4,A)')' ZONE I = ',ngx,', J = ', ngy, &
									', DATAPACKING = BLOCK, VARLOCATION=(3=CELLCENTERED, &
										   4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED, &
										   7=CELLCENTERED)'
	 WRITE(fout,'(10(E15.7,2X))') Xout, Yout, Uout, Vout, Pout, Vor, REAL(iBlankOut)
 END IF	 
 CLOSE(fout)
 
 WRITE(*,'(A)') ' DONE!'
 PRINT*,''
 
 DEALLOCATE (Xout, Yout, Uout, Vout, Pout, Vor, Phi_ex)
 IF (iIntlBody .EQ. 1)DEALLOCATE(iBlankOut)
 
 END SUBROUTINE out_data
 !----------------------------------------------------------------------------!
  
 
 SUBROUTINE out_marker_data
 
 USE grid_params
 USE global
 USE flo_vars
 USE GCM
 IMPLICIT NONE
 
 CHARACTER (LEN = 50) :: fname_out
 
 REAL(KIND = 8) :: scalefactor
 
 INTEGER, PARAMETER :: fout = 102
 INTEGER :: i, j
 
 WRITE(fname_out,'("marker_out.",I6.6,".dat")') iterG
 OPEN(fout, FILE = TRIM(fname_out))
 WRITE(fout,'(A)') 'VARIABLES = "X","Y","P"'
 
 scalefactor = MAXVAL(ABS(BCVal))
 
 DO i = 1, nElemMarker
	WRITE(fout,'(3(E14.7,2X))')xBodyMarker(i), yBodyMarker(i), &
	 pBodyMarker(i)/(0.5*rho*scalefactor**2)
 END DO
 
 CLOSE(fout)
 
 END SUBROUTINE out_marker_data
 !----------------------------------------------------------------------------!
  
 
 SUBROUTINE out_restart_data
 
 USE global
 USE grid_params, ONLY:time
 USE flo_vars
 IMPLICIT NONE
 
 INTEGER, PARAMETER :: fout = 101
 
 CHARACTER*50 :: fname_res
 
 WRITE(*,'(A)', ADVANCE = 'NO') ' WRITING RESTART FILES ...'
 
 WRITE(fname_res,'(A,I1,A)') 'restart.',iToggle+1,'.dat'
 
 !Ifort
 !OPEN(fout, FILE = TRIM(fname_res), FORM = 'BINARY')
 !Gfortran
 OPEN(fout, FILE = TRIM(fname_res), ACCESS = 'STREAM')
 
 WRITE(fout)iterG
 WRITE(fout)time
 WRITE(fout)velx
 WRITE(fout)vely
 WRITE(fout)pr
 
 CLOSE(fout)
 
 IF (iToggle .EQ. 0) THEN
	iToggle = 1
 ELSE
	iToggle = 0
 END IF
 
 WRITE(*,'(A)') ' DONE!'
 PRINT*,''
 
 END SUBROUTINE out_restart_data
  !----------------------------------------------------------------------------!
  
 
 SUBROUTINE read_restart_data
 
 USE global
 USE grid_params, ONLY:time
 USE flo_vars
 USE GCM, ONLY: iIntlBody
 IMPLICIT NONE
 
 INTEGER, PARAMETER :: fin = 101
 
 CHARACTER*50 :: fname_res
 
 WRITE(*,'(A)', ADVANCE = 'NO') ' READING RESTART DATA ...' 
 
 WRITE(fname_res,'(A,I1,A)')'restart.', nread, '.dat'
 !Ifort
 !OPEN(fin, FILE = TRIM(fname_res), FORM = 'BINARY', ACTION = 'READ')
 !Gfortran
 OPEN(fin, FILE = TRIM(fname_res), ACCESS = 'STREAM', ACTION = 'READ')
 
 READ(fin)iterG
 READ(fin)time
 READ(fin)velx
 READ(fin)vely
 READ(fin)pr
 
 CLOSE(fin)
 
 WRITE(*,'(A)') ' DONE!'
 PRINT*,''
 
 IF (iIntlBody .EQ. 1) CALL create_GCM
 
 END SUBROUTINE read_restart_data
 !----------------------------------------------------------------------------!
  
 
 SUBROUTINE compute_out_vars
 
 USE grid_params
 USE global
 USE flo_vars
 USE GCM
 IMPLICIT NONE
 
 REAL (KIND = 8) :: dist1, dist2, dist
 REAL (KIND = 8) :: scalefactor
 REAL (KIND = 8) :: vne, vse, vnw, vsw, uen, uwn, ues, uws
 
 INTEGER :: i, j, ki, kj
 
 DO j = 1, ngy
	DO i = 1, ngx
		Xout(i,j) = x(i+1,j+1); Yout(i,j) = y(i+1,j+1)
	END DO
 END DO 
 
 !IF (iIntlBody .EQ. 1) gmapOut = gmap(2:nci-1,2:ncj-1)
 IF (iIntlBody .EQ. 1) iBlankOut = iBlank(2:nci-1,2:ncj-1)
 
 !Transfer all vars to cell-centers
 DO i = 2, nci-1
	DO j = 2, ncj-1
		ki = i-1; kj = j-1
		!Compute primary vars
		Uout(ki,kj) = 0.5*(velx(i+1,j)+velx(i,j))
		Vout(ki,kj) = 0.5*(vely(i,j+1)+vely(i,j))
		Pout(ki,kj) = pr(i,j)
		!Calc vorticity
		vne = vely(i,j+1) + (x(i+1,j+1)-vx(i,j+1))* &
		(vely(i+1,j+1)-vely(i,j+1))/(vx(i+1,j+1)-vx(i,j+1))
		vse = vely(i,j) + (x(i+1,j)-vx(i,j))* &
		(vely(i+1,j)-vely(i,j))/(vx(i+1,j)-vx(i,j))
		vnw = vely(i-1,j+1) + (x(i,j+1)-vx(i-1,j+1))* &
		(vely(i,j+1)-vely(i-1,j+1))/(vx(i,j+1)-vx(i-1,j+1))
		vsw = vely(i-1,j) + (x(i,j)-vx(i-1,j))* &
		(vely(i,j)-vely(i-1,j))/(vx(i,j)-vx(i-1,j))
		
		uen = velx(i+1,j) + (y(i+1,j+1)-uy(i+1,j))* &
		(velx(i+1,j+1)-velx(i+1,j))/(uy(i+1,j+1)-uy(i+1,j))
		ues = velx(i+1,j-1) + (y(i+1,j)-uy(i+1,j-1))* &
		(velx(i+1,j)-velx(i+1,j-1))/(uy(i+1,j)-uy(i+1,j-1))
		uwn = velx(i,j) + (y(i,j+1)-uy(i,j))* &
		(velx(i,j+1)-velx(i,j))/(uy(i,j+1)-uy(i,j))
		uws = velx(i,j-1) + (y(i,j)-uy(i,j-1))* &
		(velx(i,j)-velx(i,j-1))/(uy(i,j)-uy(i,j-1))
		
		Vor(ki,kj) = 0.5*(((vne+vse)-(vnw+vsw))/(ux(i+1,j)-ux(i,j)) - &
						  ((uen+uwn)-(ues+uws))/(vy(i,j+1)-vy(i,j)))
		
		!~Vor(ki,kj) = 0.5*( (vely(i+1,j+1)+vely(i+1,j)) -(vely(i-1,j+1)+vely(i-1,j)) )/&
		!~	(0.5*( (vx(i+1,j+1)+vx(i+1,j)) -(vx(i-1,j+1)+vx(i-1,j)) )) - &
		!~	0.5*( (velx(i+1,j+1)+velx(i,j+1)) -(velx(i+1,j-1)+velx(i,j-1)) )/&
		!~	(0.5*( (uy(i+1,j+1)+uy(i,j+1)) -(uy(i+1,j-1)+uy(i,j-1)) ))
		IF (Var_ex .EQ. 1) THEN
			Phi_ex(ki,kj) = -1.*exp(-2.*time)*(COS(xc(i,j))*SIN(yc(i,j)))
		ELSEIF (Var_ex .EQ. 2) THEN
			Phi_ex(ki,kj) = exp(-2.*time)*(SIN(xc(i,j))*COS(yc(i,j)))
		ELSEIF (Var_ex .EQ. 3) THEN
			Phi_ex(ki,kj) = -1./4.*exp(-4.*time)*(COS(2.*xc(i,j))+COS(2.*yc(i,j)))
		END IF
		
		IF (iIntlBody .EQ. 1) THEN
			Uout(ki,kj) = REAL(1-iBlankOut(ki,kj))*Uout(ki,kj)
			Vout(ki,kj) = REAL(1-iBlankOut(ki,kj))*Vout(ki,kj)
			Pout(ki,kj) = REAL(1-iBlankOut(ki,kj))*Pout(ki,kj)
			Vor(ki,kj) = REAL(1-iBlankOut(ki,kj))*Vor(ki,kj)
		END IF
	END DO
 END DO
 
 !Normalize the quantities if necessary
 IF (iProblem .EQ. DrivenCavity) THEN
	scalefactor = MAXVAL(ABS(BCVal))
	Xout = Xout/(Xmax-Xmin)
	Yout = Yout/(Ymax-Ymin)
	Uout = Uout/scalefactor
	Vout = Vout/scalefactor
 END IF
 
 IF (iProblem .EQ. Cylinder) THEN
	scalefactor = MAXVAL(ABS(BCVal))
	Uout = Uout/scalefactor
	Vout = Vout/scalefactor
	Pout = Pout/(0.5*rho*scalefactor**2)
	Vor = Vor/scalefactor
 END IF
 
 END SUBROUTINE compute_out_vars
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE close_mon_files
 
 USE global, ONLY: monBase, nmon, iCDCL, CDCLbase
 IMPLICIT NONE
 
 INTEGER :: ind
 
 DO ind = 1, nmon
	CLOSE(monBase+ind)
 END DO
 
 IF (iCDCL .EQ. 1) CLOSE(CDCLbase)
 
 END SUBROUTINE close_mon_files
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE free_flo_vars
 
 USE flo_vars
 IMPLICIT NONE
 
 DEALLOCATE(velx, vely, velxH, velyH, pr)
 DEALLOCATE(Fc_u, Fc_v, Fv_u, Fv_v, Fp_u, Fp_v)
 DEALLOCATE(Gi_u, Gi_v)
 
 END SUBROUTINE free_flo_vars
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE free_grid_params
 
 USE grid_params
 IMPLICIT NONE
 
 DEALLOCATE(xgm, ygm, px, py, ux, uy, vx, vy)
 DEALLOCATE(x, y, xc, yc)
 DEALLOCATE(dxue, dxuw, dyun, dyus)
 DEALLOCATE(dxve, dxvw, dyvn, dyvs)
 
 END SUBROUTINE free_grid_params
 !----------------------------------------------------------------------------!
 
 
 SUBROUTINE disp_simulation_message
 
 USE grid_params
 USE global
 USE flo_vars
 IMPLICIT NONE
 
 PRINT*,'!------------------------------------------------------!'
 PRINT*,'               2D NAVIER-STOKES SOLVER                  '
 PRINT*,'!------------------------------------------------------!'
 PRINT*,''
 
 IF (iRestart .EQ. 0) THEN
	WRITE(*,'(A)', ADVANCE = 'NO')'  NEW SIMULATION RUN FOR '
 ELSE
	WRITE(*,'(A)', ADVANCE = 'NO')'  CALCULATION RESTART FOR '
 END IF
 
 IF (iProblem .EQ. TaylorGreen) THEN
	WRITE(*,'(A)')'TAYLOR-GREEN VORTEX PROBLEM'
 ELSEIF (iProblem .EQ. DrivenCavity) THEN
	WRITE(*,'(A)')'SHEAR-DRIVEN CAVITY FLOW PROBLEM'
	Re_num = MAXVAL(ABS(BCVal))*MAX(Xmax-Xmin,Ymax,Ymin)/nu
	WRITE(*,'(A,F7.2)')' REYNOLDS NUMBER: ', Re_Num
 ELSEIF (iProblem .EQ. Cylinder) THEN
	WRITE(*,'(A)')'FLOW PAST CYLINDER PROBLEM'
 END IF
 
 PRINT*,'!------------------------------------------------------!'
 PRINT*,''
 
 END SUBROUTINE disp_simulation_message
 !----------------------------------------------------------------------------!
 
 
  
