 !Module Section
 MODULE grid_params
 
 	REAL (KIND = 8) :: Xmin, Xmax, Ymin, Ymax, XminG, XmaxG, YminG, YmaxG
	REAL (KIND = 8) :: dx, dy
	REAL (KIND = 8) :: Tmin, Tmax, time, dt, dt2
	REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: xgm, ygm, px, py, ux, uy, vx, vy
	REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: x, y, xc, yc
	REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: dxue, dxuw, dyun, dyus
	REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: dxve, dxvw, dyvn, dyvs
	
	INTEGER :: gen_opt
	INTEGER :: ngx, ngy, ngi, ngj, ncx, ncy, nci, ncj
	INTEGER :: nodes, elem, ncc, nintcc
 
 END MODULE grid_params
 !----------------------------------------------------------------------------!
 
 
 MODULE global
 
	INTEGER :: iProblem, Sol_Method, Mode, BoundCond(4)
	INTEGER :: time_scheme

	INTEGER, PARAMETER ::  WALL     = 1, &
						   INLET	= 2, &
						   OUTLET	= 3, &
						   PERIODIC = 4, &
						   CONVECT  = 5, &
						   SYMM 	= 6
	
	INTEGER, PARAMETER :: TaylorGreen  = 1, &
						  DrivenCavity = 2, &
						  Cylinder	   = 3
						   
	INTEGER, PARAMETER ::  RK3 		= 1, &
						   AdamBash = 2
	INTEGER :: iter, iterG, itermax, ntec
	
	!Restart parameters
	INTEGER :: iRestart, nRead, itoggle
	
	!Other parameters
	INTEGER :: iCDCL, iVerbose
	
	REAL(KIND = 8), PARAMETER :: pi = 3.14159265359
	REAL(KIND = 8), PARAMETER :: nu = 1., rho = 1.
	REAL(KIND = 8) :: Re_num
	REAL(KIND = 8) :: BCVal(4)
	REAL(KIND = 8) :: freq, omega
	
	!MasterLog parameter
	INTEGER, PARAMETER :: fmasterlog = 801
	
	!Probe parameters
	INTEGER :: nmon
	INTEGER, PARAMETER :: monBase = 700, monBaseEx = 800, CDCLbase = 900
	INTEGER, ALLOCATABLE, DIMENSION(:) :: imon, jmon
	REAL, ALLOCATABLE, DIMENSION(:) :: xmon, ymon, zmon

 END MODULE global
 !----------------------------------------------------------------------------!
 
 
 MODULE flo_vars
 
	REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: velx, vely, pr
	REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: velxH, velyH
	
	!Variables for fluxes
	REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: Fc_u, Fc_v
	REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: Fv_u, Fv_v
	REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: Fp_u, Fp_v
	REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: Gi_u, Gi_v
	
	!Output variables
	REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: Xout, Yout, Uout, &
													Vout, Pout, Vor, &
													Phi_ex

	REAL(KIND = 8) :: uinf, vinf
	
	INTEGER :: Var_ex
	
 END MODULE flo_vars
  !----------------------------------------------------------------------------!
 
 
 MODULE cal_params
	
	!Parameters/ Vars associated with Runge-Kutta (RK3) time integration
	REAL(KIND = 8), PARAMETER  :: dtFactor(3) = (/ 1./3., 15./16., 8./15. /)
	REAL(KIND = 8), PARAMETER  :: GiFactor(3) = (/ 0., -5./9., -153./128. /)
	REAL(KIND = 8), PARAMETER  :: dtFactorCor(3) = (/ 1./3., 5./12., 1./4. /)
	
	!Computation time parameters
	REAL(KIND = 8) :: tstart, tend, sec
	INTEGER :: minute
	
 END MODULE cal_params 
 !----------------------------------------------------------------------------!
 
 
  MODULE Poisson_cal_params
 
 INTEGER, PARAMETER ::  JACOBI  = 1, &
 					    SOR     = 2, &
					    ConGrad = 3
 
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: U, Uold, Uex, f
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: aP, aE, aW, aN, aS
 REAL(KIND = 8) :: residG, residMon, resid0
 
 !Data structures related to Conjugate Gradient method
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: A_CG, b_CG, resid_CG, &
												 phi_CG, phiO_CG, Ar_CG

 REAL(KIND = 8) :: alpha
 
 LOGICAL :: converged
 
 END MODULE Poisson_cal_params
 !----------------------------------------------------------------------------!
 
 
 MODULE GCM
 
 INTEGER :: iIntlBody
 INTEGER :: nBodyMarker, nElemMarker
 INTEGER, PARAMETER :: nvarsGC = 3
 
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:) :: xBodyMarker, yBodyMarker, pBodyMarker
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:) :: AreaBodyMarker
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:) :: xElemCent, yElemCent
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:) :: normElX, normElY
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:) :: xBI, yBI, prBI, xIP, yIP, prIP
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:) :: xBIu, yBIu, xIPu, yIPu, uIP
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:) :: xBIv, yBIv, xIPv, yIPv, vIP
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:) :: coeffsnb
 REAL(KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: LSFMarker, LSFIP, LSFIPu, LSFIPv
 
 INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ElemMarkerConnex
 INTEGER, ALLOCATABLE, DIMENSION (:,:) :: gmap, gmapU, gmapV, gmapOut 
 INTEGER, ALLOCATABLE, DIMENSION (:,:) :: GCp, GCu, GCv
 INTEGER, ALLOCATABLE, DIMENSION (:,:) :: iBlank, iBlanku, iBlankv, iBlankOut
 INTEGER, ALLOCATABLE, DIMENSION (:,:) :: nb4v1, nb4v1u, nb4v1v
 INTEGER, ALLOCATABLE, DIMENSION (:,:) :: nbMarker
 INTEGER :: nGC, nGCu, nGCv
 
 END MODULE GCM 
 !----------------------------------------------------------------------------!  
