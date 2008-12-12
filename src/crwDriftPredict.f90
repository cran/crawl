SUBROUTINE crwdrift_predict( tau2y, tau2x, sig2, b, bd, sig2d, delta, x, y, loctype, &
                             stay, ay, ax, Py, Px, lonadj, N, lly, llx, &
                             predy, predx, vary, varx)
 !! LOADED VALUES/ARRAYS !!
  INTEGER           N
  INTEGER           loctype(N), stay(N)
  DOUBLE PRECISION  lly, llx, tau2y(N), tau2x(N)
  DOUBLE PRECISION  sig2(N), b(N), bd(N), sig2d, delta(N), x(N), y(N), lonadj(N)
  DOUBLE PRECISION  ax(3,1), ay(3,1), Py(3,3), Px(3,3)
  DOUBLE PRECISION  predy(N,3), predx(N,3), vary(3,3,N), varx(3,3,N)

 !! DERIVED ARRAYS !!
  DOUBLE PRECISION vy(N), vx(N), Fy(N), Fx(N)
  DOUBLE PRECISION Z(1,3)
  DOUBLE PRECISION Ky(3,1), Kx(3,1)
  DOUBLE PRECISION Qy(3,3), Qx(3,3), T(3,3)
  DOUBLE PRECISION PyArr(3,3,N+1), PxArr(3,3,N+1), LyArr(3,3,N), LxArr(3,3,N)
  DOUBLE PRECISION ayArr(3,1,N+1), axArr(3,1,N+1)
  DOUBLE PRECISION ry(3,1), Ny(3,3), rx(3,1), Nx(3,3)

 !! INTIAL CONDITIONS !!
  ayArr(:,:,1) = ay
  axArr(:,:,1) = ax
  PyArr(:,:,1) = Py
  PxArr(:,:,1) = Px
  Z            = RESHAPE((/1.0, 0.0, 0.0/),(/1, 3/))
  Qy           = 0.0
  Qx           = 0.0
  T            = 0.0
  T(1,1)       = 1.0
  vy           = 0.0
  vx           = 0.0
  Fy           = tau2y
  Fx           = tau2x
  ry           = 0.0
  rx           = 0.0
  Ny           = 0.0
  Nx           = 0.0

 !! BEGIN FILTER LOOP !!
  DO i=1,N

   !! GENERATE Q AND T MATRICES !!
    IF(stay(i)==1) THEN
      Qy = 0.0
      Qx = 0.0
      T(1,2) = 0.0
      T(2,2) = 0.0
      T(1,3) = 0.0
      T(3,3) = 0.0
    ELSE
      Qy(1,1) = (sig2(i)/(b(i)*b(i)))*(delta(i) - 2*(1-exp(-b(i)*delta(i)))/b(i)  + (1-exp(-2*b(i)*delta(i)))/(2*b(i))) &
                  + (sig2d/(bd(i)*bd(i)))*(delta(i) - 2*(1-exp(-bd(i)*delta(i)))/bd(i) + (1-exp(-2*bd(i)*delta(i)))/(2*bd(i)))
      Qy(2,1) = (sig2(i)/(b(i)*b(i)))*((1-2*exp(-b(i)*delta(i))+exp(-2*b(i)*delta(i)))/2)
      Qy(3,1) = (sig2d/(bd(i)*bd(i)))*((1-2*exp(-bd(i)*delta(i))+exp(-2*bd(i)*delta(i)))/2)
      Qy(1,2) = Qy(2,1)
      Qy(2,2) = (sig2(i)/b(i))*((1-exp(-2*b(i)*delta(i)))/2)
      Qy(1,3) = Qy(3,1)
      Qy(3,3) = (sig2d/bd(i))*((1-exp(-2*bd(i)*delta(i)))/2)
      Qx = Qy/(lonadj(i)*lonadj(i))
      T(1,2) = (1-exp(-b(i)*delta(i)))/b(i)
      T(2,2) = exp(-b(i)*delta(i))
      T(1,3) = (1 - exp(-bd(i)*delta(i)))/bd(i)
      T(3,3) = exp(-bd(i)*delta(i))
    END IF

   !! GENERAL KF FILTER !!
    Fy(i) = PyArr(1,1,i) + tau2y(i)
    Fx(i) = PxArr(1,1,i) + tau2x(i)/(lonadj(i)*lonadj(i))
    IF(loctype(i)==1 .OR. Fy(i)==0.0) THEN
      ayArr(:,:,i+1) = MATMUL(T,ayArr(:,:,i))
      PyArr(:,:,i+1) = MATMUL(MATMUL(T,PyArr(:,:,i)),TRANSPOSE(T)) + Qy
      LyArr(:,:,i) = T
    ELSE
      vy(i) = y(i)-ayArr(1,1,i)
      lly = lly - (log(Fy(i)) + vy(i)*vy(i)/Fy(i))/2
      Ky = MATMUL(MATMUL(T,PyArr(:,:,i)),TRANSPOSE(Z))/Fy(i)
      LyArr(:,:,i) = T - MATMUL(Ky,Z)
      ayArr(:,:,i+1) = MATMUL(T,ayArr(:,:,i)) + Ky*vy(i)
      PyArr(:,:,i+1) = MATMUL(MATMUL(T,PyArr(:,:,i)),TRANSPOSE(LyArr(:,:,i))) + Qy
    END IF
    IF(loctype(i)==1 .OR. Fx(i)==0.0) THEN          
      axArr(:,:,i+1) = MATMUL(T,axArr(:,:,i))      
      PxArr(:,:,i+1) = MATMUL(MATMUL(T,PxArr(:,:,i)),TRANSPOSE(T)) + Qx      
      LxArr(:,:,i) = T
    ELSE
      vx(i) = x(i)-axArr(1,1,i)
      llx = llx - (log(Fx(i)) + vx(i)*vx(i)/Fx(i))/2     
      Kx = MATMUL(MATMUL(T,PxArr(:,:,i)),TRANSPOSE(Z))/Fx(i)      
      LxArr(:,:,i) = T - MATMUL(Kx,Z)
      axArr(:,:,i+1) = MATMUL(T,axArr(:,:,i)) + Kx*vx(i)      
      PxArr(:,:,i+1) = MATMUL(MATMUL(T,PxArr(:,:,i)),TRANSPOSE(LxArr(:,:,i))) + Qx
    END IF
  END DO

 !! BEGIN SMOOTHING LOOP!!
  DO j=N,1,-1
    IF(loctype(j)==1 .OR. Fy(j)==0.0) THEN
      ry = MATMUL(TRANSPOSE(LyArr(:,:,j)),ry)
      Ny = MATMUL(MATMUL(TRANSPOSE(LyArr(:,:,j)), Ny), LyArr(:,:,j))
    ELSE
      ry = TRANSPOSE(Z)*vy(j)/Fy(j) + MATMUL(TRANSPOSE(LyArr(:,:,j)),ry)
      Ny = MATMUL(TRANSPOSE(Z),Z)/Fy(j) + MATMUL(MATMUL(TRANSPOSE(LyArr(:,:,j)), Ny), LyArr(:,:,j))
    END IF
    IF(loctype(j)==1 .OR. Fx(j)==0.0) THEN
      rx = MATMUL(TRANSPOSE(LxArr(:,:,j)),rx)
      Nx = MATMUL(MATMUL(TRANSPOSE(LxArr(:,:,j)), Nx), LxArr(:,:,j))
    ELSE      
      rx = TRANSPOSE(Z)*vx(j)/Fx(j) + MATMUL(TRANSPOSE(LxArr(:,:,j)),rx)      
      Nx = MATMUL(TRANSPOSE(Z),Z)/Fx(j) + MATMUL(MATMUL(TRANSPOSE(LxArr(:,:,j)), Nx), LxArr(:,:,j))
    END IF
    predy(j,:) = RESHAPE(ayArr(:,:,j) + MATMUL(PyArr(:,:,j), ry), (/3/))
    predx(j,:) = RESHAPE(axArr(:,:,j) + MATMUL(PxArr(:,:,j), rx), (/3/))
    vary(:,:,j) = PyArr(:,:,j) - MATMUL(MATMUL(PyArr(:,:,j),Ny), PyArr(:,:,j))
    varx(:,:,j) = PxArr(:,:,j) - MATMUL(MATMUL(PxArr(:,:,j),Nx), PxArr(:,:,j))
  END DO
END SUBROUTINE crwdrift_predict
