SUBROUTINE crw_predict( tau2y, tau2x, sig2, b, bg, sig2g, delta, x, y, loctype, &
                             stay, ay, ax, Py, Px, lonadj, N, lly, llx, &
                             predy, predx, vary, varx)
 !! LOADED VALUES/ARRAYS !!
  INTEGER           N
  INTEGER           loctype(N), stay(N)
  DOUBLE PRECISION  lly, llx, tau2y(N), tau2x(N)
  DOUBLE PRECISION  sig2(N), b(N), bg(N), sig2g(N), delta(N), x(N), y(N), lonadj(N)
  DOUBLE PRECISION  ax(2,1), ay(2,1), Py(2,2), Px(2,2)
  DOUBLE PRECISION  predy(N,2), predx(N,2), vary(2,2,N), varx(2,2,N)

 !! DERIVED ARRAYS !!
  DOUBLE PRECISION vy(N), vx(N), Fy(N), Fx(N)
  DOUBLE PRECISION Z(1,2)
  DOUBLE PRECISION Ky(2,1), Kx(2,1)
  DOUBLE PRECISION Qy(2,2), Qx(2,2), T(2,2)
  DOUBLE PRECISION PyArr(2,2,N+1), PxArr(2,2,N+1), LyArr(2,2,N), LxArr(2,2,N)
  DOUBLE PRECISION ayArr(2,1,N+1), axArr(2,1,N+1)
  DOUBLE PRECISION ry(2,1), Ny(2,2), rx(2,1), Nx(2,2)

 !! INTIAL CONDITIONS !!
  ayArr(:,:,1) = ay
  axArr(:,:,1) = ax
  PyArr(:,:,1) = Py
  PxArr(:,:,1) = Px
  Z            = RESHAPE((/1.0, 0.0/),(/1, 2/))
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
    ELSE
      Qy(1,1) = (sig2(i)/(b(i)*b(i)))*(delta(i) - 2*(1-exp(-b(i)*delta(i)))/b(i)  + (1-exp(-2*b(i)*delta(i)))/(2*b(i)))
      Qy(2,1) = (sig2(i)/(b(i)*b(i)))*((1-2*exp(-b(i)*delta(i))+exp(-2*b(i)*delta(i)))/2)
      Qy(1,2) = Qy(2,1)
      Qy(2,2) = (sig2(i)/b(i))*((1-exp(-2*b(i)*delta(i)))/2)
      Qx = Qy/(lonadj(i)*lonadj(i))
      T(1,2) = (1-exp(-b(i)*delta(i)))/b(i)
      T(2,2) = exp(-b(i)*delta(i))
    END IF

   !! GENERAL KF FILTER !!
    IF(loctype(i)==1) THEN
      ayArr(:,:,i+1) = MATMUL(T,ayArr(:,:,i))
      axArr(:,:,i+1) = MATMUL(T,axArr(:,:,i))
      PyArr(:,:,i+1) = MATMUL(MATMUL(T,PyArr(:,:,i)),TRANSPOSE(T)) + Qy
      PxArr(:,:,i+1) = MATMUL(MATMUL(T,PxArr(:,:,i)),TRANSPOSE(T)) + Qx
      LyArr(:,:,i) = T
      LxArr(:,:,i) = T
    ELSE
      vy(i) = y(i)-ayArr(1,1,i)
      vx(i) = x(i)-axArr(1,1,i)
      Fy(i) = PyArr(1,1,i) + tau2y(i)
      Fx(i) = PxArr(1,1,i) + tau2x(i)/(lonadj(i)*lonadj(i))
      Ky = MATMUL(MATMUL(T,PyArr(:,:,i)),TRANSPOSE(Z))/Fy(i)
      Kx = MATMUL(MATMUL(T,PxArr(:,:,i)),TRANSPOSE(Z))/Fx(i)
      LyArr(:,:,i) = T - MATMUL(Ky,Z)
      LxArr(:,:,i) = T - MATMUL(Kx,Z)
      ayArr(:,:,i+1) = MATMUL(T,ayArr(:,:,i)) + Ky*vy(i)
      axArr(:,:,i+1) = MATMUL(T,axArr(:,:,i)) + Kx*vx(i)
      PyArr(:,:,i+1) = MATMUL(MATMUL(T,PyArr(:,:,i)),TRANSPOSE(LyArr(:,:,i))) + Qy
      PxArr(:,:,i+1) = MATMUL(MATMUL(T,PxArr(:,:,i)),TRANSPOSE(LxArr(:,:,i))) + Qx
      lly = lly - (log(Fy(i)) + vy(i)*vy(i)/Fy(i))/2
      llx = llx - (log(Fx(i)) + vx(i)*vx(i)/Fx(i))/2
    END IF
  END DO

 !! BEGIN SMOOTHING LOOP!!
  DO j=N,1,-1
    IF(loctype(j)==1) THEN
      ry = MATMUL(TRANSPOSE(LyArr(:,:,j)),ry)
      rx = MATMUL(TRANSPOSE(LxArr(:,:,j)),rx)
      Ny = MATMUL(MATMUL(TRANSPOSE(LyArr(:,:,j)), Ny), LyArr(:,:,j))
      Nx = MATMUL(MATMUL(TRANSPOSE(LxArr(:,:,j)), Nx), LxArr(:,:,j))
    ELSE
      ry = TRANSPOSE(Z)*vy(j)/Fy(j) + MATMUL(TRANSPOSE(LyArr(:,:,j)),ry)
      rx = TRANSPOSE(Z)*vx(j)/Fx(j) + MATMUL(TRANSPOSE(LxArr(:,:,j)),rx)
      Ny = MATMUL(TRANSPOSE(Z),Z)/Fy(j) + MATMUL(MATMUL(TRANSPOSE(LyArr(:,:,j)), Ny), LyArr(:,:,j))
      Nx = MATMUL(TRANSPOSE(Z),Z)/Fx(j) + MATMUL(MATMUL(TRANSPOSE(LxArr(:,:,j)), Nx), LxArr(:,:,j))
    END IF
    predy(j,:) = RESHAPE(ayArr(:,:,j) + MATMUL(PyArr(:,:,j), ry), (/2/))
    predx(j,:) = RESHAPE(axArr(:,:,j) + MATMUL(PxArr(:,:,j), rx), (/2/))
    vary(:,:,j) = PyArr(:,:,j) - MATMUL(MATMUL(PyArr(:,:,j),Ny), PyArr(:,:,j))
    varx(:,:,j) = PxArr(:,:,j) - MATMUL(MATMUL(PxArr(:,:,j),Nx), PxArr(:,:,j))
  END DO
END SUBROUTINE crw_predict
