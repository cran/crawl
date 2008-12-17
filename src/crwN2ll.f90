SUBROUTINE crwn2ll(tau2y, tau2x, sig2, b, bd, sig2d, delta, x, y, loctype, stay, &
                    ay, ax, Py, Px, lonadj, N, lly, llx)
  INTEGER N
  INTEGER loctype(N), stay(N)
  DOUBLE PRECISION lly, llx, tau2y(N), tau2x(N)
  DOUBLE PRECISION sig2(N), b(N), bd, sig2d, delta(N), x(N), y(N), lonadj(N)
  DOUBLE PRECISION ay(2,1), ax(2,1), Py(2,2), Px(2,2)
  DOUBLE PRECISION vy, vx, Fy, Fx
  DOUBLE PRECISION Z(1,2)
  DOUBLE PRECISION Ky(2,1), Kx(2,1)
  DOUBLE PRECISION Qy(2,2), Qx(2,2), T(2,2), Lx(2,2), Ly(2,2)

  !! INITIAL VALUES !!
  Z = RESHAPE((/1.0, 0.0/),(/1, 2/))
  T = 0.0
  T(1,1) = 1.0
  Qy = 0.0
  Qx = 0.0

 !! BEGIN FILTER LOOP !!
  DO i=1,N

   !! GENERATE Q AND T MATRICES !!
    IF(stay(i)==1) THEN
      Qy = 0.0
      Qx = 0.0
      T(1,2) = 0.0
      T(2,2) = 0.0
    ELSE
      Qy(1,1) = (sig2(i)/(b(i)*b(i)))*(delta(i) - 2*(1-exp(-b(i)*delta(i)))/b(i) + (1-exp(-2*b(i)*delta(i)))/(2*b(i)))
      Qy(2,1) = (sig2(i)/(b(i)*b(i)))*(1-2*exp(-b(i)*delta(i))+exp(-2*b(i)*delta(i)))/2
      Qy(1,2) = Qy(2,1)
      Qy(2,2) = (sig2(i)/b(i))*(1-exp(-2*b(i)*delta(i)))/2
      Qx = Qy/(lonadj(i)*lonadj(i))
      T(1,2) = (1-exp(-b(i)*delta(i)))/b(i)
      T(2,2) = exp(-b(i)*delta(i))
    END IF

   !! GENERAL KF LOOP !!
    IF(loctype(i)==1) THEN
      ay = MATMUL(T,ay)
      ax = MATMUL(T,ax)
      Py = MATMUL(MATMUL(T,Py),TRANSPOSE(T)) + Qy
      Px = MATMUL(MATMUL(T,Px),TRANSPOSE(T)) + Qx
    ELSE
      vy = y(i)-ay(1,1)
      vx = x(i)-ax(1,1)
      Fy = Py(1,1) + tau2y(i)
      Fx = Px(1,1) + tau2x(i)/(lonadj(i)*lonadj(i))
      Ky = MATMUL(MATMUL(T,Py),TRANSPOSE(Z))/Fy
      Kx = MATMUL(MATMUL(T,Px),TRANSPOSE(Z))/Fx
      Ly = T - MATMUL(Ky,Z)
      Lx = T - MATMUL(Kx,Z)
      ay = MATMUL(T,ay) + Ky*vy
      ax = MATMUL(T,ax) + Kx*vx
      Py = MATMUL(MATMUL(T,Py),TRANSPOSE(Ly)) + Qy
      Px = MATMUL(MATMUL(T,Px),TRANSPOSE(Lx)) + Qx
      lly = lly - (log(Fy) + vy*vy/Fy)/2
      llx = llx - (log(Fx) + vx*vx/Fx)/2
    END IF
  END DO
END SUBROUTINE crwn2ll
