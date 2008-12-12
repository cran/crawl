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
    Fy = Py(1,1) + tau2y(i)
    Fx = Px(1,1) + tau2x(i)/(lonadj(i)*lonadj(i))
    IF(loctype(i)==1 .OR. Fy==0.0) THEN
      ay = MATMUL(T,ay)
      Py = MATMUL(MATMUL(T,Py),TRANSPOSE(T)) + Qy
    ELSE
      vy = y(i)-ay(1,1)
      lly = lly - (log(Fy) + vy*vy/Fy)/2
      Ky = MATMUL(MATMUL(T,Py),TRANSPOSE(Z))/Fy
      Ly = T - MATMUL(Ky,Z)
      ay = MATMUL(T,ay) + Ky*vy
      Py = MATMUL(MATMUL(T,Py),TRANSPOSE(Ly)) + Qy
    END IF
    IF(loctype(i)==1 .OR. Fx==0.0) THEN      
      ax = MATMUL(T,ax)
      Px = MATMUL(MATMUL(T,Px),TRANSPOSE(T)) + Qx
    ELSE
      vx = x(i)-ax(1,1)
      llx = llx - (log(Fx) + vx*vx/Fx)/2
      Kx = MATMUL(MATMUL(T,Px),TRANSPOSE(Z))/Fx
      Lx = T - MATMUL(Kx,Z)      
      ax = MATMUL(T,ax) + Kx*vx      
      Px = MATMUL(MATMUL(T,Px),TRANSPOSE(Lx)) + Qx
    END IF
  END DO
END SUBROUTINE crwn2ll
