   
        PROGRAM STPDC
        IMPLICIT NONE
        
        INTEGER,PARAMETER::                NN = 3,ITR = 20
       	INTEGER::                          I,J,K,L,M,P,Q,t
        REAL::                             IDSD(NN,NN),NEBF(NN),X0(NN),X(NN),FF,XX(NN)
                          

                            !f(x,y,z) = (x-1)^3 + (y-2)^3 + (z-3)^3
                               t = 0
                               x0(1) = 100.0
                               x0(2) = 100.0
                               x0(3) = 100.0
               
                          25   DO k = 1,NN
                                   NEBF(k) = (3*(x0(k)-k)**2)
                               end do
                               
                              
                               DO l = 1,NN
                               DO M = 1,NN
                               IDSD(l,m) = 0.0
                               if (l .eq. m) then
                                     IDSD(l,m) = IDSD(l,m) + (1/(6*(x0(l)-l)))
                               else if (l .ne. m) then
                                     IDSD(l,m) = 0.0
                               end if
                               !print*,"IDSD",l,m,IDSD(L,M)
                               END DO
                               END DO
                               
                               DO p = 1,NN
                                  x(p) = 0.0
                                  xx(p) = 0.0
                               DO l = 1,NN
                                xx(p) = xx(p)  - (IDSD(p,l)*NEBF(l))
                               END DO
                                x(p) = x(p) + (x0(p) + xx(p))
                               print*,"x(",p,")=" ,x(p)
                               END DO

                               FF = (x(1)-1)**2+(x(2)-2)**2 + (x(3)-3)**2
                               print*, t,"FF=",FF 
                               
                               DO K = 1,NN
                                x0(k) = x(k)
                               END DO

                               t = t+1
                               if (t .lt. itr) then
                               go to 25
                               else if (t.eq. itr) then
                               go to 26
                               end if 
                              
                       

                       26 END PROGRAM STPDC !************************************************