program fern
implicit none
integer ,parameter::ITR = 1000
real*8 :: x(ITR,2),A(2,2),B(2,2),C(2,2),D(2,2),AD(4,2),ct,r1
integer:: i,j,k,l,m,n,p,q,u,t




                  

                   x(1,1) = 0.d0
                   x(1,2) = 0.d0
        
                   A(1,1) = 0.d0
                   A(1,2) = 0.d0
                   A(2,1) = 0.d0
                   A(2,2) = 0.16d0


                   B(1,1) = 0.85d0
                   B(1,2) = 0.04d0
                   B(2,1) = -0.04d0
                   B(2,2) = 0.85d0


                   C(1,1) = 0.20d0
                   C(1,2) = -0.26d0
                   C(2,1) = 0.23d0
                   C(2,2) = 0.22d0

                   D(1,1) = -0.15d0
                   D(1,2) = 0.28d0
                   D(2,1) = 0.26d0
                   D(2,2) = 0.24d0

                   AD(2,1) =0.d0
                   AD(2,2) =1.6d0
                   AD(3,1) =0.d0
                   AD(3,2) =1.6d0
                   AD(4,1) =0.d0
                   AD(4,2) =0.44d0
 
                    t = 2

     25         CALL RANDOM_NUMBER(harvest = r1)
                        ct = 100*(r1)

                
                 If (ct  .lt. 1.d0) then
                      
                      do p = 1,2
                      x(t,p) = 0.d0
                      do q = 1,2
                           x(t,p) = x(t,p) + A(p,q)*x(t-1,q)
                      end do
                      end do
                 else if ( ct .gt. 1.d0 .and. ct .lt. 86.d0) then
                      do p = 1,2
                      x(t,p) = 0.d0
                      do q = 1,2
                           x(t,p) = x(t,p) + B(p,q)*x(t-1,q)
                      end do
                      end do
                       
                      do p = 1,2
                           x(t,p) = x(t,p)  + AD(2,p)
                      end do
                 else if ( ct .gt. 86.d0 .and. ct .lt. 93.d0) then
                      do p = 1,2
                      x(t,p) = 0.d0
                      do q = 1,2
                           x(t,p) = x(t,p) + C(p,q)*x(t-1,q)
                      end do
                      end do
                       
                      do p = 1,2
                           x(t,p) = x(t,p)  + AD(3,p)
                      end do
                 else if ( ct .gt. 93.d0 .and. ct .lt. 100.d0) then
                      do p = 1,2
                      x(2,p) = 0.d0
                      do q = 1,2
                           x(t,p) = x(t,p) + D(p,q)*x(t-1,q)
                      end do
                      end do
                       
                      do p = 1,2
                           x(t,p) = x(t,p)  + AD(4,p)
                      end do
                 end if


                      write(300,*) x(t,1), x(t,2)

                      t = t+1

                      if ( t .lt. ITR+1) then
                      go to 25
                      else if (t .eq. ITR+1) then
                      go to 28
                      end if



           28         end program fern





















