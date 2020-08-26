program mandelbrotset
implicit none
integer ,parameter::ITR = 1000,trial = 1000000000
real*8 :: c1,c2,abz,x,y,pi,r1,r2,nrmz
real*8, parameter:: mlt = 2.d0
integer:: i,j,k,l,m,n,p,q,u,t
complex*16 :: z,c
 pi = dacos(-1.d0)

       u = 1

        26        t = 1
                  z = dcmplx(0.d0,0.d0)
                  call random_number(harvest = r1)
    	          call random_number(harvest = r2)
                  c1 = 5*(2*r1 - 1)
                  c2 = 5*(2*r2 - 1)
                  c = dcmplx(c1,c2)

        25        z = z*z + c
                  nrmz = z*z
                  abz = dabs(nrmz)

      t = t+1

                  if (abz .lt. MLT) then
                     if (t .lt. ITR + 1) then
                     go to 25 
                     else if ( t .eq. ITR + 1) then
                     x = REAL(REAL(c))
                     y = REAL(AIMAG(c))
                     write(100,*) x,y
                     u = u +1
                         if (u .lt. trial + 1) then
                         go to 26
                         else if (u .eq. trial + 1) then
                         go to 28
                         end if   
                     end if
                  else if (abz .gt. MLT) then
                      u = u + 1
                         if (u .lt. trial + 1) then
                         go to 26
                         else if (u .gt. trial + 1) then
                         go to 28
                         end if
                  end if
                  


     28          end program mandelbrotset




















