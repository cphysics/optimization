


       !------------------------------------------------------------------
       !This program calculates the expectation value of kinetic energy
       !of an idealgas using Metropolis algorithm.
       !-----------------------------------------------------------------


        PROGRAM  idealgas
        IMPLICIT NONE
        INTEGER,PARAMETER::                   N = 100,itr = 1000,maxv = 1000
       	Integer::	                          i,j,k,l,p,q,t,r,a,success
        real*8,dimension(N)::                 v,u,KE
        real*8::                              SS,SSN,m,beta,Temp,av_E,DSS,kb
        real,dimension(ITR)::                 rn,r1,r2



        Temp = 300.0
        kb = 1.d0
        beta = 1/(T*kb)
        m = 3.84
        success = 0





                            do k = 1,N
                            v(k) = maxv
                            u(k) = v(k)
                            end do
                            do k = 1,N
                            KE(k) = 0.5*m*(v(k)**2)
                            end do
                            SS = 0.d0
                            do k = 1,N
                            SS = SS +  0.5*m*(v(k)**2)
                            end do


          t = 1

     98     call random_number(harvest=rn(t))
            a = int(rn(t)*N)
            if (a .eq. 0) then
            a = a+1
            end if

            !updating v

            call random_number(harvest=r1(t))
            v(a) = maxv*r1(t)
            KE(a) = 0.5*m*(v(a)**2)
            SSN = 0.d0
            do k = 1,N
            SSN = SSN +  KE(k)
            end do

            !difference SSN - SS
            DSS = SSN - SS
            if (DSS .lt. 0) then
            success = success+1
            av_E = SSN/(N)
            write(200,*) success,av_E
            write(333,*) success,sqrt(av_E*2*m)
            !store new configuration for next cycle
            do k = 1,N
            u(k) = v(k)
            end do
            SS = SSN
            else if (DSS .gt. 0) then
            call random_number(harvest = r2(t))
            if (r2(t) .gt. dexp(-DSS)) then
            success = success+1
            av_E = SSN/(N)
            write(200,*) success,av_E
            write(333,*) success,sqrt(av_E*2*m)
            !store new configuration for next cycle
            do k = 1,N
            u(k) = v(k)
            end do
            SS = SSN
            else if (r2(t) .lt. dexp(-(beta*DSS))) then
            v(a) = u(a)
            end if
            end if

     t = t+1


                        if (t .lt. itr+1) then
                        go to 98
                        else if (t .eq. itr+1) then
                        print*,success
                        end if

                        end program idealgas