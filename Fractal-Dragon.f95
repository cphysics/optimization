program coch
implicit none
integer ,parameter::ITR = 14
integer, parameter :: NITR = (2**(ITR+1) + 1)
real*8 :: x(NITR,2),A(2,2)
integer:: i,j,k,l,m,n,p,q,u,t,nbr,mbr




                   

                   x(1,1) = 0
                   x(1,2) = 0

                   x(2,1) = 1
                   x(2,2) = 0
        
                   A(1,1) = 0.d0
                   A(1,2) = 1.d0
                   A(2,1) = -1.d0
                   A(2,2) = 0.d0
                   
                   Do p = 1,2
                   x(3,p) = 0.d0
                   Do q = 1,2
                   x(3,p) = x(3,p) + A(p,q)*x(2,q)
                   end do
                   end do
                  

                   do k = 1,3
                   do l = 1,2
                   x(k,l) = x(k,l) - x(3,l)
                   end do
                   end do
     
                   do k = 4,5
                   Do p = 1,2
                   x(k,p) = 0.d0
                   Do q = 1,2
                       if ( k .eq. 4) then
                        l = 2
                       else if (k .eq. 5) then
                        l = 1
                       end if
                      x(k,p) = x(k,p) + A(p,q)*x(l,q)
                   end do
                   end do
                   end do


                   do k = 1,5
                   do l = 1,2
                   x(k,l) = x(k,l) - x(5,l)
                   end do
                   end do

                   Do k = 6,9
                   Do p = 1,2
                   x(k,p) = 0.d0
                   Do q = 1,2
                       if ( k .eq. 6) then
                        l = 4
                       else if (k .eq. 7) then
                        l = 3
                        else if ( k .eq. 8) then
                        l = 2
                       else if (k .eq. 9) then
                        l = 1
                       end if
                   x(k,p) = x(k,p) + A(p,q)*x(l,q)
                   end do
                   end do
                   end do


                   t = 3
               
      25           nbr = (2)**(t) + 1

                   do k = 1,nbr
                   do l = 1,2
                   x(k,l) = x(k,l) - x(nbr,l)
                   end do
                   end do              

                   mbr = 2*nbr - 1

                   do k = nbr+1, mbr
                   Do p = 1,2
                   x(k,p) = 0.d0
                   Do q = 1,2

                        do i = 1,nbr
                        if ( k .eq. nbr+i) then
                        l = nbr - i
                        end if
                        end do

                        x(k,p) = x(k,p) + A(p,q)*x(l,q)

                   end do
                   end do
                   end do

                   do k = 1,mbr
                      write(200,*) x(k,1),x(k,2)
                   end do
             
                   t = t+1

                   if (t .lt. ITR) THEN
                   go to 25
                   ELSE IF (t .gt. ITR) then
                   go to 28
                   end if



        28         end program coch





















