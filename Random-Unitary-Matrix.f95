  
        PROGRAM STPDC
        IMPLICIT NONE

        Integer,parameter::	           NN = 4
        Integer::	                   i,j,k,l,p,q,t,u,s  
   	real(KIND=8)::                     T1,T2,T3,phii,xi,theta,pi
        Complex*16::                       SU2(2,2),SUN(NN,NN),SUNM(NN,NN),SUNP(NN,NN),WW(NN,NN)
        complex*16::                       ai,a_a,b_b,II(NN,NN),WC(NN,NN),UC(NN,NN),III(NN,NN)
       
        
        

       
        pi = dacos(-1.d0)
        ai = dcmplx(0.d0,1.d0)


                                      s = 1
                       26             t = s+1
                              
                                            
                                       
                       25              call random_number(harvest=T1)
                                       xi = (pi*(2*T1-1))
				       call random_number(harvest=T2)
				       theta =( 0.5*pi*(T2))
				       call random_number(harvest=T3)
				       phii = (pi*(2*T3-1))
                                       a_a = dcos(theta)*(cdexp(ai*phii))
				       b_b = dsin(theta)*(cdexp(ai*xi))
                                                   SU2(1,1) = a_a
				                   SU2(1,2) = b_b
                                                   SU2(2,2) = dconjg(SU2(1,1))
				                   SU2(2,1) = -dconjg(SU2(1,2))
                                               
                  

       
                            
                            
                             Do p = 1,NN
                             DO q = 1,NN
                                IF (p .EQ. q) THEN
                                SUN(p,q) = DCMPLX(1.d0,0.d0)
                                ELSE IF (p .NE. q) THEN
                                SUN(p,q) = DCMPLX(0.d0,0.d0)
                                end if
                             end do
                             end do


                     
                          SUN(s,s) =     SU2(1,1)
                          SUN(s,t) =   SU2(1,2)
                          SUN(t,s) =   SU2(2,1)
                          SUN(t,t) = SU2(2,2)
                          
                                                           DO p =1,NN
                                                           Do q = 1,NN
                                                           UC(p,q) = dconjg(SUN(q,p))
                                                           End do
                                                           End do 

                                                           II = matmul(UC,SUN)

                                                           DO p =1,NN
                                                           Do q = 1,NN
                                                           !write(501,*)"II", p,q, II(p,q)
                                                           End do
                                                           End do 
                                                           Write(501,*) "space"
                        


                          IF (s .eq.1 .and. t .eq. 2) then
                             SUNM = SUN
                          else 
                             SUNP = MATMUL(SUNM,SUN)
                             SUNM = SUNP
                          end if

                              t = t+1
                         
                          if (t .lt. NN+1) then
                              go to 25
                          else if (t .eq. NN+1) then
                               s = s+1
                               if (s .lt. NN) THEN
                               go to 26
                               ELSE if (s .eq. NN) then
                               WW = SUNM

                                                           DO p =1,NN
                                                           Do q = 1,NN
                                                           write(120,*)  WW(p,q)
                                                           End do
                                                           End do 

                                                           DO p =1,NN
                                                           Do q = 1,NN
                                                           WC(p,q) = dconjg(WW(q,p))
                                                           End do
                                                           End do 

                                                           III = matmul(WC,WW)

                                                           DO p =1,NN
                                                           Do q = 1,NN
                                                           !write(*,*)"III", p,q, III(p,q)
                                                           End do
                                                           End do 



                               GO TO 27
                               end if
                          end if

           

                                                   
               27       END PROGRAM STPDC 
         
