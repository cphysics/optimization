   
        



        PROGRAM ISING
        IMPLICIT NONE

        Integer,parameter::	          N = 100 ,ITR = 100       
        Real,PARAMETER::                  mu = 0.2d0,J = 2.d0,q = 0.2d0,kk = 0.00000143,pn = 0.5d0
        Integer::	                  i,j,k,l,p,q,t,u,s,Na,Nb,latt(N),lb
   	real(KIND = 8)::                  r1,r2,r3,r4,r5,mu,DEL_E ,prb           
        real(KIND = 8)::                  L,M,MM

                      
                      Na = 0
                      Nb = 0

                   


                      t = 1

              25      Call random_number(harvest = r1)

                      if (r1 .lt. 0.5d0) then
                       latt(t) = -1
                       Na = Na + 1
                      else
                       latt(t) = 1
                       Nb = Nb + 1
                      end if
                      
                      if (t .lt. N) then
                      t = t+1
                      go to 25
                      else if( t .eq. N) then
                      go to 26
                      end if

            26         LB = 0
                      Do k = 1,N
                      LB = LB + latt(k)
                      LB = LB/N
                      end do
                      
                      !BM = (q*J*LB*)/mu
                      
                     L = (2.d0*(Nb/Na) - 1)   
                     M = (Nb - Na)*mu      
                     MM = lat*mu*L     
                     !Print*,M,MM
                 
                    
                     E(1) = - 0.5d0*(q*J*N*L**2) - (mu*BA*N*LB)

                     
                              tt = 2
                   
            27             latt(tt) = latt(tt)*(-1)
                            if (latt(tt) .eq. 1) then
                                Na = Na - 1
                                Nb = Nb + 1
                            else if(latt(tt) .eq. -1) then
                                Nb = Nb - 1
                                Na = Na + 1
                            end if

                    
                                LB = 0
                                    Do k = 1,N
                     		    LB = LB + latt(tt)
                     		    LB = LB/N
                      		    end do
                      
                      		   ! BM = (q*J*LB*)/mu
                      
                        	    L = (2.d0*(Nb/Na) - 1)   
                                    M = (Nb - Na)*mu 
                                         
                                    MM = lat*mu*L  
                         PRINT*,MM,
                                     tt = tt +1    
                     
                    
                     E(tt) =  - 0.5d0*(q*J*N*L**2) - (mu*BA*N*LB)
                     DEL_E =  (E(tt) - E(tt-1))/(kk*tmp)
                     prb = dexp(-DEL_E/(KK*tmp))
                     if (prb .lt. pn) then
                     latt(tt) = (-1)*(latt(tt))
                     end if

           
                     
                     
                    













                                                   
                        END PROGRAM ISING
    

         
