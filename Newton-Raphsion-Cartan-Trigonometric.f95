   
        



        PROGRAM STPDC
        IMPLICIT NONE
        
        INTEGER,PARAMETER::                NN = 3,ITR = 30
       	INTEGER::                          I,J,K,L,M,P,Q,t
        REAL::                             IDSD(NN,NN),NEBF(NN),X0(NN),X(NN),FF,XX(NN),F(NN)
        REAL::                             DF(NN,NN),DDF(NN,NN,NN),DFT(NN),DDFT(NN,NN)
        REAL,DIMENSION(NN,NN)::            HM,THM,TDDF,NA
         
       	Integer::	                   INFO,LWORK,IPIV(NN),NT(NN)          
        Integer,parameter::	           LDAA=NN-1, LWMAX = 10000
        REAL::                             pi, RWORK( 3*NN-2 ),AAA(NN-1,NN-1)
        
        REAL,DIMENSION(NN)::               NX,NNX,RTX,new_X,TDF,C(NN,NN)
       
        REAL::                             PDDF(NN-1,NN-1),IDDF(NN-1,NN-1),PDF(NN-1),III(2,2)
       
        complex::                           ai
        REAL::                              WORK(LWMAX),A(LDAA,NN)
        CHARACTER*1::                       UPLO
                            
                          

                            !f(x,y,z) = sin(x-x0)^2 * sin(y-x0)^2 * sin(z-x0)^2
                                     !+ sin(x-y0)^2 * sin(y-y0)^2 * sin(z-z0)^2
                                     !+ sin(x-z0)^2 * sin(y-z0)^2 * sin(z-z0)^2
                                     !constraint:x+y+z=0


                               t = 0
                               x0(1) = 0.15
                               x0(2) = 0.30
                               x0(3) = -0.45
                               x(1) = 0.18
                               x(2) = 0.36
                               x(3) = -0.54
                               do i = 1,NN
                               print*,t,"x0 =",x0(i)
                               end do
                               do i = 1,NN
                               print*,t,"x =",x(i)
                               end do
                         25    DO i = 1,NN
                                  F(i) = 1.0
                               DO j = 1,NN
                                  F(i) = F(i)*(sin(x(j) - x0(i)))**2
                               end do
                               end do

                              Do i = 1,NN
                                  FF = 0.0
                                  FF= FF + F(i)
                              END DO
                              WRITE(91,*)t,"FF=",FF
                              Do i = 1,NN
                              DO P = 1,NN
                                  C(i,p) =(1/(tan(x(p) - x0(i))))
                              END DO
                              END DO

                              Do i = 1,NN
                              DO P = 1,NN
                                  DF(i,p) = 0.0
                                  DF(i,p) = DF(i,p) + (2.0*F(i)*C(i,p))
                              END DO
                              END DO
                              Do p = 1,NN
                                  DFT(p) = 0.0
                                  DO i = 1,NN
                                  DFT(p) = DFT(p) + DF(i,p)
                                  END DO
                                  WRITE(92,*) t,"DFT",p,DFT(p)
                              END DO
                           
                              
                              Do i = 1,NN
                              DO p = 1,NN
                              DO q = 1,NN
                                  DDF(i,p,q) =0.0
                                  if (p .eq.q) then
                                  DDF(i,p,q) = DDF(i,p,q) + ((2.0*DF(i,q)*C(i,p)) &
                                                           &- (2.0*F(i)*(1+C(i,p)**2)))
                                  else if (p.ne.q) then
                                  DDF(i,p,q) = DDF(i,p,q) + (2.0*DF(i,q)*C(i,p))
                                  end if
                              END DO
                              END DO
                              END DO
                             
                              DO p = 1,NN
                              DO q = 1,NN
                                  DDFT(p,q) =0.0
                                  DO i = 1,NN                                  
                                  DDFT(p,q) = DDFT(p,q) + DDF(i,p,q)
                                  END DO
                                  !print*,p,q,DDFT(p,q)
                              END DO
                              END DO
         !CARTAN'S METHOD-------------------------------------
         DO P = 1,NN
         DO q = 1,NN
    	 IF (p .lt. NN)then
         if (q .eq. p) then
         HM(p,q) = (1.d0/dsqrt(dfloat(p*(p+1))))
         ELSE If (q .lt. p) then 
         HM(p,q)  = (1.d0/dsqrt(dfloat(p*(p+1))))
         else If ( q .EQ. p+1) then
         HM(p,q) = ((0.d0 - dfloat(p))/dsqrt(dfloat(p*(p+1))))
         else IF (q .gt. p+1) then
         HM(p,q) = (0.d0)
         end if
         ELSE IF (P .EQ. NN) THEN
         HM(p,q) = (1.d0/dsqrt(dfloat(p*(p+1))))
      	 END IF
     	 End do
     	 END DO


    
    								Do  p = 1,NN
   								TDF(p) = 0.d0
    								Do l= 1,NN
     								TDF(p) = TDF(p) + HM(p,l)*DFT(l)
    								end do
    								end do
    

       		DO k= 1,NN-1
       		PDF(k)  = TDF(k)
       		END DO


   
     								DO p = 1,NN
     								DO q = 1,NN
    							        THM(p,q) = HM(q,p)
     								end do 
     								end do
    
     								DO p = 1,NN
     								DO q = 1,NN
       								NA(p,q) = 0.d0
      							        do l = 1,NN
       								NA(p,q) =NA(p,q) + (DDFT(p,l)*THM(l,q))
      							        end do
   							        end do 
         							end do

    								DO p = 1,NN
    								DO q = 1,NN
       								TDDF(p,q) = 0.d0
       								do l = 1,NN
       								TDDF(p,q) =TDDF(p,q) + (HM(p,l)*NA(l,q))
       								end do
    								end do 
    								end do

    		DO p = 1,NN-1
    		Do q = 1,NN-1
    		PDDF(p,q) = TDDF(p,q)
   	        AAA(p,q) = TDDF(p,q)
                !print*,"AAA",AAA(p,q)
   	        END DO
                END DO


                        CALL SGETRF( NN-1,NN-1,PDDF,LDAA,IPIV,INFO)
                        If (info .eq.0) then
    	                LWORK = -1
     		        CALL SGETRI( NN-1, PDDF, LDAA, IPIV, WORK,LWORK, INFO )
      		        LWORK = min( LWMAX, INT( WORK( 1 ) ) )
     		        CALL SGETRI( NN-1, PDDF, LDAA, IPIV, WORK, LWORK, INFO )
     		        if (info .ne. 0) then
     		        PRINT*, 'Matrix inversion failed!'
     		        end if
     		        end if


     		DO l = 1,NN-1
     		Do m = 1,NN-1
     		IDDF(l,m) = PDDF(l,m)
                !print*,"IDDF",IDDF(l,m)
     		End do
     		End do
                III = matmul(IDDF,AAA)
                
                DO l = 1,NN-1
     		Do m = 1,NN-1
     	        !WRITE(*,*)"Identity",III(l,m)
     		End do
     		End do
               

                

     

     		DO p = 1,NN-1
     		NX(p)  = 0.0
     		DO i= 1,NN-1
     		NX(p)  = NX(p) + IDDF(p,i)*PDF(i)
     		end do
     		end do
   
                   		Do p = 1,NN
                   		if (p .lt.NN) then 
                   		NNX(p) = NX(p)
                   		else if (p .eq. NN) then
                   		NNX(p) = 0.0
                   		END IF
                   		end do
                  
      
    
                							DO l = 1,NN
     									RTX(l)  = 0.0
     									DO i= 1,NN
     									RTX(l)  = RTX(l) + THM(l,i)*NNX(i)
     									end do
                                                                        !print*,t,l,rtx(l)
     									end do
               
                		do p = 1,NN
                		new_X(P) = 0.0
                		new_X(p) = new_X(p) + (X(p) - RTX(p))
                		end do

                             ai = (0.0,1.0)
                             do p = 1,NN
                             X(p) = -ai*clog(cexp( ai*new_X(p)))
                             end do
                             t = t+1
                             print*,t,"x", x(1),x(2),x(3)
                             
                              
                               if (t .lt. itr) then
                               go to 25
                               else if (t.eq. itr) then
                               go to 26
                               end if 
                              
                       

                       26 END PROGRAM STPDC !************************************************

















