   
        
        !this program checks the taylor series of lambda(k) in summer research 2013(flux-problem)


        PROGRAM STPDC
        IMPLICIT NONE
        
        INTEGER,PARAMETER::                NN = 3,ITR = 1
       	Integer::	                   m,mm,i,j,k,l,p,q,t,u,INFO,LWORK,IPIV(NN),NT(NN)          
        Integer,parameter::	           LDA = NN ,LDAA=NN-1, LWMAX = 10000
        REAL(KIND=8)::                     S,grad_s,pi, RWORK( 3*NN-2 ),trace,MDIF
        REAL(KIND=8),DIMENSION(NN)::       C,NNEBS,new_phi,mthet,n_thet,RTPHI,NNPHI,W,phi0,phi1
        REAL(KIND=8),DIMENSION(NN)::       xx,yy,thet,phi,SS,NBPHI,NPHI,NEBS,n_phi
        REAL(KIND=8),DIMENSION(NN,NN)::    HM,THM,NA,NNA,DC,DSD,dif_thet,dtheta,DS,CC,dft,CKK,WWT,IW
        REAL(KIND=8),DIMENSION(NN,NN,NN):: DDS,DDS1,DDS2,DDS3,DDS4,DDC,ddtheta,DCC
        REAL(KIND=8)::                     NDSD(NN-1,NN-1),IDSD(NN-1,NN-1),NNS(NN-1),PDF(6),AAA(2,2),III(2,2)
        complex*16::                       ai,bi,lambda(NN),NR(NN),lam1(NN),lam0(NN),lam11(NN),lam111(NN)
        Complex*16,dimension(NN,NN)::      DD,WW,UU,DEL,dlamb,JPT,AC,AK,Bk,CK,AA,AKK,BKK,DGZ1,DGZ2
        Complex*16,DIMENSION(NN,NN,NN)::   ddlamb1,ddlamb2,ddlamb,U_D,COF
        Complex*16,DIMENSION(NN,NN,NN,NN)::U_DD
        complex*16::                       WORK(LWMAX),x(NN),A(LDA,NN),n_A(NN,NN),n_lambda(NN),lambdda(NN)
        CHARACTER*1::                      UPLO

        
        

       
                        t = 1
                        pi = dacos(-1.d0)
                        ai = dcmplx(0.d0,1.d0)
                        bi = dcmplx(0.d0,1.d0)
                        
     			OPEN (unit = 110, file = 'fort.110')
     			DO l = 1,NN
     			READ(110,*)phi(l)
     			!PRINT*,"phi(",l,")=",phi(l)
     			End do
     			CLOSE(unit = 110)
       
     			
     			Write(555,*) t,phi(1),phi(2),phi(3)
                                
                  
		25	     Do p = 1,NN
                             DO q = 1,NN
                                DD(p,q) = DCMPLX(0.d0,0.d0)
                                IF (p .EQ. q) THEN
                                DD(p,q) = DD(p,q) + cdexp(-ai*phi(p))
                                ELSE IF (p .NE. q) THEN
                                DD(p,q) = DD(p,q) + DCMPLX(0.d0,0.d0)
                                end if
                             end do
                             end do

                       
                        OPEN (unit = 120, file = 'fort.120')
     			DO l = 1,NN
     			DO m = 1,NN
     			READ(120,*) WW(l,m)
     			!PRINT*,l,m,WW(l,m)
     			End do
     			END DO
     			CLOSE(unit = 120)

                      !UU = MATMUL(DD,WW)
                         Do p = 1,NN
                         DO q = 1,NN
                            UU(p,q) = dcmplx(0.d0,0.d0)
                         Do l = 1,NN
                            UU(p,q) = UU(p,q) +(DD(p,l)*(WW(l,q)))
                         END DO
                            !print*,"UU",P,Q,UU(p,q)
                         END DO
                         END DO

                       Do p = 1,NN
                       DO q = 1,NN
                           A(p,q) =dcmplx(0.d0,0.d0)
                           A(p,q) = A(p,q) + (UU(p,q) + dCONJG(UU(q,p)))*0.5
                       END DO
                       END DO 

                         AA = A
                    
                        !Lapack-for eigen value
                        LWORK = -1
                        CALL ZHEEV( "V", "L", NN, A, LDA, W, WORK, LWORK, RWORK, INFO )
                        LWORK = min( LWMAX, INT( WORK( 1 ) ) )
                        CALL ZHEEV( "V", "L", NN, A, LDA, W, WORK, LWORK, RWORK, INFO )
                   
                     IF( INFO.GT.0 ) THEN
                     WRITE(*,*)'The algorithm failed to compute eigenvalues.'
                     else if(INFO.EQ.0) THEN
                     Do p = 1,NN
                     !WRITE(*,*)t, p,"th eigen value => ",w(p)
                     !write(*,*)t,p,"th-eigen-vector=>"
                     Do q = 1,NN
                     !write(*,*) t,"A(",q,p,")=>",A(q,p)
                     End do
                     end do
                     END IF
                     Do p=1,NN
                     Do q =1,NN
                        A(q,p) = A(Q,P)
                     end do
                     end do

                     !check--------------------------------------------------
                     do p = 1,NN
                     do q = 1,NN
                     AC(p,q) = dconjg(A(q,p))
                     end do
                     end do
                     do p = 1,NN
                     do q = 1,NN
                     JPT(p,q) = 0.D0
                     DO l= 1,NN
                     JPT(p,q) = JPT(p,q) + (AC(P,l)*A(l,q))
                      END DO
                      !PRINT*,"jpt",t,JPT(P,Q)
                      END DO
                      END DO

                     

                         Do k = 1,NN
                         lambda(k) = dcmplx(0.d0,0.d0)
                         Do p = 1,NN
                         Do q = 1,NN
                         lambda(k) = lambda(k) + (dconjg(A(p,k))*UU(p,q)*A(q,k))
                         End do
                         End do
                        !print*,"lambda",lambda(k)
                         !print*,t,(dconjg(lambda(k))*lambda(k))
                         End do
                      
                        
                                If (t .eq. itr) then
                                do p = 1,NN
                                do q = 1,NN
                                WRITE(320,*) A(p,q)
                                END DO
                                END DO
                                END IF
        
!thett--------------------------------------------------------------------              
                            
                     Do k = 1,NN
                     thet(k) = ai*cdlog(lambda(k))
                     end do
                                  Do k = 1,NN
                                  !WRITE(207,*) t,"THET(",K,")=",THET(K)
                                  thet(k) = -ai*cdlog(cdexp(ai*thet(k)))
                                  !WRITE(207,*)"II", t,"THET(",K,")=",THET(K)
                                   end do
                                  !print*,"sumthett",thet(1)+thet(2)+thet(3)

                                  !Check------------------------------------------------------------
                                    DO k = 1,NN
                                   lambdda(k) = cdexp(-ai*thet(k))
                                   !print*,"check",lambdda(k)
                                    End do
   

 !CALL UDAS(t,NN,DEL,UU,U_D,ai)
                        do p= 1,NN
                        do q = 1,NN
                           if (p .eq.q) then
                           DEL(p,q) =  dcmplx(1.d0,0.d0)
                           else if(p .ne. q) then
                           DEL(p,q) =  dcmplx(0.d0,0.d0)
                           end if
                       end do
                       end do 

                       ai = dcmplx(0.d0,1.d0)
                       DO l = 1,NN
                       Do p= 1,NN
                       DO q = 1,NN
                        U_D(l,p,q) = dcmplx(0.d0,0.d0)
                        U_D(l,p,q) =  U_D(l,p,q) - (ai)*DEL(p,l)*UU(p,q)
                        !write(*,*) l,p,q,U_D(l,p,q)
                             end Do
                             end do
                             End do
                          
 !CALL UDDAS(t,NN,DEL,UU,U_DD)
                       DO l = 1,NN
                       DO m = 1,NN
                       Do p= 1,NN
                       DO q = 1,NN
                        U_DD(l,m,p,q) = dcmplx(0.d0,0.d0)
                        U_DD(l,m,p,q) =  U_DD(l,m,p,q) - DEL(p,l)*DEL(p,m)*UU(p,q)
                        !write(*,*) l,m,p,q,U_DD(l,m,p,q)
                             end Do
                             end do
                             End do
                             END DO
 ! CALL LLAMB(t,NN,dlamb,A,U_D)
                       Do k = 1, NN
                       Do l = 1, NN
                       dlamb(k,l) = dcmplx(0.d0,0.d0)
                       Do p = 1,NN
                       Do q =  1,NN
                       dlamb(k,l) = dlamb(k,l) + dconjg(A(p,k))*U_D(l,p,q)*A(q,k)
                       End do
                       end do
                       !print*,k,l,"dlamb",dlamb(k,l)
                       end do
                       end do
                   

  !CALL COFF(t,NN,A,U_D,lambda,COF)
                      Do l = 1,NN
                      DO k = 1,NN
                      DO q = 1,NN
                       COF(l,k,q) = DCMPLX(0.D0,0.D0)
                         IF (k .ne. q) then
                           DO i = 1,NN
                           DO j = 1,NN
                            COF(l,k,q) =  COF(l,k,q) +&
                             & (((dconjg(A(i,q)))*U_D(l,i,j)*A(j,k))/(lambda(k)-lambda(q)))
                           END do
                           end do
                        end if
                       END do
                       END DO
                       end do

                       !write(*,*) COF(1,1,3),-dconjg(COF(1,3,1))

 !CALL GLAMBDA(t,NN,lambda,ddlamb1,ddlamb2,ddlamb,A,U_DD,COF)
                        !dlamb1----------------------------------------------------------             
                       DO k = 1,NN
                       Do l = 1,NN
                       Do m = 1,NN
                       ddlamb1(k,l,m) = dcmplx(0.d0,0.d0)
                           Do p = 1,NN
                           Do q = 1,NN
                           ddlamb1(k,l,m) =  ddlamb1(k,l,m) + (dconjg(A(p,k))*U_DD(l,m,p,q)*A(q,k))
                           END DO
                           END DO

                       END DO
                       END DO
                       END DO 
                  !dlamb2------------------------------------------------------------------
                       DO k = 1,NN
                       Do l = 1,NN
                       Do m = 1,NN
                       ddlamb2(k,l,m) = dcmplx(0.d0,0.d0)
                          Do p = 1,NN
                           If (p.ne.k) then
                            ddlamb2(k,l,m) =  ddlamb2(k,l,m) + ((COF(l,p,k)*COF(m,k,p)&
                                           & + COF(m,p,k)*COF(l,k,p))*(lambda(k)-lambda(p)))
                           end if
                          END DO
                       END DO
                       END DO
                       END DO 
                !dlamb--------------------------------------------------------------------------
                       DO k = 1,NN
                       Do l = 1,NN
                       Do m = 1,NN
                          ddlamb(k,l,m) = dcmplx(0.d0,0.d0)
                          ddlamb(k,l,m) =  ddlamb(k,l,m) + (ddlamb1(k,l,m)+ddlamb2(k,l,m))
                          !write(*,*) k,l,m,"=", ddlamb(k,l,m)
                       END DO
                       END DO
                       END DO 

                !Hello! check taylor series--------------------------------------------------------------------------------------

                       if (t .EQ. 1) then
                         DO k = 1,NN
                           lam0(k) = lambda(k)
                           phi0(k) = phi(k)
                           !print*,t,k,phi0(k),lambda(k),lam0(k)
                         end do
                       end if


                       if (t .EQ. 2) then
                         DO k = 1,NN
                           lam1(k) = lambda(k)
                           phi1(k) = phi(k)
                           !print*,t,k,phi1(k),lambda(k),lam1(k)
                         end do
                      
                          do k = 1,NN
                          lam11(k) = 0.d0
                          do l = 1,NN
                          lam11(k) = lam11(k) +(dlamb(k,l)*(phi1(l)-phi0(l)))
                          end do
                          end do
                          do k = 1,NN
                          lam11(k) = lam11(k) + lam0(k)
                          end do
                      
                      
                       do k = 1,NN
                       !print*,k,"lam1=",lam1(k),"lam11=",lam11(k)
                       end do

                     
                          do k = 1,NN
                          lam111(k) = 0.d0
                          do l = 1,NN
                          do m = 1,NN
                          lam111(k) = lam111(k) + ( (0.5)*ddlamb(k,l,m)*(phi1(l)-phi0(l))*(phi1(m)-phi0(m)))
                          end do
                          end do
                          end do
                          do k = 1,NN
                          lam111(k) = lam111(k) + lam11(k)
                          end do
                       



                       do k = 1,NN
                       !print*,k,"lam1=",lam1(k),"lam111=",lam111(k)
                       write(71,*)lam0(k),lam1(k)
                       end do
                       !print*,t

                       end if
                       !end check-------------------------------------------------   



 !CALL GTHETA(t,NN,bi,dtheta,ddtheta,dlamb,ddlamb,lambda)
                     !dtheta-----------------------------------------------------    
                       bi = dcmplx(0.d0,1.d0)
                       DO k = 1,NN
                       Do l = 1,NN
                       dtheta(k,l) = 0.d0
                       dtheta(k,l) = dtheta(k,l) +  bi*dconjg(lambda(k))*dlamb(k,l)
                       write(*,*)k,l, "dtheta=",dtheta(k,l)
                       end do
                       end do                     
                  !ddtheta--------------------------------------------------------     
                       DO k = 1,NN
                       Do l = 1,NN
                       DO m = 1,NN
                       ddtheta(k,l,m) = dcmplx(0.d0,0.d0)
                       ddtheta(k,l,m) = ddtheta(k,l,m) -(bi*((dconjg(lambda(k)))**2)*(dlamb(k,m)*dlamb(k,l))) &
                                         & + (bi*(dconjg(lambda(k))*ddlamb(k,l,m)))
                       write(*,*)k,l,m,"ddtheta=",ddtheta(k,l,m)
                       end do
                       end do
                       END DO 


                              
                                t = t+1

                        phi(1) = 0.13d0
                        phi(2) = 0.26d0
                        phi(3) = -0.39d0
                       

                        if (t .lt.ITR+1 ) then
                        GO TO 25
                        else if (t .eq. ITR+1) then
                        go to 28
                        end if

                        28 END PROGRAM STPDC !************************************************




        
          
           
     
      
                     

            
                   

   


   
