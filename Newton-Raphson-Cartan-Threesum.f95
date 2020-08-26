   
        !this program assums the action:
        ! S = (thet(1) - phi(1))**2 + (thet(2) - phi(2))**2 + (thet(3) - phi(3))**2

        ! -0.76918774851248162     
        !2.5652925507454345     
        !-1.7961048022329535     


        PROGRAM STPDC
        IMPLICIT NONE
        
        INTEGER,PARAMETER::                NN = 3,ITR = 10
       	Integer::	                   m,mm,i,j,k,l,p,q,t,u,INFO,LWORK,IPIV(NN),NT(NN)          
        Integer,parameter::	           LDA = NN ,LDAA=NN-1, LWMAX = 10000
        REAL(KIND=8)::                     S,grad_s,pi, RWORK( 3*NN-2 ),trace,MDIF
        REAL(KIND=8),DIMENSION(NN)::       CC,NNEBS,new_phi,mthet,n_thet,RTPHI,NNPHI,W,NR
        REAL(KIND=8),DIMENSION(NN)::       xx,yy,thet,phi,SS,NBPHI,NPHI,NEBS,n_phi
        REAL(KIND=8),DIMENSION(NN,NN)::    HM,THM,NA,NNA,DSD,dif_thet,DS,dCC,dft,CKK,dtheta,HMTHM
        REAL(KIND=8),DIMENSION(NN,NN,NN):: DDS,DDS1,DDS2,DDS3,DDS4,DDC,ddtheta
        REAL(KIND=8)::                     NDSD(NN-1,NN-1),IDSD(NN-1,NN-1),NNS(NN-1),PDF(6),AAA(2,2),III(2,2)
        complex*16::                       ai,bi,lambda(NN)
        Complex*16,dimension(NN,NN)::      DD,WW,UU,DEL,dlamb,JPT,AC,AK,Bk,CK,AA,AKK,BKK,DGZ1,DGZ2
        Complex*16,DIMENSION(NN,NN,NN)::   ddlamb1,ddlamb2,ddlamb,U_D,COF
        Complex*16,DIMENSION(NN,NN,NN,NN)::U_DD
        complex*16::                       WORK(LWMAX),x(NN),A(LDA,NN),n_A(NN,NN),n_lambda(NN)
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
                                
                  
			     Do p = 1,NN
                             DO q = 1,NN
                                DD(p,q) = DCMPLX(0.d0,0.d0)
                                IF (p .EQ. q) THEN
                                DD(p,q) = DD(p,q) + cdexp(-ai*phi(p))
                                ELSE IF (p .NE. q) THEN
                                DD(p,q) = DD(p,q) + DCMPLX(0.d0,0.d0)
                                end if
                             end do
                             end do

                       
                 25     OPEN (unit = 120, file = 'fort.120')
     			DO l = 1,NN
     			DO m = 1,NN
     			READ(120,*) WW(l,m)
     			!PRINT*,l,m,WW(l,m)
     			End do
     			END DO
     			CLOSE(unit = 120)

                        UU = MATMUL(DD,WW)

                       
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
                     do p = 1,NN
                     do q = 1,NN
                     A(q,p) = A(q,p)
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
                      !PRINT*,"jpt",JPT(P,Q)
                      END DO
                      END DO


                         Do k = 1,NN
                         lambda(k) = dcmplx(0.d0,0.d0)
                         Do p = 1,NN
                         Do q = 1,NN
                         lambda(k) = lambda(k) + (dconjg(A(p,k))*UU(p,q)*A(q,k))
                         End do
                         End do
                         !print*,(dconjg(lambda(k))*lambda(k))
                         End do
                      
                        
                              
                      
                              ! CHECH-DIAGONALIZE--------------------------------------

                              DGZ1 = MATMUL(AA,A)
                              DGZ2 = MATMUL(AC,DGZ1)
                              do p = 1,NN
                              do q = 1,NN
                              !WRITE(*,*)"DGZ",p,q,DGZ2(p,q)
                              END DO
                              END DO

                              !Check(AA)-----------------------------------------------
                               Do p = 1,NN
                               Do q = 1,nn
                                 AKK(q,p) = DCMPLX(0.d0,0.d0)
                               END DO
                               END DO
                               
                               Do k = 1,NN
                               DO p = 1,NN
                               DO q = 1,NN
                                 AKK(p,k) = AKK(p,k) + AA(p,q)*A(q,k)
                               END DO
                               END DO
                               END DO
                               
                               Do p = 1,NN
                               Do q = 1,NN
                                 BKK(q,p) = W(p)*A(q,p)
                               END DO
                               END DO 
                               
                               Do p = 1,NN
                               Do q = 1,NN
                                 CKK(q,p) =  AKK(q,p)/A(q,p)
                               END DO
                               END DO
                                

                               Do p = 1,NN
                               DO q =  1,NN
                                !PRINT*,"W(P)",W(p),"CKK",CKK(q,p)
                               End do
                               End do

                               Do p = 1,NN
                               DO q =  1,NN
                               !PRINT*,"AKK",AKK(p,q),"BKK",BKK(p,q)
                               End do
                               End do
                              !Check(UU)-----------------------------------------
                               
                               Do p = 1,NN
                               Do q = 1,nn
                                 AK(q,p) = DCMPLX(0.d0,0.d0)
                               END DO
                               END DO
                               
                               Do k = 1,NN
                               DO p = 1,NN
                               DO q = 1,NN
                                 AK(p,k) = AK(p,k) + UU(p,q)*A(q,k)
                               END DO
                               END DO
                               END DO
                               
                               Do p = 1,NN
                               Do q = 1,nn
                                 BK(q,p) = lambda(p)*A(q,p)
                               END DO
                               END DO 
                               
                               Do p = 1,NN
                               Do q = 1,nn
                                 CK(q,p) =  AK(q,p)/A(q,p)
                               END DO
                               END DO 

                               Do p = 1,NN
                               DO q =  1,NN
                               !PRINT*,"lambda",lambda(p),"CK",CK(q,p)
                               End do
                               End do

                               Do p = 1,NN
                               DO q =  1,NN
                               !PRINT*,"AK",AK(p,q),"BK",BK(p,q)
                               End do
                               End do
                        !----------------------------------------------------------------
    !CALL THETT(NN,XX,YY,THET,LAMBDA,t)
                     Do k = 1,NN
                     thet(k) = ai*cdlog(lambda(k))
                     end do
                                  Do k = 1,NN
                                  !WRITE(207,*) t,"THET(",K,")=",THET(K)
                                  thet(k) = -ai*cdlog(cdexp(ai*thet(k)))
                                  !WRITE(206,*)"II", t,"THET(",K,")=",THET(K)
                                   end do
                                  !print*,"sumthett",thet(1)+thet(2)+thet(3)
    !Call REARR(NN,dft,PDF,NT,MDIF,n_lambda,n_thet,lambda,thet,n_A,A,t,phi)
                 Do l = 1,NN
                 Do m = 1,NN
                 dft(l,m) = 0.d0
                 dft(l,m) = dft(l,m) + dabs(phi(l) - thet(m))
                 !print*,t,phi(l),thet(m),dft(l,m)
                 WRITE(206,*) t,"difthet",l,m, dft(l,m)
                 end do
                 end do

                 PDF(1) =   dsqrt((dft(1,1))**2 + (dft(2,2))**2 + (dft(3,3))**2)
                 PDF(2) =   dsqrt((dft(1,1))**2 + (dft(2,3))**2 + (dft(3,2))**2)  
                 PDF(3) =   dsqrt((dft(1,2))**2 + (dft(2,3))**2 + (dft(3,1))**2)
                 PDF(4) =   dsqrt((dft(1,2))**2 + (dft(2,1))**2 + (dft(3,3))**2)
                 PDF(5) =   dsqrt((dft(1,3))**2 + (dft(2,1))**2 + (dft(3,2))**2)
                 PDF(6) =   dsqrt((dft(1,3))**2 + (dft(2,2))**2 + (dft(3,1))**2)
 
                  
                   MDIF = MIN( PDF(1),PDF(2),PDF(3),PDF(4),PDF(5),PDF(6))
                   
                   if (MDIF .eq. PDF(1)) then
                    NT(1) = 1
                    NT(2) = 2
                    NT(3) = 3
                   else if ( MDIF .eq. PDF(2) ) then
                    NT(1) = 1
                    NT(2) = 3
                    NT(3) = 2
                   else if ( MDIF .eq. PDF(3)) then
                    NT(1) = 2
                    NT(2) = 3
                    NT(3) = 1
                  else if ( MDIF .eq. PDF(4)) then
                    NT(1) = 2
                    NT(2) = 1
                    NT(3) = 3
                 else if ( MDIF .eq. PDF(5)) then
                    NT(1) = 3
                    NT(2) = 1
                    NT(3) = 2
                 else if ( MDIF .eq. PDF(6)) then
                    NT(1) = 3
                    NT(2) = 2
                    NT(3) = 1
                 end if
                 !print*,t,PDF(1),PDF(2),PDF(3),PDF(4),PDF(5),PDF(6)
                 !print*,t,MDIF,NT(1),NT(2),NT(3)


                 do l = 1,NN
                    n_lambda(l) = lambda(NT(l))
                 end do
                 do l = 1,NN
                    n_thet(l) = thet(NT(l))
                 end do
                 DO l = 1,NN
                 DO m = 1,NN
                    n_A(m,l) = A(m,NT(l))
                 end do
                 END DO

                 !newly arranged thet and eigen vectors------------
                 do l = 1,NN
                 lambda(l) = n_lambda(l)
                 end do
                 do l = 1,NN
                 thet(l) = n_thet(l)
                 !write(207,*) t,NT(l), "thet-rearr",l,thet(l)
                 !print*,t,thet(l)
                 end do
                 write(207,*) thet(1),thet(2),thet(3)
                 do l = 1,NN
                 do m = 1,NN
                 A(l,m) = n_A(l,m)
                 !print*,l,m,"mat-A",A(l,m)
                 end do
                 end do



!record values--------------------------------------------------

                           If (t .eq. itr) then
                                do p = 1,NN
                                do q = 1,NN
                                WRITE(320,*) A(p,q)
                                END DO
                                END DO

                                do p = 1,NN
                                WRITE(420,*) thet(p)
                                END DO

                               do k = 1,NN
                               print*,"thet=",thet(k),"phi",phi(k)
                               end do

                           END IF

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
                        U_D(l,p,q) =  U_D(l,p,q) -( (ai)*DEL(p,l)*UU(p,q))
                        !write(109,*) l,p,q,U_D(l,p,q)
                             end Do
                             end do
                             End do
                          
 !CALL UDDAS(t,NN,DEL,UU,U_DD)
                       DO l = 1,NN
                       DO m = 1,NN
                       Do p= 1,NN
                       DO q = 1,NN
                        U_DD(l,m,p,q) = dcmplx(0.d0,0.d0)
                        U_DD(l,m,p,q) =  U_DD(l,m,p,q) -( DEL(p,l)*DEL(p,m)*UU(p,q))
                        !write(200,*) l,m,p,q,U_DD(l,m,p,q)
                             end Do
                             end do
                             End do
                             END DO
 ! CALL LLAMB(t,NN,dlamb,A,U_D)
                       Do l = 1, NN
                       Do k = 1, NN
                       dlamb(k,l) = dcmplx(0.d0,0.d0)
                       Do p = 1,NN
                       Do q =  1,NN
                       dlamb(k,l) = dlamb(k,l) + (dconjg(A(p,k))*U_D(l,p,q)*A(q,k))
                       End do
                       end do
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
                       !write(*,*)"cof",l,k,q, COF(l,k,q)
                       END do
                       END DO
                       end do
 !CALL GLAMBDA(t,NN,lambda,ddlamb1,ddlamb2,ddlamb,A,U_DD,COF)
                        !dlamb1----------------------------------------            
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
                  !dlamb2-----------------------------------------------
                       DO k = 1,NN
                       Do l = 1,NN
                       Do m = 1,NN
                       ddlamb2(k,l,m) = dcmplx(0.d0,0.d0)
                          Do p = 1,NN
                           If (p.ne.k) then
                            ddlamb2(k,l,m) =  ddlamb2(k,l,m) + ((COF(l,p,k)*COF(m,k,p)&
                                           & + COF(m,p,k)*COF(l,k,p))*(lambda(p)-lambda(k)))
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


 !CALL GTHETA(t,NN,bi,dtheta,ddtheta,dlamb,ddlamb,lambda)
                     !dtheta-----------------------------------------------------    
                       bi = dcmplx(0.d0,1.d0)
                       DO k = 1,NN
                       Do l = 1,NN
                       dtheta(k,l) = 0.d0
                       dtheta(k,l) = dtheta(k,l) +  bi*dconjg(lambda(k))*dlamb(k,l)
                       !write(205,*)k,l, "dtheta=",dtheta(k,l)
                       end do
                       end do                     
                  !ddtheta--------------------------------------------------------     
                       DO k = 1,NN
                       Do l = 1,NN
                       DO m = 1,NN
                       ddtheta(k,l,m) = 0.d0
                       ddtheta(k,l,m) = ddtheta(k,l,m) -(bi*((dconjg(lambda(k)))**2)*(dlamb(k,m)*dlamb(k,l))) &
                                         & + (bi*(dconjg(lambda(k))*ddlamb(k,l,m)))
                       !write(205,*)k,l,m,"ddtheta=",ddtheta(k,l,m)
                       end do
                       end do
                       END DO 



                   




  !ACTION-PART!!
  !CALL NEBLA_S(t,NN,cc,ss,ds,thet,phi,dtheta,DEL)
                              !cc---------------------------------------------------
                                      DO i = 1,NN
                                      cc(i) = 0.d0
                                      cc(i) = cc(i) + (1/(thet(i)-phi(i)))
                                      END DO 
                              !ss-----------------------------------------------------       
                                      DO i = 1,NN
                                      ss(i) = 1.d0
                                      ss(i) = ss(i)*(thet(i)-phi(i))**2
                                      END DO 
                              !ds---------------------------------------------------------     
                                      Do i = 1,NN
                                      Do l = 1,NN
                                      ds(i,l) = 0.d0
                                      ds(i,l) = ds(i,l) + (2.d0*(cc(i)*ss(i)*(dtheta(i,l) - DEL(i,l))))
                                      end do 
                                      end do
                                     
 
                           
                        
                  !DCC----------------------------------------------------------     
                         DO i = 1,NN
                         DO m = 1,NN
                           dcc(i,m) = 0.d0
                           dcc(i,m) = dcc(i,m) - ((cc(i))**2)*(dtheta(i,m)-DEL(i,m))
                         END DO
                         END DO
                
              !(DDS1)------------------------------------------------
                         DO i = 1,NN
                         DO l = 1,NN
                         DO m = 1,NN
                           DDS1(i,l,m) = 0.d0
                           DDS1(i,l,m) = DDS1(i,l,m) + (2.d0*ds(i,m)*cc(i)*(dtheta(i,l)))
                         END DO
                         END DO
                         END DO
              !(DDS2)----------------------------------------------------
                         DO i = 1,NN
                         DO l = 1,NN
                         DO m = 1,NN
                           DDS2(i,l,m) = 0.d0
                           DDS2(i,l,m) = DDS2(i,l,m) + (2.d0*SS(i)*DCC(i,m)*(dtheta(i,l)))
                         END DO
                         END DO
                         END DO
             !(DDS3)-----------------------------------------------------------------
                         DO i = 1,NN
                         DO l = 1,NN
                         DO m = 1,NN
                           DDS3(i,l,m) = 0.d0
                           DDS3(i,l,m) = DDS3(i,l,m) + (2.d0*ss(i)*CC(i)*(ddtheta(i,l,m)))
                         END DO
                         END DO
                         END DO
             !(DDS4)---------------------------------------------------------------------
                         DO i = 1,NN
                         DO l = 1,NN
                         DO m = 1,NN
                           DDS4(i,l,m) = 0.d0
                           DDS4(i,l,m) = DDS4(i,l,m) +((2.d0*DS(i,m)*cc(i)*DEL(i,l)) + (2.d0*ss(i)*dcc(i,m)*DEL(i,l)))
                         END DO
                         END DO
                         END DO
           !(DDS)-------------------------------------------------------------------------
                         DO i = 1,NN
                         DO l = 1,NN
                         DO m = 1,NN
                         DDS(i,l,m) = 0.d0 
                         DDS(i,l,m) = DDS(i,l,m) + (DDS1(i,l,m)+ DDS2(i,l,m)+ DDS3(i,l,m) - DDS4(i,l,m))
                         END DO
                         END DO
                         END DO
            !H-matrix-----------------------------------------------------------------------
                        DO p = 1,NN
                        DO q = 1,NN
                        DSD(p,q) = 0.d0
                          Do i = 1,NN
                          DSD(p,q) = DSD(p,q) + DDS(i,p,q)
                          end do
                        !WRITE(*,*) "DSD(",p,q,")=", DSD(p,q)
                        end do
                        end do
!CALL SNGREDS(NN,S,ss,NEBS,ds,grad_s,t)
                     !S-----------------------------------------------
                     S = 0.d0
                     Do i = 1,NN
                     S = S + ss(i)
                     END DO
                     !WRITE(*,*) t, "S=",S 
                     Write(111,*) t,S
             !NEB-S----------------------------------
                     Do k = 1,NN
                     NEBS(k) = 0.d0
                     Do l = 1,NN
                     NEBS(k) = NEBS(k) + ds(l,k)
                     END DO
                     !WRITE(*,*) "NEBS(",K,")=",NEBS(K) 
                     !WRITE(401,*) NEBS(K)
                     END DO
                     grad_s = 0.d0
                     DO P = 1,NN
                     grad_s = grad_s + NEBS(p)**2
                     end do
                     grad_s = dsqrt(grad_s)
                     !print*,"grad_s=", t,grad_s
                     write(112,*) t,grad_s
 !CARTAN'S METHOD OF TRANSFORMATION
 !CALL BNC(t,NN,DSD,NEBS,NDSD,NNS,HM,NNEBS,THM,NA,NNA,AAA)
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
        !PRINT*, p,q,HM(P,Q)*sqrt(float(p*(p+1)))
        !PRINT*, p,q,"hm=",HM(P,Q)
     End do
     END DO


     !NEW-NEBS-----------------------------------------------------------------
    Do  p = 1,NN
    NNEBS(p) = 0.d0
    Do l= 1,NN
     NNEBS(p) = NNEBS(p) + HM(p,l)*NEBS(l)
    end do
    end do
    
!A-----------!REGAIN NEBS=C---------------------------
       DO k= 1,NN-1
        NNS(k)  = NNEBS(k)
       !PRINT*,"projected-NNS(",k,")",NNS(k)
       END DO


    !NEW-NA-NNA-------------------------------------------------------------------------
     DO p = 1,NN
     DO q = 1,NN
     THM(p,q) = HM(q,p)
     end do 
     end do


   HMTHM = matmul(HM,THM)

   DO p = 1,NN
   DO q = 1,NN
   !print*,p,q,"hmthm",HMTHM(p,q)
   end do
   end do

   !NA= MATMUL(A,HM)---------
    
     DO p = 1,NN
     DO q = 1,NN
       NA(p,q) = 0.d0
       do l = 1,NN
       NA(p,q) =NA(p,q) + (DSD(p,l)*THM(l,q))
       !print*,"NA=",p,q,NA(p,q)
       end do
    end do 
    end do
    
    !A = MATMUL(THM,NA)------------------------------------------------------------
   
    DO p = 1,NN
    DO q = 1,NN
       NNA(p,q) = 0.d0
       do l = 1,NN
       NNA(p,q) =NNA(p,q) + (HM(p,l)*NA(l,q))
       end do
    !print*,"projected-NEW-A=",p,q,NNA(p,q)
    !write(50,*)p,q, NNA(p,q)
    end do 
    end do
    


!B----------------! Regain matrix- DSD=B----------------------------------------------


   
    DO p = 1,NN-1
    Do q = 1,NN-1
    NDSD(p,q) = NNA(p,q)
   ! print*,"projected-NDSD(",p,q,")",NDSD(p,q)
    END DO
    END DO
    DO p = 1,NN-1
    Do q = 1,NN-1
    AAA(p,q) = NNA(p,q)
    END DO
    END DO
                       !lapack inversion of DSD
                        CALL DGETRF( NN-1,NN-1,NDSD,LDAA,IPIV,INFO)
                        If (info .eq.0) then
    	                LWORK = -1
     		        CALL DGETRI( NN-1, NDSD, LDAA, IPIV, WORK,LWORK, INFO )
      		        LWORK = min( LWMAX, INT( WORK( 1 ) ) )
     		        CALL DGETRI( NN-1, NDSD, LDAA, IPIV, WORK, LWORK, INFO )
     		        if (info .ne. 0) then
     		        PRINT*, 'Matrix inversion failed!'
     		        end if
     		        end if
                        !REVERSE TRANSFORM AND NEW-PHI
                        
!CALL FNLPHI (t,NN,IDSD,NDSD,NBPHI,NNPHI,NNS,RTPHI,THM,new_phi,phi,DD,trace,ITR,III,AAA)
                        DO l = 1,NN-1
     		Do m = 1,NN-1
     		IDSD(l,m) = NDSD(l,m)
     		!WRITE(*,*)"inversed-IDSD(",l,m,")", IDSD(l,m)
     		End do
     		End do
                III = matmul(IDSD,AAA)
                !CHECK INVERSE
                DO l = 1,NN-1
     		Do m = 1,NN-1
     	        !WRITE(*,*)"Identity",III(l,m)
     		End do
     		End do
               

                

     !NEB-PHI--------------------------------------------------------

     		DO p = 1,NN-1
     		NBPHI(p)  = 0.d0
     		DO i= 1,NN-1
     		NBPHI(p)  = NBPHI(p) + IDSD(p,i)*NNS(i)
     		end do
                !print*,"prod-inverseH*nebphi-NBPHI(",p,")",NBPHI(p)
     		end do
    !Final new-phi---------------------------------------------------
                Do p = 1,NN
                   if (p .lt.NN) then 
                   NNPHI(p) = NBPHI(p)
                   else if (p .eq. NN) then
                   NNPHI(p) = 0.d0
                   END IF
                   !print*,"NNPHI",NNPHI(p)
                end do
                  !PRINT*,"new", NNPHI(1)+ NNPHI(2)+ NNPHI(3)
      
     !RTPHI =  matmul(THM,NNPHI)
                DO l = 1,NN
     		RTPHI(l)  = 0.d0
     		DO i= 1,NN
     		RTPHI(l)  = RTPHI(l) + THM(l,i)*NNPHI(i)
     		end do
               !print*,"reverse-transformed-final-neb-PHI(",l,")",RTPHI(l)
     		end do
                !PRINT*,"new", RTPHI(1)+ RTPHI(2)+ RTPHI(3)
                do p = 1,NN
                new_phi(p) = 0.d0
                new_phi(p) = new_phi(p) + (phi(p) - RTPHI(p))
                end do
                !print*,"1new_phi",(new_phi(1)+new_phi(2)+new_phi(3))

                do p = 1,NN
                new_phi(p) = -ai*cdlog(cdexp(ai*new_phi(p)))
                end do
                !print*,"2new_phi",new_phi(1)+new_phi(2)+new_phi(3)
                

     !redefine new set of phi----------------------------------------
                             do p = 1,NN
                             phi(p) = 0.d0
                             phi(p) = phi(p) + new_phi(p)
                             write(206,*) t,"phi=",phi(p)
                             end do
                             !print*,phi(1),phi(2),phi(3)
                            

			     Do p = 1,NN
                             DO q = 1,NN
                                DD(p,q) = DCMPLX(0.d0,0.d0)
                                IF (p .EQ. q) THEN
                                DD(p,q) = DD(p,q) + cdexp(-ai*new_phi(p))
                                ELSE IF (p .NE. q) THEN
                                DD(p,q) = DD(p,q) + DCMPLX(0.d0,0.d0)
                                end if
                                !WRITE(103,*) DD(p,q)
                             end do
                             end do
                 
                                
                               

                                Write(555,*) t, phi(1),phi(2),phi(3)
                                trace =  phi(1) + phi(2) + phi(3)
                                print*,t,"trace=",trace
                         t  = t+1
                        if (t .lt.ITR+1) then
                        GO TO 25
                        else if (t .eq. ITR+1) then
                        go to 28
                        end if

                        28 END PROGRAM STPDC !**************