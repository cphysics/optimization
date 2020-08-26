   
        
       !f(x,y,z) = sin(x-x0)^2 * sin(y-x0)^2 * sin(z-x0)^2
                                     !+ sin(x-y0)^2 * sin(y-y0)^2 * sin(z-z0)^2
                                     !+ sin(x-z0)^2 * sin(y-z0)^2 * sin(z-z0)^2
                                     !constraint:x+y+z=0



                           



        PROGRAM STPDC
        IMPLICIT NONE
        
        INTEGER,PARAMETER::                NN = 3,ITR = 166
       	Integer::	                   m,mm,i,j,k,l,p,q,t,u,INFO,LWORK,IPIV(NN),NT(NN)          
        Integer,parameter::	           LDA = NN ,LDAA=NN-1, LWMAX = 10000
        REAL(KIND=8)::                     S,grad_s,pi, RWORK( 3*NN-2 ),trace,MDIF
        REAL(KIND=8),DIMENSION(NN)::       C,NNEBS,new_phi,mthet,n_thet,RTPHI,NNPHI,W,NR
        REAL(KIND=8),DIMENSION(NN)::       xx,yy,thet,phi,SS,NBPHI,NPHI,NEBS,n_phi
        REAL(KIND=8),DIMENSION(NN,NN)::    HM,THM,NA,NNA,DC,DSD,dif_thet,dtheta,DS,CC,dft,CKK,DEL
        REAL(KIND=8),DIMENSION(NN,NN,NN):: DDS,DDS1,DDS2,DDS3,DDS4,DDC,ddtheta,DCC
        REAL(KIND=8)::                     NDSD(NN-1,NN-1),IDSD(NN-1,NN-1),NNS(NN-1),PDF(6),AAA(2,2),III(2,2)
        complex*16:: ai,bi,WORK(LWMAX)

        
        

       
                        t = 0
                        pi = dacos(-1.d0)
                        ai = dcmplx(0.d0,1.d0)
                        bi = dcmplx(0.d0,1.d0)
                        
     			  


                               t = 0
                               phi(1) = 0.15
                               phi(2) = 0.30
                               phi(3) = -0.45


                        25     thet(1) = (phi(1)**3)+(phi(2)**2) + phi(3)
                               thet(2) = -(phi(2)**3)+(phi(3)**2) + phi(1)
                               thet(3) = -(phi(3)**3)-(phi(1)**2) + phi(2)

                               do i = 1,NN
                               !print*,t,"phi =",phi(i)
                               end do
                               do i = 1,NN
                               !print*,t,"thet=",thet(i)
                               end do
                        do p= 1,NN
                        do q = 1,NN
                           if (p .eq.q) then
                           DEL(p,q) = 1.d0
                           else if(p .ne. q) then
                           DEL(p,q) =  0.d0
                           end if
                       end do
                       end do 
         


 !CALL GTHETA(t,NN,bi,dtheta,ddtheta,dlamb,ddlamb,lambda)
                     !dtheta-----------------------------------------------------    
                       dtheta(1,1) =3*phi(1)**2
                       dtheta(1,2) =2*phi(2)
                       dtheta(1,3) =1.0
                       dtheta(2,1) =1.0
                       dtheta(2,2) =-3*phi(2)**2
                       dtheta(2,3) =2*phi(3)
                       dtheta(3,1) =-2*phi(1)
                       dtheta(3,2) =1.0
                       dtheta(3,3) =-3*phi(3)**2
                                          
                  !ddtheta--------------------------------------------------------     
                       DO k = 1,NN
                       Do l = 1,NN
                       DO m = 1,NN
                       if (k .eq. 1)  then 
                          if (p .eq. 1 .and. q .eq. 1) then
                           ddtheta(k,l,m) = 6*phi(1)
                          else if (p .eq. 1 .and. q .eq. 1) then  
                           ddtheta(k,l,m) = 2.0
                          else
                           ddtheta(k,l,m) = 0.0
                          end if
                      else if (k .eq. 2)  then 
                          if (p .eq. 2 .and. q .eq. 2) then
                           ddtheta(k,l,m) = -6*phi(2)
                          else if (p .eq. 3 .and. q .eq. 3) then  
                           ddtheta(k,l,m) = 2.0
                          else
                           ddtheta(k,l,m) = 0.0
                          end if
                      else if (k .eq. 3)  then 
                          if (p .eq. 1 .and. q .eq. 1) then
                           ddtheta(k,l,m) = -2.0
                          else if (p .eq. 3 .and. q .eq. 3) then  
                           ddtheta(k,l,m) = -6*phi(3)
                          else
                           ddtheta(k,l,m) = 0.0
                          end if
                     end if
                       end do
                       end do
                       END DO

  !ACTION-PART!!
  !CALL NEBLA_S(t,NN,cc,ss,ds,thet,phi,dtheta,DEL)
                          !cc---------------------------------------------------
                                      DO i = 1,NN
                                      DO k = 1,NN
                                      cc(i,k) = 0.d0
                                      cc(i,k) = cc(i,k) + (1/(dtan(0.5*(thet(k)-phi(i)))))
                                      !write(207,*)"cc(",i,k,") =", cc(i,k)
                                      End do 
                                      END DO 
                              !ss-----------------------------------------------------       
                                      DO i = 1,NN
                                      ss(i) = 1.d0
                                      DO j = 1,NN
                                      ss(i) = ss(i)*(dsin(0.5*(thet(j)-phi(i))))**2
                                      END DO
                                      !write(207,*) "ss(",i,")",ss(i)
                                      END DO 
                              !ds---------------------------------------------------------     
                                      Do i = 1,NN
                                      Do l = 1,NN
                                      ds(i,l) = 0.d0
                                      DO k= 1,NN
                                      ds(i,l) = ds(i,l) + (cc(i,k)*(dtheta(k,l) - DEL(i,l)))
                                      END DO
                                      end do 
                                      end do
                                      do i = 1,NN
                                      do l = 1,NN
                                      ds(i,l) = ss(i)*ds(i,l)
                                      !WRITE(*,*)"ds",i,l, ds(i,l)
                                      end do
                                      end do
 !CALL DDAS(t,NN,cc,dc,c,ss,ds,dds,dds1,dds2,dds3,dds4,dcc,ddtheta,dtheta,DEL,DSD)
                            !C---------------------------------------------------
                         DO i = 1,NN
                         c(i) = 0.d0
                         DO k = 1,NN
                           c(i) = c(i) + cc(i,k)
                         END DO
                         !WRITE(*,*) "C(",I,")=",C(i)
                         END DO
                  !DCC----------------------------------------------------------     
                         DO i = 1,NN
                         DO k = 1,NN
                         DO m = 1,NN
                           dcc(i,k,m) = 0.d0
                           dcc(i,k,m) = dcc(i,k,m) - (0.5)*(1+(cc(i,k))**2)*(dtheta(k,m)-DEL(i,m))
                           !WRITE(*,*) "DCC(",i,k,m,")=",DCC(i,k,m)
                         END DO
                         END DO
                         END DO
                 !DC-------------------------------------------
                         DO i = 1,NN
                         DO m = 1,NN
                           dc(i,m) = 0.d0
                           DO k = 1,NN
                           dc(i,m) = dc(i,m) + dcc(i,k,m)
                           END DO
                         !WRITE(*,*) "DC(",i,m,")=",DC(i,m)
                         END DO
                         END DO
              !(DDS1)------------------------------------------------
                         DO i = 1,NN
                         DO l = 1,NN
                         DO m = 1,NN
                           DDS1(i,l,m) = 0.d0
                           DO k = 1,NN
                           DDS1(i,l,m) = DDS1(i,l,m) + cc(i,k)*(dtheta(k,l))
                           end do
                           DDS1(i,l,m) = DS(i,m)*DDS1(i,l,m)
                         !WRITE(*,*) "DDS1", DDS1(i,l,m)
                         END DO
                         END DO
                         END DO
              !(DDS2)----------------------------------------------------
                         DO i = 1,NN
                         DO l = 1,NN
                         DO m = 1,NN
                           DDS2(i,l,m) = 0.d0
                           DO k = 1,NN
                           DDS2(i,l,m) = DDS2(i,l,m) + DCC(i,k,m)*(dtheta(k,l))
                           end do
                           DDS2(i,l,m) = SS(i)*DDS2(i,l,m)
                         !WRITE(*,*) "DDS2", DDS2(i,l,m)
                         END DO
                         END DO
                         END DO
             !(DDS3)-----------------------------------------------------------------
                         DO i = 1,NN
                         DO l = 1,NN
                         DO m = 1,NN
                           DDS3(i,l,m) = 0.d0
                           DO k = 1,NN
                           DDS3(i,l,m) = DDS3(i,l,m) + CC(i,k)*(ddtheta(k,l,m))
                           end do
                         DDS3(i,l,m) = SS(i)*DDS3(i,l,m)
                         !WRITE(*,*) "DDS3", DDS3(i,l,m)
                         END DO
                         END DO
                         END DO
             !(DDS4)---------------------------------------------------------------------
                         DO i = 1,NN
                         DO l = 1,NN
                         DO m = 1,NN
                           DDS4(i,l,m) = 0.d0
                           DDS4(i,l,m) = DDS4(i,l,m) +((DS(i,m)*c(i)*DEL(i,l)) + (ss(i)*DC(i,m)*DEL(i,l))) 
                         !WRITE(*,*) "DDS4", DDS4(i,l,m)
                         END DO
                         END DO
                         END DO
           !(DDS)-------------------------------------------------------------------------
                         DO i = 1,NN
                         DO l = 1,NN
                         DO m = 1,NN
                         DDS(i,l,m) = 0.d0 
                         DDS(i,l,m) = DDS(i,l,m) + (DDS1(i,l,m)+ DDS2(i,l,m)+ DDS3(i,l,m) - DDS4(i,l,m))
                         !WRITE(399,*) "DDS(",i,l,m,")=", DDS(i,l,m)
                         WRITE(400,*) i,l,m, DDS(i,l,m)
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
                        !write(400,*) DSD(p,q)
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
                        t = t+1
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
                             phi(p) =  new_phi(p)
                             write(206,*) t,"phi=",phi(p)
                             end do
                             !print*,phi(1),phi(2),phi(3)
                            
                                Write(*,*) t, phi(1),phi(2),phi(3)
                                trace =  phi(1) + phi(2) + phi(3)
                                !print*,t,"trace=",trace
                        if (t .lt.ITR+1) then
                        GO TO 25
                        else if (t .eq. ITR+1) then
                        go to 28
                        end if

                        28 END PROGRAM STPDC !************************************************




        
          
           
     
      
                     

            
                   

   


   
                         
                     
















