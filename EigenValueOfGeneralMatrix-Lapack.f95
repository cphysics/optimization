   
        



        PROGRAM STPDC
        IMPLICIT NONE
        
        INTEGER,PARAMETER::                NN = 4
       	Integer::	                   i,j,k,l,p,q,INFO,LWORK,NT(NN)          
        Integer,parameter::	           LDA = NN , LWMAX = 10000
        REAL(KIND=8)::                     RWORK( 3*NN-2 ),NR(NN),W(NN)
        complex*16::                       WORK(LWMAX),A(LDA,NN),AC(NN,NN),JPT(NN,NN)
        CHARACTER*1::                      UPLO

                          



                          !http://www.nag.com/lapack-ex/node79.html

                         
                   A(1,1:4) =   (/(1.d0,0.d0), (2.d0,-1.d0), (3.d0,-1.d0),  (4.d0,-1.d0)/)
                   A(2,1:4) =   (/(2.d0,1.d0), (2.d0, 0.d0), (3.d0,-2.d0),  (4.d0,-2.d0)/)
                   A(3,1:4) =   (/ (3.d0,1.d0), (3.d0, 2.d0), (3.d0, 0.d0),  (4.d0,-3.d0)/)
                   A(4,1:4) =   (/ (4.d0,1.d0), (4.d0, 2.d0), (4.d0, 3.d0),  (4.d0, 0.d0)/)
                                  
                   WRITE(*,*)     " ! Eigenvalues" 
                   WRITE(*,*)      "-4.2443,     -0.6886 ,     1.1412,      13.7916"

                   WRITE(*,*)        "Eigenvectors"
                   WRITE(*,*) 		"        1       2       3       4 "
 		   WRITE(*,*)           " 1  -0.3839 -0.3975  0.3746 -0.3309"
   		   WRITE(*,*)               "-0.2941  0.5105  0.2414  0.1986"
                       
 		   WRITE(*,*)           " 2  -0.4512  0.3953 -0.2895 -0.3728"
     		   WRITE(*,*)                "0.1102 -0.3238  0.4917  0.2419"
                     
 		   WRITE(*,*)           " 3   0.0263 -0.4309 -0.3768 -0.4870"
     		   WRITE(*,*)               " 0.4857  0.0383 -0.3994  0.1938"
 
                   WRITE(*,*)           " 4   0.5602  0.3648  0.4175 -0.6155"
    		   WRITE(*,*)               "-0.0000  0.0000  0.0000  0.0000 "

                    
                         
                      
                        !LAPACK
                        LWORK = -1
                        CALL ZHEEV( "V", "L", NN, A, LDA, W, WORK, LWORK, RWORK, INFO )
                        LWORK = min( LWMAX, INT( WORK( 1 ) ) )
                        CALL ZHEEV( "V", "L", NN, A, LDA, W, WORK, LWORK, RWORK, INFO )
                      IF( INFO.GT.0 ) THEN
                     WRITE(*,*)'The algorithm failed to compute eigenvalues.'
                     else if(INFO.EQ.0) THEN
                    
                     DO p= 1,NN
                     NR(p) = 0.d0
                     DO q = 1,NN
                     NR(P) = NR(P) + (dconjg(A(q,p))*A(q,p))
                     END DO
                     NR(p) = DSQRT(NR(p))
                     PRINT *,"NR", NR(p)
                     END DO
                     
                     DO p= 1,NN
                     DO q = 1,NN
                     A(q,p) = A(q,p)/NR(p)
                     END DO
                     END DO
                   
                     Do p = 1,NN
                     WRITE(*,*) p,"th eigen value => ",w(p)
                     write(*,*) p,"th-eigen-vector=>"
                     Do q = 1,NN
                     write(*,*) "A(",q,p,")=>",A(q,p)
                     End do
                     end do
                     END IF

           !check--------------------------------------------------
                     do p = 1,NN
                     do q = 1,NN
                     AC(p,q) = dconjg(A(Q,P))
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
                       STOP




                        28 END PROGRAM STPDC !************************************************





    
    
          
