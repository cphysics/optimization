   
        !!!!!!!!!!!!!!!!!!!!!!!
        ! This program calculates the eigen value of a 3x3 hermitian matrix using lapack               
        !    Date -2013-November -2 
        !                                                                           
        !!!!!!!!!!!!!!!!!!!!!!!

       

        PROGRAM STPDC
        IMPLICIT NONE
        
        INTEGER,PARAMETER::                NN = 3
        real,parameter::                   LLM = 0.0000000001
       	Integer::	                   m,i,j,k,l,p,q,t,ut,INFO,LWORK         
        Integer,parameter::	           LDA = NN ,LDAA=NN-1, LWMAX = 10000
        REAL(KIND=8)::                     RWORK( 3*NN-2 )
        REAL(KIND=8),DIMENSION(NN)::       W
        complex*16::                       ai,bi,pi
        complex*16::                       WORK(LWMAX),A(LDA,NN)
        CHARACTER*1::                      UPLO

        
                    

              
                        pi = dacos(-1.d0)
                        ai = dcmplx(0.d0,1.d0)
                        bi = dcmplx(0.d0,1.d0)


                      A(1,1) = dcmplx(0,2)
                      A(1,2) = dcmplx(2,3)
                      A(1,3) = dcmplx(3,4)
                      A(2,1) = dcmplx(2,-3)
                      A(2,2) = dcmplx(0,3)
                      A(2,3) = dcmplx(2,4)
                      A(3,1) = dcmplx(3,-4)
                      A(3,2) = dcmplx(2,-4)
                      A(3,3) = dcmplx(0,4)
                       


                      
                       
                        !Lapack-for eigen value
                        LWORK = -1
                        CALL ZHEEV( "V", "L", NN, A, LDA, W, WORK, LWORK, RWORK, INFO )
                        LWORK = min( LWMAX, INT( WORK( 1 ) ) )
                        CALL ZHEEV( "V", "L", NN, A, LDA, W, WORK, LWORK, RWORK, INFO )
                   
                     IF( INFO.GT.0 ) THEN
                     WRITE(*,*)'The algorithm failed to compute eigenvalues.'
                     else if(INFO.EQ.0) THEN
                     Do p = 1,NN
                     WRITE(*,*)t, p,"th eigen value => ",w(p)
                     write(*,*)t,p,"th-eigen-vector=>"
                     Do q = 1,NN
                     write(*,*) t,"A(",q,p,")=>",A(q,p)
                     End do
                     end do
                     END IF
                     do p = 1,NN
                     do q = 1,NN
                     A(q,p) = A(q,p)
                     end do
                     end do

                        
        



                     END PROGRAM STPDC !************************************************