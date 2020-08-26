program UNTM
implicit none
integer ,parameter::N = 10, D = 3,M=1
integer::i,j,k,l,p,q
real*8,dimension(N)::t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18
complex*16,dimension(M,D)::X,Y,Z,U,V,W,UU,VV,WW
Real*8,dimension(M):: NX,NNX,NY,NNY,NZ,NNZ,NU,NNU,NV,NNV,NW,NNW,YU,ZU,ZV
COMPLEX*16::TMU(M,D,D),MAT_U(M,D,D),II(M,D,D),DUU(M)



		Do i = 1,M

		call random_number(harvest = t1(i))
		call random_number(harvest = t2(i))
		call random_number(harvest = t3(i))
		call random_number(harvest = t4(i))
		call random_number(harvest = t5(i))
		call random_number(harvest = t6(i))
		call random_number(harvest = t7(i))
		call random_number(harvest = t8(i))
		call random_number(harvest = t9(i))
		call random_number(harvest = t10(i))
		call random_number(harvest = t11(i))
		call random_number(harvest = t12(i))
		call random_number(harvest = t13(i))
		call random_number(harvest = t14(i))
		call random_number(harvest = t16(i))
		call random_number(harvest = t17(i))
		call random_number(harvest = t18(i))

                   !first set of vectors------------------------
                  X(i,1) = dcmplx(t1(i),t2(i))
                  X(i,2) = dcmplx(t3(i),t4(i))
                  X(i,3) = dcmplx(t5(i),t6(i))
                    do k = 1,D
                    write(1,*) X(i,k)
                    end do
                  Y(i,1) = dcmplx(t7(i),t8(i))
                  Y(i,2) = dcmplx(t9(i),t10(i))
                  Y(i,3) = dcmplx(t11(i),t12(i))
                    do k = 1,D
                    write(2,*) Y(i,k)
                    end do
                  Z(i,1) = dcmplx(t13(i),t14(i))
                  Z(i,2) = dcmplx(t15(i),t16(i))
                  Z(i,3) = dcmplx(t17(i),t18(i))
                    do k = 1,D
                    write(3,*) Z(i,k)
                    end do
                  
                    !normalization-------------------------------            
                    NX(i) =  dot_product(X(i,1:D),X(i,1:D))
                    NNX(i)= dsqrt(NX(i))
                    NY(i) =  dot_product(Y(i,1:D),Y(i,1:D))
                    NNY(i)= dsqrt(NY(i))
                    NZ(i) =  dot_product(Z(i,1:D),Z(i,1:D))
                    NNZ(i)= dsqrt(NZ(i))

                    !NEW SET OF VECTORS------------
                               
                    !VECTOR- U-------------------------------------------------
                    DO K = 1,D
                    U(i,k) = X(i,k)
                    end do
                    NU(i) = dot_product(U(i,1:D),U(i,1:D))

                    NNU(i) = dsqrt(NU(i))

                    do k = 1,D
                    UU(i,k) = U(i,k)/NNU(i)
                    END DO

                    

                    !VECTOR -V----------------------------------------------------
                    YU(i) = dot_product(Y(i,1:D),U(i,1:D))/NU(i)
                    Do k = 1,D
                    V(i,k) = dcmplx(0.d0,0.d0)
                    V(i,k) = V(i,k) + (Y(i,k) - (YU(i)*U(i,k)))
                    end do
                    
                    NV(i) = dot_product(V(i,1:D),V(i,1:D))
                    NNV(i) = dsqrt(NV(i))
                    Do k = 1,D
                    VV(i,k) = V(i,k)/NNV(i)
                    end do
                   
                   
                    !VECTOR -W-------------------------------------------------------
                    ZU(i) = dot_product(Z(i,1:D),U(i,1:D))/NU(i)
                    ZV(i) = dot_product(Z(i,1:D),V(i,1:D))/NV(i)
                    Do k = 1,D
                    W(i,k) = dcmplx(0.d0,0.d0)
                    W(i,k) = W(i,k) + (Z(i,k) - (ZU(i)*U(i,k))&
                                  & - (ZV(i)*V(i,k)))
                    END DO
                   
                    NW(i) = dot_product(W(i,1:D),W(i,1:D))
                    NNW(i) = dsqrt(NW(i))
                    DO k = 1, D
                    WW(i,k) = W(i,k)/NNW(i)
                    end do
                   

                    !MATRIX -U-------------------------------------
                    Do q = 1,D
                    MAT_U(i,q,1) = UU(i,q)
                    MAT_U(i,q,2) = VV(i,q)
                    MAT_U(i,q,3) = WW(i,q)
                    end do
                    !print U
                   Do p = 1,D
                   DO q = 1,D
                     WRITE(320,*)Mat_U(i,p,q)
                   end do
                   end do
                   !GETTING TRANSPOSE--------
                   Do p = 1,D
                   DO q = 1,D
                     TMU(i,p,q) = DCONJG(MAT_U(i,q,p))
                   end do
                   end do
                   !Unitarity-------------------------
                   !Matmul(TMU,MAT_U)
                   Do p = 1,D
                   DO q = 1,D
                     II(i,p,q) = dcmplx(0.d0,0.d0)
                   Do l = 1,D
                     II(i,p,q) = II(i,p,q) + (TMU(i,p,l)*MAT_U(i,l,q))
                   end do
                   end do
                   end do
                  !print II
                   Do p = 1,D
                   DO q = 1,D
                     PRINT*,II(i,p,q)
                   end do
                   end do

                 END DO


end program UNTM



