program homework
  implicit none 
  
  double precision, dimension(:,:), allocatable :: S, Z, X, Sbar, Sbp, Tt, Sp, h, g, Dpq, R,fock, C
  double precision, dimension(:),   allocatable :: A, diag, off_diag, d_sorted, alpha,coeff,  e02, E2, E3
  double precision, dimension(:,:,:,:), allocatable :: pqrs, sum2
  integer,          dimension(:),   allocatable :: d_id
  double precision :: Summ, prod, e0, E, E_old, J11, EI, EIkoop
  double precision, dimension(:,:), allocatable :: coeff2, FRS, SRF, coeff3
  double precision, dimension(:,:,:), allocatable ::fock2
  integer, parameter :: N = 2
  integer :: i, j, k, l, iter		
 
  ! ############## RESULTS FROM SNOW AND BILLS to compare #######################################

  allocate(coeff2(N*10,N), fock2(N*10,N,N), e02(N*10), E2(10*N))
  allocate(coeff3(N*10,N), E3(N*10))
  coeff3(1,1)=1        ; coeff3(1,2)=0
  coeff3(2,1)=0.809249 ; coeff3(2,2)=0.21906   
  coeff3(3,1)=0.847034 ; coeff3(3,2)=0.176952
  coeff3(4,1)=0.839638 ; coeff3(4,2)=0.185241
  coeff3(5,1)=0.84191  ; coeff3(5,2)=0.183615
  coeff3(6,1)=0.840806 ; coeff3(6,2)=0.183934
  coeff3(7,1)=0.840862 ; coeff3(7,2)=0.183871
  coeff3(8,1)=0.840851 ; coeff3(8,2)=0.183884
  coeff3(9,1)=0.840853 ; coeff3(9,2)=0.183881
  E3(1)=-2.83308; E3(2)=-2.86061; E3(3)=-2.86163; E3(4)=-2.86167
  E3(5)=-2.86167; E3(6)=-2.86167; E3(7)=-2.86167; E3(8)=-2.86167
  E3(9)=-2.86167
  coeff2=0.0d0; fock2=0.0d0; e02=0.0d0; E2=0.0d0
  

  allocate( A(N), S(N,N), Z(N,N), X(N,N), Sbar(N,N), Sbp(N,N), Tt(N,N), Sp(N,N),h(N,N),R(N,N), alpha(N),pqrs(N,N,N,N))
  S = 0.0d0 ; Z = 0.0d0 ; X = 0.0d0; Sbar=0.0d0; Sbp=0.0d0
  Allocate(Dpq(N,N), C(N,N), fock(N,N), coeff(N))
  !
  allocate( diag(N), off_diag(N), d_sorted(N), FRS(N,N), SRF(N,N))
  diag = 0.0d0 ; off_diag = 0.0d0 ; d_sorted = 0.0d0; 
  !
  allocate( d_id(N) )

  A(1)=1.45d0 ; A(2)=2.9d0 

  ! ######### the overlap matrix before diagolanization ################

  do i=1, N 
    do j=1,N
       S(i,j)=(2*sqrt(A(i)*A(j))/(A(i)+A(j)))**3
    end do
       S(i,i)=1
  end do   
       
  write(*,*)
  write(*,*) 'Overlap matrix before diagonalization:'
  write(*,*)
  do i = 1, N
     write(*,'(100f16.8)') (S(i,j), j = 1, N)
  end do
 
  Z=S



!######### Diagonalizing S and COMPUTING S-1/2   #######################
  call tred2(Z, N, N, diag, off_diag)
  call tqli(diag, off_diag, N, N, Z)
  call bubble_sort_mat( diag, N, d_sorted, d_id, Z)
   X = Z

   do i=1, N 
     Sbar(i,i)=diag(i)
  end do   
   
  !
  write(*,*) 'sorted eigenvalues:'
  !
  do i = 1, N
     write(*,'(f16.8,X,I4)') diag(i), d_id(i)
  end do
  !
  write(*,*)
  write(*,*) 'sorted eigenvectors (column format):'
  !
  open(unit=1, file='eigenvectors.txt', status='unknown')
  do j = 1, N
     write(*,'(100f16.8)') (X(j,i), i = 1, N)
     write(1,'(100f16.8)') (X(j,i), i = 1, N)
  end do
  close(1)


  write(*,*)
  write(*,*) 'diagonal matrix with eigenvalues S̄ :'
  write(*,*)
  do i = 1, N
     write(*,'(100f16.8)') (Sbar(i,j), j = 1, N)
  end do


 !################## Sbar power -1/2 #################@
  do i=1 ,N
    Sbp(i,i)=1.0d0/sqrt(abs(Sbar(i,i)))
  end do

  Tt=transpose(X)
  
  Sp=matmul(X,matmul(Sbp,Tt))
  
  write(*,*)
  write(*,*) 'S**(-1/2) :'
  write(*,*)
  do i = 1, N
     write(*,'(100f16.8)') (Sp(i,j), j = 1, N)
  end do


! ########### hpq and (pq/rs) ###############
  do i=1, N
     do j=1, N
        h(i,j)=4.0d0*((sqrt(A(i)*A(j))/(A(i)+A(j)))**3)*((A(i)*A(j))-2.0d0*((A(i)+A(j))))
     end do
  end do
  write(*,*)
  write(*,*) 'hpq Matrix:'
  write(*,*)
  do i = 1, N
     write(*,'(100f16.8)') (h(i,j), j = 1, N)
  end do

!######## alpha as array ###########
  alpha=A
  do i=1, N
     do j=1, N
       do k=1, N 
         do l=1,N 
            summ=alpha(i)+alpha(j)+alpha(k)+alpha(l)
            prod=alpha(i)*alpha(j)*alpha(k)*alpha(l)
            pqrs(i,j,k,l)= -((32.0d0 * prod**1.5) / ((alpha(i)+ alpha(j))**3 * summ**2)) & 
                  + ((32.0d0 * prod**1.5) / ((alpha(i) + alpha(j))**3 * (alpha(k) + alpha(l))**2)) &
                  - ((32.0d0 * prod**1.5) / ((alpha(i) + alpha(j))**2 * (summ**3)))
         end do 
      end do    
     end do  
  end do

  write(*,*)
  write(*,*) '(pq/rs) :'
  write(*,*)
  do i = 1, N
    do j=1, N
      do k=1,N
         do l=1, N
            write(*,*) i, j, k, l
             write(*,'(100f16.8)') pqrs(i,j,k,l)
  
         end do
      end do
    end do
  end do 


!#############################@ EXERCISE 2      ##################################@
  E=0
  E_old=2.0d0
  coeff=0.0d0
  coeff(1)=1.0D0
  allocate(sum2(N,N,N,N),g(N,N))
  iter=0
  do while (abs(E-E_old)>10e-8)       !!!!!!! START OF LOOP ##################@ 
    E_old=E
    iter=iter+1
    write(*,*) '##########################################################################################'
    write(*,*) 'RESULTS OF ITERATION NUMBER :      ', iter
    sum2=0.0d0
    g=0.0d0
    
    

    write(*,*)
        write(*,*) 'coefficients '
    do i=1, N
        write(*,'(100f16.8)') coeff(i)
    end do
    write(*,*)   



   

    do i=1, N
      do j=1, N
        do k=1, N
          do l=1,N
              g(i,j)= g(i,j)+ (coeff(k) * coeff(l) *pqrs(i,j,k,l))         
          end do
        end do
      end do
    end do
    

 
    write(*,*) 'Fock Matrix '
    do i=1, N
      do j=1,N
        fock(i,j)=h(i,j)+g(i,j)
      end do
      write(*,'(100f16.8)') (fock(i,j), j=1,N)
    end do

    if (iter==1) then
       J11=g(1,1)
    end if 

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N
      coeff2(iter,i)=coeff(i)
     end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    

! ######################### EIGEN VECTOR of F ##############################@
    Z = 0.0d0 ; X = 0.0d0; Sbar=0.0d0;
    diag = 0.0d0 ; off_diag = 0.0d0 ; d_sorted = 0.0d0
    Z=matmul(Sp,matmul(fock,Sp))
    call tred2(Z, N, N, diag, off_diag)
    call tqli(diag, off_diag, N, N, Z)
    call bubble_sort_mat( diag, N, d_sorted, d_id, Z)
    X = Z
    e0=minval(diag)
  
    write(*,*) 'sorted eigenvalues of fock matrix:'
  
    do i = 1, N
       write(*,'(f16.8,X,I4)') diag(i), d_id(i)
    end do
  
    write(*,*)
    write(*,*) 'sorted eigenvectors of fock matrix, C_dash (column format):'
  
    open(unit=1, file='eigenvectors.txt', status='unknown')
    do j = 1, N
      write(*,'(100f16.8)') (X(j,i), i = 1, N)
      write(1,'(100f16.8)') (X(j,i), i = 1, N)
    end do
    close(1)

 !  Finding C= S**-0.5 * C'
    C=matmul(Sp,X)
  
    write(*,*) 'lowest eigen value:   ', e0
!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N
      coeff(i)=C(i,1)
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!
    summ=0.0d0
 !  Desnity matrix D=2R, R=C(i,1)*C(j,1)
    do i = 1, N
      do j= 1, N
         R(i,j) = coeff(i)*coeff(j)
         summ= summ +((R(i,j)*2)*S(i,j))
      end do 
    end do
    R=2*R
! printing density matrix
    write(*,*) 'Desnity matrix' 
    do i=1,N  
      write(*,'(100f16.8)') (R(i,j), j=1,N)    
    end do
  
    write(*,'(A, F8.2,A, I3)') ' equation (8)= ' , summ, '  and should be equal to    N=' ,N

    
    
    E=0.0d0
    do i=1,N
      do j=1,N
       E= E + coeff(i)*coeff(j)*(2*h(i,j)+ g(i,j))  
      end do
    end do
   
 
    
   
    write(*,*)
    write(*,*) 'E value is'
    write(*,'(100f16.8)') E
    if (abs(E)>10.0d0) then
       write(*,*) ' error in convergence '
       exit
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N
      e02(iter)=e0
      E2(iter)=E
      do j=1,N
        fock2(iter,i,j)=fock(i,j)
      end do
     end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  
  end do        !!!!!!!! END WHILE LOOP. ####################################
 
  ! #####################@ printing all results ##################
  write(*,*) 'iter    c1             c2             F11             F12             F22             min(e)             E'  
  
  do i=1,iter
     write(*,'(I3,100f16.8)') i,coeff2(i,1), coeff2(i,2),  fock2(i,1,1), fock2(i,1,2), fock2(i,2,2),e02(i), E2(i)
  end do 
  

  ! ##############@ The deviation % from Snow and Bills data  ###############

  write(*,*)
  write(*,*) "Results compared to Snow and Bills' data"
  write(*,*) 'iter    c1 deviation%         c2 deviation%         E deviation%' 
  do i=1,iter
     coeff3(i,1)= abs((coeff2(i,1)-coeff3(i,1))*100/coeff2(i,1))
     coeff3(i,2)=abs((coeff2(i,2)-coeff3(i,2))*100/coeff2(i,2))
     E3(i)= abs((E2(i)-E3(i))*100/E2(i))
     write(*,'(I3,100f16.8)') i,coeff3(i,1), coeff3(i,2), E3(i)
  end do 

  E=2*h(1,1)+J11
  write(*,*)
  write(*,*) '          EXERCISE A-3(a)            '
  write(*,*) ' E using (E=2I1 + J11) =  ', E

  write(*,*)
  write(*,*) 'Density matrix EXERCISE A-3(b)'
  do i=1,N
    do j=1,N
       R(i,j)=coeff(i)*coeff(j)
    end do
  end do

   
  do i=1,N
      write(*,'(100f16.8)') (R(j,i), j = 1, N)
  end do


  write(*,*)
  write(*,*) ' expressing R=RSR  '
  Tt=matmul(R,matmul(S,R))
  do i=1,N
      write(*,'(100f16.8)') (Tt(j,i), j = 1, N)
  end do

  write(*,*)
  write(*,*) 'the commutation relation is expressed as: FRS = SRF'
  FRS=matmul(Fock,matmul(R,S))
  SRF=matmul(S,matmul(R,Fock))
  write(*,*) ' FRS = '
  do i=1,N
      write(*,'(100f16.8)') (FRS(j,i), j = 1, N)
  end do
  write(*,*) ' SRF = '
  do i=1,N
      write(*,'(100f16.8)') (SRF(j,i), j = 1, N)
  end do


  E=2.0d0*e02(1) - J11
  write(*,*) 
  write(*,*) 'the total Energy using E=2e1-J11 is ', E
  EI=0.0d0
  do i=1,N
    do j=1,N
    EI= EI + coeff(i)*coeff(j)*(h(i,j))  
    end do
  end do

  write(*,*) 
  write(*,*) 'total energy compared to experimental results have deviation of %', abs((-2.904-E)*100/2.904)
  EI=EI-E

  EIkoop= -e0
  write(*,*)
  write(*,'(A,F8.4)') 'first ionization energy of He  ', EI
  write(*,'(A,F8.4)') 'EI using koopman-theorem ', EIkoop
  write(*,'(A,F8.2)')  'difference % compared to experimental EI=0.904 =', abs((EI-0.904)*100/0.904d0)
  write(*,'(A,F8.2)')  'difference % compared to Koopmans EI',  abs((EI-EIkoop)*100/EIkoop)
END program homework 

   



!##################### SUBROUTINES ######################

SUBROUTINE TQLI(D,E,N,NP,Z)

  IMPLICIT NONE

  INTEGER                             :: N, NP
  DOUBLE PRECISION, DIMENSION(NP, NP) :: Z
  DOUBLE PRECISION, DIMENSION(NP)     :: D, E

  INTEGER          :: I, ITER, K, L, M
  DOUBLE PRECISION :: B, C, DD, F, G, P, R, S

  IF ( N .GT. 1 ) THEN
     DO 11 I = 2, N
        E(I-1) = E(I)
11   END DO
     E(N) = 0.0D0
     DO 15 L = 1, N
        ITER = 0
1       DO 12 M = L, N-1
           DD = ABS(D(M))+ABS(D(M+1))
           IF ( ABS(E(M))+DD .EQ. DD ) GO TO 2
12      END DO
        M=N
2       IF ( M .NE. L) THEN
           IF ( ITER .EQ. 30 ) EXIT
           !IF(ITER.EQ.30) EXIT 'too many iterations'
           ITER = ITER + 1
           G = (D(L+1)-D(L))/(2.0D0*E(L))
           R = SQRT(G**2+1.0D0)
           G = D(M)-D(L)+E(L)/(G+SIGN(R,G))
           S = 1.0D0
           C = 1.0D0
           P = 0.0D0
           DO 14 I = M-1, L, -1
              F = S*E(I)
              B = C*E(I)
              IF ( ABS(F) .GE. ABS(G) ) THEN
                 C = G/F
                 R = SQRT(C**2 + 1.0D0)
                 E(I+1) = F*R
                 S = 1.0D0/R
                 C = C*S
              ELSE
                 S = F/G
                 R = SQRT(S**2 + 1.0D0)
                 E(I+1) = G*R
                 C = 1.0D0/R
                 S = S * C
              ENDIF
              G     = D(I+1) - P
              R     = (D(I)-G)*S + 2.0D0*C*B
              P     = S*R
              D(I+1)= G + P
              G     = C*R - B
              DO 13 K = 1, N
                 F        = Z(K,I+1)
                 Z(K,I+1) = S*Z(K,I) + C*F
                 Z(K,I  ) = C*Z(K,I) - S*F
13            END DO
14         END DO
           D(L)=D(L)-P
           E(L)=G
           E(M)=0.0D0
           GO TO 1
        ENDIF
15   END DO
  ENDIF

  RETURN
END SUBROUTINE TQLI




SUBROUTINE TRED2(A,N,NP,D,E)

  IMPLICIT NONE

  INTEGER                         :: N, NP
  REAL(KIND=8), DIMENSION(NP, NP) :: A
  REAL(KIND=8), DIMENSION(NP)     :: D, E

  INTEGER      :: I, J, K, L
  REAL(KIND=8) :: F, G, H, HH, SCALE

  IF(N.GT.1)THEN
     DO 18 I=N,2,-1
        L=I-1
        H=0.
        SCALE=0.
        IF(L.GT.1)THEN
           DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11         END DO
           IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
           ELSE
              DO 12 K=1,L
                 A(I,K)=A(I,K)/SCALE
                 H=H+A(I,K)**2
12            END DO
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.
              DO 15 J=1,L
                 A(J,I)=A(I,J)/H
                 G=0.
                 DO 13 K=1,J
                    G=G+A(J,K)*A(I,K)
13               END DO
                 IF(L.GT.J)THEN
                    DO 14 K=J+1,L
                       G=G+A(K,J)*A(I,K)
14                  END DO
                 ENDIF
                 E(J)=G/H
                 F=F+E(J)*A(I,J)
15            END DO
              HH=F/(H+H)
              DO 17 J=1,L
                 F=A(I,J)
                 G=E(J)-HH*F
                 E(J)=G
                 DO 16 K=1,J
                    A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16               END DO
17            END DO
           ENDIF
        ELSE
           E(I)=A(I,L)
        ENDIF
        D(I)=H
18   END DO
  ENDIF
  D(1)=0.
  E(1)=0.
  DO 23 I=1,N
     L=I-1
     IF(D(I).NE.0.)THEN
        DO 21 J=1,L
           G=0.
           DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19         END DO
           DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20         END DO
21      END DO
     ENDIF
     D(I)=A(I,I)
     A(I,I)=1.
     IF(L.GE.1)THEN
        DO 22 J=1,L
           A(I,J)=0.
           A(J,I)=0.
22      END DO
     ENDIF
23 END DO

  RETURN
END SUBROUTINE TRED2




subroutine bubble_sort_mat(u_inp, n, u_out, id, mat)
  !
  implicit none
  !
  integer, intent(in)   :: n  ! input variable
  double precision, dimension(n),   intent(in)   :: u_inp  ! input variable to sort
  double precision, dimension(n,n), intent(inout):: mat    ! matrix of eigenvector to sort
  double precision, dimension(n), intent(out)    :: u_out  ! output variable sorted
  integer,          dimension(n), intent(out)    :: id     ! output new indexation
  !
  double precision, dimension(n)   :: u ! working variable
  double precision, dimension(n,n) :: tmp_mat   ! working matrix
  !
  double precision :: tmp
  integer :: i, j, tmp_id
  logical :: swap

  ! initialisation
  u    = u_inp
  swap = .true.
  do i = 1, n
     id(i) = i
  end do
  !
  do while ( swap .eqv. .true. )
     !
     swap = .false.
     !
     loop_do: do i = 1, n - 1
        !
        if ( u(i) > u(i+1) ) then
           !
           ! swap!
           tmp    = u(i  ) ; tmp_id  = id(i)
           u(i  ) = u(i+1) ; id(i)   = id(i+1)
           u(i+1) = tmp    ; id(i+1) = tmp_id
           swap   = .true.
        end if
        !
     end do loop_do

  end do
  !
  u_out = u
  !
  ! Sort eigenstates from mat accordingly using id
  !
  do i = 1, N
     j = id(i)
     tmp_mat(:,i) = mat(:,j)
  end do
  !
  mat = tmp_mat
  !
end subroutine bubble_sort_mat


