!  MatrixMultipl.f90 
!
!  FUNCTIONS:
!  MatrixMultipl - Entry point of console application.
!

module procedures
  implicit none

contains

    !Взяла у старших
    
    subroutine multiplication_matrix_vector(matrix, vector, result_vector)
        integer :: i, j, length_m, length_v
        real, dimension(:,:), intent(in) :: matrix
        real, dimension(:), intent(in) :: vector
        real, dimension(:), intent(out) :: result_vector
        
        length_m=size(matrix,1)
        length_v=size(vector,1)
        
        do i=1,length_m
            result_vector(i)=0 
            do j=1,length_v
                result_vector(i)=result_vector(i)+matrix(i,j)*vector(j)
            end do
        end do
    end subroutine multiplication_matrix_vector  
        
    !Сама
    subroutine find_norma(vector, n, l, norma)
        integer, intent(in):: n
        integer :: i, j
        real, intent(in) :: l
        real, dimension(:), intent(in) :: vector
        real :: summ, stepen
        real, intent(out) :: norma
        
        summ=0;
        
        do i=1,n
            summ=summ+vector(i)**l
        end do
        
        norma=summ**(1/l)    
        
        return
    end subroutine find_norma
    
    subroutine AIG2R(A,N)
    ! Обращение вещественной матрицы методом Жордана
    ! с выбором ведущего элемента по всей матрице
    ! http://num-anal.srcc.msu.ru/lib_na/cat/ai/aig2r.htm

    real:: A(:,:)
    integer, intent(in):: N
    real, allocatable:: S(:,:)
    integer F,G,O,P,Q,R,T,U,BB,BD,J,I,M,K,L
    real W,BA,BJ,BN,BQ,BE,X,Y,Z
    real, parameter:: ZERO=0.0
!
    if (N==1) then
        A(1,1)=1./A(1,1); return
    end if
    allocate(S(N,2))
    F=1; G=1
1   continue    
    U=G; W=ZERO
    DO 3 J=G,N
        BB=F; BA=ZERO
        if (F.GT.N) goto 22
        DO 2 I=F,N
            BJ=abs(A(I,J))
            if (BA.GE.BJ) goto 2
            BA=BJ; BB=I
2       continue
22      if (W.GT.BA) goto 3
        W=BA; T=BB; U=J
3   continue
    do I=1,N,1
        BN=A(I,U); A(I,U)=A(I,F); A(I,F)=BN
    end do
    do J=1,N,1
        BQ=A(T,J); A(T,J)=A(F,J); A(F,J)=BQ
    end do
    S(F,2)=U; S(F,1)=T
    J=G
    X=1./A(F,F)
    do I=1,N,1
        A(I,F)=-(A(I,F)*X)
    end do
    A(F,F)=X
    O=0
7   O=O+1
    if (O.EQ.F) goto 7
    if (N.LT.O) goto 10
    BE=A(F,O)
    BD=F-1
    if (BD.EQ.0) goto 23
    do P=1,BD,1
        A(P,O)=A(P,O)+A(P,F)*BE
    end do
23  BD=F+1
    if (BD.GT.N) goto 7
    do P=BD,N,1
        A(P,O)=A(P,O)+A(P,F)*BE
    end do    
    goto 7
10  M=J+1
    if (M.GT.N) goto 24
    do O=M,N,1
        A(F,O)=A(F,O)*A(F,G)
    end do
24  M=J-1
    if (M.EQ.0) goto 25
    do O=1,M,1
        A(F,O)=A(F,O)*A(F,G)
    end do
25  F=F+1; G=G+1
    if (F.GT.N) goto 13
    if (G.LE.N) goto 1
13  J=1
14  Q=ISIGN(1,1-N)
    K=N-Q
15  K=K+Q
    I=S(K,2)
    Y=A(I,J); A(I,J)=A(K,J); A(K,J)=Y
    if (K.NE.1) goto 15
    J=J+1
    if (J.LE.N) goto 14
    I=1
16  R=ISIGN(1,1-N)
    L=N-R
17  L=L+R
    J=S(L,1)
    Z=A(I,J); A(I,J)=A(I,L); A(I,L)=Z
    if (L.NE.1) goto 17
    I=I+1
    if (I.LE.N) goto 16
    deallocate(S)
    return
end subroutine AIG2R

end module

!****************************************************************************
!
!  PROGRAM: MatrixMultipl
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program MatrixMultipl
    use procedures
	use omp_lib
	

    implicit none
    

    ! Variables

    integer :: n, i, j, k

    real, allocatable :: A(:,:), B(:,:), D(:,:), bb(:), x(:), x0(:), xsol(:), &
    vector_difference(:), Bx(:)
    
    real :: summ, norma, e, l
    

    real :: start_time, finish_time
    
    write(*,*) "Enter n"
    read(*,*) n
    
    allocate (A(n,n))
    allocate (B(n,n))
    allocate (D(n,n))
    allocate (bb(n))
    allocate (x(n))
    allocate (x0(n))
    allocate (xsol(n))
    allocate (vector_difference(n))
    allocate (Bx(n))
    
    call cpu_time(start_time)
    
    call random_seed
    call omp_set_nested(.true.)
    !Здаём матрицу А

    CALL RANDOM_NUMBER(A)
    
    print *, 'A = '
    write(*,*) ((A(i,j),j=1,n),i=1,n)
    
    !В конспекте написано, что нужно выечть 0.5
    
    A=A-0.5
    
    summ=0
    
    !Чтобы матрица А наверняка была не вырожденная
    !считаем сумму всех элементов строки (по модулю)
    !и ставим вместо диагонального элемента
    
    do i=1,n
    	    do j=1,n
    	        if (i.ne.j) then
    	            summ=summ+abs(A(i,j))
    	        end if
    	    end do
    	    A(i,i)=summ
    	    summ=0;
    end do
    
    !Рандомим xsol(x_solution - наше решение)
    
    CALL RANDOM_NUMBER(xsol)
    
    !Получаем b (bb; A*xsol=b)
    !НЕ ВЫЗЫВАЕТСЯ
    call multiplication_matrix_vector(A,xsol,bb)
    
    !Забываем про xsol...
    
    !Рандомим первое решение x0
    
    CALL RANDOM_NUMBER(x0)
    
    !Вводим матрицу D (вроде как это главная диагональ матрицы А)
    
    do i=1,n
    	    do j=1,n
    	        if (i==j) then
    	            D(i,j)=A(i,j)
    	        else
    	            D(i,j)=0
    	        endif    
    	    end do
    end do
    
    !Получаем матрицу B (D-B=A)
    
    B=D-A
    
    !Нахождение обратной матрицы D^(-1)
    !НЕ ВЫЗЫВАЕТСЯ
    call AIG2R(D, n)
    
    !Порядок нормы
    l=2
    
    !Точность
    e=0.01
    
    !Начальная норма
    norma=1
    
    do while (norma>=e)
        call multiplication_matrix_vector(B,x0,Bx)
        Bx=Bx+bb
        call multiplication_matrix_vector(D,Bx,x)
        vector_difference=x-x0
        call find_norma(vector_difference, n, l, norma)
        x0=x
    end do
    
    call cpu_time(finish_time)

    print *, 'Calculation time = ', finish_time - start_time
    
    print *, 'x = '
    write(*,*) (x(i),i=1,n)
    
    print *, 'xsol = '
    write(*,*) (xsol(i),i=1,n)

    write(*,*) "Press num button..."
    read(*,*) n
    

    end program MatrixMultipl
    
    
   
