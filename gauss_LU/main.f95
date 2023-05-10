PROGRAM gauss_ppal

    USE mod_clreal
    IMPLICIT NONE
    INTEGER::n,info,i
    REAL(clreal), ALLOCATABLE :: A(:,:),B(:),U(:)
    REAL(clreal)::deter
    real :: start, finish
	
	call cpu_time (start)
    

    print*, "======GAUSS LU======"

    read*, info

    read*,n
    if ( info==1 ) then
        print*,"Tamaño: ",n
    end if
    
	ALLOCATE (A(n,n), b(n), u(n))
	
    CALL lecmat(A,n)
    if ( info==1 ) then
        print*, "Matriz A:"
        CALL prinmat(A,n)
    end if
    
    read*,b
    if ( info==1 ) then
        print*, "Vector de términos independientes:", b
    end if

    CALL gauss_LU(n,A)

    print*, "Matriz modificada LU:"
    CALL prinmat(A,n)

    deter=1
    do i = 1, n
        deter=deter*A(i,i)
    end do

    CALL descenso_L(n,A,b,u)
    CALL remonte(n,A,u,b)
    PRINT*, "Soulución:", b
    print*, "Determinante: ",deter

    DEALLOCATE(A,b,u)
    
    call cpu_time (finish)
	
	print*, 'Tiempo total de calculo: ',finish - start
END PROGRAM gauss_ppal