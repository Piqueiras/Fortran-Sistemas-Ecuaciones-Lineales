PROGRAM main
    
    USE mod_clreal
    IMPLICIT NONE

    REAL(clreal), DIMENSION(:,:), ALLOCATABLE :: A
    REAL(clreal), DIMENSION(:), ALLOCATABLE :: b,u
    INTEGER :: n,info
    REAL(clreal) :: eps
    LOGICAL::dominante

    print*, "======JACOBI======"

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

    read*,eps
    print*, "Epsilon: ", eps

    if ( .not. dominante(A,n) ) then
        stop "La matriz no es apta para Jacobi"
    end if

    CALL jacobi(n,A,b,u,eps)

    print*, "Resultado: ",u

END PROGRAM main
