PROGRAM main
    
    USE mod_clreal
    IMPLICIT NONE

    REAL(clreal), DIMENSION(:,:), ALLOCATABLE :: A
    REAL(clreal), DIMENSION(:), ALLOCATABLE :: b,u
    INTEGER :: n
    REAL(clreal) :: eps

    read*, n
    print*, "Tamaño: ",n

    ALLOCATE(A(n,n),b(n),u(n))

    CALL lecmat(A,n)
    print*, "Matriz A: "
    CALL prinmat(A,n)

    read*, b
    print*, "Vector b: ",b

    read*,eps
    print*, "Epsilon: ", eps

    CALL jacobi(n,A,b,u,eps)

    print*, "Resultado: ",u

END PROGRAM main
