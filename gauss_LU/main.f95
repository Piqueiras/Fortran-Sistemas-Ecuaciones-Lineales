PROGRAM gauss_ppal

    USE mod_clreal
    IMPLICIT NONE
    INTEGER::n
    REAL(clreal), ALLOCATABLE :: A(:,:),B(:),U(:)

    print*, "======GAUSS LU======"

    read*,n
    print*,"Tamaño: ",n

    ALLOCATE (A(n,n),b(n),u(n))

    CALL lecmat(A,n)
    print*, "Matriz A:"
    CALL prinmat(A,n)

    read*,b
    print*, "Vector de términos independientes:", b

    CALL gauss_LU(n,A)

    print*, "Matriz modificada LU:"
    CALL prinmat(A,n)

    CALL descenso_L(n,A,b,u)
    CALL remonte(n,A,u,b)
    PRINT*, b

    DEALLOCATE(A,b,u)

END PROGRAM gauss_ppal