PROGRAM gauss_ppal
    
    USE mod_clreal
    IMPLICIT NONE
    INTEGER::n
    REAL(clreal), ALLOCATABLE :: A(:,:),AA(:,:),B(:),BB(:),U(:),R(:)
    REAL(clreal)::deter

    print*, "======GAUSS LU======"

    read*,n
    print*,"Tamaño: ",n

    ALLOCATE (A(n,n),AA(n,n),b(n),bb(n),u(n),r(n))

    CALL lecmat(A,n)
    print*, "Matriz A:"
    CALL prinmat(A,n)

    read*,b
    print*, "Vector de términos independientes:", b

    CALL gauss(n,A,b,u,deter)
    PRINT*, deter
    PRINT*, u

    DEALLOCATE(A,b,u)

END PROGRAM gauss_ppal