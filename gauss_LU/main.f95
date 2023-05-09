PROGRAM gauss_ppal
    
    USE mod_clreal
    IMPLICIT NONE
    INTEGER::n
    REAL(clreal), ALLOCATABLE :: A(:,:),AA(:,:),B(:),BB(:),U(:),R(:)
    REAL(clreal)::deter

    READ*, n
    ALLOCATE (A(n,n),AA(n,n),b(n),bb(n),u(n),r(n))

    CALL lecmat(A,n)
    CALL prinmat(A,n)
    READ*, b
    PRINT*, b

    CALL gauss(n,A,b,u,deter)
    PRINT*, deter
    PRINT*, u

    DEALLOCATE(A,b,u)

END PROGRAM gauss_ppal