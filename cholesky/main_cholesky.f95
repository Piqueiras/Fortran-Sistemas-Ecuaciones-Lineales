PROGRAM main
    
    USE mod_clreal
    IMPLICIT NONE

    INTEGER::n
    INTEGER,ALLOCATABLE::ip(:)
    REAL(clreal),ALLOCATABLE::A(:,:),b(:),u(:),v(:)

    read*,n
    print*,"Tamaño: ",n

    ALLOCATE(A(n,n),b(n),u(n),v(n))

    CALL lecmat(A,n)
    print*, "Matriz A:"
    CALL prinmat(A,n)

    read*,b
    print*, "Vector de términos independientes:", b

    CALL cholesky(n,A)

    print*, "Nueva matriz reducida:"
    CALL prinmat(A,n)

    CALL descenso(n,A,b,u)
    CALL remonte(n,TRANSPOSE(A),u,v)
    
    print*, "Solución: ",v

    DEALLOCATE(a,b,u,v)
END PROGRAM main