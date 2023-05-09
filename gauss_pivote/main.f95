PROGRAM main
    
    USE mod_clreal
    IMPLICIT NONE

    INTEGER::n
    REAL(clreal),ALLOCATABLE::A(:,:),b(:),u(:),r(:),ip(:)
    REAL(clreal)::deter

    print*, "======GAUSS PIVOTE======"

    read*,n
    print*,"Tamaño: ",n

    ALLOCATE(A(n,n),b(n),u(n),r(n),ip(n))

    CALL lecmat(A,n)
    print*, "Matriz A:"
    CALL prinmat(A,n)

    read*,b
    print*, "Vector de términos independientes:", b

    CALL gauss_pivote(n,A,b,deter,ip)
    print*, "Determinante: ",deter

    print*, "Nueva matriz reducida:"
    CALL prinmat(A,n)
    print*, "Nuevo vector independiente:", b

    CALL remonte_permutado(n,A,b,u,ip)

    print*, "Solución: ",u
END PROGRAM main