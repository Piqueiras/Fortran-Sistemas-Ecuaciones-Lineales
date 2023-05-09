PROGRAM main
    
    USE mod_clreal
    IMPLICIT NONE

    INTEGER::n
    REAL(clreal),ALLOCATABLE::A(:,:),b(:),u(:),r(:),Acopia(:,:),Bcopia(:)
    REAL(clreal)::deter

    read*,n
    print*,"Tamaño: ",n

    ALLOCATE(A(n,n),b(n),u(n),r(n),Acopia(n,n),Bcopia(n))

    CALL lecmat(A,n)
    print*, "Matriz A:"
    CALL prinmat(A,n)

    Acopia=A

    read*,b
    print*, "Vector de términos independientes:", b

    Bcopia=b

    CALL householder(n,A,b,deter)
    print*, "Determinante: ",deter

    print*, "Nueva matriz reducida:"
    CALL prinmat(A,n)
    print*, "Nuevo vector independiente:", b

    CALL remonte(n,A,b,u)

    print*, "Solución: ",u



END PROGRAM main