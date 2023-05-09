PROGRAM main
    
    USE mod_clreal
    IMPLICIT NONE

    INTEGER::n
    INTEGER,ALLOCATABLE::ip(:)
    REAL(clreal),ALLOCATABLE::A(:,:),b(:),u(:),r(:)
    REAL(clreal)::deter

    print*, "======GAUSS LU PIVOTE======"

    read*,n
    print*,"Tamaño: ",n

    ALLOCATE(A(n,n),b(n),u(n),r(n),ip(n))

    CALL lecmat(A,n)
    print*, "Matriz A:"
    CALL prinmat(A,n)

    read*,b
    print*, "Vector de términos independientes:", b

    CALL gauss_PALU(n,A,deter,ip)
    print*, "Determinante: ",deter

    print*, "Nueva matriz reducida:"
    CALL prinmat(A,n)

    print*, "Permutacion:",ip

    CALL descenso_L_permutado(n,A,b,u,ip)
    CALL remonte_permutado(n,A,u,b,ip)
   

    print*, "Solución: ",b, u
END PROGRAM main