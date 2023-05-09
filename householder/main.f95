PROGRAM main
    
    USE mod_clreal
    IMPLICIT NONE

    INTEGER::n,info
    REAL(clreal),ALLOCATABLE::A(:,:),b(:),u(:),r(:)
    REAL(clreal)::deter

    print*, "======HOUSEHOLDER======"

    read*, info

    read*,n
    if ( info==1 ) then
        print*,"Tamaño: ",n
    end if
    
	ALLOCATE (A(n,n), b(n), u(n), r(n))
	
    CALL lecmat(A,n)
    if ( info==1 ) then
        print*, "Matriz A:"
        CALL prinmat(A,n)
    end if
    
    read*,b
    if ( info==1 ) then
        print*, "Vector de términos independientes:", b
    end if

    CALL householder(n,A,b,deter)
    

    print*, "Nueva matriz reducida:"
    CALL prinmat(A,n)
    print*, "Nuevo vector independiente:", b

    CALL remonte(n,A,b,u)

    print*, "Solución: ",u
    print*, "Determinante: ",deter
END PROGRAM main