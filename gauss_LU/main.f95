PROGRAM gauss_ppal

    USE mod_clreal
    IMPLICIT NONE
    INTEGER::n,info
    REAL(clreal), ALLOCATABLE :: A(:,:),B(:),U(:)

    print*, "======GAUSS LU======"

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

    CALL gauss_LU(n,A)

    print*, "Matriz modificada LU:"
    CALL prinmat(A,n)

    CALL descenso_L(n,A,b,u)
    CALL remonte(n,A,u,b)
    PRINT*, "Soulución:", b

    DEALLOCATE(A,b,u)

END PROGRAM gauss_ppal