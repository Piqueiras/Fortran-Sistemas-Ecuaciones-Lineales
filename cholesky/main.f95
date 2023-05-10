PROGRAM main
    
    USE mod_clreal
    IMPLICIT NONE

    INTEGER::n,info
    INTEGER,ALLOCATABLE::ip(:)
    REAL(clreal),ALLOCATABLE::A(:,:),b(:),u(:),v(:)
    real :: start, finish
	
	call cpu_time (start)
    

    print*, "======CHOLESKY======"

    read*, info

    read*,n
    if ( info==1 ) then
        print*,"Tamaño: ",n
    end if
    
	ALLOCATE (A(n,n), b(n), u(n), v(n))
	
    CALL lecmat(A,n)
    if ( info==1 ) then
        print*, "Matriz A:"
        CALL prinmat(A,n)
    end if
    
    read*,b
    if ( info==1 ) then
        print*, "Vector de términos independientes:", b
    end if

    CALL cholesky(n,A)

    print*, "Nueva matriz reducida:"
    CALL prinmat(A,n)

    CALL descenso(n,A,b,u)
    CALL remonte(n,TRANSPOSE(A),u,v)
    
    print*, "Solución: ",v

    DEALLOCATE(a,b,u,v)
    call cpu_time (finish)
	
	print*, 'Tiempo total de calculo: ',finish - start
END PROGRAM main