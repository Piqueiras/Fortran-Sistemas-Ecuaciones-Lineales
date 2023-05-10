PROGRAM main
    
    USE mod_clreal
    IMPLICIT NONE

    INTEGER::n,info,i,j
    REAL(clreal),ALLOCATABLE::A(:,:),b(:),u(:),r(:)
    REAL(clreal)::deter
    INTEGER,ALLOCATABLE::ip(:)
    real :: start, finish
	
	call cpu_time (start)
    

    print*, "======GAUSS PIVOTE======"

    read*, info

    read*,n
    if ( info==1 ) then
        print*,"Tamaño: ",n
    end if
    
	ALLOCATE (A(n,n), b(n), u(n), r(n), ip(n))
	
    CALL lecmat(A,n)
    if ( info==1 ) then
        print*, "Matriz A:"
        CALL prinmat(A,n)
    end if

    read*,b
    if ( info==1 ) then
        print*, "Vector de términos independientes:", b
    end if

    CALL gauss_pivote(n,A,b,deter,ip)

    print*, 'Matriz pivotada y reducida:'
    do i=1,n
        print*, (0._clreal,j=1,i-1),a(ip(i),i:n) 
    end do

    print*, 'Vector independiente permutado y modificado:', b(ip(1:n))

    print*, "Permutacion: ",ip

    CALL remonte_permutado(n,A,b,u,ip)

    print*, "Solución: ",u
    print*, "Determinante: ",deter
    
    DEALLOCATE(A,b,u,r,ip)
    
    call cpu_time (finish)
	
	print*, 'Tiempo total de calculo: ',finish - start
END PROGRAM main