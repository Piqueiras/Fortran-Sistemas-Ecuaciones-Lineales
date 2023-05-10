PROGRAM main
    
    USE mod_clreal
    IMPLICIT NONE

    INTEGER::n,info,i,j
    INTEGER,ALLOCATABLE::ip(:)
    REAL(clreal),ALLOCATABLE::A(:,:),b(:),u(:),r(:)
    REAL(clreal)::deter

    print*, "======GAUSS LU PIVOTE======"

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

    CALL gauss_PALU(n,A,deter,ip)
    

    print*, 'Matriz pivotada y factorizada:'
    do i=1,n
        print*, a(ip(i),1:n) 
    end do

    print*, "Permutacion:",ip

    CALL descenso_L_permutado(n,A,b,u,ip)
    CALL remonte_permutado_LU(n,A,u,b,ip)
   

    print*, "Solución: ",b
    print*, "Determinante: ",deter
END PROGRAM main