program main
	use mod_clreal
	implicit none
	integer :: n, info
	real (kind = clreal) :: deter
	real (kind = clreal), allocatable :: a(:,:), aa(:,:), b(:), x(:), y(:), r(:)

    print*, "======DOOLITTLE======"
	
    read*, info

    read*,n
    if ( info==1 ) then
        print*,"Tamaño: ",n
    end if
	
	ALLOCATE(a(n,n), b(n), x(n), y(n), r(n))
	
    CALL lecmat(A,n)
    if ( info==1 ) then
        print*, "Matriz A:"
        CALL prinmat(A,n)
    end if
    
    read*,b
    if ( info==1 ) then
        print*, "Vector de términos independientes:", b
    end if

	call doolittle_LU (n, a, deter)
	call descenso_L(n, a, b, y)
	call remonte (n, a, y, x)

    print*, "Factorizacion: "
    CALL prinmat(A,n)

    print*, "Solución: ", x

    print*, "Determinante: ", deter
	
end program main