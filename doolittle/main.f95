program main
	use mod_clreal
	implicit none
	integer :: n	
	real (kind = clreal) :: deter
	real (kind = clreal), allocatable :: a(:,:), aa(:,:), b(:), x(:), y(:), r(:)

    print*, "======DOOLITTLE======"
	
    read*,n
    print*,"Tamaño: ",n
	
	allocate (a(n,n), b(n), x(n), y(n), r(n))
	
    CALL lecmat(A,n)
    print*, "Matriz A:"
    CALL prinmat(A,n)

    read*,b
    print*, "Vector de términos independientes:", b

	call doolittle_LU (n, a, deter)
	call descenso_L(n, a, b, y)
	call remonte (n, a, y, x)

    print*, "Factorizacion: "
    CALL prinmat(A,n)

    print*, "Solucion: ", x

    print*, "Determinante: ", deter
	
end program main