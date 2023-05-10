subroutine jacobi(n,A,b,u,eps)
    use mod_clreal
    implicit none
    integer,intent(in)::n
    real(clreal),intent(out)::u(n)
    real(clreal),intent(in)::A(n,n),b(n),eps
    integer :: i,j,k
    real(clreal)::uold(n),error
    error = 69420.0
    
    k=0
    do while ( error>eps .and. k<50)
        uold=u
        do i = 1, n
            u(i)=(b(i)-sum(a(i,1:i-1)*uold(1:i-1))-sum(a(i,i+1:n)*uold(i+1:n)))/a(i,i)
        end do
    
        error = SQRT(DOT_PRODUCT(ABS(u - uold), ABS(u - uold)))
        k=k+1
        print*, "Iteracion ",k,": ",uold

        
    end do

end subroutine jacobi

logical function dominante(A, n)
  implicit none
  integer, intent(in) :: n
  real, dimension(n,n), intent(in) :: A
  integer :: i, j
  real :: sum

  dominante = .true.  ! Asumimos que si es

  do j = 1, n
    sum = 0.0
    do i = 1, n
      sum = sum + abs(A(i,j))
    end do
    if (sum - abs(A(j,j)) >= abs(A(j,j))) then
      dominante = .false.  ! A no es dominante
      return
    end if
  end do

end function dominante
