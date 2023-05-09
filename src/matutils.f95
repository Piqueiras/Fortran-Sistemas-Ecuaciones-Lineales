subroutine lecmat(a,n)
    !Subrutina para leer matriz, porque Fortran tiene daño cerebral e indexa en 1
    use mod_clreal
    implicit none
    integer::i
    integer,INTENT(IN)::n
    real(clreal),INTENT(OUT)::a(n,n)
    do i = 1, n
        read*, a(i,:)
    end do
end subroutine lecmat

subroutine prinmat(a,n)
    !Subrutina para escribir matriz, porque Fortran tiene daño cerebral y guarda por columnas
    use mod_clreal
    implicit none
    integer::i
    integer,INTENT(IN)::n
    real(clreal),INTENT(IN)::a(n,n)
    do i = 1, n
        print*, a(i,:)
    end do
end subroutine prinmat

subroutine printrig(n, U, D, L)
    use mod_clreal
    implicit none
    integer, intent(in) :: n
    real(clreal), intent(in) :: U(n-1), D(n), L(n-1)
    integer :: i, j

    do i = 1, n
        do j = 1, n
            if (i == j) then
            write(*,"(F20.16)", advance="no") D(i)
            else if (j == i-1) then
            write(*,"(F20.16)", advance="no") L(i-1)
            else if (j == i+1) then
            write(*,"(F20.16)", advance="no") U(i)
            else
            write(*,"(F20.16)", advance="no") 0.0
            end if
        end do
        write(*,*)
    end do
end subroutine printrig
