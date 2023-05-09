subroutine cholesky(n,A)
    !Subrutina que calcula la factorización A=BBt para matrices simétricas
    use mod_clreal
    implicit none
    integer,intent(in)::n
    real(clreal),intent(inout)::A(n,n)
    integer :: i,j,k

    do j=1,n
        !Esquina
        if ( a(j,j) <= 0 ) then
            stop "Matriz no definida positiva (o no es simétrica)"
        end if
        a(j,j)=a(j,j)-DOT_PRODUCT(a(j,1:j-1),a(j,1:j-1))
        a(j,j)=SQRT(a(j,j))
        !Resto columna
        do i=j+1,n
            a(i,j)=a(i,j)-DOT_PRODUCT(a(i,1:j-1),a(j,1:j-1))
            a(i,j)=a(i,j)/a(j,j)
            a(j,i)=0
        enddo
    enddo

end subroutine cholesky

subroutine cholesky_tridiagonal(n,d,s)
    !Cholesky modificado para una tridiagonal, de la cual solo necesitamos una diagonal secundaria por ser simetrica
    use mod_clreal
    implicit none
    integer :: i
    integer,intent(in)::n
    real(clreal),intent(inout)::d(n),s(n-1)

    d(1)=SQRT(d(1))
    do i = 1, n-1
        s(i)=s(i)/d(i)
        d(i+1)=SQRT(d(i+1)-s(i)*s(i))
    end do
    
end subroutine cholesky_tridiagonal