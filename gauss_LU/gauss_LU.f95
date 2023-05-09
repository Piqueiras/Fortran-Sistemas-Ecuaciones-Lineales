subroutine gauss_LU(n,A)
    !Subrutina anterior modificada para simplemente calcular la LU, que encaja en la matriz original
    use mod_clreal
    implicit none
    integer,intent(in)::n
    real(clreal),intent(inout)::A(n,n)
    integer :: i,j,k
    real(clreal) :: z, eps
    eps = 1.e-12

    do k=1,n-1
        if ( abs(a(k,k)) < eps ) then
            stop "Elemento diagonal nulo"
        end if
        do i=k+1,n !Iteramos por filas
            z=a(i,k)/a(k,k) !z va a ser el factor de multiplicación       
            A(i,k+1:n)=A(i,k+1:n)-z*A(k,k+1:n) !Así se hace más facil la transformación
            a(i,k)=z !Guardamos el factor de multiplicación. Por si acaso. Total, la otra opción era ponerlo a 0.
        enddo
    enddo
    if(abs(a(n,n)) < eps ) then
        stop "Elemento diagonal nulo"
    end if
end subroutine gauss_LU