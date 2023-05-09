subroutine householder(n,A,b,deter)
    !El método de householder se aprovecha de las propiedades de la perpendicularidad y las reflexiones 
    use mod_clreal
    implicit none

    integer,intent(in)::n
    real(clreal),intent(out)::deter
    real(clreal),intent(inout)::A(n,n),b(n)

    integer :: i,j,k
    real(clreal) :: z, eps, s2, alfa, beta, p ,q

    deter = 1

    do k=1,n-1
        s2=sum(a(k:n,k)*a(k:n,k))
        if ( abs(s2) < eps ) then
            stop "Vector nulo: matriz singular"
        else if ( sum(a(k+1:n,k)*a(k+1:n,k)) < eps) then
            print*, "No se necesita eliminacion en ", k
            deter=deter*a(k,k)
            cycle
        else
            alfa=sign(sqrt(s2),a(k,k))
            beta=1.0/(alfa*(alfa+a(k,k)))
            a(k,k)=a(k,k)+alfa
        end if

        !Modificacion de columnas
        do j = k+1, n
            q=sum(a(k:n,k)*a(k:n,j))
            p=beta*q
            a(k:n,j)=a(k:n,j)-p*a(k:n,k)
        end do

        !Modificacion termino independiente

        q=sum(a(k:n,k)*b(k:n))
        p=beta*q
        b(k:n)=b(k:n)-p*a(k:n,k)

        !Guardamos -alfa

        a(k,k)=-alfa

        !Actualizar determinante

        deter=-deter*a(k,k)
    enddo

    deter=deter*a(n,n)
    
end subroutine householder