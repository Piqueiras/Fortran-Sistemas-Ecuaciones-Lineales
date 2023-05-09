subroutine gauss(n,A,b,deter)
    !Subrutina que emplea reducción de Gauss para convertir la matriz en triangular superior
    !En el resto de elementos se guardan los factores de multiplicado, así que la matriz resultante se va a tener que calcular siempre como triangular superior
    !Además se calcula el determinante gracias a los elementos diagonales
    !Nótese que, por tanto, A queda descompuesta en su factorización LU, de la forma U+L-I (con I la identidad puesto que la diagonal de L se pierde, pero son todos 1s)
    use mod_clreal
    implicit none
    integer,intent(in)::n
    real(clreal),intent(out)::deter
    real(clreal),intent(inout)::A(n,n),b(n)

    integer :: i,j,k
    real(clreal) :: z, eps

    eps = 1.e-12
    deter = 1.0

    do k=1,n-1
        if ( abs(a(k,k)) < eps ) then
            stop "Elemento diagonal nulo"
        end if
        do i=k+1,n !Iteramos por filas
            z=a(i,k)/a(k,k) !z va a ser el factor de multiplicación

            !do j=k+1,n !Iteramos por columnas
            !    a(i,j)=a(i,j)-z*a(k,j) !Se trata de la trasformación Fila i = Fila i - z * Fila k
            !enddo
            
            A(i,k+1:n)=A(i,k+1:n)-z*A(k,k+1:n) !Así se hace más facil la transformación
            b(i)=b(i)-z*b(k) !También tenemos que hacer lo mismo con el vector de términos independientes
            a(i,k)=z !Guardamos el factor de multiplicación. Por si acaso. Total, la otra opción era ponerlo a 0.
        enddo
        deter = deter * a(k,k)
    enddo
    if(abs(a(n,n)) < eps ) then
        stop "Elemento diagonal nulo"
    end if
    deter = deter * a(n,n)
end subroutine gauss

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

subroutine gauss_pivote(n,A,b,deter,ip)
    !Para evitar errores de redondeo vamos a usar el metodo del pivote parcial
    !Consiste basicamente en mover arriba la fila cuyo pivote tenga mayor valor absoluto
    use mod_clreal
    implicit none
    integer::i,j,k,cont,ipi,ipk,ipiv
    real(clreal)::piv,factor,eps
    integer, intent(in)::n
    real(clreal), intent(inout)::a(n,n),b(n)
    real(clreal), intent(out)::deter
    integer,dimension(n), intent(out)::ip

    deter=1
    cont=0
    eps=1.e-12

    ip=(/(i,i=1,n)/) !Constructor de arrays (1,2,...,n)

    do k=1,n-1
        piv=a(ip(k),k)
        ipiv=k

        do i=k+1,n !Buscamos el mayor pivote
            if(abs(piv)<abs(a(ip(i),k))) then
                piv=a(ip(i),k)
                ipiv=i
            end if
        end do

        if(abs(piv)<eps) then
            stop "Pivote nulo (aqui si aseguramos matriz singular)"
        end if

        if(ipiv/=k) then
            ipk=ip(ipiv)
            ip(ipiv)=ip(k)
            ip(k)=ipk
            cont=cont+1
        else
            ipk=ip(k)
        end if
        deter=deter*piv

        !Eliminacion
        do i=k+1,n
            ipi=ip(i)
            factor=a(ipi,k)/piv
            do j=k+1,n
                a(ipi,j)=a(ipi,j)-factor*a(ipk,j)
            end do
            b(ipi)=b(ipi)-factor*b(ipk)
        end do
    end do

    if(abs(a(n,n))<eps) then
    stop "Elemento final nulo"
    end if

    deter=deter*a(ip(n),n)*(-1)**cont
end subroutine gauss_pivote

subroutine gauss_PALU(n,A,P)
    use mod_clreal
    implicit none
    integer,intent(in)::n
    real(clreal),intent(inout)::A(n,n)
    real(clreal),intent(out)::P(n,n)
end subroutine gauss_PALU