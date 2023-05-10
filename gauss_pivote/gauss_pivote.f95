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