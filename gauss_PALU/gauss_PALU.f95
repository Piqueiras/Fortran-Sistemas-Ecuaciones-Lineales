subroutine gauss_PALU(n,A,deter,ip)
    !Es como el metodo anterior de A=LU pero mezclado con el pivote parcial
    use mod_clreal
    implicit none
    integer::i,j,k,cont,ipi,ipk,ipiv
    real(clreal)::piv,factor,eps
    integer, intent(in)::n
    real(clreal), intent(inout)::a(n,n)
    real(clreal), intent(out)::deter
    integer, intent(out)::ip(n)

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
        end do
    end do

    if(abs(a(n,n))<eps) then
    stop "Elemento final nulo"
    end if

    deter=deter*a(ip(n),n)*(-1)**cont
end subroutine gauss_PALU