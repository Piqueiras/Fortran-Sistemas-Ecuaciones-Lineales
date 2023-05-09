subroutine doolittle_LU(n,a,deter)
    use mod_clreal
    implicit none
    integer,intent(in)::n
    real(clreal),intent(inout)::a(n,n)
    real(clreal), intent(out):: deter
    integer:: i,j,k
    real(clreal):: piv, l, eps

    eps = 1.e-12
    deter = a(1,1)

    if(abs(a(1,1))<eps) then
        stop "Pivote nulo"
    end if

    do j=2,n
        a(j,1)=a(j,1)/a(1,1)
    end do

    do i=2,n
        !fila i de U
        do j=i,n
            do k=1,i-1
                a(i,j)=a(i,j)-a(i,k)*a(k,j)
            end do
            !a(i,j)=a(i,j)-sum(a(i,1:i-1)*a(1:i-1,j))
        end do

        !comprobación de que el k-esimo pivote no es nulo
        if(abs(a(i,i))<eps) then
            print*, 'Pivote nulo en la etapa: ',i
            stop 
        end if

        !columna i de L
        do j=i+1,n
            do k=1,i-1
                a(j,i)=a(j,i)-a(j,k)*a(k,i)
            end do
            !a(j,i)=a(j,i)-sum(a(j,1:i-1)*a(1:i-1,i))

            a(j,i)=a(j,i)/a(i,i)
        end do

        deter=deter*a(i,i)

    end do

end subroutine doolittle_LU