subroutine residuo(m,n,a,b,u,r)
    !Programa que calcula el residuo r = Au - b, con el objetivo de ver si se acerca a 0
    use mod_clreal
    implicit none
    integer,intent(in):: m,n
    real(kind=clreal),intent(in)::a(m,n),u(n),b(m)
    real(kind=clreal),intent(out)::r(m)
    integer::j

    r = -b
    do j=1,n
        r=r+a(:,j)*u(j)
    end do
end