SUBROUTINE gauss(n,a,b,u,deter)

    USE mod_clreal
    IMPLICIT NONE
    INTEGER::i,j,k
    REAL(clreal)::z
    INTEGER, INTENT(IN)::n
    REAL(clreal),INTENT(OUT)::deter
    REAL(clreal),INTENT(INOUT)::a(n,n)
    REAL(clreal),INTENT(INOUT)::b(n)
    REAL(clreal),INTENT(OUT)::u(n)
    
    !Eliminacion
    do k=1,n-1
    !comprobacion de que o
    !k-esimo elemento diagonal non e nulo
        if(abs(a(k,k))<1.e-12) then
            print*,"Elemento diagonal nulo"
            print*,"na etapa: ",k
            stop
        end if
        !eliminacion
        do i=k+1,n
            z=a(i,k)/a(k,k)
            a(i,k)=0. !Posta a cero da parte inferior da matriz - non obrigatoriodo j=k+1,n
            do j=k+1,n
            a(i,j)=a(i,j)-z*a(k,j)
            end do
        b(i)=b(i)-z*b(k)
        end do
    end do
    !comprobacion de que o
    !ultimo elemento diagonal non e nulo
    if(abs(a(n,n))<1.e-12) then
        print*,"Elemento diagonal nulo"
        print*,"na etapa: ",n
        stop
    end if
    !C´alculo del determinante
    deter=1.
    do i=1,n
        deter=deter*a(i,i)
    end do

END SUBROUTINE gauss