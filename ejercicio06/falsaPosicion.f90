program falsa_posicion
    implicit none
    real (kind=8) xl,xu,epsilon_s,f,fl,fu
    real (kind=8) xanterior,xr,epsilon_a
    integer n, max_iter,j

    open(unit=99,file='falsa_posicion.txt',status='replace',err=23)
    ! Intervalo que encierra la raiz    
    xl=8.d0
    xu=9.d0
    !decimales de presición 
    n=4
    !Tolerancia (ec de 1a unidad 2)
    epsilon_s=0.5d0*10.d0**(2-n) ![%]

    !Número máximo de iteraciones
    max_iter=100 
    !inicializamos el valor anterior
    xanterior=0.d0
    77 format(I3,3x,F16.10,3x,F16.10,3x,F16.10,3x,F16.10)
    20 format(2x,A1,12x,A2,17x,A2,17x,A2,13x,A9)
    write(99,20) 'j','hl','hu','hr','epsilon_a'
    write(*,20) 'j','hl','hu','hr','epsilon_a'
    do j=1,max_iter,1
        ! calculo de fl y fu para poder evaluar el cambio de signo
        fl=f(xl)
        fu=f(xu)
        ! aproximación de la raíz xr
        xr=xu-(f(xu)*(xl-xu)/(f(xl)-f(xu)))
        epsilon_a=abs(100.d0*(xr-xanterior)/xr)

        write(*,77)j,xl,xu,xr,epsilon_a
        write(99,77)j,xl,xu,xr,epsilon_a
        !se realizaran las pruebas (a) y (b) del "Algoritmo parael método de disección"
        if (fl*f(xr)<0) then
            xu=xr 
            else
                xl=xr
        end if 
        ! se almacena el último valor calculado de xr como xanterior, 
        !para obtener progresivamente el error relativo aproximado
        xanterior=xr
        !secuencia de escape en caso de que epsilo a<epsilon
        if (epsilon_a < epsilon_s)exit

    end do

    close(unit=99,status='keep',err=23)    
23 end program falsa_posicion
!***********************************!
function f(h)
    implicit none
    real (kind=8) f,h,pi,v,R
    v = 30
    R = 3
    PI=2.d0*DASIN(1.d0)
    f=v-(PI*h**2)*(((3*R)-h)/3)
end function