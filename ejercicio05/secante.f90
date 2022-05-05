program secante
    implicit none
    real (kind=8) xi,f,Df,ximas1,epsilon_s,epsilon_a
    integer itermax,n,i
    !archivo donde se almacenarán as iteraciones
    open(unit=99,file='secante.txt',status='replace',err=23) 
    
    !valor inicial
    xi=-2.d0
    !iteraciones máximas
    itermax=100
    !núnero de cifras significativas en el resultado
    n=3
    epsilon_s=0.5d0*10.d0**(2-n)

    32 format(i3,4x,F11.8,4x,f11.8,4x,f11.8)
    20 format(2x,A1,9x,A3,11x,A5,7x,A9)
    write(99,20) 'i','x_i','x_i+1','epsilon_a'
    write(*,20) 'i','x_i','x_i+1','epsilon_a'

    do i=1, itermax
        ximas1=xi-f(xi)/Df(xi)
        epsilon_a=abs(100*(ximas1-xi)/ximas1)

        write(*,32) i,xi,ximas1,epsilon_a
        write(99,32) i,xi,ximas1,epsilon_a

        !condición de escape en caso de convergencia anticipada
        if(epsilon_a<epsilon_s)exit
        xi=ximas1
    end do 
    !cierre de archivo
    close(unit=99,status='keep',err=23)
23 end program
!******f(x)+*******
function f(x)
    implicit none
    real(kind=8) f,x
    f=-12.d0-(21.d0*x)+(18.d0*x**2)-(2.4d0*x**3)
end function
!******Df(x)+*******
function Df(x)
    implicit none
    real(kind=8) Df,x,h,f
    h=0.0001d0
    Df=(f(x+h)-f(x-h))/(2.d0*h)
end function