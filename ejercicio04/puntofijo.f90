program Punto_Fijo 
    implicit none
    real(kind=8) g,epsilon_a,anterior,x_k,x_k1
    integer i,max_iter

    open(unit=99,file='punto_fijo.txt',status='replace',err=23)
    max_iter=100

    anterior = 1.d0
    x_k = 0.5d0

    77 format(I3,3x,F16.10,3x,F16.10,3x,F14.10)
    20 format(2x,A1,12x,A3,15x,A5,10x,A9)
    write(99,20) 'i','x_k','x_k+1','epsilon_a'
    write(*,20) 'i','x_k','x_k+1','epsilon_a'

    do i=1,max_iter,1
        x_k1 = g(x_k)
        epsilon_a=abs(100.d0*(x_k1-anterior)/x_k1)

        if (epsilon_a < 0.01)exit

        write(99,77) i,x_k,x_k1,epsilon_a
        write(*,77) i,x_k,x_k1,epsilon_a

        anterior = x_k1
        x_k = x_k1
        
    end do

    close(unit=99,status='keep',err=23)
23 end program Punto_Fijo 
function g(x)
    implicit none
    real(kind=8) g,x
    g=sin(x**0.5)
end function
