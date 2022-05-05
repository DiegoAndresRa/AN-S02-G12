program Newton_Rapson_V_Var
    implicit none
    real(kind=8) xi,yi,h,ximas1,yimas1,epsilon_s
    real(kind=8) dudx,dudy,dvdx,dvdy,det_jacobiano,u,v,epsilon_ax,epsilon_ay
    integer n,iter_max,i
    open(unit=13,file='Newton_Raphson.txt',status='replace',err=23)
    !vaores iniciales i=0
    xi=-0.65d0
    yi=1.4d0
    !incremento 'h' para las derivadas
    h=0.0001d0
    !cifras significativas
    n=8
    !tolerancia del error
    epsilon_s=0.5d0*10.d0**(2-n)
    !iteraciones máximas
    iter_max=100
    20 format(I3,2x,6(F21.16,2x))
    24 format(A3,13x,A2,21x,A2,20x,A6,17x,A6,14x,A10,14x,A10)
    write(*,24) 'i','xi','yi','ximas1','yimas1','epsilon_ax','epsilon_ay'
    write(13,24) 'i','xi','yi','ximas1','yimas1','epsilon_ax','epsilon_ay'
    do i=1,iter_max,1
        !Cálcilo de la raóz i
        ximas1=xi-(u(xi,yi)*dvdy(xi,yi,h)-v(xi,yi)*dudy(xi,yi,h))/det_jacobiano(xi,yi,h)
        yimas1=yi-(v(xi,yi)*dudx(xi,yi,h)-u(xi,yi)*dvdx(xi,yi,h))/det_jacobiano(xi,yi,h)

        !Cálculo de los errores relativos aproximados
        epsilon_ax=abs(100.d0*(ximas1-xi)/ximas1)
        epsilon_ay=abs(100.d0*(yimas1-yi)/yimas1)

        write(*,20) i,xi,yi,ximas1,yimas1,epsilon_ax,epsilon_ay
        write(13,20) i,xi,yi,ximas1,yimas1,epsilon_ax,epsilon_ay
        ! Condición de escape
        if(epsilon_ax < epsilon_s .and. epsilon_ay < epsilon_s) exit
        ! En caso de no converger, se reitilizan los valores recien calculados
        xi=ximas1
        yi=yimas1
    end do 
    close(unit=13,status='keep',err=23)
23 end program

!*********  funciones  *********!
! u y v
function u(x,y)
    implicit none
    real(kind=8) u,x,y
    u=y-x**2-1.d0
end function

function v(x,y)    
    implicit none
    real(kind=8) v,x,y
    v=y-2.d0*cos(x)
end function


function dudx(x,y,h)
    implicit none
    real(kind=8) dudx,u,x,y,h
    dudx=(-u(x+2.d0*h,y)+8.d0*u(x+h,y)-8.d0*u(x-h,y)+u(x-2.d0*h,y))/(12.d0*h)
end function

function dudy(x,y,h)
    implicit none
    real(kind=8) dudy,u,x,y,h
    dudy=(-u(x,y+2.d0*h)+8.d0*u(x,y+h)-8.d0*u(x,y-h)+u(x,y-2.d0*h))/(12.d0*h)
end function

function dvdx(x,y,h)
    implicit none
    real(kind=8) dvdx,v,x,y,h
    dvdx=(-v(x+2.d0*h,y)+8.d0*v(x+h,y)-8.d0*v(x-h,y)+v(x-2.d0*h,y))/(12.d0*h)
end function

function dvdy(x,y,h)
    implicit none
    real(kind=8) dvdy,v,x,y,h
    dvdy=(-v(x,y+2.d0*h)+8.d0*v(x,y+h)-8.d0*v(x,y-h)+v(x,y-2.d0*h))/(12.d0*h)
end function

!Jacobiano
function det_jacobiano(x,y,h)
    implicit none
    real(kind=8) dudx,dudy,dvdx,dvdy,x,y,h,det_jacobiano
    det_jacobiano=dudx(x,y,h)*dvdy(x,y,h)-dudy(x,y,h)*dvdx(x,y,h)
end function