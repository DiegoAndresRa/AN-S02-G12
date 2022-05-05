PROGRAM metodo_de_Bairstow
    IMPLICIT NONE
    REAL (KIND=8)::a(0:4),b(0:4),c(0:4)
    REAL (KIND=8) dr,ds,epsilon_s,epsilon_ar,epsilon_as,s,s0,r,r0
    REAL (KIND=8)::Re(0:4),Im(0:4)  !Parte real e imaginaria de las respectivas raíces localizadas
    REAL (KIND=8) det,r1,i1,r2,i2
    INTEGER itermax,j,jj,i,k,contador,n,nn
    OPEN(UNIT=20,FILE='incisoC.txt',STATUS='replace',err=500)
    !***********Valores iniciales que se usan para obtener el
    !***********polinomio f(x)=x^2-rx-s, cuyas raíces nos permitiran localizar raíces
    r0=-0.5D0
    s0=-0.5D0
    r=r0
    s=s0
    !***********Declaración del polinomio y grado
    a(4)=1.d0;a(3)=-2.d0;a(2)=6.d0;a(1)=-2.d0;a(0)=5
    n=4      !Grado del polinomio del cual se determinarán las raíces
    nn=n
    !**********************************************************************
    itermax=100     !úmero de iteraciones máximas y tolerancia
    epsilon_s=0.0000000005D0  !0.000000000005D0  ![%] Tolerancia del error
    epsilon_as=0.D0 ! Inicialización del error aproximado de "s"
    epsilon_ar=0.D0 ! Inicialización del error aproximado de "r"
    contador=0      ! Inicialización del contador para saber en que iteración converge a la solución
    !*****Nota: El DO exterior se emplea para repetir el ciclo de iteraciones 
    !*****     una y otra vez hasta que se alcance la tolerancia deseada
    !**************************************************************************************
    write(*,300) 'k','Re','Im'
    write(20,300) 'k','Re','Im'
    300 format(A1,17x,A2,23x,A2)
    DO WHILE (n.GT.2 .AND. contador.LT.itermax)  !***Ciclo principal [1]***Calculo de las "n-2" raíces descompuestas en Raiz(n)=Re(n)+Im(n)*i
            !Se calculan los coeficientes b(n) y c(n) necesarios para el método de Bairstow  
        contador=0
            DO j=1,itermax,1 !****Ciclo interno [2]***Inicio del cálculo de "r" y "s"  
             contador=contador+1
                b(n)=a(n)
                b(n-1)=a(n-1)+r*b(n)
                c(n)=b(n)
                c(n-1)=b(n-1)+r*c(n)
                DO jj=n-2,0,-1
                 b(jj)=a(jj)+r*b(jj+1)+s*b(jj+2)
                 c(jj)=b(jj)+r*c(jj+1)+s*c(jj+2)
                ENDDO   !Fin del calculo de los coeficientes b(n) y c(n) 

            !****Calculo del determinante del sistema (55) de los apuntes 
            !**[7.34 y 7.35 de de Métodos numéricos para ingenieros, 7ª edición]
                det=c(2)*c(2)-c(1)*c(3)
            !***Obtención de r y s
                IF (det.NE.0) THEN
                    dr=(-b(1)*c(2)+b(0)*c(3))/det
                    ds=(-b(0)*c(2)+b(1)*c(1))/det
                    r=r+dr
                    s=s+ds
                    IF (r.NE.0) THEN
                        epsilon_ar=ABS(dr/r)*100.D0
                    ENDIF       
                    IF (s.NE.0) THEN
                        epsilon_as=ABS(ds/s)*100.D0
                    ENDIF  
                ELSE
                !si la aproximación de r y s nos produce un sistema inconsistente, se aumenta en 
                ! una unidad para enmendar el problema y así, el sistema sea consistente
                        r=r+1.D0
                        s=s+1.d0
                        contador=0 !se reinicia el contador
                ENDIF
        !*****Condición de convergencia o de máximo número de iteraciones
                IF ((epsilon_s>ABS(epsilon_as) .AND. epsilon_s>ABS(epsilon_ar)).OR.CONTADOR.GT.itermax) EXIT
            END DO           !****CICLO INTERNO [2]***Fin del cálculo de "r" y "s"  
      
    !*Cálculo de la Raíz(n) 
        CALL CUADRATICA(r,s,r1,i1,r2,i2)
        Re(n)=r1;Im(n)=i1;Re(n-1)=r2;Im(n-1)=i2
    !El cociente que resulta (en términos de b(n)) se renombra ahora como a(n), y este reduce su grado en 2
    !para volver insertarse al método y así determinar todas sus raíces
        n=n-2
        DO i=0,n
            a(i)=b(i+2)
        ENDDO 
        contador=contador
    ENDDO!***Ciclo principal [1]***Calculo de las "n-2" raíces descompuestas en Raiz(n)=Re(n)+Im(n)*i
             !**************************

    IF (contador<itermax) THEN
        IF (n.EQ.2) THEN
            r=-a(1)/a(2)
            s=-a(0)/a(2)
            CALL CUADRATICA(r,s,r1,i1,r2,i2)
            Re(n)=r1
            Im(n)=i1
            Re(n-1)=r2
            Im(n-1)=i2
        ELSE 
            Re(n)=-a(0)/a(1)
            Im(n)=0.D0
        END IF
    ELSE
        contador=1
    ENDIF
    !Impresión de las raíces en archivo de datos
    30 FORMAT(I3,2f25.16)
    DO k=nn,1,-1
        WRITE(20,30) k,Re(k),Im(k) 
       WRITE(*,30) k,Re(k),Im(k) 
    ENDDO
    CLOSE(UNIT=20,STATUS='keep',ERR=500)        
        
500 END PROGRAM metodo_de_Bairstow

!***********Subrutina para calcular las raíces de x^2-rx-s=0
SUBROUTINE CUADRATICA(r,s,r1,i1,r2,i2)
    IMPLICIT NONE
    REAL (KIND=8) r,s,r1,i1,r2,i2,disc
    disc=r**2+4.D0*s
    IF (disc.GT.0) THEN
            r1=(r+SQRT(disc))/2.D0
            r2=(r-SQRT(disc))/2.D0
            i1=0.D0
            i2=0.D0
        ELSE 
            r1=r/2
            r2=r1
            i1=SQRT(ABS(disc))/2.D0
            i2=-i1
    ENDIF  
ENDSUBROUTINE 