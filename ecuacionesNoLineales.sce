//Recibe una funcion continua en un intervalo y dos puntos del intervalo
//Convergencia lineal asegurada

function x= Biseccion(f,a,b)
    if f(a)*f(b) > 0 then
        x = %nan;
    else
        c = (a+b)/2;
        while b-c > 0.001
            if f(b)*f(c) <= 0 then
                a = c;
            else
                b = c;
            end
         c = (a+b)/2;
     end
     x=c
    end
endfunction


//Recibe una funcion y dos puntos de un intervalo
//Convergencia mas rapida que lineal
//Convergencia no asegurada
//Tratar de que F'(a) != 0

function x = Secante(f,a,b)
    eps = 0.0000001
    actual = b - f(b)*(b-a)/(f(b)-f(a));
    anterior = b;
    while abs(f(actual)-f(anterior)) >= eps
        siguiente = actual - f(actual) * (actual-anterior)/(f(actual)-f(anterior))
        anterior = actual
        actual = siguiente
    end
    x = actual
endfunction

function x=PuntoFijo(a,c)
    actual = a
    while abs(actual+sqrt(5)) > 0.001 then
        actual = actual+c*((actual^2)-5)
    end
    x = actual
endfunction

function x = LongitudDeOnda(h,t,aproxl)
    actual = t^2/(2*%pi) * 9.8 * tanh((8*%pi)/aproxl)
        disp(actual)
    while abs(actual-aproxl) > 0.1 then
        aproxl = actual
        actual = t^2/(2*%pi) * 9.8 * tanh((8*%pi)/aproxl)
        disp(actual)
    end
    x = actual
endfunction

//Recibe una funcion y un punto del espacio vectorial
//El punto debe ser una estimacion cercana de la raiz
//Convergencia cuadratica
//Convergencia no asegurada
//Asegurarse que F'(a) != 0 

function y = MetodoNewton(f,a)
    anterior = a
    jacobiana = numderivative(f,a)
    y = a - (inv(jacobiana))*f(a)
    while norm(f(y)-f(anterior)) > 10^(-12)
        jacobiana = numderivative(f,anterior)
        anterior = y
        y = anterior - (inv(jacobiana))*f(anterior)
    end
endfunction

//Recibe una funcion y dos puntos del intervalo
//Convergencia asegurada

function x = RegulaFalsi(f,a,b)
    eps = 0.0000001
    c = b - f(b)*(b-a)/(f(b)-f(a))
    while f(c) > eps
        if(f(a)*f(c) < 0) then
            b = c
        else
            a = c
        c = b - f(b)*(b-a)/(f(b)-f(a))
        end
    end
    x = c
endfunction
