function y = Lk(x,k)
    [Xn,Xm] = size(x)
    r = [x(1:k-1) x(k+1:Xm)]
    p = poly(r,"x","roots")
    pk = horner(p,x(k))
    y = p / pk
endfunction

function z = interpolLagrange(x,y)
    [Xn,Xm] = size(x)
    pol = 0
    for k = 1:Xm
        pol = pol + (Lk(x,k)*y(k))
    end
    z = pol
endfunction

function z = diferenciaDividida(x,y)
    [Xn,Xm] = size(x)
    if Xm == 1 then
        z = y(1)
    else
        z =((diferenciaDividida(x(2:Xm),y(2:Xm))-diferenciaDividida(x(1:Xm-1),y(1:Xm-1))) / (x(Xm)-x(1)))
    end
endfunction

function z = interpolNewton(x,y)
    [Xm,Xn] = size(x)
    z= 0
    for k = 1:Xn
        pol = 1
        for j = 1:k-1
            pol = poly(x(1:j),"x","roots")
        end
        disp(pol)
        dif = diferenciaDividida(x(1:k),y(1:k))
        
        
       disp(dif)
       dif = dif*pol
        z = z + dif
    end
endfunction

function z = minimosCuadrados(x,y,n)
    //Ejecutar con sistemasDeEcuaciones.sce
    [Xm,Xn] = size(x)
    y = y'
    A = zeros(Xn,n)
    for i = 1:Xn
        for j = 1: n
            A (i,j) = x(i) ^ (j-1)
        end
    end
    nuevo = A' * A
    disp(nuevo)
    a = inversa(nuevo)*A'
    disp(a)
        a = a*y
    pol = poly(a, "x", "coeff")
    z = pol
endfunction

// Recibe un numero n
// Devuelve el polinomio de chebyshev de ese grado con sus raíces
function [y,x] = Chebyshev(n)
    t(1) = 1
    t(2) = poly([0],"x","r")
    for i = 3:n+1
        t(i) = poly([0 2], "x", "coeff")*t(i-1)-t(i-2) 
    end
    y = t(n+1)
    x = roots(y)
endfunction

function y = NodosChebyshev(n,a,b)
    [pol,r] = Chebyshev(n)
    for i = 1 : n
        y(i) = ((b+a) + r(i) * (b - a))/2
    end
endfunction