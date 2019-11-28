function y = trapecio(f,x0,x1)
    // pol1 = poly([x1],"x",r)
    //pol1 = pol1 / (x0-x1)
    //pol1 = pol1 * f(x0)
    //pol0 = poly([x0],"x",r)
    //pol0 = pol0 / (x1-x0)
    //pol0 = pol0 * f(x1)
    //y = pol1 + pol0
    y = (x1-x0)*(f(x0)+f(x1))/2
endfunction

function y = trapecioComp(f,a,b,n)
    h = (b-a)/n
    contador = f(a)*1/2 + f(b)*1/2
    for j = 1:n-1
        contador = contador + f(a+j*h) 
    end
    y = h*contador
endfunction

function y = trapecioCompDoble(f,xi,a,b,n)
    h = (b-a)/n
    contador = f(xi,a)*1/2 + f(xi,b)*1/2
    for j = 1:n-1
        contador = contador + f(xi,a+j*h) 
    end
    y = h*contador
endfunction

function y = simpson(f,a,b)
    h = (b-a)/2
    c = a+h
    contador = h/3
    funciones = f(a)+4*f(c)+f(b)
    y = contador*funciones
endfunction

function y = simpsonComp(f,a,b,n)
    if modulo(n,2) <> 0 then
        disp("N no es par")
        y = %nan
        return
    end
    h = (b-a)/n
    y = f(a)+f(b)
    for i = 1:n-1
        if modulo(i,2) == 0 then
            y = y+2*f(a+i*h)
        else
            y = y+4*f(a+i*h)
        end
    end
    y = y*h/3
endfunction

function y = simpsonCompDoble(f,xi,a,b,n)
    if modulo(n,2) <> 0 then
        disp("N no es par")
        y = %nan
        return
    end
    h = (b-a)/n
    y = f(xi,a)+f(xi,b)
    for i = 1:n-1
        if modulo(i,2) == 0 then
            y = y+2*f(xi,a+i*h)
        else
            y = y+4*f(xi,a+i*h)
        end
    end
    y = y*h/3
endfunction

function y = polinomio(f,x0,x1)
    pol1 = poly([x1],"x",r)
    pol1 = pol1 / (x0-x1)
    pol1 = pol1 * f(x0)
    pol0 = poly([x0],"x",r)
    pol0 = pol0 / (x1-x0)
    pol0 = pol0 * f(x1)
    y = pol1 + pol0
endfunction

function y = TrapecioBi(f,a,b,c,d,n,m)
    // a b contantes externas
    // c,d funciones de segunda integral
    // n intervalor entre a y b
    // m intervalos en c(x) d(x)
    hx = (b-a)/n
    for i = 0:n
        xi = a+i*hx
        G(i+1) = trapecioCompDoble(f,xi,c(xi),d(xi),m)
    end
    contador = G(1)*1/2 + G(n+1)*1/2
    for j = 1:n-1
        contador = contador + G(j+1) 
    end
    y = hx*contador
endfunction


function y = SimpsonCompBi(f,a,b,c,d,n,m)
    // a b contantes externas
    // c,d funciones de segunda integral
    // n intervalor entre a y b
    // m intervalos en c(x) d(x)
    if modulo(n,2) <> 0 then
        disp("N no es par")
        y = %nan
        return
    end
    hx = (b-a)/n
    for i = 0:n
        xi = a+i*hx
        G(i+1) = simpsonCompDoble(f,xi,c(xi),d(xi),m)
    end
    y = G(1)+G(n+1)
    for i = 1:n-1
        if modulo(i,2) == 0 then
            y = y+2*G(i+1)
        else
            y = y+4*G(i+1)
        end
    end
    y = y*hx/3
endfunction

function x = rangosX(y)
    c = 1-y^2
    c = c ^ 1/2
    x(1) = c + 1
    x(2) = -c + 1
endfunction

function z = d(y)
    z = rangosX(y)(1)
endfunction

function z = c(y)
    z = rangosX(y)(2)
endfunction

