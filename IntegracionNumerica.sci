// Metodo del trapcio Simple
// Tiene mucho error de aproximacion
function y = trapecio(f,x0,x1)
    y = (x1-x0)*(f(x0)+f(x1))/2
endfunction

// Metodo del trapecio compuesto para n subdivisiones
function y = trapecioComp(f,a,b,n)
    h = (b-a)/n
    contador = f(a)*1/2 + f(b)*1/2
    for j = 1:n-1
        contador = contador + f(a+j*h) 
    end
    y = h*contador
endfunction

// Trapecio para f de 2 variables con x fijo
function y = trapecioCompDoblex(f,xi,a,b,n)
    h = (b-a)/n
    contador = f(xi,a)*1/2 + f(xi,b)*1/2
    for j = 1:n-1
        contador = contador + f(xi,a+j*h) 
    end
    y = h*contador
endfunction

// Trapecio para f de 2 variables con y fijo
function y = trapecioCompDobley(f,yi,a,b,n)
    h = (b-a)/n
    contador = f(a,yi)*1/2 + f(b,yi)*1/2
    for j = 1:n-1
        contador = contador + f(a+j*h,yi) 
    end
    y = h*contador
endfunction

// Metodo de simpson simple
function y = simpson(f,a,b)
    h = (b-a)/2
    c = a+h
    contador = h/3
    funciones = f(a)+4*f(c)+f(b)
    y = contador*funciones
endfunction

// Metodo de simpson compuesto para n subdivisiones
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

// Metdodo de simpson para f de dos variables con x fija
function y = simpsonCompDoblex(f,xi,a,b,n)
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

// Metdodo de simpson para f de dos variables con y fija
function y = simpsonCompDobley(f,yi,a,b,n)
    if modulo(n,2) <> 0 then
        disp("N no es par")
        y = %nan
        return
    end
    h = (b-a)/n
    y = f(a,yi)+f(b,yi)
    for i = 1:n-1
        if modulo(i,2) == 0 then
            y = y+2*f(a+i*h,yi)
        else
            y = y+4*f(a+i*h,yi)
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

// Metodo del Trapecio para f de 2 variables
// Con dydx (Primero se integra 'y' luego 'x')
function y = TrapecioBi1(f,a,b,c,d,n,m)
    // a,b valores de integral exterior
    // c,d valores de integral interior (Se reciben como funcion, si se quiere constante ingresar funciones constantes)
    // n intervalos entre a y b
    // m intervalos en c(x) d(x)
    hx = (b-a)/n
    for i = 0:n
        xi = a+i*hx
        G(i+1) = trapecioCompDoblex(f,xi,c(xi),d(xi),m)
    end
    contador = G(1)*1/2 + G(n+1)*1/2
    for j = 1:n-1
        contador = contador + G(j+1) 
    end
    y = hx*contador
endfunction

// Metodo del Trapecio para f de 2 variables
// Con dxdy (Primero se integra 'x' luego 'y')
function y = TrapecioBi2(f,a,b,c,d,n,m)
    // a,b valores de integral exterior
    // c,d valores de integral interior (Se reciben como funcion, si se quiere constante ingresar funciones constantes)
    // n intervalos entre a y b
    // m intervalos en c(y) d(y)
    hy = (b-a)/n
    for i = 0:n
        yi = a+i*hx
        G(i+1) = trapecioCompDobley(f,yi,c(yi),d(yi),m)
    end
    contador = G(1)*1/2 + G(n+1)*1/2
    for j = 1:n-1
        contador = contador + G(j+1) 
    end
    y = hy*contador
endfunction

// Metodo de Simpson para f de 2 variables
// Con dydx (Primero se integra 'y' luego 'x')
function y = SimpsonCompBi(f,a,b,c,d,n,m)
    // a,b valores de integral exterior
    // c,d valores de integral interior (Se reciben como funcion, si se quiere constante ingresar funciones constantes)
    // n intervalos entre a y b
    // m intervalos en c(x) d(x)
    // m y n deben ser par
    if modulo(n,2) <> 0 then
        disp("N no es par")
        y = %nan
        return
    end
    hx = (b-a)/n
    for i = 0:n
        xi = a+i*hx
        G(i+1) = simpsonCompDoblex(f,xi,c(xi),d(xi),m)
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

// Metodo de Simpson para f de 2 variables
// Con dxdy (Primero se integra 'x' luego 'y')
function y = SimpsonCompBi(f,a,b,c,d,n,m)
    // a,b valores de integral exterior
    // c,d valores de integral interior (Se reciben como funcion, si se quiere constante ingresar funciones constantes)
    // n intervalos entre a y b
    // m intervalos en c(y) d(y)
    // m y n deben ser par
    if modulo(n,2) <> 0 then
        disp("N no es par")
        y = %nan
        return
    end
    hy = (b-a)/n
    for i = 0:n
        yi = a+i*hx
        G(i+1) = simpsonCompDobley(f,yi,c(xi),d(xi),m)
    end
    y = G(1)+G(n+1)
    for i = 1:n-1
        if modulo(i,2) == 0 then
            y = y+2*G(i+1)
        else
            y = y+4*G(i+1)
        end
    end
    y = y*hy/3
endfunction
