function x = metodoRobusto(p)
    //Solucion de ecuaciones cuadraticas con discriminante positivo
    c = coeff(p,0)
    b = coeff(p,1)
    a = coeff(p,2)
    if((b*b-4*a*c)>0) then
        if(b < 0) then
            x(1) = (2*c)/(-b + sqrt(b^2 - 4*a*c))
            x(2) = (-b + sqrt(b^2 - 4*a*c)) / (2*a)
        elseif (b > 0) then
            x(1) = (-b-sqrt(b^2-4*a*c))/(2*a)
            x(2) = (2*c)/(-b - sqrt(b^2 - 4*a*c))
        else
            x(1) = (sqrt(-c*a))
            x(2) = -(sqrt(-c*a))
        end
    else
        x(1) = %nan
        x(2) = %nan
    end
endfunction

function x = evaluarHorner(p,a)
    //Evaluacion eficiente utilizando el metodo de Horner
    coefAnterior = coeff(p,degree(p))
    coefActual = 0
    for i = (degree(p)-1):-1:0
        coefActual = coeff(p,i)
        coefAnterior = (coefActual + a*coefAnterior)
    end
    x = coefAnterior
endfunction

function x = derivarYEvaluarHorner(p,a)
    //Evalua un polinomio en un punto y devuelve el valor y su derivada
    coefAnterior = coeff(p,degree(p))
    coefActual = 0
    coeffDerivado = 0
    for i = (degree(p)-1):-1:0
        coeffDerivado = coeffDerivado + coefAnterior * (a**(i))
        coefActual = coeff(p,i)
        coefAnterior = (coefActual + a*coefAnterior)
    end
    x(1) = coefAnterior
    x(2) = coeffDerivado
endfunction

function e = errorAbsoluto(x,y)
    e(1) = abs(x(1)-y(1))
    e(2) = abs(x(2)-y(2))
endfunction

function e = errorRelativo(x,y)
    e(1) = abs(x(1)-y(1))/abs(x(1))
    e(2) = abs(x(2)-y(2))/abs(x(2))
endfunction

function x = derivar(f,v,n,h)
    //Recibe una funcion, un valor, un orden y un paso
    //Consume MUCHA memoria
    if (n ==0) then
        x = f(v)
    else
    x = (derivar(f,v+h,n-1,h)-derivar(f,v,n-1,h))/h
    end
endfunction

function x = derivarIterativo(f,v,n,h)
    //Recibe una funcion, un valor, un orden y un paso
    //Metodo iterativo
    deff("y=D0F(x)", "y="+f)
    for i = 1:1:n-1
        deff("y=D"+string(i)+"F(x)", "y=(D"+string(i-1)+"F(x+h)-D"+string(i-1)+"F(x))/h")
    end
    deff ("y=DnF(x)", "y=(D"+string(n-1)+"F(x+h)-D"+string(n-1)+"F(x))/h")
    x = DnF(v)
endfunction

//-*- Epsilon
e = 0.0001
//-*-
function x = taylor(f,n,a,v,e)
    x=f(a)
    for i = 1:n
        x = x + (derivar(f,a,i,e)*(v-a))/factorial(i)
    end
endfunction
