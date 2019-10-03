//Recibe una funcion, un valor, un orden y un paso
//Consume MUCHA memoria
function x = derivar(f,v,n,h)
    if (n ==0) then
        x = f(v)
    else
    x = (derivar(f,v+h,n-1,h)-derivar(f,v,n-1,h))/h
    end
endfunction

//Recibe una funcion, un valor, un orden y un paso
//Metodo iterativo
function x = derivarIterativo(f,v,n,h)
    deff("y=D0F(x)", "y="+f);
    for i = 1:1:n-1
        deff("y=D"+string(i)+"F(x)", "y=(D"+string(i-1)+"F(x+h)-D"+string(i-1)+"F(x))/h");
    end
    deff ("y=DnF(x)", "y=(D"+string(n-1)+"F(x+h)-D"+string(n-1)+"F(x))/h");
    x = DnF(v);
endfunction

//Evaluacion eficiente utilizando el metodo de Horner
function x = evaluarHorner(p,a)
    coefAnterior = coeff(p,degree(p));
    coefActual = 0;
    for i = (degree(p)-1):-1:0
        coefActual = coeff(p,i);
        coefAnterior = (coefActual + a*coefAnterior);
    end
    x = coefAnterior
endfunction

//Evalua upolinomio en un punto y devuelve el valor y su derivada

function x = derivarYEvaluarHorner(p,a)
    coefAnterior = coeff(p,degree(p));
    coefActual = 0;
    coeffDerivado = 0;
    for i = (degree(p)-1):-1:0
        coeffDerivado = coeffDerivado + coefAnterior * (a**(i));
        coefActual = coeff(p,i);
        coefAnterior = (coefActual + a*coefAnterior);
    end
    x(1) = coefAnterior;
    x(2) = coeffDerivado;
endfunction

function x = taylor(f,n,a,v)
    x=f(a);
    for i = 1:n
        x = x + (derivar(f,a,i,0.0001)*(v-a))/factorial(i);
    end
endfunction
