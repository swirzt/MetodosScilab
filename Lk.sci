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
