function c = circ(r,x,y)
    rect = [x-r,y-r,x+r,y+r]
    p = 0
    xgrid(1)
    plot2d(p,rect=rect)
    xarc(x-r,y+r,r*2,r*2,0,360*64)
    c = 0
endfunction

function z = Gersch(A)
    [nA,mA] = size(A)
    //Obtenemos los centros y radios
    for i = 1:nA
        centro(i) = A(i,i)
        radio(i) = 0
        posibleradio(1) = 0
        posibleradio(2) = 0
        for j = 1:nA
            if j <>i then
                posibleradio(1) = posibleradio(1) + abs(A(i,j))
            end
        end
        for j = 1:nA
            if j <>i then
                posibleradio(2) = posibleradio(2) + abs(A(j,i))
            end
            radio(i) = min(posibleradio)
        end
    end
    disp(centro)
    //Dibujamos los circulos
    centroMenor = min(centro)
    centroMayor = max(centro)
    radioMayor = max(radio)
    rect = [centroMenor-radioMayor,-radioMayor,centroMayor+radioMayor, radioMayor]
    //plot2d(0,rect=rect)
    for i = 1:nA
        circ(radio(i),centro(i),0)
    end
    z=0
endfunction

function [x,p] = metodoPotencia(A,z)
    if A == A' then
    z0 = z 
    w1 = A * z0
    z1 = w1 / norm(w1,'inf')
    [wk,k] = max(w1)
    if wk == 0 then // en caso que el maximo factor en w sea 0
        [s,n] = size(w1)
        if k <> n then
            k = k+1
            wk = w1(k)
        else
            k = 1
            wk = w1(k)
        end
    end
    l1 = wk / z0(k)
    w2 = A * z1
    z2 = w2 / norm(w2,'inf')
    [wk,k] = max(w2)
    if wk == 0 then // en caso que el maximo factor en w sea 0
        [s,n] = size(w2)
        if k <> n then
            k = k+1
            wk = w2(k)
        else
            k = 1
            wk = w2(k)
        end
    end
    l2 = wk / z1(k)
    w3 = A * z2
    z3 = w3 / norm(w3,'inf')
    p=3
    [wk,k] = max(w3)
    if wk == 0 then // en caso que el maximo factor en w sea 0
        [s,n] = size(w3)
        if k <> n then
            k = k+1
            wk = w3(k)
        else
            k = 1
            wk = w3(k)
        end
    end
    l3 = wk / z2(k)
    r = (l3 - l2)/(l2-l1)
    while(r/(1-r)*(l3-l2) > 1e-12)
        z2 = z3
        w3 = A*z3
        z3 = w3 / norm(w3,'inf')
        l1 = l2
        l2 = l3
        [wk,k] = max(w3)
        l3 = wk / z2(k)
        r = (l3 - l2)/(l2-l1)
        p = p+1
    end
    x = l3
else 
   x = %nan
end
endfunction
