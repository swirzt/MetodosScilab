function c=circ(r,x,y)
    rect = [x-r,y-r,x+r,y+r]
    p = 0
    xgrid(1)
    plot2d(p,"031",rect=rect)
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
            if j <>i
                posibleradio(1) = posibleradio(1) + abs(A(i,j))
            end
        end
        for j = 1:nA
           if j <>i
               posibleradio(2) = posibleradio(2) + abs(A(j,i))
           end
         radio(i) = min(posibleradio)
    end
    disp(centro)
    //Dibujamos los circulos
    centroMenor = min(centro)
    centroMayor=max(centro)
    radioMayor = max(radio)
    rect = [centroMenor-radioMayor,-radioMayor,centroMayor+radioMayor, radioMayor]
    //plot2d(0,rect=rect)
    for i = 1:nA
        circ(radio(i),centro(i),0)
    end
    z=0
endfunction
