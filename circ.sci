//Arreglar lo de centrar y adecuar el radio, por ahora lo usamos pasando como parametro 2*r
function c=circ(r,x,y)
    rect = [x-r,y-r,x+r,y+r]
    p = 0
    plot2d(p,rect=rect)
    xarc(x-r/2,y+r/2,r,r,0,360*64)
    c = 0
endfunction

function z = Gersch(A)
    [nA,mA] = size(A)
    //Obtenemos los centros y radios
    for i = 1:nA
        centro(i) = A(i,i)
        radio(i) = 0
        for j = 1:nA
            if j <>i
                radio(i) = radio(i) + abs(A(i,j))
            end
        end
    end
    disp(centro)
    //Dibujamos los circulos
    centroMenor = min(centro)
    centroMayor=max(centro)
    radioMayor = max(radio)
    rect = [centroMenor-radioMayor,-radioMayor,centroMayor+radioMayor, radioMayor]
    //plot2d(0,rect=rect)
    for i = 1:nA
        circ(2*radio(i),centro(i),0)    
    end
    z=0
endfunction
