function y = jacobi(A,b,x,eps)
    [nA,mA] = size(A)
    I = eye(nA,mA)
    for k=1:nA-1
        [v,i]=max(abs(A(k:nA,k)))
            kpivot = k-1+i
            temp = A(kpivot,:); A(kpivot,:) = A(k,:); A(k,:) = temp
            temp = b(kpivot,:); b(kpivot,:) = b(k,:); b(k,:) = temp
            temp = x(kpivot,:); x(kpivot,:) = x(k,:); x(k,:) = temp
        end
    N = diag(diag(A))
    inversaN = inversa(N)
    Norma = I-inversaN*A
    n(1) = norm(Norma, 1)
    n(2) = norm(Norma, 'inf')
    n(3) = norm(Norma, 'fro')
    n(4) = norm(Norma)
    if min(n) >= 1 then
        if max(abs(spec(A))) > = 1
            disp("La solucion no converge para todo punto inicial") 
        end
    end
    contador = 0
        while(norm(A*x-b) > eps)
            xAnt = x
            for i = 1:nA
                x(i) = b(i)
                suma = 0
                for j = 1:nA
                    if j <> i then
                        suma = suma + A(i,j)*xAnt(j)
                    end
                end
                x(i) = (x(i) - suma)/ A(i,i)
            end 
            contador = contador+1
        end
    y = x
    disp(contador)
endfunction

function y = gauseidel(A,b,x,eps)
    [nA,mA] = size(A)
    if diagonalDominante(A) == 0 then
        I = eye(nA,mA)
        for k=1:nA-1
            [v,i]=max(abs(A(k:nA,k)))
                kpivot = k-1+i
                temp = A(kpivot,:); A(kpivot,:) = A(k,:); A(k,:) = temp
                temp = b(kpivot,:); b(kpivot,:) = b(k,:); b(k,:) = temp
                temp = x(kpivot,:); x(kpivot,:) = x(k,:); x(k,:) = temp
            end
        N = A
        for i = 1:nA-1
            for j = i+1:nA
                N(i,j) = 0
            end
        end
        inversaN = inversa(N)
        Norma = I-inversaN*A
        n(1) = norm(Norma, 1)
        n(2) = norm(Norma, 'inf')
        n(3) = norm(Norma, 'fro')
        n(4) = norm(Norma)
        if min(n) >= 1 then
            if max(abs(spec(A))) > = 1
                disp("La solucion no converge para todo punto inicial") 
            end
        end
    end
    contador = 0
    while(norm(A*x-b) > eps)
         for i = 1:nA
            for i = 1:nA
                x(i) = b(i)
                suma = 0
                for j = 1:nA
                    if j <> i then
                        suma = suma + A(i,j)*x(j)
                    end
                end
                x(i) = (x(i) - suma)/ A(i,i)
            end 
         end
         contador = contador+1
     end
     y=x
     disp(contador)
endfunction

function x =diagonalDominante(A)
    [nA,mA] = size(A)
    for i = 1:nA
        suma = 0
        for j = 1:mA
            if j <> i then
                suma = suma + abs(A(i,j))
            end
        end
        if suma >= abs(A(i,i))
            x = 0
            return
        end 
    end
    x = 1
endfunction

function x = inversa(A)
    [nA,mA] = size(A)
    if nA<>mA then
        error('La matriz A debe ser cuadrada')
        abort
    end
    i = eye(nA,nA)
    x = gausselim(A,i)
endfunction

function [x,a] = gausselim(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dadas las matrices de coeficientes A y b.
// La función implementa el método de Eliminación Gaussiana sin pivoteo.  

[nA,mA] = size(A) 
[nb,mb] = size(b)

if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada')
    abort
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b')
    abort
end

a = [A b] // Matriz aumentada
// Eliminación progresiva
for i = 1:(nA-1)
    for j = (i+1):nA
        mjk = a(j,i)/a(i,i)
        a(j,i)=0
        a(j,(i+1):(nA+mb)) = a(j,(i+1):(nA+mb)) - mjk*a(i,(i+1):(nA+mb))
    end
end
for k = 1: mb
    x(nA,k) = a(nA,nA+k)/a(nA,nA) 
    for i = (nA-1):-1:1
        suma = 0
        for j = i+1:nA 
            suma = suma + a(i,j)*x(j,k)

        end
        x(i,k) = (a(i,nA+k) - suma) / a(i,i)
    end
    end

endfunction
