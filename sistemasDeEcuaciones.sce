// Funciones de scilab para solucionar problemas de sistemas de ecuaciones

// TriangularSuperior : Matriz Vector -> Vector
// Recibe matriz A triangular superior sin 0 en su diagonal
function x = TriangularSuperior(A,b)

    if nA<>mA then
    error('Error - La matriz A debe ser cuadrada')
    abort
    elseif mA<>nb then
    error('Error - dimensiones incompatibles entre A y b')
    abort
    end

    t = size(A)
    t = t(1)
    x(t) = b(t) / A(t,t) 
    for i = (t-1):-1:1
        suma = 0
        for j = i+1:t 
            suma = suma + A(i,j) * x(j)
        end
        x(i) = (b(i) - suma) / A(i,i)
    end
endfunction

// TriangularInferior : Matriz Vector -> Vector
// Recibe matriz A triangular inferior sin 0 en su diagonal
function x = TriangularInferior(A,b)

    if nA<>mA then
    error('Error - La matriz A debe ser cuadrada')
    abort
    elseif mA<>nb then
    error('Error - dimensiones incompatibles entre A y b')
    abort
    end

    x(1) = b(1) / A(1,1) 
    tamano = size(A)
    for i = 2:(tamano(1))
        suma = 0
        for j = 1:i-1 
            suma = suma + A(i,j) * x(j)
        end
        x(i) = (b(i) - suma) / A(i,i)
    end
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
contador = 0
// Eliminación progresiva
for i = 1:(nA-1)
    for j = (i+1):nA
        mjk = a(j,i)/a(i,i)
        a(j,i)=0
        a(j,(i+1):(nA+mb)) = a(j,(i+1):(nA+mb)) - mjk*a(i,(i+1):(nA+mb))
        contador = contador + 1
    end
end

// Sustitución regresiva
disp(a)
for k = 1: mb
    x(nA,k) = a(nA,nA+k)/a(nA,nA) 
    contador = contador +1
    for i = (nA-1):-1:1
        suma = 0
        for j = i+1:nA 
            suma = suma + a(i,j)*x(j,k)
            contador = contador + 1
        end
        x(i,k) = (a(i,nA+k) - suma) / a(i,i)
        contador = contador + 1
    end
    end
    disp(contador)
endfunction

function x = inversa(A)
    [nA,mA] = size(A)
    if nA<>Ma then
        error('La matriz A debe ser cuadrada')
        abort
    end
    i = eye(nA,nA)
    x = gausselim(A,i)
endfunction

function x = determinante(a)
    [nA,mA] = size(a)
    if nA<>mA then
    error('determinante - La matriz A debe ser cuadrada')
    abort
    end
    x=1
    for i = 1:(nA)
        for j = (i+1):nA
            mjk = a(j,i)/a(i,i)
            a(j,i)=0
            a(j,(i+1):nA) = a(j,(i+1):nA) - mjk*a(i,(i+1):nA)
        end
        x = x * a(i,i)
    end
endfunction

function [x,a,P] = gausselimPP(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana con pivoteo parcial.

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
n = nA    // Tamaño de la matriz
P = eye(n,n)
// Eliminación progresiva con pivoteo parcial
for k=1:n-1
    [v,i]=max(abs(a(k:n,k)))
    kpivot = k-1+i
    temp = a(kpivot,:); a(kpivot,:) = a(k,:); a(k,:) = temp
    temp = P(kpivot,:); P(kpivot,:) = P(k,:); P(k,:) = temp
    for i=k+1:n
        for j=k+1:n+1
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k)
        end
        for j=1:k        // no hace falta para calcular la solución x
            a(i,j) = 0  // no hace falta para calcular la solución x
        end              // no hace falta para calcular la solución x
    end
end

// Sustitución regresiva
x(n) = a(n,n+1)/a(n,n)
for i = n-1:-1:1
    sumk = 0
    for k=i+1:n
        sumk = sumk + a(i,k)*x(k)
    end
    x(i) = (a(i,n+1)-sumk)/a(i,i)
end
endfunction

// FALTA HACER MATRIZ TRIDIAGONAL

function [L,U,P] = FactorizacionLU(A)
    U = A
    [nA,mA] = size(A)
    disp(nA)
    disp(mA)
    L = eye(nA,nA)
    P = eye(nA,nA)
    for k = 1:mA-1
        maxI = k
        for i = k+1:nA
            if abs(U(i,k))> abs(U(maxI,k))
                maxI = i
            end
            end
            Aux = U(k, k:mA)
            U(k,k:mA) = U(maxI, k:mA)
            U(maxI, k:mA) = Aux
            Aux = L(k, 1:k-1)
            L(k, 1:k-1) = L(maxI,1:k-1)
            L(maxI,1:k-1) = Aux
            Aux = P(k,:)
            P(k,:) = P(maxI,:)
            P(maxI,:) = Aux
            for j = k+1:mA
                L(j,k) = U(j,k)/U(k,k)
                U(j, k:mA) = U(j, k:mA) - (L(j,k)*U(k,k:mA))
            end
         
    end
endfunction

function x = ResolverLUP(L,U,P,b)
    B = P*b
    y = TriangularInferior(L,B)
    x = TriangularSuperior(U,y)
endfunction

function [U,ind] = cholesky(A,b)
// Factorización de Cholesky.
// Trabaja únicamente con la parte triangular superior.
//
// ind = 1  si se obtuvo la factorización de Cholesky.
//     = 0  si A no es definida positiva
//
//******************
eps = 1.0e-8
//******************

n = size(A,1)
U = zeros(n,n)
A = [A b]
for k = 1:n
    t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
    if t <= eps then
        printf('Matriz no definida positiva.\n')
        ind = 0
        return
    end
    U(k,k) = sqrt(t)
    for j = k+1:n
        U(k,j) = ( A(k,j) - U(1:k-1,k)'*U(1:k-1,j) )/U(k,k)
    end
end
ind = 1

endfunction
