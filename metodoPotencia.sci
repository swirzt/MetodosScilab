function [x,p] = metodoPotencia(A,z)
    if A == A' then
    z0 = z 
    w1 = A * z0
    z1 = w1 / norm(w1,'inf')
    [wk,k] = max(w1)
    l1 = wk / z0(k)
    w2 = A * z1
    z2 = w2 / norm(w2,'inf')
    [wk,k] = max(w2)
    l2 = wk / z1(k)
    w3 = A * z2
    z3 = w3 / norm(w3,'inf')
    p=3
    [wk,k] = max(w3)
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
