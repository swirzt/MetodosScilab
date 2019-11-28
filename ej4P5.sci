Ns = [500 1000]
exps = [-6 -12];
 
for N = Ns
  for exp = exps
    eps = 10^(exp);
    printf("N = %d, eps = 1e%d \n", N, exp);
 
 
    A = 8*eye(N,N) + 2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1) + diag(ones(N-3,1),3) + diag(ones(N-3,1),-3);
    b = ones(N,1);
 
    tic();
    [x, a, P] = gausselimPP(A, b)
    t = toc();
 
    printf("resolver_gauss_pp: %f s \n", t);
 
    tic();
    x = gauseidelItera(A, b, zeros(N, 1), eps);
    t = toc();
    printf("resolver_gauss_seidel: %f s \n", t);
    // disp(x);
  end
end
