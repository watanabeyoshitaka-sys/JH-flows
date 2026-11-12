function [uh,b] = Newton_Raphson(K,h,d1,d2,d3,uh,b,sp_node,D,L,TOL,MAXIT)
  H = h*ones(K,1);
  for k = 1:MAXIT
    u = [uh;b];
    G2 = matrix_G2(K,h,uh,sp_node);
    f = [(D-d1*L-d2*G2)*uh + H*b; H'*uh-d3];
    G = [D-d1*L-2*d2*G2, H; H', 0];
    u = [uh;b];
    w = u - G\f;
    uh = w(1:K); b = w(K+1);
    t = norm(f,inf)/norm(w,inf);
    if (norm(u,inf)<1e-2), t = norm(f,inf); end
    fprintf('step:%d',k);
    fprintf('  relative (or absolute) error: %d\n',t);
    if t < TOL
      disp('Newton-Raphson iteration converges numerically.')
      % save("G.mat","G")
      return
    end
    if (k==MAXIT) || (t>=1e3) || any(isnan(w))
      error('Newton-Raphson iteration does not converge.')
    end
  end