function [u,beta] = FN_Int(K,h,d1,d2,d3,uh,b,el_node,sp_node,D,L,MAXIT)
  uh = intval(uh); H = h*ones(K,1);
  G2 = matrix_G2_intval(K,h,uh,sp_node);     
  r1 = (-D+d1*L+d2*G2)*uh - H*intval(b);
  r2 = d3 - H'*uh;
  fprintf('  inf-norm of res1: %d\n',sup(max(mag(r1))));
  fprintf('      norm of res2: %d\n',sup(mag(r2)));
  G = [D-d1*L-2*d2*G2, H; H', 0];
  unit=infsup(-1,1);  d0=1d-2;  C=h/intval('pi'); st=sqrt(intval('2'));           
  u = unit*eps*ones(K+1,1);  beta0 = eps;                 
  for k = 1:MAXIT
  fprintf('== verification step:%d\n',k);
    u = (1+d0*unit)*u;  beta = (1+d0)*beta0;  % inflation  
    G2 = matrix_G2_intval(K,h,u(1:K),sp_node);
    w0 = d2*G2*u(1:K);
    wh = uh + u(1:K);
    fprintf('diam of wh: %d\n',max(diam(wh)));
    fprintf('diam of u: %d\n',max(diam(u)));
    w1 = C*beta*local_L2_norm(K,h,d1,d2,wh);
    fprintf('diam of w0: %d\n',max(diam(w0)));
    fprintf('diam of w1: %d\n',max(diam(w1)));
    w = [r1+w0+unit*(w1+abs(d2)*C/st*intval(beta)^2); ...
         r2+unit*C*beta*sqrt((2*h)/3)];
    fprintf('diam of w: %d\n',max(diam(w)));
    w =  G\w;
    fprintf('diam of w: %d\n',max(diam(w)));
    % ---- infinite dimensional part ---- %
    t0 = L2_norm_of_f(K,h,d1,d2,wh,b,el_node);  
    t1 = C*beta*inf_norm(K,d1,d2,wh,sp_node);
    t2 = abs(d2)*C*(intval(beta)^2)/st;  
    t3 = b+u(K+1)+d3-h*sum(wh)+2*unit*C*beta;
    fprintf('t0: %d\n',sup(t0));
    fprintf('t1: %d\n',sup(t1));
    fprintf('t2: %d\n',sup(t2));
    fprintf('t3: %d\n',sup(t3));
    beta0 = sup(C*( sqrt((t0 + t1 + t2)^2 + t3^2) ));
    fprintf('  inflated beta: %d\n',beta);
    fprintf('      next beta: %d\n',beta0);
    if beta > 1,  error('Give Up!');  end
    res = sum(in(w,u));
    if beta0 <= beta,  res = res + 1;  end 
    fprintf('contraction rate: %7.3f\n',res/(K+2));
    if res == K+2
      disp('Verification has been completed!')
      break;
    end
    u = w; 
  end
 
