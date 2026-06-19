function G2 = matrix_G2_intval(K,h,u,sp_node)
  w1 = intval(zeros(K,1)); w2 = intval(zeros(K-1,1));
  G2 = diag(w1) + diag(w2,1) + diag(w2,-1);
  s = h/12; b0 = intval(zeros(3,1));  
  E = [ 1 1 0; 1 6 1; 0 1 1];
  for k = 1:K
    m = sp_node(k,:); t = find(m>0);
    b = b0; b(t) = u(m(t));
    v = [b(1) b(2) b(3)]; U = E*v';
    G2(k,m(t)) = U(t)';
  end
  G2 = s*G2;