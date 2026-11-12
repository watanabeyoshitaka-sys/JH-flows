function s = inf_norm(K,d1,d2,u,sp_node)
  w = intval(zeros(K,1));  c0 = intval(zeros(3,1)); 
  for k=1:K
    m = sp_node(k,:); t = find(m>0); c = c0; c(t) = u(m(t));
    w(k) = max(mag(d1+2*d2*c));
  end
  s = intval(max(mag(w)));