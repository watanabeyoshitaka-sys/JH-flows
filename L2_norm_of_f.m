function d = L2_norm_of_f(K,h,d1,d2,u,b,el_node)
  I = intval(zeros(4,1)); s = I; c0 = I(1:2);
  for k=1:K+1
    m = el_node(k,:);  t = find(m>0); 
    c = c0;  c(t) = u(m(t));
    I(1) = c(1) + c(2);
    for J = 2:4
      I(J) = c(1)*I(J-1) + c(2)^J;
    end
    s = s + I;
  end
  p = [-d1*b (d1^2-2*d2*b)/3 d1*d2/2 (d2^2)/5]*s;
  d = sqrt(intval(mag(2*b^2+h*p)));