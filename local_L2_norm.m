function w = local_L2_norm(K,h,d1,d2,wh)
w = intval(zeros(K,1));
for I=1:K
  if (I==1); w0=0; else w0=wh(I-1); end
  w1=wh(I);
  if (I==K); w2=0; else w2=wh(I+1); end
  w(I) = sqrt(intval(mag(h*( -45*(d1^2) + (4*d2*w0+5*d1)^2 + (4*d2*w2+5*d1)^2 + ...
         3*(8*d2*w1+5*d1)^2 + 48*(d2^2)*w1*(w0+w2))/120)));
end