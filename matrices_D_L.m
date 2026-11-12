function [D,L] = matrices_D_L(K,h)
  w = ones(K-1,1);
  D = (2*diag(ones(K,1)) - diag(w,1) - diag(w,-1))/h;
  L = (4*diag(ones(K,1)) + diag(w,1) + diag(w,-1))*h/6;