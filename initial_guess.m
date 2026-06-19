function [uh,b] = initial_guess(N,R,flag)
  if flag == 1
    r = linspace(-1,1,N+1); x = r(2:N)+1;
    uh = -100*mid(R)*sin(2*x*pi)';
    %uh = 100*mid(R)*sin(x*pi)';
    %uh = mid(R)*sin(x*pi/2)';
    %uh = -99*cos(x*pi/6)';
    if (R==0), b=1; else, b=mid(R); end  % case of R=0 should be avoided!
  else
    load("app.mat")
    np = size(uh,1)+1; xp = -1+2/np:2/np:1-2/np; x = -1+2/N:2/N:1-2/N;
    uh = interp1(xp,uh,x,'pchip'); uh=uh';
  end

  