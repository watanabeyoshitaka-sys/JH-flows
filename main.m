format long;                 
intvalinit('DisplayInfsup'); 
intvalinit('SharpIVmult');   
Pi = intval('pi'); one = intval('1'); % some constants
a = 4*Pi/16;           % alpha
R = -3*one; % Reynolds number
%R = -4*one;
gn = 2;  % initial guess number
eq = 3;
switch eq 
  case 1
    d1 = 4*a^2; d2 = a*R; d3 = 2;  sc=mid(R/a);
  case 2
    d1 = 4*a^2; d2 = 2*a; d3 = R;  sc=mid(2/a);
  otherwise
    d1 = 4*a^2; d2 = a^2; d3 = (2*R)/a; sc=1;
end
N = 2000;            % partition number of [-1,1]
h = intval('2')/N;  % width of the partition
K = N-1;            % system dimension of Sh
MAXIT = 50;         % maximum iteration number
TOL = N^2*eps;      % convergence criterion

disp(['      N = ' num2str(N)'']) 
disp(['      R = ' num2str(mid(R))'']) 
disp(['      a = ' num2str(mid(a))'']) 
[el_node,sp_node] = generate_local_node(N,K);
[uh,b] = initial_guess(N,R,gn);  

[D,L] = matrices_D_L(K,h);
[uh,b] = Newton_Raphson(K,mid(h),mid(d1),mid(d2),mid(d3),uh,b,sp_node,mid(D),mid(L),TOL,MAXIT);
plot_approximate_solution(N,mid(h),mid(R),mid(a),uh,eq)
save("app.mat","uh","b")
disp(['maximum of uh = ' num2str(max(mag(uh)))''])
disp(['translated f(0) = ' num2str(mid(uh(N/2)*sc))''])
[u,beta] = FN_Int(K,h,d1,d2,d3,uh,b,el_node,sp_node,D,L,MAXIT);