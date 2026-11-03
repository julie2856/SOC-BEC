function [x,k,dx,dk,U]=grid1D(L,N,k1)

%%% note different from using grid 1D function
% this assumes L=pi/k1 and that the functions are expanded in the
% planewaves:
% k1 +n*2*k1, n=0,\pm1,\pm2,.. [Eq. (3) Yun Li PRL (2013)]

x=linspace(-0.5*L,0.5*L,N+1);
x=x(1:N);
dx=x(2)-x(1);
K=2*pi/dx;

k=linspace(-0.5*K,0.5*K,N+1);
k=k(1:N);dk=k(2)-k(1);


%% for odd grids, shift to ensure x=0, k=0
if(round(N/2)~=N/2)
    k=k+0.5*dk;
    x=x+0.5*dx;
end
if(nargin==2)
    k1=pi/L;
end
k=k-k1;

U=exp(-1i*k(:)*x)/sqrt(N);

 

 