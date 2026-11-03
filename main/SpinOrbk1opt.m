function [x,psi0,E0,mu,k1,output]=SpinOrbk1opt(G1,G2,k0,Omega,M,k1in,dk,ksteps,psiin)
% function to vary k1 in the determination of a spin-orbit coupled ground
% state
%
% usage: 
% [x,psi,E,mu,k1,output]=SpinOrbk1opt(G1,G2,k0,Omega,M,k1in,dk,ksteps)
% k1in - the initial value of k1
% dk  - the initial dk used to find the minimum
% ksteps - the number of steps to take 

%% solution at initial k1 point
if(nargin==9)
    [x,psi0,Ei,mu,output]=SpinOrb_gflow(G1,G2,k0,Omega,k1in,M,'psiin',psiin,'disp',1e4);
else
    [x,psi0,Ei,mu,output]=SpinOrb_gflow(G1,G2,k0,Omega,k1in,M,'disp',1e4);
end
EpN0=Ei*k1in/(2*pi); % energy per particle is E/L
k1=k1in;
disp(['init val (EpN,k)=(',num2str(EpN0,8),',',num2str(k1,4),')'])
for pcks=1:ksteps
    k1v= k1+[-1,1]*dk;
    [x,psitmp,Ev(1),mu,output]=SpinOrb_gflow(G1,G2,k0,Omega,k1v(1),M,'psiin',psi0,'disp',1e4);
    [x,psitmp,Ev(2),mu,output]=SpinOrb_gflow(G1,G2,k0,Omega,k1v(2),M,'psiin',psi0,'disp',1e4);
    EpNv=Ev.*k1v/(2*pi);
    dEpNk = (EpNv(2)- EpNv(1))/(2*dk);
    d2EpNk = (EpNv(2)-2*EpN0+EpNv(1))/(dk^2);
    deltak=-dEpNk/d2EpNk;
    k1=k1+deltak;
    [x,psi0,E0,mu,output]=SpinOrb_gflow(G1,G2,k0,Omega,k1,M,'psiin',psi0,'disp',1e4);
    EpN0=E0*k1/(2*pi);
    if(abs(deltak)<2*dk);dk=dk/4;end
    disp(['Result of k1 optimisaiton, step#',num2str(pcks)])
    disp([' step #',num2str(pcks),' (EpN,k)=(',num2str(EpN0,8),',',num2str(k1,4),')'])





end




