clear;
%% This script reproduce Fig.5.1 of Martone's PhD Thesis.
% It compares the variational solution to the GP solution for a stripe
% phase
Omega=1;
k0=1; 
G1=0.3;
G2=0.08; 
k1in=1;
dk=0.1;
ksteps=2;
M=20;

%% optimise a GP solution and vary k1
[x,psi,E,mu,k1,output]=SpinOrbk1opt(G1,G2,k0,Omega,M,k1in,dk,ksteps);
  
%% make better state!
[x,psi,E,mu,output]=SpinOrb_gflow(G1,G2,k0,Omega,k1,M,'resid',1e-11,'disp',50,'psiin',psi,'dt',10e-2);
 
%% BdG
qxv=linspace(1e-3,2*k1-1e-3,64);Nbands=12;
[ev,u1v,u2v,v1v,v2v] = SpinOrb_BdG(qxv,Nbands,Omega,G1,G2,psi,mu,k1);

 
plot(qxv/k1,ev,'-x')
  
 ylim([0,1.8])
 xlim([0,2])
 set(gca,'fontsize',12,'ytick',[0.3:0.2:0.7],'xtick',[-2:2])
 xlabel('$q_x/k_1$','FontSize',16,'interpreter','latex')
 ylabel('$\omega_{\nu q}$','FontSize',16,'interpreter','latex')



 