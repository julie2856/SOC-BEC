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
M=100;

%% optimise a GP solution and vary k1
[x,psi,E,mu,k1,output]=SpinOrbk1opt(G1,G2,k0,Omega,M,k1in,dk,ksteps);
 
[~,kx,~,~,UGP]=grid1D(pi/k1,M,k1);
%% find the variational solution
[k1var,EpNvar,beta,other]=SpinOrb_var(k0,G1,G2,Omega);
  
[xv,kv,~,~,Uv]=grid1D(pi/k1var,M,k1var);
c=sqrt((1+(k1var/k0))/2);s=sqrt((1-(k1var/k0))/2);    
psivara= (other.C1*c*exp(1i*k1var*xv)+other.C2*s*exp(-1i*k1var*xv));
psivarb= (-other.C1*s*exp(1i*k1var*xv)-other.C2*c*exp(-1i*k1var*xv));


%% iterpolate
% final grid  
[xf,~,~,~,Uf]=grid1D(2*2.2*pi/k0,1000,k1); % match Martone 
% take solutions to k-space and then back to final grids
phia=UGP*psi(1:M);
phib=UGP*psi(1+M:end);
phivara=Uv*psivara(:);
phivarb=Uv*psivarb(:);
UfGP=(1/sqrt(M))*exp(-1i*kx'*xf);
Ufv=(1/sqrt(M))*exp(-1i*kv'*xf);
psifa=UfGP'*phia;
psifb=UfGP'*phib;
psivarfa=Ufv'*phivara;
psivarfb=Ufv'*phivarb; 
  

 
 plot(xf*k0/pi,abs(psifa).^2,'b-','LineWidth',1);hold on
 plot(xf*k0/pi,abs(psivarfa).^2,'r--','LineWidth',1);hold on
 plot(xf*k0/pi,abs(psifb).^2,'b.','LineWidth',1);hold on
 plot(xf*k0/pi,abs(psivarfb).^2,'r.','LineWidth',1);hold on
 

 ylim([0.25,0.75] )
 xlim([-1,1]*2.2)
 set(gca,'fontsize',12,'ytick',[0.3:0.2:0.7],'xtick',[-2:2])
 xlabel('$k_0x/\pi$','FontSize',16,'interpreter','latex')
 ylabel('$n_{\uparrow,\downarrow}/\bar{n}$]','FontSize',16,'interpreter','latex')



 