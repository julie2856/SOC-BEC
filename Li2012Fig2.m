clear;
%% Fig. 2 of Li et al., PRL 108 225301 (2012)
 
 
k0=1; 
gamma=0.0012; 
G1c=k0^2/(4*gamma*(1+gamma));
G2c=gamma*G1c;
  
  
n=linspace(1e-2,1.5,100);
Omegav=linspace(0,2.2*k0^2,100);   

OmegaI_II=zeros(1,length(n));
OmegaII_III=zeros(1,length(n));
OmegaI_III=zeros(1,length(n));

k1=zeros(length(Omegav),length(n));
sigmaz=zeros(length(Omegav),length(n));
nc=1;
for j=1:length(n)
    G1=n(j)*G1c;
    G2=n(j)*G2c;
    
    % phase boundaries
    OmegaI_II(j)=2*sqrt( (k0^2+G1)*(k0^2-2*G2) * 2*G2/(G1+2*G2));
    OmegaII_III(j)=2*(k0^2-2*G2);
    OmegaI_III(j)=2*(k0^2+G1)-2*sqrt( (k0^2+G1)*G1 );
     
    for l=1:length(Omegav)    
        %% I changed this =SpinOrb_var2d->=SpinOrb_var2d1
        [k1(l,j),EpN,beta(l,j),other]=SpinOrb_var(k0,G1,G2,Omegav(l));
        sigmaz(l,j)=other.sigmaz;
        EI(l,j)=other.Ph1.E;
        EIII(l,j)=other.Ph3.E;
        k1I(l,j)=other.Ph1.k1;
        k1III(l,j)=other.Ph3.k1;
    end
end

figure('Position',[300,100,700,700])

subplot(2,1,1)
    imagesc(n./nc,Omegav./k0^2,sigmaz);hold on
    plot(n(n<=nc)/nc,OmegaI_II(n<=nc)/k0^2,'w','LineWidth',2)
    plot(n(n<=nc)/nc,OmegaII_III(n<=nc)/k0^2,'w','LineWidth',2)
    plot(n(n>=nc)/nc,OmegaI_III(n>=nc)/k0^2,'w','LineWidth',2)
        set(gca,'fontsize',16,'YTick',0:1:2,'YDir','normal')
        ylabel('$\Omega/k_0^2$','interpreter','latex','fontsize',16)
        xlabel('$n/n^{(c)}$','interpreter','latex','fontsize',16)
        axis([0,1.5,0,2.2]); v=axis;
        text(0.85*(v(2)-v(1))+v(1), 0.85*(v(4)-v(3))+v(3),'(a)','FontSize',16,'Color','w')
    colormap jet
    colorbar('northoutside')   

subplot(2,1,2)
    imagesc(n./nc,Omegav./k0^2,k1./k0);hold on
    plot(n(n<=nc)/nc,OmegaI_II(n<=nc)/k0^2,'w','LineWidth',2)
    plot(n(n<=nc)/nc,OmegaII_III(n<=nc)/k0^2,'w','LineWidth',2)
    plot(n(n>=nc)/nc,OmegaI_III(n>=nc)/k0^2,'w','LineWidth',2)
        set(gca,'fontsize',16,'YTick',0:1:2,'YDir','normal')
        ylabel('$\Omega/k_0^2$','interpreter','latex','fontsize',16)
        xlabel('$n/n^{(c)}$','interpreter','latex','fontsize',16)
        axis([0,1.5,0,2.2]); v=axis;
        text(0.85*(v(2)-v(1))+v(1), 0.85*(v(4)-v(3))+v(3),'(b)','FontSize',16,'Color','w')
 
toc

Npts = 256;  

% try to make the Li colour map 
cmap_points = [
    0.5 0 0;  % dark red
    1 0 0;  % red
    1 1 0   % yellow
    ]; 
x = linspace(1, size(cmap_points,1), Npts); 
colormap(interp1(1:size(cmap_points,1), cmap_points, x)); 