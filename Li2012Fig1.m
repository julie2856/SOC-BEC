oclear;
%% Fig. 1 of Li et al., PRL 108 225301 (2012)
% slices through the phase diagram showing properties calculated with the
% variational solution

Omegav=linspace(0,2.4,256);  
q=8;    %plot every qth point

k0=1; 
G1=0.2;
G2v=[0.05 0.16];
nv=[1 1 5];

figure('Position',[100,100,700,700])


%%
for i=[1 2] 
G2=G2v(i);

labels={'a','e','b','f','c','g','d','h'};
titles={'$n<n^{(c)}$','$n>n^{(c)}$'};

    for pc=1:length(Omegav)
        [k1(pc),E(pc),beta,other]=SpinOrb_var(k0,G1,G2,Omegav(pc));
            sigmax(pc)=other.sigmax;
            sigmaz(pc)=other.sigmaz;

           k1I(pc)  =other.Ph1.k1;
           k1II(pc) =other.Ph2.k1;
           k1III(pc)=other.Ph3.k1;
           EI(pc)   =other.Ph1.E;
           EII(pc)  =other.Ph2.E;
           EIII(pc) =other.Ph3.E;
           sigmaxI(pc)=other.Ph1.sigmax;
           sigmaxII(pc)=other.Ph2.sigmax;
           sigmaxIII(pc)=other.Ph3.sigmax;
           sigmazI(pc)=other.Ph1.sigmaz;
           sigmazII(pc)=other.Ph2.sigmaz;
           sigmazIII(pc)=other.Ph3.sigmaz;
        
    end

      
        

    subplot(4,2,i)
        plot(Omegav(1:q:end),k1(1:q:end)/k0,'ko','markersize',5,'LineWidth',1);hold on
        plot(Omegav,k1I/k0,'r--','linewidth',1);
        plot(Omegav,k1II/k0,'b--','linewidth',1); 
        plot(Omegav,k1III/k0,'g-','linewidth',1);
            set(gca,'fontsize',12)
            ylabel('$k_1/k_0$','interpreter','latex','fontsize',14)
            axis([Omegav(1),Omegav(end),0,1]);v=axis;
            text(0.05*(v(2)-v(1))+v(1),0.2*(v(4)-v(3))+v(3),['(' labels{i} ')'],'FontSize',12)
            title([titles{i}],'interpreter','latex','fontsize',14)
    
    subplot(4,2,i+2)
        plot(Omegav(1:q:end),E(1:q:end)/k0^2,'ko','markersize',5,'linewidth',1);hold on
        plot(Omegav,EI/k0^2,'r--','linewidth',1);
        plot(Omegav,EII/k0^2,'b--','linewidth',1); 
        plot(Omegav,EIII/k0^2,'g-','linewidth',1);
            set(gca,'fontsize',12)
            ylabel('$E/Nk_0^2$','interpreter','latex','fontsize',14)
            axis([Omegav(1),Omegav(end),-0.5,0.5]);v=axis;
            text(0.05*(v(2)-v(1))+v(1),0.2*(v(4)-v(3))+v(3),['(' labels{i+2} ')'],'FontSize',12)


    subplot(4,2,i+4)
        plot(Omegav(1:q:end),sigmax(1:q:end),'ko','markersize',5,'linewidth',1);hold on
        plot(Omegav,sigmaxI,'r--','linewidth',1);
        plot(Omegav,sigmaxII,'b--','linewidth',1); 
        plot(Omegav,sigmaxIII,'g-','linewidth',1); 
            set(gca,'fontsize',12)
            ylabel('$\langle \sigma_x\rangle$','interpreter','latex','fontsize',14)
            axis([Omegav(1),Omegav(end),-1,0]);v=axis;
            text(0.05*(v(2)-v(1))+v(1),0.2*(v(4)-v(3))+v(3),['(' labels{i+4} ')'],'FontSize',12)

    subplot(4,2,i+6)
        plot(Omegav(1:q:end),abs(sigmaz(1:q:end)),'ko','markersize',5,'linewidth',1);hold on
        plot(Omegav,abs(sigmazI),'r--','linewidth',1);
        plot(Omegav,abs(sigmazII),'b--','linewidth',1); 
        plot(Omegav,abs(sigmazIII),'g-','linewidth',1); 
            set(gca,'fontsize',12)
            xlabel('$\Omega/k_0^2$','interpreter','latex','fontsize',14)
            ylabel('$|\langle \sigma_z\rangle|$','interpreter','latex','fontsize',14)
            axis([Omegav(1),Omegav(end),0,1]);v=axis;
            text(0.05*(v(2)-v(1))+v(1),0.2*(v(4)-v(3))+v(3),['(' labels{i+6} ')'],'FontSize',12)
end
legend('gs','PhI','PhII','PhIII')
