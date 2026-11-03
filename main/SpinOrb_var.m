function [k1,EpN,beta,other] = SpinOrb_var(k0,G1,G2,Omega)
%% SpinOrb_var - is a function to solve for the variational state of a spin-orbit 
% coupled BEC following the treatment presented in Li PRL 2012  
%
% usage:
% [k1,EpN,beta,other] = SpinOrb_var(k0,G1,G2,Omega)
%
% here k1, EpN, and beta are properties of the energy-per-particle
% minimising state. THis could be:
% Phase I  (Ph1): stripe/spin mixed;
% Phase II (Ph2): separated phase;
% Phase III(Ph3): single min/zero mtm;
% other is a structured array containing the other parameters of interest
% and the properties of the other phases


Ffn =@(beta)((k0^2-2*G2)+4*(G1+2*G2)*beta);
EpNfn =@(beta)(-Omega^2./(8*Ffn(beta))+G1+G2*(1-4*beta));

Sxf =@(k1)(-sqrt(k0^2-k1^2)/k0);
Szf =@(k1,C1,C2)(k1/k0*(abs(C1)^2-abs(C2)^2));  
% in general minimising states are boundary values, only test beta=0,1/4;
betav=[1/4,0,0]; 
      
%% phase I:
betaI=1/4;C1I=1/sqrt(2);C2I=1/sqrt(2);
k1I=k0*sqrt(1-Omega^2/(4*Ffn(betaI)^2));
EpNI=(-Omega^2./(8*Ffn(betaI))+G1+G2*(1-4*betaI));
validI=(1-Omega^2/(4*Ffn(betaI)^2))>0;
%% phase II:
betaII=0;C1II=1;C2II=0;
k1II=k0*sqrt(1-Omega^2/(4*Ffn(betaII)^2));
EpNII=(-Omega^2./(8*Ffn(betaII))+G1+G2*(1-4*betaII));
validII=(1-Omega^2/(4*Ffn(betaII)^2))>0;
%% phase III:
betaIII=0;C1III=1;C2III=0;
k1III=0;
EpNIII=k0^2/2-Omega/2+G1; % NOTE! 
validIII=true;
 
validv= [validI,validII,validIII]; % only check valid solutions
EpNv=[EpNI,EpNII,EpNIII];EpNv=EpNv(validv);
k1v=[k1I,k1II,k1III];k1v=k1v(validv);
[EpN,indm]=min(EpNv);
beta=betav(indm);
k1=k1v(indm);
 
C1=sqrt((1+sqrt(1-4*beta))/2); 
C2=sqrt((1-sqrt(1-4*beta))/2);
 
other.theta=acos(k1/k0)/2;
other.sigmax=Sxf(k1);
other.sigmaz=Szf(k1,C1,C2);
other.C1=C1;
other.C2=C2;

if(nargout>3) % other outputs
    other.Ph1=struct('C1',NaN,'C2',NaN,'k1',NaN,'E',NaN,'sigmax',NaN,'sigmaz',NaN);
    other.Ph2=struct('C1',NaN,'C2',NaN,'k1',NaN,'E',NaN,'sigmax',NaN,'sigmaz',NaN);
    other.Ph3=struct('C1',NaN,'C2',NaN,'k1',NaN,'E',NaN,'sigmax',NaN,'sigmaz',NaN);
 
 
% Phase I
betaa=1/4;
     if(1-Omega^2/(4*Ffn(betaa)^2)>0)
        C1=1/sqrt(2);C2=1/sqrt(2);
        k1a=real(k0*sqrt(1-Omega^2/(4*Ffn(betaa)^2)));
       % k1a=real(k0*sqrt(1-Omega^2/(4*(k0^2+G1)^2)));
        EI=EpNfn(betaa);
        other.Ph1=struct('C1',C1,'C2',C1,...
            'k1',k1a,'E',EI, ...
            'sigmax',Sxf(k1a),'sigmaz',Szf(k1a,C1,C2)); 
     else
        EI=NaN;
     end

% Phase II 
betaa=0;
    if ((1-Omega^2/(4*Ffn(betaa)^2))>0)
        beta=0;C1=1;C2=0; 
        %k1a=real(k0*sqrt(1-Omega^2/(4*Ffn(betaa)^2)));
        k1a=real(k0*sqrt(1-Omega^2/(4*(k0^2-2*G2)^2)));
        EII=EpNfn(betaa);
        other.Ph2=struct('C1',C1,'C2',C1,...
            'k1',k1a,'E',EII, ...
            'sigmax',Sxf(k1a),'sigmaz',Szf(k1a,C1,C2));   
    else
        EII=NaN;
    end

% Phase III
betaa=0;
    if Omega>2*Ffn(betaa) 
        C1=1;C2=0;     
        k1a=0;  
        EIII=k0^2/2-Omega/2+G1;
        other.Ph3=struct('C1',C1,'C2',C1,...
            'k1',k1a,'E',EIII, ...
            'sigmax',Sxf(k1a),'sigmaz',Szf(k1a,C1,C2)); 
    else
        EIII=NaN;
    end

    if(min([EI,EII,EIII])<EpN-1e-9)
        beta
        Omega
        EpN
        [EI,EII,EIII]
        error(' ')
    end

end
