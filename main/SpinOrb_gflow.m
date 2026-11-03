function [x,psi,E,mu,output]=SpinOrb_gflow(G1,G2,k0,Omega,k1,M,varargin)
% SpinOrb_gflow - a function to find GP solutions to a spin-orbit
% coupled 1D uniform BEC.
%
% usage: (1 - use variational initial state)
% [x,psi,E,mu,output]=SpinOrb_gflow3(g,gab,k0,Omega,k1,M)
%
% usage: (2 - supply initial state)
% [x,psi,E,mu,output]=SpinOrb_gflow3(g,gab,k0,Omega,k1,M,'psiin',psi)
p = inputParser;
p.CaseSensitive = true;
p.KeepUnmatched = true;
%% system parameters
addRequired(p,'G1',@isnumeric);   
addRequired(p,'G2',@isnumeric); 
addRequired(p,'k0',@isnumeric); 
addRequired(p,'Omega',@isnumeric);
addRequired(p,'k1',@isnumeric);
addRequired(p,'M',@isnumeric);


%% inputs,grid/potential control
addParameter(p,'psiin',nan,@isnumeric);  
 

%% algorithm control 
addParameter(p,'dt',1e-2,@isnumeric);
addParameter(p,'maxsteps',20e3,@isnumeric); 
addParameter(p,'resid',1e-7,@isnumeric);  


%% feedback/output control 
addParameter(p,'disp',500,@isnumeric);
addParameter(p,'graphics',true,@islogical);
 

parse(p,G1,G2,k0,Omega,k1,M,varargin{:}) 
solnp=p.Results;
dt=solnp.dt;  
tol=solnp.resid;
displ=solnp.disp;      
maxsteps=solnp.maxsteps;
% conversions
% G1=(g+gab)/4; G2=(g-gab)/4;
g  =2*(G1+G2);
gab=2*(G1-G2);

%% care with grids
% the stripe phase fields can be expanded as
% psi ~\sum_n c_n*exp(1i*(1+2*n)k1*x),
% this exp(-1i*k1*x)\psi is periodic with period L=pi/k1, 
% i.e. the relevant reciprocal lattice vector is 2*k1!
% grid1D constructs appropriate grids.

L=pi/k1;
[x,k,dx,~,U]=grid1D(L,M,k1);
D1=U'*diag(1i*k)*U;
quad = @(f)(sum(f(:)*dx));  
renorm = @(psi)(psi/sqrt(quad(abs(psi).^2))*sqrt(L));% renorm so that the average value of |psi|^2=1
if(isnan(solnp.psiin)) % no supplied solution : us variational
    % constructing ansatz Eq.(3) Li2012 
    [k1,~,~,other] = SpinOrb_var(k0,G1,G2,Omega);

    c=cos(other.theta);s=sin(other.theta);    
    psia= (other.C1*c*exp(1i*k1*x)+other.C2*s*exp(-1i*k1*x));
    psib= (-other.C1*s*exp(1i*k1*x)-other.C2*c*exp(-1i*k1*x)); 
    psi=[psia(:);psib(:)];psii=psi;
else
    psiin=solnp.psiin;
    Min=length(psiin)/2;
    psia=interp1(linspace(x(1),x(end),Min),psiin(1:Min),x,'linear','extrap');
    psib=interp1(linspace(x(1),x(end),Min),psiin(1+Min:end),x,'linear','extrap');
    psi=[psia(:);psib(:)];psii=psi;
end



psi=psi(:);psi=renorm(psi); 
i1=1:M;i2=M+1:2*M;

% pauli matrices 
sigmaz=[eye(M),zeros(M);zeros(M),-eye(M)];
sigmax=[zeros(M),eye(M);eye(M),zeros(M)];  
%momentum operator, p_x in block diagonal form
px=blkdiag(-1i*D1,-1i*D1);

%single particle hamiltonian h_0, Eq. (2)
h0=0.5*((px-k0*sigmaz)^2 + Omega*sigmax);

%non-linear part of overall operator (L)
ga=g*abs(psi(i1)).^2 + gab*abs(psi(i2)).^2; 
gb=gab*abs(psi(i1)).^2 + g*abs(psi(i2)).^2; 
gblock = blkdiag(diag(ga),diag(gb));

%overall operator similar to L in Eq.(68) in the comp. physics notes
LGP=h0+gblock;

%overall Hamiltonian - from Eq.(1) Li2012
H=h0+0.5*gblock;

%% energy=psi'*H*psi - Eq.(1) Li2012 
E =real(psi(:)'*H*psi(:)*dx);Ei=E;
 

%% chemical potential - Eq.(61) comp. notes  
mu  = real(psi(:)'*LGP*psi(:)*dx)/L;
 

%% residual 
resid=max(abs(LGP*psi(:)-mu*psi(:)));
if(isnan(resid))

    disp('')
end
keepiter=true;pc=0;
while keepiter
    pc=pc+1;
     
    A=eye(2*M)+dt*LGP;    %backwards euler
    psi=A\psi;    
    psi=renorm(psi); 

    ga=g*abs(psi(i1)).^2 + gab*abs(psi(i2)).^2; 
    gb=gab*abs(psi(i1)).^2 + g*abs(psi(i2)).^2; 
    gblock = blkdiag(diag(ga),diag(gb));
    LGP=h0+gblock;
    H=h0+0.5*gblock;
    E =real(psi(:)'*H*psi(:)*dx); 
    mu=real(psi(:)'*LGP*psi(:)*dx)/L;
  
    resid=max(abs(LGP*psi(:)-mu*psi(:)));
     

     if (pc/displ==round(pc/displ))
        disp(['step#',num2str(pc),', E=',num2str(E,7),', resid=',num2str(resid,3)])
        %[pc,K,resid]
        figure(88)
        plot(x,real(psi(i1)),'b-');hold on
        plot(x,imag(psi(i1)),'r-'); 
        plot(x,real(psi(i2)),'b--'); 
        plot(x,imag(psi(i2)),'r--'); 

        plot(x,real(psii(i1)),'b:'); 
        plot(x,imag(psii(i1)),'r:'); 
        plot(x,real(psii(i2)),'b.'); 
        plot(x,imag(psii(i2)),'r.');hold off
        disp('')
     end




     %% test
    if(pc>maxsteps)
         keepiter=false;
    elseif(resid<tol)
         keepiter=false;
    end



 
end

output.resid=resid;
output.sigmax=real(psi(:)'*sigmax*psi(:)*dx)/L;
output.sigmaz=real(psi(:)'*sigmaz*psi(:)*dx)/L;
 

 