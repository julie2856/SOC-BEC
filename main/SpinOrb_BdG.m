function [ev,u1v,u2v,v1v,v2v] = SpinOrb_BdG(qxv,Nbands,Omega,G1,G2,psiv,mu,k1)

 
 M=length(psiv)/2;
[x,k,dx,dk,U]=grid1D(pi/k1,M,k1);
 
 %[x,k,dx,dk,D1,D2,U]=grid(pi/k1,M); %??

  

%% operators
quad=@(f)(sum(f(:))*dx) ;
%T{1} =  (U'*diag(0.5*(k-1).^2)*U); 
%T{2} =  (U'*diag(0.5*(k+1).^2)*U);    
 
% Martone 2021 notation for interactions (Gdd,Gss)
Gdd=2*G1;
Gss=2*G2;
g=(Gdd+Gss);
g12=(Gdd-Gss);
psi1=psiv(1:M);cpsi1=conj(psi1);
psi2=psiv(1+M:end);cpsi2=conj(psi2); 
n1=abs(psi1).^2;n2=abs(psi2).^2;
n=n1+n2;fz=n1-n2;
hD = Gdd*[diag(n+n1),diag(psi1.*cpsi2);diag(psi2.*cpsi1),diag(n+n2)] + ...
     Gss*[diag(fz+n1),-diag(psi1.*cpsi2);-diag(cpsi1.*psi2),diag(-fz+n2)];

hC = Gdd*[diag(psi1.^2),diag(psi1.*psi2);diag(psi2.*psi1),diag(psi2.^2)] + ...
     Gss*[diag(psi1.^2),-diag(psi1.*psi2);-diag(psi2.*psi1),diag(psi2.^2)];
%Tp2= U'*diag(k.^2)*U + 2i*U'*diag(1i*k)*U+eye(M);
Tp2=U'*diag((k-1).^2)*U;
%Tm2= U'*diag(k.^2)*U - 2i*U'*diag(1i*k)*U+eye(M);  
Tm2=U'*diag((k+1).^2)*U;
hSO=0.5*[Tp2, Omega*eye(M); Omega*eye(M),Tm2];
chSO=0.5*[Tm2 Omega*eye(M); Omega*eye(M),Tp2]; 
BdGH = [hSO-mu*eye(2*M)+hD,hC;-conj(hC),-chSO-conj(-mu*eye(2*M)+hD)];

tmp=BdGH*[psi1;psi2;-conj(psi1);-conj(psi2)];
[max(abs(tmp(1:M))),max(abs(tmp(M+1:2*M))),max(abs(tmp(2*M+1:3*M))),max(abs(tmp(3*M+1:4*M)))]

 
for pc=1:length(qxv)
    qx=qxv(pc);    
    %Tp2= U'*diag(k.^2)*U + 2i*U'*diag(1i*k)*U+eye(N);
    %Tm2= U'*diag(k.^2)*U - 2i*U'*diag(1i*k)*U+eye(N);
    hSO=0.5*[U'*diag((k+qx-1).^2)*U, Omega*eye(M); Omega*eye(M),U'*diag((k+qx+1).^2)*U];
    chSO=0.5*[U'*diag((k+qx+1).^2)*U, Omega*eye(M); Omega*eye(M),U'*diag((k+qx-1).^2)*U];
    BdGH = [hSO-mu*eye(2*M)+hD,hC;-conj(hC),-chSO-conj(-mu*eye(2*M)+hD)];
    [W,D]=eig(BdGH);
    evf=diag(D);
    diag(D);
    evr=real(evf);
    [evr,ind]=sort(evr);evf=evf(ind);
    W=W(:,ind);
    indpos=find(evr>-1e-9);
  %  indpos=1:4*N;
    evf=evf(indpos);evr=evr(indpos);W=W(:,indpos);
    %evf=evf(N+1:end);evr=evr(N+1:end);W=W(:,N+1:end);
    
    for pcj=1:Nbands 
            if(abs(imag(evf))>1e-4)
                disp(['warning eigenvalue q=',num2str(qx),', band ',num2str(pcj),' has significant imaginary part'])
            end
           u1t=W(1:M,pcj);u2t=W(M+1:2*M,pcj);
           v1t=W(2*M+1:3*M,pcj);v2t=W(3*M+1:4*M,pcj);
           normuv(pc,pcj) = quad(abs(u1t).^2+abs(u2t).^2-abs(v1t).^2-abs(v2t).^2);
       
            u1v{pc,pcj}=u1t/sqrt(normuv(pc,pcj));
            v1v{pc,pcj}=v1t/sqrt(normuv(pc,pcj));
            u2v{pc,pcj}=u2t/sqrt(normuv(pc,pcj));
            v2v{pc,pcj}=v2t/sqrt(normuv(pc,pcj));

            ev(pc,pcj)=evr(pcj);
             
    end
    
end
fprintf('\n')  
 


 