clear all
close all
%alpha-Zr
deg=pi/180;
%hex
c=5.149;
a=3.231;
%system {10-11}<11-20>
dfpln=[1 0 -1 1;...
       1 -1 0 1;...
       0 1 -1 1;...
       -1 1 0 1;...
       -1 0 1 1;...
       0 -1 1 1];

dfdir=[1 -2 1 0;...
       -1 -1 2 0;...
       -2 1 1 0;...
       -1 -1 2 0;...
       1 -2 1 0;...
       2 -1 -1 0];

sum(dfpln.*dfdir,2)

X=[1; 0; 0];
Y=[0; 1; 0];
Z=[0; 0; 1];

PLN=[2/sqrt(3) 1/sqrt(3) 0 0;...
    0 1 0 0;...
    0 0 0 a/c];

DIR=[sqrt(3) sqrt(3)/2 0 0;...
    0 1.5 0 0;...
    0 0 0 c/a];

dfpln=dfpln';
dfpln=PLN*dfpln;
% dfpln=dfpln1;

dfdir=dfdir';
dfdir=DIR*dfdir;
% dfdir=dfdir1;

dfpln
Y
sum(dfpln.*Y)
pln(:,1)=sum(dfpln.*X)./sqrt(sum(dfpln.^2,1));
pln(:,2)=sum(dfpln.*Y)./sqrt(sum(dfpln.^2,1));
pln(:,3)=sum(dfpln.*Z)./sqrt(sum(dfpln.^2,1)); 
dir(:,1)=sum(dfdir.*X)./sqrt(sum(dfdir.^2,1));
dir(:,2)=sum(dfdir.*Y)./sqrt(sum(dfdir.^2,1));
dir(:,3)=sum(dfdir.*Z)./sqrt(sum(dfdir.^2,1));

dir;
sum(pln(1,:),2)
step=5;
nn=90/step+1;
psi=linspace(0,0.5*pi,nn); phi=linspace(0,0.5*pi,nn);
[psi,phi]=meshgrid(psi,phi);
xx=tan(0.5*psi).*cos(phi);
yy=tan(0.5*psi).*sin(phi);
tdx=sin(psi).*cos(phi);
tdy=sin(psi).*sin(phi); 
tdz=cos(psi);
length(dfpln)
tau=zeros(nn,nn);
for j=1:length(dfpln)
    chi=tdx.*pln(j,1)+tdy.*pln(j,2)+tdz.*pln(j,3);
    chi=chi./sqrt(tdx.^2+tdy.^2+tdz.^2)/sqrt(sum(pln(j,:).^2,2)); 
    lambda=tdx.*dir(j,1)+tdy.*dir(j,2)+tdz.*dir(j,3); 
    lambda=lambda./sqrt(tdx.^2+tdy.^2+tdz.^2)/sqrt(sum(dir(j,:).^2,2)); 
    tau1=abs(chi.*lambda);
    for m = find (tau1>tau)
        tau(m)=tau1(m);
    end
end
tau;
max(max(tau));

contourf(xx,yy,tau)
colormap winter
title('Фактор Шмида')
