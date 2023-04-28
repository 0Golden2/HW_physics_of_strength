clear all
close all
deg=pi/ 180;

epsilon = [0 0.2 0.6 1.2 1.8 2.4 3.0]'
sigma = [0 2 4.5 6 6.5 7 7;...
       0 6 9.5 10.5 11 11 11;...
       0 15 24 27 27.5 28 28]';
ksi = [3e-5 3e-4 3e-3];


color = ["black";"red";"blue"];
xx = 0: 0.01:3;

for i = 1:length(ksi)
    ss(:,i)=spline(epsilon, sigma(:,i),xx);
    plot(epsilon, sigma(:,i), "o", xx,ss(:,i), "color", color(i))
    hold on
end

objfun=@(x)sum ((x(1) * xx(:).^x(2)*ksi(1)^x(3) - ss(:,1)).^2) +    ...
sum((x(1) * xx(:).^x(2)*ksi(2)^x(3)-ss(:,2)).^2)+...
sum((x(1) * xx(:).^x(2)*ksi(3)^x(3)-ss(:,3)).^2);

x0 = [500,0.05,1];
objfun(x0)
[x, fval,exitflag, output] = fminunc(objfun,x0);

x
figure
for i=1:length (ksi)
    for j=1:length(xx)
        ss(j,i) = YS(xx(j), ksi(i), x(1), x(2), x(3));
    end
    plot (epsilon, sigma(:,i),'o',xx,ss(:,i),'color',color(i))
    hold on
end

function f = YS(x,ksi,K,n,m)
    f =K*x^n*ksi^m;
end