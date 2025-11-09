for i=1:length(Tt)-1
A=@(x)-2.*x(:,1).*x(:,3).*sin(Tt(i));
B=@(x)2.*x(:,3).*(x(:,4)-x(:,1).*cos(Tt(i)));
C=@(x)(x(:,2)).^2-(x(:,1)).^2-(x(:,3)).^2-(x(:,4)).^2+2.*x(:,1).*x(:,4).*cos(Tt(i));
G=@(x)alpha.*pi./180-abs(76.*pi./180+60.*pi./180.*sin((3./4).*(Tt(i)-95.5.*pi./180))-2.*atan((A(x)+sqrt(A(x).^2+B(x).^2-C(x).^2))./(B(x)+C(x))));
[P1(i,:),  beta1(i)]=AFOSM_solu(mu,diag(sigma.^2), G);
clear G A B C
A=@(x)-2.*x(:,1).*x(:,3).*sin(Tt(i+1));
B=@(x)2.*x(:,3).*(x(:,4)-x(:,1).*cos(Tt(i+1)));
C=@(x)(x(:,2)).^2-(x(:,1)).^2-(x(:,3)).^2-(x(:,4)).^2+2.*x(:,1).*x(:,4).*cos(Tt(i+1));
G=@(x)alpha.*pi./180-abs(76.*pi./180+60.*pi./180.*sin((3./4).*(Tt(i+1)-95.5.*pi./180))-2.*atan((A(x)+sqrt(A(x).^2+B(x).^2-C(x).^2))./(B(x)+C(x))));
[P2(i,:), beta2(i)]=AFOSM_solu(mu,diag(sigma.^2),G);

alphat(i)=-sum((((P1(i,:)-mu)./sigma)./beta1(i)).*(((P2(i,:)-mu)./sigma))./beta2(i));
v(i)=mvncdf([beta1(i),-beta2(i)],[0,0],[1,alphat(i);alphat(i),1])./h; %跨越率
end
pf=1-(1-normcdf(-beta1(1))).*exp(-sum(v).*h) %时变失效概率