%算例1
clear all
O=0;%计数器
G=@(x)x(:,1).^2.*x(:,2)-5.*x(:,1).*x(:,3)+(x(:,2)+1).*(x(:,3).^2)-20;%时变功能函数
%随机输入变量的分布参数
mu=[3.5,3.5];
sigma=[0.3,0.3];
fx=@(x)normpdf(x(:,1),mu(1),sigma(1)).*normpdf(x(:,2),mu(2),sigma(2));
%时间区间
t=[0,5];
n=length(mu);
p=sobolset(n,'Skip',10000);
Nx=8192;%随机变量样本数
N0=5;%外层初始训练样本点个数
P1=p(1:Nx,:);
for i=1:n
    PP(:,i)=norminv(P1(:,i),mu(i),sigma(i));
end

P0=PP(1:N0,:);
Pt=0:0.01:5;
theta=0.01.*ones(1,n+1);
lob=1e-5.*ones(1,n+1);
upb=20.*ones(1,n+1);

for i=1:N0
	O=O+1;
	PPP=ones(length(Pt),1)*P0(i,:);
	PPP(:,n+1)=Pt';
	Pt0=PPP([1,200,400,500],:);
	ypt0=G(Pt0);
	%内层 Kriging 代理模型训练过程
	for j=1:1000
		if j==1
		    x=Pt0;
		    y=ypt0;
		    dmodelt=dacefit(x,y,@regpoly0,@corrgauss,theta,lob,upb);
		    [ugt,sigmagt]=predictor(PPP,dmodelt);
		    prxit=(ones(length(ugt),1).*min(ugt)-ugt).*normcdf((ones(length(ugt),1).*min(ugt)-ugt)./sqrt(sigmagt))+sqrt(sigmagt).*normpdf((ones(length(ugt),1).*min(ugt)-ugt)./sqrt(sigmagt));
		    [PDt(j),It]=max(prxit);
		else
		    x(length(x(:,1))+1,:) = PPP(It,:);
		    y(length(y)+1,1)	  = G(PPP(It,:));
		    
		    dmodelt=dacefit(x,y,@regpoly0, @corrgauss,theta,lob,upb);
		    [ugt,sigmagt]=predictor(PPP,dmodelt);

		    prxit=(ones(length(ugt),1).*min(ugt)-ugt)    .*normcdf((ones(length(ugt),1).*min(ugt)-ugt)./sqrt(sigmagt))+sqrt(sigmagt).*normpdf((ones(length(ugt),1).*min(ugt)-ugt)./sqrt(sigmagt));
		    
		    [PDt(j),It]=max(prxit);
		end
		Crt(j)=max(prxit);
		if Crt(j)<=0.01.*abs(min(ugt)) %内层训练停止准则
		break
		end
	clear ugt sigmagt prxit
	end
	NCALL(O)=length(y);
	Gmin(i)=min(ugt);
	clear PDt It Crt x y ugt sigmagt
end

x=P0;
y=Gmin';
thetax=0.01.*ones(1,n);
lobx=1e-5.*ones(1,n);
upbx=20.*ones(1,n);
%外层Kriging代理模型训练过程
for i=1:1000
    if i==1
        dmodel=dacefit(x,y,@regpoly0,@corrgauss,thetax,lobx,upbx);
        [ug,sigmag]=predictor(PP,dmodel);
        prxi=normcdf((abs(ug))./sqrt(sigmag)).*fx(PP).*sigmag;
        [PD(i),I]=max(prxi);
    else
        x(length(x(:,1))+1,:)=PP(I,:);
        O=O+1;
        PPP=ones(length(Pt),1)*PP(I,:);
        PPP(:,n+1)=Pt';
        Pt0=PPP([1,200,400,500],:);
        ypt0=G(Pt0);
        %嵌套内层Kriging代理模型训练过程
        for j=1:1000
            if j==1
                xx=Pt0;
                yy=ypt0;
                dmodelt=dacefit(xx,yy,@regpoly0,@corrgauss,theta,lob,upb);
                [ugt,sigmagt]=predictor(PPP,dmodelt);
                prxit=(ones(length(ugt),1).*min(ugt)-ugt).*normcdf((ones(length(ugt),1).*min(ugt)-ugt)./sqrt(sigmagt))+sqrt(sigmagt).*normpdf((ones(length(ugt),1).*min(ugt)-ugt)./sqrt(sigmagt));
                [PDt(j),It]=max(prxit);
            else
            xx(length(xx(:,1))+1,:)=PPP(It,:);
            yy(length(yy)+1,1)=G(PPP(It,:));
            dmodelt=dacefit(xx,yy,@regpoly0, @corrgauss,theta,lob,upb);
            [ugt,sigmagt]=predictor(PPP,dmodelt);
            prxit=(ones(length(ugt),1).*min(ugt)-ugt).*normcdf((ones(length(ugt),1).*min(ugt)-ugt)./sqrt(sigmagt))+sqrt(sigmagt).*normpdf((ones(length(ugt),1).*min(ugt)-ugt)./sqrt(sigmagt));
            [PDt(j),It]=max(prxit);
            end
            Crt(j)=max(prxit);
            if Crt(j)<=0.01.*abs(min(ugt))%内层训练停止准则
            break
            end
        clear ugt sigmagt prxit
        end
        NCALL(O)=length(yy);              %真实模型调用统计
        y(length(y)+1,1)=min(ugt);
        clear PDt It Crt xx yy ugt sigmagt
        dmodel=dacefit(x,y,@regpoly0, @corrgauss,thetax,lobx,upbx);
        [ug,sigmag]=predictor(PP,dmodel);
        
        prxi = normcdf((abs(ug))./sqrt(sigmag)).*fx(PP).*sigmag;
        [PD(i),I]=max(prxi);
    end
    
    Cr(i)=mean(normcdf((abs(ug))./sqrt(sigmag)));
    Cr(i)
        if Cr(i)>=0.9999        %外层训练停止准则 -> 
        break
        end
    clear ug sigmag
end

pf=length(find(ug<=0))./length(ug)%时变失效概率
cov_pf=sqrt((1-pf)./(Nx-1)./pf)%变异系数