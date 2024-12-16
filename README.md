# kalman-filter
n=500;
t=0:1:n;
v=0.01;
x_theo=t*v;
x0=x_theo+1;
r=rand(1,n+1)-0.5;
x=x0+0.5*r;
xn=kalman_s(x,x_theo,v);
plot(t,x0,'b');
hold on;
plot(t,x,'g');
hold on;
plot(t,xn,'r');
x_theo(1,:)=0.011*ones(1,n+1);
x_theo(2,:)=x_theo(1,:).*t;
x0=x_theo;
x0(1,:)=x0(1,:)+0.01;
x0(2,:)=x0(2,:)+1;
r1=rand(1,n+1)-0.5;
r2=rand(1,n+1)-0.5;
hold off;
x_theo(1,:)=0.011*ones(1,n+1);
x_theo(2,:)=x_theo(1,:).*t;
x=x0;
x(1,:)=x0(1,:)+0.005*r1;
x(2,:)=x0(2,:)+0.07*r2;
xn=kalman(x,x_theo);
subplot(211)
plot(t,x0(1,:),'b');
hold on;
plot(t,x(1,:),'g');
hold on;
plot(t,xn(1,:),'r');
hold on;
title('V');
subplot(212)
plot(t,x0(2,:),'b');
hold on;
plot(t,x(2,:),'g');
hold on;
plot(t,xn(2,:),'r');
hold on;
title('X');

function xn=kalman_s(x,x_theo,v)
    r=x-x_theo;
    xn=zeros(1,length(x));
    xn(1)=x_theo(1);
    sigma2=zeros(1,length(x));
    sigma2(1)=7;
    for i=2:1:length(x)
        xnp=xn(i-1)+v;
        q=std(r(1:i));
        sigma2p=sigma2(i-1)+q;
        kn=sigma2p/(sigma2p+sigma2(1)^2);
        xn(i)=xnp+kn*(x(i)-xnp);
        sigma2(i)=(1-kn)*sigma2p;
    end
end

function xn=kalman_v(x,x_theo)
    r=x-x_theo;
    xn=zeros(2,length(x));
    xn(:,1)=x_theo(:,1);
    p=zeros(2,length(x(1,:)));
    p(:,1)=[7;7];
    F=[1,0;1,1];
    for i=2:1:length(x(1,:))
        xnp=F*xn(:,i-1);
        q(1,1)=std(r(1,1:i));
        q(2,1)=std(r(2,1:i));
        pp=p(:,i-1)+q;
        kn=pp/(pp+p(:,1));
        xn(:,i)=xnp+kn*(x(2,i)-xnp);
        p(:,i)=([1,0;0,1]-kn)*pp;
    end
end
