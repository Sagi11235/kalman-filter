# kalman-filter
n=500;
t=0:1:n;
v=0.01;
x_theo=t*v;
x0=x_theo+1;
r=rand(1,n+1)-0.5;
x=x0+0.5*r;
xn=kalman(x,x_theo,v);
plot(t,x0,'b');
hold on;
plot(t,x,'g');
hold on;
plot(t,xn,'r');

function xn=kalman(x,x_theo,v)
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
