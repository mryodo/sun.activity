N_grid=4000;

periods=5;
L=periods*pi;
W=2*pi;

dw=0.1;
T=2*pi/W;


t=linspace(0, L, N_grid);
t_interp=linspace(0, L, N_grid);
h=t(2)-t(1);

p=round(T/h)

k0=0.25*ones(1, N_grid);
for i=3*p:N_grid-p
    k0(i)=0.3;
end

init=asin(2*dw/k0(1));

[t2, theta]=ode45(@ (t, x) myode(t,x,t_interp,k0,dw), t, init);
x0=sin(W*t);
y0=sin(W*t+theta');

%figure;
%hold on;
%grid on;
%plot(t, x0);
%plot(t, y0);

th=round(T/(2*h));

C0=zeros(N_grid, 1);
for i=1:N_grid
   if ((t(i)-T/2>=0) &  (t(i)+T/2<=L))
       D1=x0(i-th:i+th);
       D2=y0(i-th:i+th);
       C0(i)=corr(D1',D2');
   else
       C0(i)=NaN;
   end
end

phi0=acos(C0);

%figure;
%hold on;
%plot(t, theta);
%plot(t, phi0);

k_hat=2*dw*ones(1, N_grid)./sin(phi0);


figure;
hold on;
plot(t, k0);
plot(t, k_hat);




function dydt = myode(t,y,t_interp, k,dw)
    k = interp1(t_interp,k,t);
    dydt = 2*dw-k.*sin(y);
end
