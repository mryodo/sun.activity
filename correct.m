N_grid=1000;
L=5*pi;

W=pi;
dw=0.1;


t=linspace(0, L, N_grid);
t_interp=linspace(0, L, N_grid);
h=t(2)-t(1);


k0=0.25*ones(1, N_grid);
for i=300:350
    k0(i)=0.7;
end

init=asin(2*dw/k0(1));

[t2, theta]=ode45(@ (t, x) myode(t,x,t_interp,k0,dw), t, init);
x0=sin(W*t);
y0=sin(W*t+theta');

figure;
hold on;
grid on;
plot(t, x0);
plot(t, y0);

figure;
hold on;
plot(t, theta);
plot(t, k0);







function dydt = myode(t,y,t_interp, k,dw)
    k = interp1(t_interp,k,t);
    dydt = 2*dw-k.*sin(y);
end
