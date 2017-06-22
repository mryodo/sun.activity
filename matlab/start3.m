N_grid=1000;
L=2*pi;

t=linspace(0, L, N_grid);
t_interp=linspace(0, L, N_grid);
C_0=zeros(1, length(t));

dw=0.1;

C_0=0.5*ones(1, length(t));
for i=300:350
    C_0(i)=0.9;
end

phi_0=acos(C_0);
figure;
hold on;
grid on;
dd=diff(phi_0)/h;
d=[dd(1), diff(phi_0)/h];

plot(t, d, 'b');
plot(t, 2*dw*ones(1, length(t)), 'r--');
plot(t, -2*dw*ones(1, length(t)), 'r--');

%size(t)
%size(diff(phi_0)
k=(2*dw*ones(1, length(t))-[diff(phi_0)/h, 0])./sin(phi_0);
% 
%C_0=@(t) 1;
%phi_0=@(t) acos(C_0(t));
%sym t
%p=@(t) diff(phi_0, t)
%k=(2*dw-diff(phi_0(t)))/sin(phi_0(t));
%k
%plot(t, k(t));

[t2, theta]=ode45(@ (t, x) myode(t,x,t_interp,k,dw), t, phi_0(1));
%plot(t2, theta);
%size(theta)

C_0_hat=cos(theta);

figure;
hold on;
grid on;
plot(t2, C_0, 'r');
plot(t2, C_0_hat, 'b');

RC=(t(2)-t(1))/(t(end)-t(1))*sum((C_0-C_0_hat').*(C_0-C_0_hat'))
RC_0=RC/std(C_0)


function dydt = myode(t,y,t_interp, k,dw)
    k = interp1(t_interp,k,t);
    dydt = 2*dw-k.*sin(y);
end