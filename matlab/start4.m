N_grid=1000;
%L=2*pi;
L=10;

t=linspace(0, L, N_grid);
t_interp=linspace(0, L, N_grid);
C_0=zeros(1, length(t));

dw=0.1;

C_0=0.5*ones(1, length(t));
for i=300:1000
    C_0(i)=0.9;
end

phi_0=acos(C_0);
k_0=2*dw*ones(1, length(t))./sin(phi_0);
figure;
hold on;
grid on;
title('Check for slow changing for reconstructing function');
xlabel('t');
ylabel('\partial\phi/\partial t');
dd=diff(phi_0)/h;
d=[dd(1), dd];

plot(t, d, 'b');
plot(t, 2*dw*ones(1, length(t)), 'r--');
plot(t, -2*dw*ones(1, length(t)), 'r--');

%size(t)
%size(diff(phi_0)
k=(2*dw*ones(1, length(t))-d)./sin(phi_0);
% 
%C_0=@(t) 1;
%phi_0=@(t) acos(C_0(t));
%sym t
%p=@(t) diff(phi_0, t)
%k=(2*dw-diff(phi_0(t)))/sin(phi_0(t));
%k
%plot(t, k(t));

[t2, theta]=ode45(@ (t, x) myode(t,x,t_interp,k,dw), t, phi_0(1));

C_0_hat=cos(theta);
ddth=diff(theta)/h;
dth=[dd(1), dd];
k_0_hat=(2*dw*ones(1, length(t))-dth)./sin(theta');

for i=1:length(t)
    if abs(k_0_hat(i))>10
        %k_0_hat(i)=NaN;
        k_0_hat(i)=k_0_hat(i-1);
    end
end

figure;
hold on;
grid on;
plot(t2, k_0, 'r', 'DisplayName', 'k_0');
plot(t2, k_0_hat, 'b', 'DisplayName', 'k_{new}');
title('Parameters before and after reconstruction');
xlabel('t');
ylabel('k(t)');
legend('show','Location','northwest');

figure;
hold on;
grid on;
plot(t2, C_0, 'r', 'DisplayName', 'C_0');
plot(t2, C_0_hat, 'b', 'DisplayName', 'C_{new}');
title('Sliding correlation before and after reconstruction');
xlabel('t');
ylabel('C(t)');
legend('show','Location','northwest');

RC=(t(2)-t(1))/(t(end)-t(1))*sum((C_0-C_0_hat').*(C_0-C_0_hat'))
RC_0=RC/std(C_0)

Rk=(t(2)-t(1))/(t(end)-t(1))*sum((k_0-k_0_hat).*(k_0-k_0_hat))
Rk_0=Rk/std(k_0)

function dydt = myode(t,y,t_interp, k,dw)
    k = interp1(t_interp,k,t);
    dydt = 2*dw-k.*sin(y);
end