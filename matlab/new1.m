N_grid=1000;
%L=2*pi;
L=10;

t=linspace(0, L, N_grid);
t_interp=linspace(0, L, N_grid);
C_0=zeros(1, length(t));

dw=0.1;

start=200;

total=0.4;

totalk=abs(2*dw*(1/sin(acos(0.9))-1/sin(acos(0.5))));

jumpC=zeros(1, N_grid-start);
jumpk=zeros(1, N_grid-start);
jump_RC=zeros(1, N_grid-start);
jump_RC_0=zeros(1, N_grid-start);
jump_Rk=zeros(1, N_grid-start);
jump_Rk_0=zeros(1, N_grid-start);

for finish=(start+1):N_grid
    C_0=0.5*ones(1, length(t));
    for i=start:finish
        C_0(i)=C_0(i)+total;
    end

    phi_0=acos(C_0);
    k_0=2*dw*ones(1, length(t))./sin(phi_0);
    dd=diff(phi_0)/h;
    d=[dd(1), dd];
   
    k=(2*dw*ones(1, length(t))-d)./sin(phi_0);
    

    [t2, theta]=ode45(@ (t, x) myode(t,x,t_interp,k,dw), t, phi_0(1));

    C_0_hat=cos(theta);
    ddth=diff(theta)/h;
    dth=[dd(1), dd];
    k_0_hat=(2*dw*ones(1, length(t))-dth)./sin(theta');
%     for j=1:length(t)
%         if abs(k_0_hat(j))>10
%             k_0_hat(j)=k_0_hat(j-1);
%         end
%     end

    jumpC(finish-start)=(C_0_hat(finish)-0.5)/total;
    jumpk(finish-start)=(k_0_hat(finish)-2*dw/sin(acos(0.5)))/totalk;
    jump_RC(finish-start)=(t(2)-t(1))/(t(end)-t(1))*sum((C_0-C_0_hat').*(C_0-C_0_hat'));
    jump_RC_0(finish-start)=RC/std(C_0);
    jump_Rk(finish-start)=(t(2)-t(1))/(t(end)-t(1))*sum((k_0-k_0_hat).*(k_0-k_0_hat));
    jump_Rk_0(finish-start)=Rk/std(k_0)';
    finish
end

% figure;
% grid on;
% hold on;
% plot(linspace(0, 0.8, N_grid-start), jumpC);
% title('Relative C_0 reconstuction jump towards intial C_0 jump');
% xlabel('Shock length/Period');
% ylabel('Relative height');
% 
% figure;
% grid on;
% hold on;
% plot(linspace(0, 0.8, N_grid-start), jump_RC);
% title('Reconstruction error on C (RC) towards intial C_0 jump');
% xlabel('Shock length / Period');
% ylabel('RC');
% 
% figure;
% grid on;
% hold on;
% plot(linspace(0, 0.8, N_grid-start), jump_RC_0);
% title('Reconstruction error on C (RC_0) towards intial C_0 jump');
% xlabel('Shock length / Period');
% ylabel('RC_0');

figure;
grid on;
hold on;
plot(linspace(0, 0.8, N_grid-start), jumpk);
title('Relative C_0 reconstuction jump towards intial K_0 jump');
xlabel('Shock length/Period');
ylabel('Relative height');

figure;
grid on;
hold on;
plot(linspace(0, 0.8, N_grid-start), jump_Rk);
title('Reconstruction error on K (Rk) towards intial K_0 jump');
xlabel('Shock length / Period');
ylabel('RC');

figure;
grid on;
hold on;
plot(linspace(0, 0.8, N_grid-start), jump_Rk_0);
title('Reconstruction error on K (Rk_0) towards intial K_0 jump');
xlabel('Shock length / Period');
ylabel('RC_0');


function dydt = myode(t,y,t_interp, k,dw)
    k = interp1(t_interp,k,t);
    dydt = 2*dw-k.*sin(y);
end