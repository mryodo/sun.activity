N_grid=4000;

periods=5;
L=periods*pi;
W=2*pi;

dw=0.1;
T=2*pi/W;


t=linspace(0, L, N_grid);
t_interp=linspace(0, L, N_grid);
h=t(2)-t(1);

p=round(T/h);

jump_k=zeros(1, N_grid-3*p);
jump_phi=zeros(1, N_grid-3*p);

jump_Rk=zeros(1, N_grid-3*p);
jump_Rk0=zeros(1, N_grid-3*p);

jump_Rphi=zeros(1, N_grid-3*p);
jump_Rphi0=zeros(1, N_grid-3*p);


for i=1:N_grid-3*p
    k0=0.25*ones(1, N_grid);
    for j=3*p:3*p+i
        k0(j)=0.5;
    end

    init=asin(2*dw/k0(1));
    [t2, theta]=ode45(@ (t, x) myode(t,x,t_interp,k0,dw), t, init);
    x0=sin(W*t);
    y0=sin(W*t+theta');

    th=round(T/(2*h));
    C0=zeros(N_grid, 1);
    for j=1:N_grid
        if ((t(j)-T/2>=0) &  (t(j)+T/2<=L))
            D1=x0(j-th:j+th);
            D2=y0(j-th:j+th);
            C0(j)=corr(D1',D2');
        else
            C0(j)=NaN;
        end
    end

    phi0=acos(C0);
    k_hat=2*dw*ones(1, N_grid)./sin(phi0);
    display(i);

    for j=1:N_grid
        if ~(isnan(C0(j)))
            start=j;
            break;
        end
    end

    for j=N_grid:-1:1
        if ~(isnan(C0(j)))
            finish=j;
            break;
        end
    end

    jump_k(i)=abs((max(k_hat(i+3*p:finish)-k0(i+3*p:finish)))/0.25);
    jump_phi(i)=abs((max(phi0(i+3*p:finish)-theta(i+3*p:finish)))/0.25);

    jump_Rk(i)=sqrt(h/L*sum((k_hat(start:finish)-k0(start:finish)).*(k_hat(start:finish)-k0(start:finish))));
    jump_Rk0(i)=jump_Rk(i)/std(k0(start:finish));

    jump_Rphi(i)=sqrt(h/L*sum((phi0(start:finish)-theta(start:finish)).*(phi0(start:finish)-theta(start:finish))));
    jump_Rphi0(i)=jump_Rphi(i)/std(theta(start:finish));

end

dlmwrite('jump_k_2.txt', jump_k);
dlmwrite('jump_phi_2.txt', jump_phi);
dlmwrite('jump_Rk_2.txt', jump_Rk);
dlmwrite('jump_Rk0_2.txt', jump_Rk0);
dlmwrite('jump_Rphi_2.txt', jump_Rphi);
dlmwrite('jump_Rphi0_2.txt', jump_Rphi0);

function dydt = myode(t,y,t_interp, k,dw)
    k = interp1(t_interp,k,t);
    dydt = 2*dw-k.*sin(y);
end

