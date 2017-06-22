t=linspace(0, 2*pi, 1000);
h=t(2)-t(1);

x=sin(t);
d=diff(x)/h;
d=[d(1), d];

figure;
hold on;
plot(t, d, 'b');
plot(t, cos(t), 'r');