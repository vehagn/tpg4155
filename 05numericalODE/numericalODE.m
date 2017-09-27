%% Problem setup
clear;
% We are looking to solve the ODE y' = f(y,t) = y - t^2 + 1 using Euler
% and compairing it with the analytical solution:
y = @(t) (t+1).*(t+1) - 0.5*exp(t); % Analytical solution
% y = @(t) (exp(3.*t).*(5.*t-1)+exp(-2.*t))./25;
t = linspace(0,2);                  % Dense t-grid to plot y(x) smoother
a = 0.5;                            % Starting value


f = @(t,y) y - t.^2 +1;
% f = @(t,y) t.*exp(3.*t)-2.*y;
n = 10;
ti = linspace(0,2,n); % Coarse t-grid for numerical method
h = ti(2) - ti(1);   

%% Initialise plots
figure(1); clf;
hold on; grid on;
plot(t ,y(t),'k','DisplayName','Analytical Exact',...
     'LineWidth',2);
legend('-DynamicLegend','Location','NorthWest');
title('Function'); xlabel('t'); ylabel('y');

figure(2); clf;
hold on; grid on;
legend('-DynamicLegend','Location','NorthWest');
title('Error'); xlabel('t'); ylabel('y');

%% Euler
ye = zeros(0,n); % Euler Approximation
ye(1) = a;

for i = 1:n-1
    ye(i+1) = ye(i) + h*f(ti(i),ye(i));
end

figure(1);
plot(ti,ye,'--o','DisplayName','Numerical Euler',...
     'LineWidth',2);

eErr = y(ti)-ye;
figure(2);
plot(ti,eErr,'--o','DisplayName','Euler','LineWidth',2);

%% Heun
yh = zeros(0,n); % Heun Approximation
yh(1) = a;

for i = 1:n-1
    k1 = f(ti(i),yh(i));
    k2 = f(ti(i+1),yh(i)+h*k1);
    yh(i+1) = yh(i) + (h/2)*(k1+k2);
end

figure(1);
plot(ti,yh,'--o','DisplayName','Heun','LineWidth',2);

hErr = y(ti)-yh;

figure(2);
plot(ti,hErr,'--o','DisplayName','Heun','LineWidth',2);

%% Midpoint
ym = zeros(0,n); % Midpint Approximation
ym(1) = a;

for i = 1:n-1
    ym(i+1) = ym(i) + h*f(ti(i)+h/2,ym(i)+(h/2)*f(ti(i),ym(i)));
end

figure(1);
plot(ti,ym,'--o','DisplayName','Midpoint','LineWidth',2);

mErr = y(ti)-ym;

figure(2);
plot(ti,mErr,'--o','DisplayName','Midpoint','LineWidth',2);

%% RK4
yr = zeros(0,n);
yr(1) = a;

for i = 1:n-1
    k1 = f(ti(i)    ,yr(i));
    k2 = f(ti(i)+h/2,yr(i)+(h/2)*k1);
    k3 = f(ti(i)+h/2,yr(i)+(h/2)*k2);
    k4 = f(ti(i)+h  ,yr(i)+    h*k3);
    yr(i+1) = yr(i) + (h/6)*(k1+2*k2+2*k3+k4);
end

figure(1);
plot(ti,yr,'--o','DisplayName','RK4','LineWidth',2);

rErr = y(ti)-yr;

figure(2);
plot(ti,rErr,'--o','DisplayName','RK4','LineWidth',2);

%% ODE45
% We use the built-in ode45() function in MATLAB with default values.
[to,yo] = ode45(f,[ti(1) ti(end)],a);

to = to';
yo = yo';
    
figure(1);
plot(to,yo,'--o','DisplayName','ode45','LineWidth',2);

oErr = y(to)-yo;

figure(2);
plot(to,oErr,'--o','DisplayName','ode45','LineWidth',2);

%% Error
fprintf('\nEndpoint errors\n');
fprintf('Euler:    %4.2e\n',abs(eErr(end)));
fprintf('Heun:     %4.2e\n',abs(hErr(end)));
fprintf('Midpoint: %4.2e\n',abs(mErr(end)));
fprintf('RK4:      %4.2e\n',abs(rErr(end)));
fprintf('ode45:    %4.2e\n',abs(oErr(end)));