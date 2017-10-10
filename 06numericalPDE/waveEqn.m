clear;
%% Initial conditions
f = @(x) sin(pi*x); % u(x,0)   = f(x)
g = @(x) 0;         % u_t(x,0) = g(x)

%% Geometry
dx = 0.10;
dt = 0.05;
c  = 2;

xEnd = 1;
tEnd = 1;

nx = xEnd/dx + 1;
nt = tEnd/dt + 1;

x = linspace(0,xEnd,nx);
t = linspace(0,tEnd,nt);

r = c*dt/dx;
r2 = r*r;

fprintf('r = %5.4f\n',r);

%% Start solving
u = zeros(nx,nt);
u(:,1) = f(x); % Initial Value

i = 2:nx-1; % inner indices to make slicing easier.

u(i,2) = (1-r2)*f(x(i))+(r2/2)*(f(x(i+1))+f(x(i-1)))+dt*g(x(i));
for j = 2:nt-1
    u(i,j+1) = 2*(1-r2)*u(i,j)+r2*(u(i+1,j)+u(i-1,j))-u(i,j-1);
end
[X,T] = meshgrid(x,t);
X = X'; T = T';

%% Analytical solution
uTrue = sin(pi*X).*cos(2*pi*T);

%% Ploting
figure(1)
surf(X,T,u); shading interp;
title('Numerical solution');
view([60,40]); xlabel('x'); ylabel('t')

figure(2)
surf(X,T,uTrue); shading interp;
title('True (analytical) solution')
view([60,40]); xlabel('x'); ylabel('t')

figure(3)
surf(X,T,(u-uTrue)); shading interp;
title('Error (u-uTrue)')
view([60,40]); xlabel('x'); ylabel('t')