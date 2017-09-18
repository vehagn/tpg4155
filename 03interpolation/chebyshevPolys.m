clear;
% Define the first two Chebyshev polynomials
T0 = @(x) 1+0*x;
T1 = @(x) x; 

n = 20; % Amount of polynomails
for i = 2:n
    % Apply recursive definition of Chebyshev polynomials
    str = sprintf('T%i = @(x) 2.*x.*T%i(x) - T%i(x);',i,i-1,i-2);
    eval(str);
end

x = linspace(-1,1,300);
figure(1); clf; hold on; grid on;
for i = 0:n
    str = sprintf('plot(x,T%i(x));',i);
    eval(str);
end
ylim([-1.1,1.1]);