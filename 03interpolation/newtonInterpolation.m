clear;

f = @(x) cos(x);
xk = [0, 1, 2, 3, 4]; % Where to evaluate our function f

x = linspace(-1,5,1000);

C = getCoeffs(xk,f);
a = diag(C); % We only need the diagonal coefficients.

N1 = @(x)  a(1) + a(2).*(x-xk(1));
N2 = @(x) N1(x) + a(3).*(x-xk(1)).*(x-xk(2));
N3 = @(x) N2(x) + a(4).*(x-xk(1)).*(x-xk(2)).*(x-xk(3));
N4 = @(x) N3(x) + a(5).*(x-xk(1)).*(x-xk(2)).*(x-xk(3)).*(x-xk(4));

figure(1); clf; hold on; grid on;
title('Newton interpolation of cos(x)');
xlabel('x'); ylabel('y');
plot(x, f(x),'LineWidth',4);
plot(x,N1(x),'--','LineWidth',2);
plot(x,N2(x),'--','LineWidth',2);
plot(x,N3(x),'--','LineWidth',2);
plot(x,N4(x),'--','LineWidth',2);
ylim([-1.1,1.1])
legend('f','N1','N2','N3','N4','Location','SouthWest');

figure(2); clf; hold on; grid on;
title('Newton interpolation error');
xlabel('x'); ylabel('y');
plot(x, f(x),'--','LineWidth',1);
plot(x,N1(x)-f(x),'-','LineWidth',2);
plot(x,N2(x)-f(x),'-','LineWidth',2);
plot(x,N3(x)-f(x),'-','LineWidth',2);
plot(x,N4(x)-f(x),'-','LineWidth',2);
ylim([-1.1,1.1])
legend('f','N1-f','N2-f','N3-f','N4-f','Location','SouthWest');

fprintf('x\t\t|N1(x)-f(x)|\t|N2(x)-f(x)|\t|N3(x)-f(x)|\t|N4(x)-f(x)|\n');
for m = 1:length(xk)
    fprintf('%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\n',...
        xk(m),...
        abs(N1(xk(m))-f(xk(m))),...
        abs(N2(xk(m))-f(xk(m))),...
        abs(N3(xk(m))-f(xk(m))),...
        abs(N4(xk(m))-f(xk(m))));
end

function fi = getCoeffs(xk,f)
% Calculate the Newton Polynomial divided differences cecoefficients. 
% We want the diagonal of the output here.
    n = length(xk);
    fi = zeros(n,n);
    fi(:,1) = f(xk);
    for j = 2:n   
        fi(j:n,j) = (fi(j:n,j-1)-fi(j-1:n-1,j-1))./(xk(j:n)'-xk(1:n-j+1)');
    end
end

