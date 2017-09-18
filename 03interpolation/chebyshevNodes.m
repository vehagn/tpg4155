clear;

f = @(x) exp(x);

cheb = @(k,n) cos(((2.*k+1).*pi)/(2.*n));

k = 4;

xk = linspace(-1,1,k); % Equidistant points
xc = cheb(0:k-1,k);    % Chebyshev points

E = getCoeffs(xk,f);
a = diag(E);
Equidist = @(x)  ...
      a(1) ...
    + a(2).*(x-xk(1)) ...
    + a(3).*(x-xk(1)).*(x-xk(2)) ...
    + a(4).*(x-xk(1)).*(x-xk(2)).*(x-xk(3));

C = getCoeffs(xc,f);
b = diag(C);
Chebyshev  = @(x)  ...
      b(1) ...
    + b(2).*(x-xc(1)) ...
    + b(3).*(x-xc(1)).*(x-xc(2)) ...
    + b(4).*(x-xc(1)).*(x-xc(2)).*(x-xc(3));

x = linspace(-1,1,1000);

figure(1); clf; hold on; grid on;
title('Equidistant vs. Cebyshev');
xlabel('x'); ylabel('y');
plot(x, f(x),'LineWidth',4);
plot(x,Equidist(x),'-.','LineWidth',2);
plot(x,Chebyshev(x),'--','LineWidth',2);
legend('exp(x)','Equidistant','Chebyshev','Location','NorthWest');

figure(2); clf; hold on; grid on;
title('Equidistant vs. Cebyshev error');
xlabel('x'); ylabel('y');
plot([-1 1], [0 0],'LineWidth',0.1); % Hack to get the same colours
he=plot(x,(Equidist(x)-f(x)),'-','LineWidth',2);
hc=plot(x,(Chebyshev(x)-f(x)),'-','LineWidth',2);
legend([he hc],'Equidistant','Chebyshev','Location','NorthWest');

[mEerr,iE] = max(abs(Equidist(x)-f(x)));
[mCerr,iC] = max(abs(Chebyshev(x)-f(x)));

fprintf('Max equidistant error: %5.4f at x=%5.4f\n',mEerr,x(iE));
fprintf('Max Chebyshev error:   %5.4f at x=%5.4f\n',mCerr,x(iC));

function fi = getCoeffs(xk,f)
% Calculate the Newton Polynomial divided differen cecoefficients. 
% We want the diagonal of the output here.
    n = length(xk);
    fi = zeros(n,n);
    fi(:,1) = f(xk);
    for j = 2:n   
        fi(j:n,j) = (fi(j:n,j-1)-fi(j-1:n-1,j-1))./(xk(j:n)'-xk(1:n-j+1)');
    end
end

