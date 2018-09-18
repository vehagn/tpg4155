clear;

f = @(x) 1./(1+25.*x.*x);

cheb = @(k,n) cos(((2.*k+1).*pi)/(2.*n));

k = 8;

xk = linspace(-1,1,k); % Equidistant points
xc = cheb(0:k-1,k);    % Chebyshev points

E = constructPoly(xk,f);
C = constructPoly(xc,f);

x = linspace(-1.0,1.0,1000);

figure(1); clf; hold on; grid on;
title('Equidistant vs. Cebyshev');
xlabel('x'); ylabel('y');
plot(x, f(x),'LineWidth',4);
plot(x,E(x),'-.','LineWidth',2);
plot(x,C(x),'--','LineWidth',2);
legend('1/(1+25xˆ2)','Equidistant','Chebyshev','Location','NorthWest');
xlim([min(x) max(x)]);

figure(2); clf; hold on; grid on;
title('Equidistant vs. Cebyshev error');
xlabel('x'); ylabel('y');
plot([-1 1], [0 0],'LineWidth',0.1); % Hack to get the same colours
he=plot(x,(E(x)-f(x)),'-','LineWidth',2);
hc=plot(x,(C(x)-f(x)),'-','LineWidth',2);
legend([he hc],'Equidistant','Chebyshev','Location','NorthWest');
xlim([min(x) max(x)]);

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

function P = constructPoly(xk,f)
% Input interpolation points and function.
    a = diag(getCoeffs(xk,f));
    k = length(xk);
    p = '';
    poly = '@(x) a(1)';
    for i = 2:k
        % Construct a string defining the polynomial
        % A bit of a hack, but it works.
        c = sprintf('a(%i)',i); 
        p = strcat(p,sprintf('.*(x-xk(%i))',i-1));
        poly = strcat(poly,'+',c,p);
    end
    poly = strcat(poly,';');
    P = eval(poly);
end
