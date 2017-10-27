function [x,p,U] = powell2(F,x0,n,eps,pSearch,xSearch)
% Powell method for finding the minimum of a function
if (nargin < 6); xSearch =     0.5; end
if (nargin < 5); pSearch = xSearch; end
if (nargin < 4); eps     =    1e-7; end
if (nargin < 3); n       =     100; end

% dim = 2;
% n = 100;
% pSearch = 1.0;
% xSearch = .1;

dim = length(x0);

x = zeros(n,dim);
x(1,:) = x0;
p = zeros(dim+1,dim);
f = zeros(dim+1,1);

U = eye(dim);

for i = 2:n
    p(1,:) = x(i-1,:);
    f(1) = F(x(i-1,:));
    for k = 1:dim
        G = @(g) F(p(k,:)+g.*U(:,k)');
        g = lineSearch(G,-pSearch,pSearch);
        p(k+1,:) = p(k,:) + g.*U(:,k)';
        f(k+1) = (f(1) - F(p(k+1,:)));
    end
    [r,ri]  = max(f(2:end));
    ri = ri +1;
    U(:,ri-1) = p(ri,:)'/norm(p(ri,:));
    
    p0 = p(1,:);
    pN = p(end,:);
    
    if ((F(2.*pN-p0) >= F(p0))) || (2*(F(p0)-2*F(pN)+F(2*pN-p0))*(F(p0)-F(pN)-r)^2 >= r*(F(p0)-F(2*pN-p0))^2)
        x(i,:) = pN;  
    else
       U(:,ri-1) = (pN-p0)'/norm(pN-p0);
       G = @(g) F(p(1,:)+g.*U(:,ri-1)');
       g = lineSearch(G,-xSearch,xSearch);
       x(i,:) = p(1,:) + g.*U(:,ri-1)';
    end
    
    if (sqrt(sum((x(i,:)-x(i-1,:)).^2)) < eps)
        break
    end
end
x = x(1:i,:);
end

