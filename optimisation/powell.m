function [x,p,U] = powell(F,x0,n,eps,pSearch,xSearch)
%POWELL Summary of this function goes here
%   Detailed explanation goes here
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

U = eye(dim);

for i = 2:n
    p(1,:) = x(i-1,:);
    for k = 1:dim
        G = @(g) F(p(k,:)+g.*U(:,k)');
        g = lineSearch(G,-pSearch,pSearch);
        p(k+1,:) = p(k,:) + g.*U(:,k)';
    end
    for j = 1:length(U)-1
       U(:,j) = U(:,j+1);
    end
    U(:,end) = (p(k+1,:)-p(1,:))';
    if norm(U(:,end)) < eps
        fprintf('U has become linearly dependent! - stopping naive Powell.\n');
        i = i-1;
        break
    end
    U(:,end) = U(:,end)/norm(U(:,end));
    
    G = @(g) F(p(1,:)+g.*U(:,end)');
    g = lineSearch(G,-xSearch,xSearch);
    
    x(i,:) = p(1,:) + g.*U(:,k)';
    
    if (sqrt(sum((x(i,:)-x(i-1,:)).^2)) < eps)
        break
    end
end
x = x(1:i,:);
end

