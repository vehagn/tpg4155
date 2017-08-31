function x = coordSearch(F,x0,n,eps,pSearch)
%POWELL Summary of this function goes here
%   Detailed explanation goes here
if (nargin < 5); pSearch =  0.5; end
if (nargin < 4); eps     = 1e-7; end
if (nargin < 3); n       =  100; end

% dim = 2;
% n = 100;
% pSearch = 1.0;
% xSearch = .1;

dim = length(x0);

x = zeros(n,dim);
x(1,:) = x0;

U = eye(dim);

for i = 2:n
    for k = 1:dim
        G = @(g) F(x(i-1,:)+g.*U(:,k)');
        g = lineSearch(G,-pSearch,pSearch);
        x(i,k) = x(i-1,k) + g;
    end
    if (sqrt(sum((x(i,:)-x(i-1,:)).^2)) < eps)
        break
    end
end
x = x(1:i,:);
end

