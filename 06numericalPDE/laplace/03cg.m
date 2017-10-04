%% Solve the Laplace equaiton using the Conjugate Gradient method.

%% Initialisation
clear;

n = 3;           % Inner points on x-axis
m = 3;           % Inner points on y-axis
eps = 1e-16;     % Error threshold
nIter = 1000;    % Maximum number of iterations
N = n*m;         % Total number of inner nodes

% Create a sparse band-matrix using spdiags
% Sparse matrices use much less memory as we don't have to store
% zeroes.
e = ones(N,1);
d = ones(N,1); d(n:n:(N-m)) = 0;
A = spdiags([e d -2*e],[-3,-1,0],N,N);
A = A + A';

b = [-25, -50,-150, 0, 0, -50, 0, 0, -25]';

%% Initial guess
p = zeros(N,1);
% Preconditioned
% p = 25*ones(N,1);

%% Begin iteration
it = 0;
r = b - A*p;
d = r;
while it < nIter
    rTr   = r'*r;
    alpha = rTr/(d'*A*d);
    p     = p + alpha*d;
    r     = r - alpha*A*d;
    beta  = r'*r/rTr; 
    d     = r + beta*d;
    it    = it + 1;
    if norm(r,inf) < eps
        break
    end
end
fprintf('Numer of iterations: %4i\n',it);

% In order to plot the result we have to reshape our solution vector
pSol = reshape(p,n,m)';
% We want to also plot the boundary conditions
x = linspace(0,0.5,5);
y = linspace(0,0.5,5);
U = zeros(n+2,m+2);
U(2:4,2:4) = pSol;        % Inner part
U(1,:)     = 200*x;       % Top part
U(:,end)   = 200*flip(y); % Right part
% surf(x,y,flipud(U)); shading('interp');
imagesc(x,y,flipud(U));
set(gca,'YDir','normal')
xlabel('x');
ylabel('y');
colorbar();