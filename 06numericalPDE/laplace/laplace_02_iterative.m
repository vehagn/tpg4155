%% Solve the Laplace equation using the Jacobi or Gauss-Seidel method

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

checkConvergence(A)

b = [-25, -50,-150, 0, 0, -50, 0, 0, -25]';

%% Initial guess
p = zeros(N,1);
% Preconditioned
% p = 25*ones(N,1);

%% Begin iteration
it = 0;
while it < nIter
%    pnew = jacobi(A,b,p);
     pnew = gaussSeidel(A,b,p);    
    if norm((pnew-p),inf) < eps
        break
    end
    p  = pnew;
    it = it + 1;
end
p  = pnew;
fprintf('Numer of iterations: %4i\n',it);

% In order to plot the result we have to reshape our solution vector
pSol = reshape(p,n,m)';
% We want to also plot the boundary conditions
x = linspace(0,0.5,5);
y = linspace(0,0.5,5);
U = zeros(n+2,m+2);
U(2:4,2:4) = pSol;        % inner part
U(1,:)     = 200*x;       % top part
U(:,end)   = 200*flip(y); % right part
% surf(x,y,flipud(U)); shading('interp');
imagesc(x,y,flipud(U));
set(gca,'YDir','normal')
xlabel('x');
ylabel('y');
colorbar();

function checkConvergence(A)
    % We can check the convergence of the algorithm by testing if A is
    % stictly diagnonally dominant.
    D = diag(A); % Creates a column vector with the diagonal of A.
    if any(abs(D) <= abs(sum(A,2)-D))
        % We sum the whole row of A, so we have to subtract the 
        % diagonal element.
        fprintf('A is not strictly diagonally dominant!\n');
        fprintf('Convergence is not guaranteed.\n');
    else
        fprintf('A is stricly diagonally dominant.\n');
        fprintf('Set sail for convergence!\n');
    end
end

function x = gaussSeidel(A,b,x)
    % We don't need to "split" into (k) and (k-1) since the new
    % values are stored over the old values. This makes for cleaner code.
    for i = 1:length(x)
        % We include j==i in the sum (which we shouldn't), 
        % so we have to subtract it.
        x(i) = x(i) - (A(i,:)*x - b(i))/A(i,i) ;
    end
end

function x = jacobi(A,b,x)
    % Jacobi iteration on component form. Written like this we use all
    % the old values stored in x before overwriting them.
    x = x - (A*x - b)./diag(A) ;
end