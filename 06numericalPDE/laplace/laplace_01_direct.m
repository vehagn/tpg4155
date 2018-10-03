%% Solve the Laplace equation by inversion

%% Initialisation 
clear;

n = 3; % inner nodes on x-axis
m = 3; % inner nodes on y-axis

N = n*m; % Total number of inner nodes

% Create a sparse band-matrix using spdiags
% Sparse matrices use much less memory as we don't have to store
% zeroes.
e = ones(N,1);
d = ones(N,1); d(n:n:(N-m)) = 0;
A = spdiags([e d -2*e],[-3,-1,0],N,N);
A = A + A';

% Convert from sparse to full matrix
Afull = full(A);

% From boundary conditions
b = [-25, -50,-150, 0, 0, -50, 0, 0, -25]';

%% Direct solve
p = A\b;

% In order to plot the result we have to reshape our solution vector
pSol = reshape(p,n,m)';
% We want to also plot the boundary conditions
x = linspace(0,0.5,n+2);
y = linspace(0,0.5,m+2);
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

