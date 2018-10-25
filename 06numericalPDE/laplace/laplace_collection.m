%% Solve the Laplace equation by inversion

%% Initialisation 
clear;

%% Geometry
n = 20; % nodes on x-axis
m = 20; % nodes on y-axis

ax = 0.0; % Start on x-axis
ay = 0.0; % Start on y-axis
bx = 0.5; % End on x-axis
by = 0.5; % End on y-axis;

%% Boundary conditions (left and bottom are zero)
gTop   = @(x,y) 200*x;
gRight = @(x,y) 200*y;

%% Create nodes
x = linspace(ax,bx,n);
y = linspace(ay,by,m);
N = (n-2)*(m-2); % Total number of inner nodes

%% Prepare matrix to hold solution
U = zeros(m,n);
U(:,end) = gRight(bx,flip(y)); % Top part
U(1,:)   = gTop(x,by);         % Right part

Udirect = U;
Uiter   = U;

%% Create A matrix
% Create a sparse band-matrix using spdiags
% Sparse matrices use much less memory as we don't have to store zeroes.
e = ones(N,1);      % For full bands
d = ones(N,1);      % For end nodes
d(n-2:n-2:N) = 0;   % Insert zeros for end-nodes
% We create just the lower portion then add the transpose for the upper
A = spdiags([e d -2*e],[-(n-2),-1,0],N,N);
A = A + A';

%% Convert from sparse to full matrix for illustration
fullA = full(A);

%% Create b vector
b = zeros(N,1);
for i = 1:n-2
    b(i) = b(i) - gTop(x(i+1),by);
end
for j = 1:m-2
    b(j*(n-2)) = b(j*(n-2)) - gRight(bx,y(end-j));
end

%% Direct solve
pDirect = A\b;
% In order to plot the result we have to reshape our solution vector
Udirect(2:m-1,2:n-1) = reshape(pDirect,n-2,m-2)'; % Insert solution

%% Iteration
pIter = zeros(N,1); % Initial guess
% pIter = 25*ones(N,1); % Precondition

eps = 1e-14;     % Error threshold (Play around with this)
w   = 1.0;       % Successive overrelaxation omega
nIter = 5000;    % Maximum number of iterations
checkConvergence(A)

%% Uncomment one of the methods (Jacobi, Gauss-Seidel, CG, or SOR)
% [pIter,it] = jacobi(A,b,pIter,nIter,eps);
% [pIter,it] = gaussSeidel(A,b,pIter,nIter,eps);
[pIter,it] = cg(A,b,pIter,nIter,eps);

Uiter(2:m-1,2:n-1) = reshape(pIter,n-2,m-2)';
% [Uiter,it] = SOR(Uiter,w,nIter,eps);

fprintf('Numer of iterations: %4i\n',it);

%% True solution
[X,Y] = meshgrid(x,y);
uTrue = 400.*X.*Y;


%% Plot
figure(1);
subplot(3,2,1); % Numerical direct solution
imagesc(x,y,flipud(Udirect));
axis('image'); colorbar();
set(gca,'YDir','normal')
title('Numerical direct');
xlabel('x'); ylabel('y');

subplot(3,2,2); % Numerical iterative solution
imagesc(x,y,flipud(Udirect));
axis('image'); colorbar();
set(gca,'YDir','normal')
title('Numerical iterative');
xlabel('x'); ylabel('y');

subplot(3,2,[3 4]); % Analytical solution
imagesc(x,y,uTrue);
axis('image'); colorbar();
set(gca,'YDir','normal')
title('Analytical');
xlabel('x'); ylabel('y');

subplot(3,2,5); % Difference direct
imagesc(x,y,uTrue-flipud(Udirect));
axis('image'); colorbar();
title('Difference direct');
set(gca,'YDir','normal')
xlabel('x'); ylabel('y');

subplot(3,2,6); % Difference direct
imagesc(x,y,uTrue-flipud(Uiter));
axis('image'); colorbar();
title('Difference iterative');
set(gca,'YDir','normal')
xlabel('x'); ylabel('y');

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

function [x,it] = gaussSeidel(A,b,x,nIter,eps)
    % We don't need to "split" into (k) and (k-1) since the new
    % values are stored over the old values. This makes for cleaner code.
    it = 0;
    xnew = zeros(size(x));
    while it < nIter
        for i = 1:length(x)
            % We include j==i in the sum (which we shouldn't), 
            % so we have to subtract it.
            xnew(i) = x(i) - (A(i,:)*x - b(i))/A(i,i) ;
        end
        if norm((xnew-x),inf) < eps
            break
        end
        x  = xnew;
        it = it + 1;
    end
end

function [x,it] = jacobi(A,b,x,nIter,eps)
    % Jacobi iteration on component form. Written like this we use all
    % the old values stored in x before overwriting them.
    it = 0;
    while it < nIter
        xnew = x - (A*x - b)./diag(A);
        if norm((xnew-x),inf) < eps
            break
        end
        x  = xnew;
        it = it + 1;
    end
end

function [x,it] = cg(A,b,x,nIter,eps)
    it = 0;
    r = b - A*x;
    d = r;
    while it < nIter
        rTr   = r'*r;
        alpha = rTr/(d'*A*d);
        x     = x + alpha*d;
        r     = r - alpha*A*d;
        beta  = r'*r/rTr; 
        d     = r + beta*d;
        it    = it + 1;
        if norm(r,inf) < eps
            break
        end
    end
end

function [U,it] = SOR(U,w,nIter,eps)
    Uk = U;
    [m,n] = size(U);
    r  = zeros(m,n);
    xi = 2:n-1;
    yj = 2:m-1;
    it = 0;
    while it < nIter                         
        r(xi,yj) = 0.25*(U(xi+1,yj  )+U(xi-1,yj  ) +...
                         U(xi  ,yj+1)+U(xi  ,yj-1) - 4*U(xi,yj));
        Uk(xi,yj) = U(xi,yj) + w*r(xi,yj);
        it = it + 1;
        if norm(U-Uk,inf) < eps
            break
        end
        U = Uk;
    end
end


