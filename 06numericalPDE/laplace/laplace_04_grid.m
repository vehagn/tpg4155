%% Solve the 2-D Laplce equation using a naive iterative approach.

%% Initialisation
clear;
n = 5; % nodes on x-axis
m = 5; % nodes on y-axis
eps = 1e-16;     % Error threshold
nIter = 1000;

x = linspace(0,0.5,n);
y = linspace(0,0.5,m);

U = zeros(n,m);

%% Insert boundary conditions
U(1,:)     = 200*x;       % top part
U(:,end)   = 200*flip(y); % right part

%% Preconditioning
% Uncomment for a better initial guess
% U(2:end-1,2:end-1) = 25;
% Uncomment for exact guess
% U(2:end-1,2:end-1) = flipud((400*x(2:end-1)'*y(2:end-1)));

%% Start iterative loop
it = 0;
Unew = U;
while it < nIter
%%  Using nested for-loops
%    for i = 2:n-1
%        for j = 2:m-1
%            Unew(i,j) = 0.25*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1));
%        end
%    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Using slicing (will often result in faster running code)
    Unew(2:end-1,2:end-1) = 0.25*(U(3:end  ,2:end-1)+U(1:end-2,2:end-1) +...
                                  U(2:end-1,3:end  )+U(2:end-1,1:end-2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
    it = it + 1;
    if norm(U-Unew,inf) < eps
        break
    end
    U = Unew;
end
fprintf('Numer of iterations: %4i\n',it);
U = Unew;

imagesc(x,y,flipud(U));
set(gca,'YDir','normal')
xlabel('x');
ylabel('y');
colorbar();