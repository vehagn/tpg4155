clear;
A = [ 4, -1, 1;
      4, -8, 1;
     -2,  1, 5];
 
b = [7; -21; 15];

x0 = [1; 2; 2];
xT = [2;4;3];

N = 10;

checkConvergence(A);

x = x0;
fprintf('Jacobi:\n');
fprintf('It: %2i\tx_k: %7.6f\ty_k: %7.6f\tz_k: %7.6f\t||x^k-x^T|| %1.4e\n',...
    0, x, sqrt(sum((x - xT).^2)));
for k = 1:N
    x = jacobi(A,b,x);
    fprintf('It: %2i\tx_k: %7.6f\ty_k: %7.6f\tz_k: %7.6f\t||x^k-x^T|| %1.4e\n',...
        k, x, sqrt(sum((x - xT).^2)));
end

x = x0; % reset starting vector
fprintf('Gauss-Seidel:\n');
fprintf('It: %2i\tx_k: %7.6f\ty_k: %7.6f\tz_k: %7.6f\t||x^k-x^T|| %1.4e\n',...
    0, x, sqrt(sum((x - xT).^2)));
for k = 1:N
    x = gaussSeidel2(A,b,x);
    fprintf('It: %2i\tx_k: %7.6f\ty_k: %7.6f\tz_k: %7.6f\t||x^k-x^T|| %1.4e\n',...
        k, x, sqrt(sum((x - xT).^2)));
end

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

function x = jacobi(A,b,x)
    % Jacobi iteration on component form. Written like this we use all
    % the old values stored in x before overwriting them.
    x = x - (A*x - b)./diag(A) ;
end

function xO = gaussSeidel1(A,b,x)
    % Naïve implementation of the Gauss-Seidel algorithm
    xO = x;
    for i = 1:length(x)
       x(i) = (b(i) - sum(A(i,1:i-1)*x(1:i-1)) - sum(A(i,i+1:end)*xO(i+1:end)))/A(i,i);
    end
    xO = x;
end

function x = gaussSeidel2(A,b,x)
    % We don't need to "split" into (k) and (k-1) since the new
    % values are stored over the old values. This makes for cleaner code.
    for i = 1:length(x)
        % We include the diagonal elements (j==i) in the sum, which we 
        % shouldn't. We therefore have to subtract the diagonal elements.
        % As an added bonus this results in cleaner code.
        x(i) = x(i) + (b(i) - A(i,:)*x)/A(i,i) ;
    end
end

