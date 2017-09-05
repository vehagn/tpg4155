clear;
A = [ 4, -1, 1;
      4, -8, 1;
     -2,  1, 5];
 
b = [7; -21; 15];

x0 = [1; 2; 2];
xT = [2;4;3];

N = 10;
x = x0;

fprintf('It: %2i\tx_k: %7.6f\ty_k: %7.6f\tz_k: %7.6f\t||x^k-x^T|| %1.4e\n',...
    0, x, sqrt(sum((x - xT).^2)));
fprintf('Gauss-Seidel:\n');
for k = 1:N
    x = gaussSeidel(A,b,x);
    fprintf('It: %2i\tx_k: %7.6f\ty_k: %7.6f\tz_k: %7.6f\t||x^k-x^T|| %1.4e\n',...
        k, x, sqrt(sum((x - xT).^2)));
end


N = 10;
x = x0;
fprintf('Jacobi:\n');
for k = 1:N
    x = jacobi(A,b,x);
    fprintf('It: %2i\tx_k: %7.6f\ty_k: %7.6f\tz_k: %7.6f\t||x^k-x^T|| %1.4e\n',...
        k, x, sqrt(sum((x - xT).^2)));
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