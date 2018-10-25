clear;
f = @(x) sin(x);

a =  0; % Start point
b = pi; % End point
n = 20; % Number of nodes

R = romberg(f,n,a,b);

fprintf('Romberg integration: %8.7f\n',R(n,n));

function R = romberg(f,n,a,b)
%Performs naive Romberg integration. 
    R = zeros(n,n); % We are not concerned about memory
    m = 1;          % Starting number of intervals
    h = (b-a)/m;    % Step length
    R(1,1) =  0.5*h*(f(a)-f(b));
    for k = 2:n
        % Calculate the first column of R.
        % This is the only part that requires function evaluations.
        h = 0.5*h; % halve the step length
        r = 0;
        for i = 1:m
           r = r + f(a+(2*i-1)*h);
        end
        R(k,1) = 0.5*R(k-1) + h*r;
        m = 2*m; % Double the number of intervals
    end
    % Calculate the higher orders using recurssion.
    for j = 2:n
        for k = j:n
            R(k,j) = R(k,j-1) + (R(k,j-1)-R(k-1,j-1))/((4^j)-1);
        end
    end       
end
