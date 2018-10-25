clear;
f = @(x) sin(x);
% f = @(x) x + x.^2 + x.^3 + x.^4;

a =     0; % Start point
b =    pi; % End point

TOL = 1e-12;
N = 20;

h  = 0.5*(b-a);
fa = f(a);
fc = f(a+h);
fb = f(b);

[I, err, n] = adaptive(f,fa,fc,fb,a,b,h,TOL,simpson(fa,fc,fb,h),0,N);
Imat = integral(f,a,b);

fprintf('The integral of f(x) from %4.3f to %4.3f\n',a,b);
fprintf('Our approximation:    %5.15e\n',I);
fprintf('MATLAB approximation: %5.15e\n',Imat);
fprintf('Difference:           %5.15e\n',abs(I-Imat));
fprintf('Target error:         %5.15e\n',TOL);
fprintf('Pessimistic error:    %5.15e\n',err);
fprintf('Deepest recusion level was: %i\n',n);

function [I, err, n] = adaptive(f,fa,fc,fb,a,b,h,TOL,I,n,N)
    n = n+1;        % Keeping track of how deep in the recursion we are.
    c = 0.5*(a+b);  % Midpoint of the current interval
    h = 0.5*h;      % Halve the stepsize for the sub-intervals
    fd = f(a+h);    % Function evaluated at midpoint of left interval
    fe = f(a+3*h);  % Function evaluated at midpoint of right interval
    % Full interval:      [a--d--c--e--b]
    % Left sub interval:  [a--d--c]
    % Right sub interval:       [c--e--b]
    Ileft  = simpson(fa,fd,fc,h); % Simpson's rule left sub-interval
    Iright = simpson(fc,fe,fb,h); % Simpson's rule right sub-interval
    err = abs(Ileft+Iright-I); % Error estimate from theory
    if n > N
        % In order to avoid infinite recursion we set a deepest level the
        % recursion is allowed to go.
        error('Deepest level exceeded');
    elseif err < 15*TOL % 15 comes from error analysis
        % Richardson extrapolation gives a correcting factor of 
        % (Ileft + Iright - I)/15 resulting in a degree of precision of 5,
        % i.e. exact for polynomials of degree 5 or less.
        I = Ileft+Iright + (Ileft + Iright - I)/15;
        return;
    else
        % If the desired accuracy is not obtained we further sub-divide the
        % interval using a recursive relation.
        [Ileft , Lerr, Li] = adaptive(f,fa,fd,fc,a,c,h,TOL/2,Ileft ,n,N);
        [Iright, Rerr, Ri] = adaptive(f,fc,fe,fb,c,b,h,TOL/2,Iright,n,N);
        I   = Ileft + Iright;   % Addming left and right integral
        err = max(Lerr,Rerr);   % Pessimistic error estimate
        n   = max(Li,Ri);       % Pick the deepest level
    end
end

function I = simpson(fa,fc,fb,h)
    I = h*(fa+4*fc+fb)/3;
end
