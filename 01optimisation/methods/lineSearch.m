function a = lineSearch(f,a,b,n,eps,type)
%LINESEARCH
if (nargin < 6); type =    0; end
if (nargin < 5); eps  = 1e-8; end
if (nargin < 4); n    = 1000; end

switch type
    case 0
        r = (sqrt(5)-1)/2;
        for i = 1:n
            c = a+(1-r)*(b-a);
            d = b-(1-r)*(b-a);
            if f(c) <= f(d)
                b = d;
            else
                a = c;
            end
            if abs(c-d) < eps
                a = (c+d)/2;
                return
            end
        end
    otherwise
        disp('Invalid line search method selected!');
end
end

