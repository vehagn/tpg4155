function [ df ] = firstDerivative(F,x,dx,type)
%Calculates the first dericate using either the midpoint difference,
%forward difference, or backward difference.
if (nargin < 4); type = 'mid'; end
if (nargin < 3); dx   =  1e-8; end

switch type
    case 'mid'
        df = (F(x-dx) - F(x+dx))/(2*norm(dx));
    case 'forw'
        df = (F(x+dx) - F(x   ))/(norm(dx));
    case 'back'
        df = (F(x)    - F(x-dx))/(norm(dx));
    otherwise
        df = 0;
        error('Derivation type not supported')
end
end

