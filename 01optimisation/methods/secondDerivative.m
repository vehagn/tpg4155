function [ ddf ] = secondDerivative(F,x,dx,dy,type)
%SECONDDERIVATIVE Summary of this function goes here
%   Detailed explanation goes here
if (nargin < 5); type = 'mid'; end
if (nargin < 4); dy   =  1e-8; end
if (nargin < 3); dx   =    dy; end

switch type
    case 'mid'
        ddf = (F(x+dx+dy)-F(x-dx+dy)-F(x+dx-dy)+F(x-dx-dy))/(4*norm(dx)*norm(dy));
    otherwise
        ddf = 0;
        error('Derivation type not supported')
end
end

