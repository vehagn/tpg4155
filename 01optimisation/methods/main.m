clear; close;
%% Preamble
% Testfunctions
f = @(x,y) x+(x-1).^2+(y-1).^2 ;
%f =@(x,y) cos(pi*x) + y.*x.*sin(pi*y) + x.*x - x + y.*y;
%f = @(x,y) cos(pi*x) + sin(pi*y);

% Initial guess
x0 = [-0.3,-1.6];

F =@(v) f(v(1),v(2));
[X,Y] = meshgrid(linspace(-2,2),linspace(-2,2));
figure(100); clf;
hold on
contourf(X,Y,f(X,Y),64,'LineColor','none');
colormap('parula');

%% Run algorithms
coord = coordSearch(F,x0);
powe1 = powell(F,x0);
powe2 = powell2(F,x0);
steep = steepDesc(F,x0);
newto = newton(F,x0);

%% Plot algorithms
hc = plot(coord(:,1),coord(:,2),'bd-');
hp = plot(powe1(:,1),powe1(:,2),'ro-');
hq = plot(powe2(:,1),powe2(:,2),'g*-');
hs = plot(steep(:,1),steep(:,2),'kh-');
hn = plot(newto(:,1),newto(:,2),'cp-');

legend([hc hp hq hs hn],...
    'Coordinate Search', 'Naive Powell', 'True Powell',...
    'Steepest Descent', "Newton's Algorithn",...
    'Location','NorthWest');

%% Output results
fprintf('Coord Search:   %+10.10f at (%+5.4f,%+5.4f) after %3i iterations.\n',F(coord(end,:)),coord(end,:),length(coord));
fprintf('Powell (naive): %+10.10f at (%+5.4f,%+5.4f) after %3i iterations.\n',F(powe1(end,:)),powe1(end,:),length(powe1));
fprintf('Powell (true):  %+10.10f at (%+5.4f,%+5.4f) after %3i iterations.\n',F(powe2(end,:)),powe2(end,:),length(powe2));
fprintf('Steep Descent:  %+10.10f at (%+5.4f,%+5.4f) after %3i iterations.\n',F(steep(end,:)),steep(end,:),length(steep));
fprintf('Newton:         %+10.10f at (%+5.4f,%+5.4f) after %3i iterations.\n',F(newto(end,:)),newto(end,:),length(newto));

