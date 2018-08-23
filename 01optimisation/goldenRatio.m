clear; close;           % Make sure that we don't have any stray variables

f = @(x) x.*x - sin(x); % Define the funciton f
x = linspace(-.2,1.2);  % Create an intervall x to plot f(x)

r = (sqrt(5)-1)/2;      % The inverse golden ratio
n = 100;                % Amount of iterations before we give up

eps = 1e-12;            % Error tolerance

% Starting values for our interval [a,b] and internal point c,d
a(1) = 0;
b(1) = 1;
alpha(1) = a(1)+(1-r)*(b(1)-a(1));
beta(1)  = b(1)-(1-r)*(b(1)-a(1));

% Evaluate the function at point c and d.
falpha = f(alpha(1)); 
fbeta  = f(beta(1));

% Plotting routine for the graph.
figure(1);
hold on
hf = plot(x,f(x),'k-');
%%%% We cheat a bit to find the true minimum point
plot(x(f(x) == min(f(x))),min(f(x)),'ko');
hmin = plot([x(f(x) == min(f(x))),x(f(x) == min(f(x)))],[-.4,0.6],'k-.');
%%%% Ignore this part as it's mostly for demonstration purposes
title('$$x^2 - \sin(x)$$','interpreter','latex');
ylabel('$$f(x)$$','interpreter','latex');
xlabel('$$x$$','interpreter','latex')

for i = 2:n
   if falpha <= fbeta
       % New endpoints
       a(i) = a(i-1);                     % Keep our a point
       b(i) = beta(i-1);                  % Our new b will be the old d
       % New internal points
       alpha(i) = a(i)+(1-r)*(b(i)-a(i)); % Calculate a new c
       beta(i) = alpha(i-1);              % Our new d point will be the old c
       % One new function evaluation
       fbeta = falpha;                    % The new f(d) is the old value of f(c)
       falpha = f(alpha(i));              % Calculate a new f(c) value based on our new point
   else
       % New endpoints
       a(i) = alpha(i-1);                 % Our new a will be the old b
       b(i) = b(i-1);                     % Keep our b point
       % New internal points
       alpha(i) = beta(i-1);              % Our new c point will be the old d
       beta(i) = b(i)-(1-r)*(b(i)-a(i));  % Calculate a new d
       % One new function evaluation
       falpha = fbeta;                    % The new f(c) is the old value of f(d)
       fbeta = f(beta(i));                % We have to calculate the new f(d)
   end
   if abs(alpha(i)-beta(i)) < eps
       % End the for loop if we've reached convergence.
       disp(['Convergence reached after ', num2str(i), ' iterations.']); 
       break
   end
end
% Plot how the interval [a,b] and internal points alpha, beta change for each
% iteration
yyaxis right
ha     = plot(a,     1:i, '>:','color',[     0    0.4470    0.7410]); 
hb     = plot(b,     1:i, '<:','color',[0.8500    0.3250    0.0980]);
halpha = plot(alpha, 1:i, '>:','color',[0.9290    0.6940    0.1250]); 
hbeta  = plot(beta,  1:i, '<:','color',[0.4940    0.1840    0.5560]);
ylabel('Iteration number');
legend([hf, ha, hb, halpha, hbeta, hmin],'f(x)','a','b','c alpha','beta',...
    'min','location','northwest'); 
xlim([-0.2,1.2]);