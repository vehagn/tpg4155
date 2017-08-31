%%%% Example 8.6
%% Preamble
clear; clc;
f =@(x,y) x.*x - 4.*x + y.*y - y -x.*y;
F =@(V) f(V(1),V(2));

V  = [0.0 ,0.0; 1.2 ,0.0; 0.0 ,0.8];  % Starting points
FV = [F(V(1,:)),F(V(2,:)),F(V(3,:))]; % Starting values
% Sort the starting points s.t. F(B) < F(G) < F(W)
[B,G,W,FB,FG,FW] = sortPoints(V,FV);

[X,Y] = meshgrid(linspace(-1,6),linspace(-2,5));
printDetails(B,G,W,FB,FG,FW,0);
plotFunc(X,Y,f,B,G,W);

i   = 0;    % Iterator
n   = 50;   % Maximum number of iterations
eps = 2e-4; % Error tolerance

%% Main loop
while (i < n) && (sum(sqrt((B-W).^2)) > eps)
    i = i + 1; % Manually increment i
    % Run Nelder-Mead algorithm
    [B,G,W,FB,FG,FW] = nelderMeadAlg(B,G,W,FB,FG,FW,F);
    % Sort vertices
    % The sorting can be done during Nelder-Mead to save computational
    % time, but for demonstrational purposes this is not done here.
    [B,G,W,FB,FG,FW] = sortPoints([B;G;W],[FB;FG;FW]);
    % Print details about our vertices
    printDetails(B,G,W,FB,FG,FW,i);
    if i == 10 % Zoom in on the plot
        [X,Y] = meshgrid(linspace(2.5,3.5),linspace(1.5,2.5)); 
    end
    if i == 20 % Zoom in even further on the plot
        [X,Y] = meshgrid(linspace(2.95,3.05),linspace(1.95,2.05)); 
    end
    plotFunc(X,Y,f,B,G,W);

end
fprintf('\nThe best point found was: (%10.10f, %10.10f) ',B(1),B(2));
fprintf(' with the function value: %10.10f\n',FB);

%% Help functions
% Functions defined inside scripts can only be seen by that script

%% Nelder-Mead Algorithm
function [B,G,W,FB,FG,FW] = nelderMeadAlg(B,G,W,FB,FG,FW,F)
    M = (B+G)/2;
    R = 2.*M - W;
    FR = F(R);
    if FR < FG
        % Reflect and possibly extend
        if FB < FR
            % Reflect
            W = R; FW = FR;
            fprintf('Reflection!\n');
        else
            % Try to expand since R is the best point
            E = 2.*R - M;
            FE = F(E);
            if FE < FB
                % Expansion is better, replace worst with expanded
                W = E; FW = FE;
                fprintf('Expansion!\n');                
            else           
                % Expansion is not better, replace worst with reflected
                W = R; FW = FR;
                fprintf('Reflection!\n');
            end
        end
    else
        % Contraction
        if FR < FW
            % Choose C2 if F(R) < F(W)
            W = R; FW = FR;
            fprintf('Reflection-');
        end
        C = (W+M)/2; FC = F(C);
        if FC < FW
            % Contract if F(C) < F(W)
            W = C; FW = FC;
            fprintf('Contraction!\n');
        else
            % Shrink towards B if F(C) > F(W)
            S = (B+W)/2; FS = F(S);
            W = S; FW = FS;
            G = M; FG = F(M);
            fprintf('Shrink!\n');
        end
    end
end

%% Sorting algorithm
function [B,G,W,FB,FG,FW] = sortPoints(V,FV)
    %Sort the points according to function value
    V1 = V(1,:); F1 = FV(1);
    V2 = V(2,:); F2 = FV(2);
    V3 = V(3,:); F3 = FV(3);
    
    if (F3 <= F1) && (F3 <= F2)
        if (F2 <= F1)
            % F3 < F2 < F1
            [B,G,W]    = deal(V3,V2,V1);
            [FB,FG,FW] = deal(F3,F2,F1);
        else
            % F3 < F1 < F2
            [B,G,W]    = deal(V3, V1, V2);
            [FB,FG,FW] = deal(F3,F1,F2);
        end
    elseif (F1 <= F2) && (F1 <= F3)
        if (F2 <= F3)
            % F1 < F2 < F3
            [B,G,W]    = deal(V1,V2,V3);
            [FB,FG,FW] = deal(F1,F2,F3);
        else
            % F1 < F3 < F2
            [B,G,W]    = deal(V1,V3,V2);
            [FB,FG,FW] = deal(F1,F3,F2);
        end
    else %(F2 < F1) && (F2 < F3)
       if (F1 <= F3)
           % F2 < F1 < F3
           [B,G,W] = deal(V2,V1,V3);
           [FB,FG,FW] = deal(F2,F1,F3);
       else
           % F2 < F3 < F1
           [B,G,W] = deal(V2,V3,V1);
           [FB,FG,FW] = deal(F2,F3,F1);
       end
    end
end

%% Printing function
function printDetails(B,G,W,FB,FG,FW,i)
    fprintf('Iteration %3i:\t',i);
    fprintf('|B: %4.2f, %4.2f : %7.6f |\t',B(1),B(2),FB);
    fprintf('|G: %4.2f, %4.2f : %7.6f |\t',G(1),G(2),FG);
    fprintf('|W: %4.2f, %4.2f : %7.6f |\n',W(1),W(2),FW);
end

%% Plotting function
function plotFunc(X,Y,f,B,G,W)
    figure(1); clf;
    hold on
    contourf(X,Y,f(X,Y),30);
    title('Nelder-Mead Algorithm')
    xlabel('x-axis');
    ylabel('y-axis');
    axis image
    colorbar();
    plot([B(1), G(1), W(1), B(1)],...
         [B(2), G(2), W(2), B(2)],...
         '-',...
         'color',[0.8500    0.3250    0.0980],...
         'LineWidth', 2.5);
    text(B(1),B(2),'B','color',[1    1    1]);
    text(G(1),G(2),'G','color',[1    1    1]);
    text(W(1),W(2),'W','color',[1    1    1]);
    pause(1);
end