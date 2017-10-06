%% Geometry
z  =  0.05; % [m] - Height of my 1-dimensional steak.
nz =   100; % Number of discrete points.
t  = 30*60; % [s] - Max simulation time.
rt =  5*60; % [s] - Rest time for steak.
dt =   0.1; % [s] - Time-step length.
% Calculate space step length and numer of time steps.
dz =  z/nz; % [m] - Discrete step length.
nt =  t/dt; % Number of timesteps.

%% Constants
alpha = 1.43e-7; % [m^2/s] - Thermal diffusivity of water.
% Calculate beta and check if it meets the requirements.
beta = (alpha*dt)/(dz*dz);
if (beta > 0.5); disp('beta is more than 0.5!'); end

%% Initial conditions
sTemp = 273 +  20;  % [K] - Initial temperature of the steak.
oTemp = 273 + 120;  % [K] - Oven temperature.
aTemp = 273 +  20;  % [K] - Ambient temperature.
fTemp = 273 +  60;  % [K] - Final minimum temperature (Medium rare).
% Initialise array to hold temperatures.
T  = sTemp*ones(1,nz);
T([1,end]) = oTemp; % The outermost gridpoints are oven temperature.

%% Cook the steak
for k = 1:nt
    T = heatTimestep(T,beta); % Do a time step.
    if ~mod(k,nt/500) % Plot every (nt/500)-th time step.
        plotTemp(T,z,k,dt,sTemp,oTemp,1);
    end
    if min(T) > fTemp % Stop simulation if final temperature is met.
        disp(['Finished cooking after '...     % Three periods lets you...
            num2str(floor(k*dt/60)) ' min '... % continue on the next...
            num2str(mod(k*dt,60)) ' s.']);     % line.
        
        break % break exits the current loop.
    end
end

%% Let the  steak rest by changing the end temperatures 
%  to ambient temperature.
disp(['Letting the steak rest for ' num2str(rt/60) ' minutes.']);
T([1,end]) = aTemp;
for k = 1:(rt/dt)
    T = heatTimestep(T,beta); % Do a time step.
    if ~mod(k,nt/500) % Plot every (nt/500)-th time step.
        plotTemp(T,z,k,dt,sTemp,oTemp,1);
    end
end
plotTemp(T,z,k,dt,sTemp,oTemp,1); % Plot one final time.

%% Function definitions
function T = heatTimestep(T,beta)
%HEATTIMESTEP Perform one step forward in time of the heat equation.
%   Forward in Time, Central in Space stencil (FTCS)
%% Using slicing
    z = 2:length(T)-1;
    T(z) = T(z) + beta*(T(z+1)-2*T(z)+T(z-1));
%% Using for loop
%     Tnew = T;
%     for i = 2:length(T)-1
%         Tnew(i) = T(i) + beta*(T(i+1)-2*T(i)+T(i-1));
%     end
%     T = Tnew
end

function plotTemp(T,z,k,dt,minTemp,maxTemp,figNum)
%PLOTTEMP Help function to plot the temperature gradient
    figure(figNum);
    % Convert the length from m to cm and the temparature from K to C.
    plot(100*linspace(0,z,length(T)),T-273);
    % sprintf() works like fprintf(), but instead of displaying the string
    % it creates an array of characters.
    title(sprintf('Time elapsed: %6.2f seconds',k*dt));
    xlabel('Cross section [cm]');
    ylabel('Temperature [C]');
    ylim([minTemp,maxTemp]-273);
    drawnow();
end