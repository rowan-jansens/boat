function boat_analysis()
% Note: This code has *almost* everything you need to generate a moment
% curve. Your task is to complete the waterline_error(d) function below.
% The waterline_error(d) function is defined at the bottom of the file;
% this code will not run until you complete it!!!

%% Setup
% Housekeeping
clf % clear the current figure
boatcolor = [0.9290 0.6940 0.1250]; % define the color of the boat
watercolor = [0 0.4470 0.7410]; % define the color of the water

% Design parameters
W = 0.4;            % width  (m)
H = 0.2;            % height (m)
L = 1;            % length (m)
n = 1;               % shape parameter (-)
dispratio = 0.25;    % ratio of boat displacement to max diplacement
CoM_offset = [0; 0]; % Offset of the center of mass, e.g. using ballast

% physical constants
wrho = 1000; % water density kg/m^3
g = 9.8;     % gravity m/s^2

%% One-time computations
% boat definition and key variables
Npts = 200;                        % number of 1D spatial points (not a design variable; don't change)
xPoints = linspace(-W/2,W/2,Npts); % set of points in the x direction (horizontal)
zPoints = linspace(0,H,Npts);      % set of points in the z direction (vertical)

[X, Z] = meshgrid(xPoints, zPoints); % create the meshgrid
P = [X(:)'; Z(:)'];                  % pack the points into a matrix
% find all the points inside the boat - a logical array
insideBoat = transpose(P(2,:) >= H*abs((2*P(1,:)/W)).^n & P(2,:) <= H);

dx = xPoints(2)-xPoints(1); % delta x
dz = zPoints(2)-zPoints(1); % delta z
dA = dx*dz;                 % define the area of each small section

boatmasses = insideBoat*wrho*dA*L; % find the water mass of each small section
maxdisp = sum(boatmasses);         % find the maximum displacement of the boat
boatdisp = dispratio*maxdisp;      % set the displacement of the boat
CoD = P*boatmasses/maxdisp;        % find the centroid of the boat
P = P - CoD;                       % center the boat on the centroid
CoD = CoD - CoD;                   % update the centroid after centering
CoM = CoD + CoM_offset;            % adjust the center of mass, e.g. for ballast

%% Looped computation
dtheta = 1; % define the angle step
R = [cosd(dtheta) -sind(dtheta);  % define rotation matrix
     sind(dtheta)  cosd(dtheta)];

j = 1; % set the counter
for theta = 0:dtheta:180 % loop over the angles
    % Bound the waterline
    dmin = min(P(2, :)); % find the minimum z coordinate of the boat
    dmax = max(P(2, :)); % find the maximum z coordinate of the boat
    % Solve for the waterline of the rotated boat
    d = fzero(@waterline_error, [dmin, dmax]); % find the waterline d

    % Analyze the boat
    underWater = (P(2,:) <= d)';
    underWaterAndInsideBoat = insideBoat & underWater;
    watermasses = underWaterAndInsideBoat*wrho*dA*L;
    watermass = sum(watermasses);
    CoB = P*watermasses./watermass;
    % Store the results
    torque(j) = boatdisp * g * (CoB(1,1) - CoM(1,1)); % find the torque (N m)
    angle(j) = theta;                                 % record the angle (deg)

%     % Uncomment these lines to display a rotated boat;
%     % useful for undertanding, but slows down the computation
%     hold off
%     scatter(P(1,insideBoat),P(2,insideBoat),[],boatcolor), axis equal, axis([-max(W,H) max(W,H) -max(W,H) max(W,H)]), hold on
%     scatter(P(1,underWaterAndInsideBoat),P(2,underWaterAndInsideBoat),[],watercolor)
%     scatter(CoM(1,1), CoM(2,1), 1000, 'r.'); % plot the COM
%     scatter(CoB(1,1), CoB(2,1), 1000, 'k.'); % plot the COB
%     drawnow

    % Rotate the system
    P = R*P;     % rotate the boat by dtheta
    CoM = R*CoM; % rotate the center of mass too
    j = j + 1;   % update the counter
end

%% Plot the torque versus the angle curve
clf
plot(angle, torque)
xlabel('heel angle (degrees)')
ylabel('Torque (N m)')
min(torque)
grid on

%% TODO: Define the bouyancy error function - should be zero when balanced
function res = waterline_error(d)
    underWater = (P(2,:) <= d)'; % test if each part of the meshgrid is under the water
    
   
    
    watermass = sum(underWater & insideBoat) .* dA * L * wrho;


    %(boatArea * rho_boat) - (waterArea * rho_water)

    % TASK: Complete this function by solving for the mass of the water
    % displaced by the boat.

    % HINT: You did *most* of the work necessary in the `boat02_centers`
    % exercise.

    % YOUR CODE HERE

    res = watermass - boatdisp; % difference between boat displacement and water displacement
end

end
