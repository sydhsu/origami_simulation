% Radial load vector generation
% ======================================================================= %
% Inputs: Node - array of initial node coordinates
%         U    - displacement with respect to origin for current time increment
%         icrm - incremental number (pseudo-time)
% Output: F - load column vector with dimensions (NDOF x # nodes) by 1
% ======================================================================= %

function [F] = RadialLoad(Node,U,icrm)
%% Compute current configuration 
Nodenw = Node;
Nodenw(:,1) = Node(:,1)+U(1:3:end);
Nodenw(:,2) = Node(:,2)+U(2:3:end);
Nodenw(:,3) = Node(:,3)+U(3:3:end);

%% Define Load here based on icrm and Nodenw
if icrm<=0 
    error('Wrong increment!'); 
else
    % Define load function
    global loadMag;
    global n_flasher;
    global N_flasher;

    % Apply loads to outer nodes only
    outerNodesIdx = (n_flasher:n_flasher:(n_flasher*N_flasher*2))';
    outerNodes = Nodenw(outerNodesIdx,:);
    outerNodes(:,3) = 0; % set all z-components to 0
    % load: magnitude = loadMag, direction = radial vector from origin
    load = outerNodes ./ vecnorm(outerNodes,2,2) * loadMag;

    % Ramp down load if it exceeds 95% of deployed radius
    global R_planar;
    norm2d = vecnorm(outerNodes(:,1:2),2,2);
    rampDown = norm2d >= R_planar*.95;
    if any(rampDown)
        fprintf("Ramping down load \n");
    end
    load(rampDown,:) = load(rampDown,:)*.1; % 10% of load

end

%% Wrap up Load info to force/displacement vector F
m = size(Node,1);
% Add nodal indices
idx = outerNodesIdx;
Load = [idx load];

F = zeros(3*m,1);
indp = Load(:,1); % nodal indices
F(3*indp-2) = Load(:,2); % x components
F(3*indp-1) = Load(:,3); % y components
F(3*indp) = Load(:,4); % z components
