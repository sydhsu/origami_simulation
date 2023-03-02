% ======================================================================= %
%               Starshade OS testbed deployment simulation
% Author: Sydney Hsu
% Morphing Space Structures Lab
% AA 290 Winter 2023
% 
% Ref: K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid    %
%      origami - An efficient computational approach.' PRSA.              %
%      K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to  % 
%      capture highly nonlinear behavior of non-rigid origami.'           %
%      Proceedings of IASS Annual Symposium 2016.                         %
%      E. T. Filipov, K. Liu, T. Tachi, M. Schenk, G. H. Paulino (2017).  %
%      'Bar and hinge models for scalable analysis of origami.'  IJSS     %
%      K. Liu, G. H. Paulino (2018). 'Highly efficient nonlinear          %
%      structural analysis of origami assemblages using the MERLIN2       %
%      software.' Origami^7.                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define geometry

close all;

tic

% Origami structure to analyze: 'simple_flasher' or 'starshade_os'
origamiStructure = 'simple_flasher';
% 1 = Displacement, 2 = Force, 3 = Adaptive load
controlsetting = 4;
% initial configuration: 'folded' or 'unfolded'
initialState = 'folded';

% Option to show axes labeled
showAxes = true;

% Access Flasher pattern gen. functions
addpath('..\crease_pattern_generator_package');

% Design Flasher
global n_flasher;
global N_flasher;

if strcmp(origamiStructure,'simple_flasher')
    A_flasher   = 400;  % radius of inner polygon
    N_flasher   = 6;   % number of sides of polygon
    h_flasher   = 10;   % thickness
    n_flasher   = 7;   % number of minor folds
    [nodes_planar, nodes_folded, edges, faces] = flasher(N_flasher, n_flasher, h_flasher, A_flasher);
elseif strcmp(origamiStructure,'starshade_os')
    A_flasher   = .734; % [m] radius of inner polygon
    N_flasher   = 14; % number of sides of polygon
    h_flasher   = 16e-3; % [m] thickness
    n_flasher = 14; % number of minor folds
    phi_flasher = 5.45; % [deg]
    [nodes_planar, nodes_folded, edges, faces] = flasherConical_v3(N_flasher, n_flasher, h_flasher, A_flasher,phi_flasher*pi/180);
else
    error('Undefined origami structure');
end

flasher_angles = foldedCreaseAngles(nodes_folded, edges, faces);
global R_planar;
R_planar  = max(sqrt(nodes_planar(1, :).^2 + nodes_planar(2, :).^2)); % deployed radius (in the x-y plane)

% Flasher visualization using pattern generation functions
% plot2dNodesEdges(nodes_planar, edges, flasher_angles);
% plotCreasePatternWithLabels(nodes_planar, edges);
% plot3dNodesEdges(nodes_folded, edges, flasher_angles);

% Re-formatting data from pattern generation function flasher() 
if strcmp(initialState,'folded')
    Node = nodes_folded';
elseif strcmp(initialState,'unfolded')
    Node = nodes_planar';
else
    error('Undefined initial state');
end

Panel = num2cell(faces,2); % TODO: quad faces

% Visualize initial configuration 
figure()
PlotOri(Node,Panel,[],'PanelColor','g')
axis equal; 
axis off;
rotate3d on;
light
% Inspect nodal index assignment
figure()
PlotOri(Node,Panel,[],'ShowNumber','on');
axis equal
rotate3d on;

%% Set up boundary conditions
m = size(Node,1); % total number of nodes

% Adaptive load magnitude
global loadMag; % used in RadialLoad.m
loadMag = 2; % [N] (+) value defined as outwards force

% Flasher boundary conditions
centralNodes = (1:n_flasher:m)'; % nodes to fix; follows numbering system from flasher()
constrainDOF = ones(length(centralNodes),3);
Supp = [centralNodes,constrainDOF]; % fixed nodes at the central polygon

%% Define material and modeling parameters
% Simulation options using the N5B8 model
if controlsetting == 1 % Displacement mode settings
    AnalyInputOpt = struct(...
        'ModelType','N5B8',...
        'MaterCalib','auto',...  
        'ModElastic', 1e3,... % [MPa] (1 GPa)
        'Poisson', 0.3,...
        'Thickness', 0.25,... % [mm]
        'LScaleFactor', 2,...
        'LoadType','Displacement',...  % Displacement load
        'DispStep',200);
elseif controlsetting == 2 % Force mode settings
    AnalyInputOpt = struct(...
        'ModelType','N5B8',...
        'MaterCalib','auto',...  
        'ModElastic', 1e3,... % [MPa] (1 GPa)
        'Poisson', 0.3,...
        'Thickness', 0.25,... % [mm]
        'LScaleFactor', 2,...
        'LoadType','Force',...  % Force load
        'InitialLoadFactor', 0.00001,...
        'MaxIcr', 100);
elseif controlsetting == 3 % Adaptive Load auto mode settings
    AnalyInputOpt = struct(...
        'ModelType','N5B8',...
        'MaterCalib','auto',...  
        'ModElastic', 1e3,... % [MPa] (1 GPa)
        'Poisson', 0.3,...
        'Thickness', 0.25,... % [mm]
        'LScaleFactor', 2,...
        'AdaptiveLoad',@RadialLoad,...
        'InitialLoadFactor', 0.00001,...
        'MaxIcr', 100); 
elseif controlsetting == 4 % Adaptive load force mode settings (overrides Load)
    C0 = 1e10; % initial tangent modulus of bar
    Abar = 3.2e-5; % cross-sectional area of bar [m^2] (32x1 mm)
    Kb = 1e1; % [Nm/rad]
    Kf = 1e-5; % [Nm/rad]
    thetaMin = 0; % [deg]
    thetaMax = 360; % [deg]

    AnalyInputOpt = struct(...
        'ModelType','N5B8',... % triangulation mode
        'MaterCalib','manual',... 
        'BarCM',@(Ex)Ogden(Ex, C0),...
        'Abar',Abar,...
        'Kb',Kb,...
        'Kf',Kf,...
        'RotSprBend', @(he,h0,Kb,L0)EnhancedLinear(he,h0,Kb,L0,thetaMin,thetaMax),...
        'RotSprFold', @(he,h0,Kf,L0)EnhancedLinear(he,h0,Kf,L0,thetaMin,thetaMax),...
        'ZeroBend','Flat',...
        'LoadType','Force',...  % Force load
        'InitialLoadFactor', 0.00001,...
        'MaxIcr', 200,...
        'AdaptiveLoad',@RadialLoad,...
        'StopCriterion',@(Node,U,icrm)DeployToRadius(Node,U,icrm));
end

% Simulation options using the N4B5 model
% AnalyInputOpt = struct(...
%     'ModelType','N4B5',...
%     'MaterCalib','manual',... 
%     'BarCM', @(Ex)Ogden(Ex, 1e4),...
%     'Abar', 2e-1,...
%     'Kb',0.3,...
%     'Kf',0.03,...
%     'RotSprBend', @(he,h0,Kb,L0)EnhancedLinear(he,h0,Kb,L0,30,330),...
%     'RotSprFold', @(he,h0,Kf,L0)EnhancedLinear(he,h0,Kf,L0,30,330),...
%     'LoadType','Displacement',...    % Displacement load
%     'DispStep', 200);

%% Perform analysis
% Assemble input data
[truss, angles, AnalyInputOpt] = PrepareData(Node,Panel,Supp,[],AnalyInputOpt);
% Specify initial deformation state
truss.U0 = zeros(3*size(truss.Node,1),1);
% Perform path-following analysis
[Uhis,Fhis] = PathAnalysis(truss,angles,AnalyInputOpt);
% Postprocess output data
Uhis = real(Uhis);
Fhis = real(Fhis);
STAT = PostProcess(Uhis,truss,angles); 

%% Visualize simulation
% This is the node it plots for. [node idx #, coord. direction]
nodeDOF = [14,-3]; 

interv = 1; % step interval
endicrm = size(Uhis,2); % endicrm = N_icrm

% Animation monitoring node-wise change
VIntensityDataInten = zeros(size(truss.Node,1),size(Uhis,2));
IntensityDataM = bsxfun(@times,STAT.bar.Sx,truss.A);
for k = 1:size(Uhis,2)
    IntensityDataIntenk = sparse(truss.Bars(:,1),truss.Bars(:,2),abs(IntensityDataM(:,k)),size(truss.Node,1),size(truss.Node,1));
    VIntensityDataInten(:,k) = sum((IntensityDataIntenk+IntensityDataIntenk'),2); 
end
% Plot node-wise change and animated load factor vs. displacement for specified node
%VisualFold(Uhis(:,1:interv:endicrm),truss,angles,Fhis(1:interv:endicrm,:),nodeDOF,'IntensityMap','Vertex','IntensityData',VIntensityDataInten)

% Plot node-wise change only
VisualFold(Uhis(:,1:interv:endicrm),truss,angles,[],[],'IntensityMap','Vertex','IntensityData',VIntensityDataInten)
if showAxes
    axis on; xlabel('x'); ylabel('y'); zlabel('z');
end 
%%
% Animation monitoring panel-wise value change
%VisualFold(Uhis(:,1:10:endicrm),truss,angles,[],[],'IntensityMap','Edge','IntensityData',STAT.bar.Sx(:,1:10:endicrm),'ShowInitial','off')
% Animation only showing the configurational change
VisualFold(Uhis(:,1:interv:endicrm),truss,angles,[],[])
% To visualize the strains in bars:
% VisualFold(Uhis(:,1:10:endicrm),truss,angles,[],[],'IntensityMap','Edge','IntensityData',STAT.bar.Sx(:,1:10:endicrm),'ShowInitial','off')

%% Plot diagrams
% Plot displacement of an individual node along the x, y, or z coordinate (not animated version of VisualFold)
dspSign = sign(nodeDOF(2)); % + or -
nodeDOFIdx = nodeDOF(1)*3-(3-abs(nodeDOF(2))); % index of specific node's specific DOF in Uhis, Fhis data structures
dspByDOF = dspSign * Uhis(nodeDOFIdx,:);

% Displacement vector magnitude of given node
% index of node x-coord: nodeDOF(1)*3-2
% x, y, z: -2, -1, -0
nodeZIdx = nodeDOF(1)*3;
nodeYIdx = nodeDOF(1)*3-1;
nodeXIdx = nodeDOF(1)*3-2;
nodeXYZUhis = Uhis(nodeXIdx:nodeZIdx,:); % x, y and z coordinate displ. history
dspMag = vecnorm(nodeXYZUhis,2,1); % magnitude of each displacement 
U_final = vecnorm(nodeXYZUhis(:,end)); % magnitude of final displacement vector
dspPct = dspMag ./ U_final; % fraction of total displacement % TODO: FIX

% Load vs displacement
% figure()
% plot(dspMag,Fhis*loadMag,'b-','linewidth',1); % in force mode: Fhis is load factor; in displ: (-) resistant forces in DOF in nodes imposed with displacement loads
% axis tight
% xlabel('Displacement [m]','fontsize',14); % (in x, y or z coord)
% ylabel('Load [N]','fontsize',14);
% legend(strcat("Node ",num2str(nodeDOF(1))))
% 
figure()
plot(dspPct,Fhis*loadMag,'b-','linewidth',1); % in force mode: Fhis is load factor; in displ: (-) resistant forces in DOF in nodes imposed with displacement loads
axis tight
xlabel('Fraction of total displacement','fontsize',14); % (in x, y or z coord)
ylabel('Load [N]','fontsize',14);
legend(strcat("Node ",num2str(nodeDOF(1))))

% Stored energy vs displacement
figure()
plot(dspPct,STAT.PE,'r-','linewidth',2);    % Red line is the total energy.
hold on                                  % Between red and cyan is the folding energy. 
plot(dspPct,STAT.bend.UB+STAT.bar.US,'c-'); % Between cyan and green is the portion of energy for bending.
plot(dspPct,STAT.bar.US,'g-');              % Below green is the stretching energy of bars.
axis tight
hold on

% Shade energy areas
fillX = [dspPct fliplr(dspPct)];
foldingEnergy = [STAT.PE, fliplr(STAT.bend.UB+STAT.bar.US)];
fill(fillX,foldingEnergy,'r','FaceAlpha',.2);
bendingEnergy = [STAT.bend.UB+STAT.bar.US, fliplr(STAT.bar.US)];
fill(fillX,bendingEnergy,'c','FaceAlpha',.2);
stretchingEnergy = [STAT.bar.US, zeros(size(STAT.bar.US))];
fill(fillX,stretchingEnergy,'g','FaceAlpha',.2);

xlabel('Fraction of total displacement','fontsize',14);
ylabel('Stored Energy [Nm]','fontsize',14);
legend('U_F', 'U_B', 'U_S');

%% Plot final configuration (no animation)
Ux = Uhis(:,end);
Nodef = truss.Node;
Nodef(:,1) = truss.Node(:,1)+Ux(1:3:end);
Nodef(:,2) = truss.Node(:,2)+Ux(2:3:end);
Nodef(:,3) = truss.Node(:,3)+Ux(3:3:end); 

% figure()
% % plot initial configuration
% PlotOri(truss.Node,angles.Panel,truss.Trigl,'FoldEdgeStyle','-','EdgeShade',0.3,'PanelColor','none');
% % plot deformed configuration
% PlotOri(Nodef,angles.Panel,truss.Trigl);
% axis equal; axis off;
% camproj('perspective')
% light
% view(117,18)
% rotate3d on

%% Export final configuration to an OBJ file
% Write2OBJ('BendedMiura5x5', Nodew, truss.Trigl, truss.Bars, angles.bend)

toc

