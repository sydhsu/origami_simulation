clear all; close all; clc;

% A   = 6e-2;
% N   = 9;
% h   = 5e-3;
% k   = 0.5;
% n   = 12;

% A   = 787.4/2;
% N   = 12;
% h   = 16;
% k   = 0.5;
% n   = 21;

A   = 800;  % radius of inner polygon
N   = 12;   % number of sides of polygon
h   = 32;   % thickness
n   = 11;   % number of minor folds

% A   = 400;
% N   = 8;
% h   = 50;
% k   = 0.5;
% n   = 7;

[nodes_planar, nodes_folded, edges, faces] = flasher(N, n, h, A);

lengths     = getEdgeLengths(nodes_planar, edges);
lengths_p   = getEdgeLengths(nodes_folded, edges);
error       = (lengths - lengths_p)./lengths_p*100;

figure();
plot(1:length(edges), error, 'o--')
xlabel('Edge index')
ylabel('Length error (percent)')

angles = foldedCreaseAngles(nodes_folded, edges, faces);

plot2dNodesEdges(nodes_planar, edges, angles);
f3 = plot3dNodesEdges(nodes_folded, edges, angles);
%plotCreasePatternWithLabels(nodes_planar, edges)

D_deployed  = 2*max(sqrt(nodes_planar(1, :).^2 + nodes_planar(2, :).^2));


figure(f3);
p1 = patch();
p1.Faces = faces;
p1.Vertices = nodes_folded';
p1.FaceAlpha = 0.9;
p1.FaceColor = [0.7 0.7 0.7];
p1.LineStyle = 'none';
axis equal
axis off
view(3)

% figure();
% p2 = patch();
% p2.Faces = faces;
% p2.Vertices = nodes_planar';
% p2.FaceAlpha = 0.2;
% p1.FaceColor = [0.7 0.7 0.7];
% axis equal
% axis off 