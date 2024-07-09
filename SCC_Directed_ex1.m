clc
close all;
clear;

g=SCCAnanyzerClass;
% Define the graph using the vectors of starting nodes and terminating nodes
s = [1 2 2 3 3 4 4 4 4 5 5 6 6 6 7 8 8 9];
t = [2 3 4 1 4 1 2 3 5 4 6 5 7 8 5 6 9 8];
G = digraph(s,t);
A1a = adjacency(G);
save('Ex1a.mat','A1a');

% Plot the graph
[n,m]=size(A1a);
XData = [1, 0, 1, 2, 3, 5, 4, 6, 7];
YData = [1.5, 1, 0, 1, 1, 1, 0, 1, 1];

h = plot(G, 'XData', XData, 'YData', YData, 'NodeColor', 'blue');


% --------
% Step 1
disp("Step 1");
[omega,disc,zeta]=g.NeighborStructure(G,n);
zeta_outcome=g.MatrixOutput(zeta(:,:,n+1))
zeta_depth = g.MatrixOutput(omega(:,:,n+1))
disp("strong connectivity=")
disp(disc) % strongly connected or not

% --------
% Step 2
disp("Step 2");
[ci]= g.InformationNumber(G,n);
eta_outcome=g.MatrixOutput(ci(:,:,end))

% --------
%Step 3
disp("Step 3");
[SC, SP]=g.findingSCC(G,n)
[Vsour, Vsink, Visol, Vmid] = g.indexSCC(G,n)

% --------
%Step 4a
disp("Step 4a");
[theta, ~] = g.SCCStructure(G, n);
theta_outcome = g.MatrixOutput(theta(:,:,n+1))

disp("Step 4b");
[indexset,thetaRdd, nuRdd] = g.SCCReducedStructure(G, n);
indexset
zeta_reduced=g.MatrixOutput(thetaRdd(:,:,end))
zeta_depth_reduced=g.MatrixOutput(nuRdd(:,:,end))

disp("Step 4c");
% Create a cell array to store node labels
nodeLabels = cell(1, numel(indexset));

% Populate the cell array with labels in the format "x(i)"
for i = 1:numel(indexset)
    nodeLabels{i} = sprintf('%d(%d)', indexset(i), i);
end
nRdd=size(thetaRdd(:,:,end),2);
xDataSs=XData(indexset);
yDataSs=YData(indexset);

% GnuRdd = reduced graph of SCCs
[nGdd, mGdd] = size(nuRdd(:, :, end));
GnuRdd = zeros(nGdd, nGdd);
for j = 1:nGdd
    for i= 1:nGdd
        if nuRdd(i, j, end) == 1
            % Condition for direct links
            GnuRdd(i, j) = nuRdd(i, j, end);
        else
            GnuRdd(i, j) = 0;
        end
    end
end

Grdd=digraph(GnuRdd');

figure()
Rddh = plot(Grdd,'XData', xDataSs, 'YData', yDataSs, 'NodeLabel', indexset);
saveas(gcf, 'red1a_plot.png');

