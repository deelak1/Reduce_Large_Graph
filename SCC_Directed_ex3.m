clc
close all;
clear;

g=SCCAnanyzerClass;
% Define the graph using the vectors of starting nodes and terminating nodes
s = [1 2 3 4 5 5 7 8 9 9 9 10 12 15 15 13 13 14 17 18 18 16 19 20];
t = [7 9 8 7 3 10 4 5 6 7 10 9 15 12 13 15 11 17 18 14 16 18 20 19];
G = digraph(s,t);
A3 = adjacency(G);
save('Ex3.mat','A3');

% Plot the graph
[n,m]=size(A3);
XData = [3, 0, 1, 4, 1, 1, 3, 2, 1, 2 5 5 5 6 5 7 7 7 8 8];
YData = [3, 3, 5, 0, 4, 0, 1, 5, 2, 3 0 4 2 3 3 0 3 2 1 2];

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
eta_outcome=g.MatrixOutput(ci(:,:,n+1))

% --------
%Step 3
disp("Step 3");
[SC, SP]=g.findingSCC(G,n)
[Vsour, Vsink, Visol, Vmid] = g.indexSCC(G,n)

% --------
%Step 4
disp("Step 4");
[indexset,thetaRdd, nuRdd] = g.SCCReducedStructure(G, n);

zeta_reduced=g.MatrixOutput(thetaRdd(:,:,end))
zeta_depth_reduced=g.MatrixOutput(nuRdd(:,:,end))

indexset



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


figure(2)
Rddh = plot(Grdd,'XData', xDataSs, 'YData', yDataSs, 'NodeLabel', indexset);
highlight(Rddh, 1, 'NodeColor', 'green');
highlight(Rddh, 3, 'NodeColor', 'green');
highlight(Rddh, 8, 'NodeColor', 'green');
highlight(Rddh, 10, 'NodeColor', 'red');
saveas(gcf, 'red3_plot.png');

