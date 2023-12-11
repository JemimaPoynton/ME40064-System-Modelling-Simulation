%% Set Boundary Conditions
% Define boundary conditions in structure BCs
BCs = struct(); 

BCs.neumann.pos = [0 0]; % boolean representing whether boundary applied at [1st node, nth node]
BCs.neumann.val = [0 1];
% could generalise so BCs can be enforced internally using BCs.neumann.x =
% 2, and feeding in the mesh

BCs.dirichlet.pos = [1 1];
BCs.dirichlet.val = [0 1];

%% Analyse positive, negative reaction term
f = 0; D = 1;
mesh = OneDimLinearMeshGen(0,1,150);
figure()

c1 = solveDiffusionReaction(mesh, D, 0, f, BCs);
plot(mesh.nvec, c1, '-', 'Color', 'black')
hold on
linewidth = 2.5;

for D = [1 2 3]
    linewidth = linewidth - 0.75/D;

    c2 = solveDiffusionReaction(mesh, D, 3, f, BCs);
    c3 = solveDiffusionReaction(mesh, D, -9, f, BCs);

    plot(mesh.nvec, c2, '--', 'Color', 'black', 'LineWidth', linewidth)
    plot(mesh.nvec, c3, '-.', 'Color', 'black', 'LineWidth', linewidth)
end

grid on
ylim([min(c3), 1.03])
ylabel('c')
xlabel('x')

legend('D = 1, lambda = 0', 'D = 1, lambda = 3', 'D = 1, lambda = -9', 'D = 2, lambda = 3', ...
    'D = 2, lambda = -9', 'D = 3, lambda = 3', 'D = 3, lambda = -9', 'Location', 'northwest')

set(gcf, 'Position',  [100, 100, 700, 350])
%% Analyse and plot positive, negative source term
lambda = 0; D = 1;
mesh = OneDimLinearMeshGen(0,1,150);

c1 = solveDiffusionReaction(mesh, D, lambda, 0, BCs);
plot(mesh.nvec, c1, '-', 'Color', 'black')
hold on
linewidth = 2.5;

for D = [1 2 3]
    linewidth = linewidth - 0.75;

    c2 = solveDiffusionReaction(mesh, D, lambda, 9, BCs);
    c3 = solveDiffusionReaction(mesh, D, lambda, -9, BCs);

    plot(mesh.nvec, c2, '--', 'Color', 'black', 'LineWidth', linewidth)
    plot(mesh.nvec, c3, '-.', 'Color', 'black', 'LineWidth', linewidth)
end

grid on
ylim([-0.72, 1.72])
ylabel('c')
xlabel('x')

legend('D = 1, f = 0', 'D = 1, f = 9', 'D = 1, f = -9', 'D = 2, f = 9', ...
    'D = 2, f = -9', 'D = 3, f = 9', 'D = 3, f = -9', 'Location', 'northwest')

set(gcf, 'Position',  [100, 100, 700, 350])