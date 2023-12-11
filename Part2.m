%% Define Material Coefficients
clear

neumann = BC([0 0], [0 -1800], 1, 'neumann');
dirichlet = BC([1 1], [70.6 0], 1, 'dirichlet');

BCs.dirichlet = dirichlet;
BCs.neumann = neumann;

%% Initialise Mesh
E = 0.00166667; % Epidermis UB
De = 0.005; % Dermin UB
B = 0.01; % Sub cutaneous UB

layers = [E, De, B];
elSizes = [16, 22, 22];

mesh = distributedMesh(layers, elSizes);

%% Initialise Time Integration Scheme
solverScheme = solverScheme(1/2, 30, 4000);

%% Create a varying field
[mesh, E_idx, De_idx, B_idx] = skinProperties(mesh, layers, [4e-6 5e-6 2e-6], [0 0.01 0.01], [0.02 0.02 0.02]);

%% Calculate
c =  solveTransientDiffReact_qb_part2(mesh, BCs, solverScheme);
c = c(:,1:2:end);

%% Choose indices
t = (0:1:solverScheme.N)*solverScheme.dt;
x = mesh.nvec;

analysis_t = [3 4 7 30];
analysis_x = De; % This will find roughly the closest

for i = 1:length(analysis_t)
    [~, tIndex(i)] = min(abs(t-analysis_t(i)));
end

[~, xIndex_lb] = min(abs(x-analysis_x));

%% Plot against x
f1 = figure();
f1s1 = subplot(3,1,1);
lineStyle = ["black-"; "black--"; "black:"; "black-."];
hold on

for i = 1:length(tIndex)
    plot(x, c(tIndex(i),:), lineStyle(i))
end

legend(['t = ' num2str(analysis_t(1))], ['t = ' num2str(analysis_t(2))], ['t = ' num2str(analysis_t(3))], ['t = ' num2str(analysis_t(4))])
grid on

ylabel('c')

f1s2 = subplot(3,1,2);
plotDiscontinuious(E_idx, De_idx, B_idx, mesh.DVec, x)
grid on
ylabel('D')
xlabel('x')
ylim([min(mesh.DVec)-1e-7 max(mesh.DVec)+1e-7])

f1s3 = subplot(3,1,3);
plotDiscontinuious(E_idx, De_idx, B_idx, mesh.betaVec, x)
ylim([min(mesh.betaVec)-0.001 max(mesh.betaVec)+0.001])

grid on
ylabel('beta')

% Formatting
f1s1.Position = [0.08 0.59 0.875 0.38];
f1s2.Position = [0.08 0.13 0.875 0.16];
f1s3.Position = [0.08 0.35 0.875 0.16];

%% Plot against t
figure()
plot(t,  c(:,xIndex_lb)', '-black')
grid on

ylabel('c(x = D,t)')
xlabel('t')


%% Caclulate Integral and Location of Sufficient Concentration 
K = 0;
[~, t30] = min(abs(t-30)); % find index of t = 30

options = optimset('TolX',0.01); 
[cDose, error] = fminbnd(@(cDose)errorFun(cDose, mesh, BCs, solverScheme, De_idx, t30),...
    50, 100, options);