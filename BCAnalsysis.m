%% Define Material Coefficients
clear

% Define a structure full of BCs for plotting
pos_n = [0 1; 0 1; 1 1; 0 0];
pos_d = [1 0; 1 0; 0 0; 1 1];

val_n = [0 0; 0 -3000; -3000 -1500; 0 0];
val_d = [30 30; 30 0; 0 0; 30 0];

for i = 1:length(pos_n)
    neumann = BC(pos_n(i,:), val_n(i,:), 0, 'neumann');
    dirichlet = BC(pos_d(i,:), val_d(i,:), 0, 'dirichlet');
    
    BCs(i).dirichlet = dirichlet; BCs(i).neumann = neumann;
end

%% Initialise Mesh
E = 0.00166667; % Epidermis UB
De = 0.005; % Dermin UB
B = 0.01; % Sub cutaneous UB

layers = [E, De, B]; % define different mesh layers
elSizes = [20, 20, 20]; % define number of elements in each layer

mesh = distributedMesh(layers, elSizes);

%% Initialise Time Integration Scheme
solver = solverScheme(1/2, 30, 400);

%% Create a varying field
[mesh, ~, ~, ~] = skinProperties(mesh, layers, [4e-6 5e-6 2e-6], [0 0.01 0.01], [0.02 0.02 0.02]);

%% Choose indices
t = (0:1:solver.N)*solver.dt; % create time vector
x = mesh.nvec;

analysis_t = [4 7 15 30]; % find indexes of these times for analysis
[xIndex, tIndex] = chooseIndices(t, x, analysis_t, 0.005);

%% Plot
figure()
hold on
grid on
lineStyle = ["black-"; "black--"; "black:"; "black-."];

for i = 1:length(BCs)
    c =  solveTransientDiffReact_qb_part2(mesh, BCs(i), solver);
    c = c(:,1:2:end);

    subplot(2,2,i)
    for n = 1:length(tIndex)
        plot(x, c(tIndex(n),:), lineStyle(n))
        hold on
    end
    ylabel('c')
    xlabel('x')
    legend(['t = ' num2str(analysis_t(1))], ['t = ' num2str(analysis_t(2))], ['t = ' num2str(analysis_t(3))], ['t = ' num2str(analysis_t(4))])
    grid on
end
