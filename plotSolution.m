%% Define Material Coefficients
clear
D = 1; lambda = 0; f = 0; el = 10;
BCs = struct(); % Defining BCs
BCs.neumann.pos = [0 0]; BCs.neumann.val = [1 0]; BCs.dirichlet.pos = [1 1]; BCs.dirichlet.val = [0 1]; BCs.dirichlet.t = 0; BCs.neumann.t = 0;

%% Initialise Mesh
mesh = OneDimLinearMeshGen(0,1,el);
mesh.DVec = D*ones(1,mesh.ngn);

%% Initialise Time Integration Scheme
solver = solverScheme(1/2, 1, 350);

%% Calculate Numerical Solution
c_num_qb = solveTransientDiffReact_qb(mesh, lambda, f, BCs, solver);

%% Determine Indices
analysis_t = [0.05 0.1 0.2 1];
x_qb = linspace(0,1,2*el + 1);
x_qb = x_qb(1:2:end); % ignore centre xi points

x = mesh.nvec; % Select and arbritrary point on the mesh
t = (0:1:solver.N)*solver.dt;

for i = 1:length(analysis_t)
    [~, tIndex(i)] = min(abs(t-analysis_t(i)));
end

%% Plot
figure()
plot(x_qb, c_num_qb(tIndex,1:2:end))
grid on

ylabel('c')
xlabel('x')
legend(['t = ' num2str(analysis_t(1))], ['t = ' num2str(analysis_t(2))], ['t = ' num2str(analysis_t(3))], ['t = ' num2str(analysis_t(4))])