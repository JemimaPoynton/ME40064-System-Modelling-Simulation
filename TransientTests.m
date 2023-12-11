%% Define Material Coefficients
clear
D = 1; lambda = 0; f = 0; el = 12; x_length = 1;

neumann = BC([0 0], [0 -1800], 1, 'neumann');
dirichlet = BC([1 1], [0 1], 1, 'dirichlet');

BCs.dirichlet = dirichlet;
BCs.neumann = neumann;

%% Initialise Mesh
mesh = OneDimLinearMeshGen(0,x_length,el);
mesh.DVec = D*ones(1,mesh.ngn);

%% Initialise Time Integration Scheme
solver = solverScheme(1/2,2,650);

%% Check that Solution Conditions are appropriate
% assert(D*intScheme.dt/(mesh.nvec(2) - mesh.nvec(1))^2 <= 1/2)
dt = (1/D)*0.5*(mesh.nvec(2) - mesh.nvec(1))^2

%% Calculate Numerical Solution
c_num_qb = solveTransientDiffReact_qb(mesh, lambda, f, BCs, solver);
c_num_lb = solveTransientDiffReact_lb(mesh, lambda, f, BCs, solver);
disp(D*solver.dt/(mesh.nvec(2) - mesh.nvec(1))^2)

%% Plot Against Analytical Solution
x = mesh.nvec; % Select and arbritrary point on the mesh
t = (0:1:solver.N)*solver.dt;

for n = 1:length(x) % !CAN THIS BE COMPRESSED INTO A VECTOR?!
    for i = 1:length(t)
        c_an(i,n) = TransientAnalyticSoln(x(n),t(i)); 
    end
end

%% Compare results for Each
analysis_t = 0.15;
analysis_x = [0.8 0.4]; % This will find roughly the closest

x_qb = linspace(0,x_length,2*el + 1);
x_qb = x_qb(1:2:end); % ignore centre xi points

[xIndex, tIndex] = chooseIndices(t, x, analysis_t, analysis_x); 

%% Plot against x
figure()
plot(x, [c_an(tIndex,:); c_num_lb(tIndex,:)]) %% error is in the intial conditions
hold on
grid on
plot(x_qb, c_num_qb(tIndex,1:2:end))

ylabel('c')
xlabel('x')
legend('Analytical Solution', 'Numerical Solution (Linear Basis)', 'Numerical Solution (Quadratic Basis)')

%% Plot against t
figure()
plot(t, [c_an(:,xIndex)'; c_num_lb(:,xIndex)'])
ylim([0 1])
hold on
grid on

c_num_qb_pl = c_num_qb(:,1:2:end); % restrict to only end of element solutions
plot(t, c_num_qb_pl(:,xIndex),'-black')
legend('Analytical Solution', 'Numerical Solution (Linear Basis)', 'Numerical Solution (Quadratic Basis)')
ylabel('C')
xlabel('t')