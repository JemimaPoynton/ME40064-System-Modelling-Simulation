%% Define Material Coefficients
clear
D = 1; lambda = 0; f = 0; el = 32;
neumann = BC([0 0], [0 1], 0, 'neumann');
dirichlet = BC([1 1], [0 1], 0, 'dirichlet');

BCs.dirichlet = dirichlet;
BCs.neumann = neumann;

%% Initialise Mesh
mesh = OneDimLinearMeshGen(0,1,el);

%% Initialise Time Integration Scheme
solver = solverScheme(1/2, 2, 1350);

%% Create a varying field
mesh.DVec(:,1:15) = linspace(1,3,15); % Set a varying D field
mesh.DVec(:,15:25) = 6;
mesh.DVec(:, 25:mesh.ngn) = 1;

%% Calculate
c =  solveTransientDiffReact_lb(mesh, lambda, f, BCs, solver);

%% Plot Against Analytical Solution
x = mesh.nvec; % Select and arbritrary point on the mesh
t = (0:1:solver.N)*solver.dt;

for n = 1:length(x) % !CAN THIS BE COMPRESSED INTO A VECTOR?!
    for i = 1:length(t)
        c_an(i,n) = TransientAnalyticSoln(x(n),t(i)); 
    end
end

%% Choose indices
analysis_t = 0.5;
analysis_x = 0.25; % This will find roughly the closest

[~, tIndex] = min(abs(t-analysis_t)); % only works if an existing point is specified
[~, xIndex_lb] = min(abs(x-analysis_x));

%% Plot against x
f1 = figure();
f1s1 = subplot(2,1,1);
plot(x, c_an(tIndex,:),'--black') %% error is in the intial conditions
grid on
hold on
plot(x, c(tIndex,:),'-black')

ylabel('c')
xlabel('x')
legend('D = 1 (constant)', 'Numerical Solution (Quadratic Basis)')

f1s2 = subplot(2,1,2);
plot(x, mesh.DVec, '-black')
grid on
ylabel('D')
xlabel('x')

% Formatting
f1s1.Position = [0.08 0.50 0.875 0.47];
f1s2.Position = [0.08 0.11 0.875 0.25];

