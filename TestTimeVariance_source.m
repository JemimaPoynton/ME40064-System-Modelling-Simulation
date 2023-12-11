%% Define Material Coefficients
clear
D = 1; lambda = 0; el = 10;
neumann = BC([0 0], [0 1], 0, 'neumann');
dirichlet = BC([1 1], [0 1], 0, 'dirichlet');

BCs.dirichlet = dirichlet;
BCs.neumann = neumann;

%% Initialise Mesh
mesh = OneDimLinearMeshGen(0,1,el);
mesh.DVec = D*ones(1,mesh.ngn);

%% Initialise Time Integration Scheme
solver = solverScheme(1/2, 2, 350);

%% Setup time varying field - realistic, varied number of points
num = solver.N;

f = zeros(2,num);

f(1,:) = linspace(0, solver.simTime, num); % time corresponding to each value of f
f(2,:) = linspace(0, 4, num);
f(2,100:200) = 3;
f(2, 200:end) = 0;
noise = 0.4 + (1.2-0.8).*rand(1,num);
f(2,:) = f(2,:).*noise; % Create noise to verify stability

c_num_qb = solveTransientDiffReact_qb(mesh, lambda, f, BCs, solver);

%% Plot Against Analytical Solution
x = mesh.nvec; % Select and arbritrary point on the mesh
t = (0:1:solver.N)*solver.dt;

for n = 1:length(x) 
    for i = 1:length(t)
        c_an(i,n) = TransientAnalyticSoln(x(n),t(i)); 
    end
end

%% Choose indices
analysis_t = 1.15;
analysis_x = 0.25; % This will find roughly the closest

x_qb = linspace(0,1,2*el + 1);
x_qb = x_qb(1:2:end); % ignore centre xi points

[~, tIndex] = min(abs(t-analysis_t)); % only works if an existing point is specified
[~, xIndex_lb] = min(abs(x-analysis_x));
[~, xIndex_qb] = min(abs(x_qb - analysis_x));  

%% Plot against x
figure()
plot(x, c_an(tIndex,:),'--black') %% error is in the intial conditions
hold on
grid on
plot(x_qb, c_num_qb(tIndex,1:2:end), '-black')

ylabel('c')
xlabel('x')
legend('Analytical Solution', 'Numerical Solution (Linear Basis)')

%% Plot against t
f1 = figure();
f1s1 = subplot(2,1,1);
plot(t, c_an(:,xIndex_lb), '--black')
hold on
c_num_qb_pl = c_num_qb(:,1:2:end); % restrict to only end of element solutions
plot(t, c_num_qb_pl(:,xIndex_qb), '-black')
hold on
grid on

legend('f = 0', 'Numerical Solution (Quadratic Basis)', 'Location', 'southeast')
ylabel('C')
xlabel('t')

f1s2 = subplot(2,1,2);
plot(f(1,:), f(2,:), '-black')
grid on
ylabel('f')
xlabel('t')

% Formatting
f1s1.Position = [0.08 0.53 0.875 0.42];
f1s2.Position = [0.08 0.11 0.875 0.30];
