%% Define Material Coefficients
clear
D = 1; lambda = 0; f = 0; el = 8;
BCs = struct(); % Defining BCs
BCs.neumann.pos = [0 0]; BCs.neumann.val = [1 0]; BCs.dirichlet.pos = [1 1]; BCs.dirichlet.val = [0 1]; BCs.dirichlet.t = 0; BCs.neumann.t = 0;
x_length = 1;

%% Initialise Mesh
mesh = OneDimLinearMeshGen(0,x_length,el);
mesh.DVec = D*ones(1,mesh.ngn);

%% Initialise Time Integration Scheme
solver = solverScheme(1,1,400);

%% Plot Against Analytical Solution
x = mesh.nvec; % Select and arbritrary point on the mesh
t = (0:1:solver.N)*solver.dt;

for n = 1:length(x) % !CAN THIS BE COMPRESSED INTO A VECTOR?!
    for i = 1:length(t)
        c_an(i,n) = TransientAnalyticSoln(x(n),t(i)); 
    end
end

%% Run crank nicholson and forward euler for the same scheme, compare accuracy and total run time
solver.theta = 1/2;
c1 = solveTransientDiffReact(mesh, D, lambda, f, BCs, solver);
% c1 = c1(:,1:2:end);

solver.theta = 1;
c2 = solveTransientDiffReact(mesh, D, lambda, f, BCs, solver);
% c2 = c2(:,1:2:end);

solver.theta = 0;
c3 = solveTransientDiffReact(mesh, D, lambda, f, BCs, solver);
% c3 = c3(:,1:2:end);

%% Compare results for Each
analysis_t = 0.95;
analysis_x = 0.75; % This will find roughly the closest

x_qb = linspace(0,x_length,2*el + 1);
x_qb = x_qb(1:2:end); % ignore centre xi points

[~, tIndex] = min(abs(t-analysis_t)); % only works if an existing point is specified
[~, xIndex] = min(abs(x-analysis_x)); 

%% Plot against x
figure()
plot(x_qb, c1(tIndex,:))
hold on
plot(x_qb, c_an(tIndex,:))
grid on
plot(x_qb, c2(tIndex,:))
plot(x_qb, c3(tIndex,:))

ylabel('c')
xlabel('x')
legend('')

%% Plot against t
figure()
plot(t(2:end), c1(2:end,xIndex),'black-')
hold on
plot(t(2:end), c_an(2:end,xIndex),'black:')

plot(t(2:end), c2(2:end,xIndex),'black--')
% 
% plot(t(2:end), c3(2:end,xIndex),'black-.')

grid on

legend('Crank Nicholson','Analytical Solution','Backward Euler','Forward Euler','location','southeast')
ylabel('C')
xlabel('t')

% xlim([0.01 0.1])
