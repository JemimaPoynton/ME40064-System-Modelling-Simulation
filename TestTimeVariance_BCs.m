%% Define Material Coefficients
clear
D = 1; lambda = 0; f = 0; el = 20;
neumann = BC([1 0], [0 1], 0, 'neumann');
dirichlet = BC([0 1], [0 1], 0, 'dirichlet');

BCs.dirichlet = dirichlet;
BCs.neumann = neumann;

%% Initialise Mesh
mesh = OneDimLinearMeshGen(0,1,el);
mesh.DVec = D*ones(1,mesh.ngn);

%% Initialise Time Integration Scheme
solverScheme = solverScheme(1/2, 2, 350);

%% Setup time varying Neumann BCs
num = solverScheme.N; % number of time steps in input data (can handle ~= N)
BC_n = linspace(0, 4, num);
BC_n(1,num/2 + 35:end) = BC_n(1,num/2 + 35);
t = linspace(0, solverScheme.simTime, num)';

BCs.neumann = setBC(BCs.neumann, 0, BC_n, t); 

%% Setup time varying Dirichlet BCs
num = solverScheme.N; % number of time steps in input data
BC_d = linspace(0, 2, num);
BC_d(1,num/2 -35:end) = BC_d(1,num/2 -35);
t = linspace(0, solverScheme.simTime, num)'; % Time of points for varying Neumann BCs

BCs.dirichlet = setBC(BCs.dirichlet, 1, BC_d, t); 

%% Introduce some noise
% Uncomment this to see response under noise (stable)
% noise = 0.8 + (1.2-0.8).*rand(1,num);
% BCs.dirichlet.val(:,2) = BC_d'.*noise';
% BCs.neumann.val(:,1) = BC_n'.*noise';

%% Solve
c_num_qb = solveTransientDiffReact_qb(mesh, lambda, f, BCs, solverScheme);

%% Plot Against Analytical Solution
x = mesh.nvec; % Select and arbritrary point on the mesh
t = (0:1:solverScheme.N)*solverScheme.dt;

for n = 1:length(x) % !CAN THIS BE COMPRESSED INTO A VECTOR?!
    for i = 1:length(t)
        c_an(i,n) = TransientAnalyticSoln(x(n),t(i)); 
    end
end

%% Choose indices
analysis_t = [0.2 0.6 0.8 1.2 2];
analysis_x_0 = 0; % This will find roughly the closest
analysis_x_1 = 1;

x_qb = linspace(0,1,2*el + 1);
x_qb = x_qb(1:2:end); % ignore centre xi points

for i = 1:length(analysis_t)
    [~, tIndex(i)] = min(abs(t-analysis_t(i)));
end

[~, xIndex_lb_0] = min(abs(x-analysis_x_0));
[~, xIndex_qb_0] =  min(abs(x_qb - analysis_x_0));
[~, xIndex_qb_1] =  min(abs(x_qb - analysis_x_1));

%% Plot against x
figure()
plot(x_qb, c_num_qb(tIndex,1:2:end))
grid on

ylabel('c')
xlabel('x')
legend(['t = ' num2str(analysis_t(1))], ['t = ' num2str(analysis_t(2))], ['t = ' num2str(analysis_t(3))], ['t = ' num2str(analysis_t(4))], ['t = ' num2str(analysis_t(5))])

%% Plot against t
figure()
subplot(2,1,1)
hold on
c_num_qb_pl = c_num_qb(:,1:2:end); % restrict to only end of element solutions
plot(t, c_num_qb_pl(:,xIndex_qb_1), '--black')
hold on
c_num_qb_pl = c_num_qb(:,1:2:end); % restrict to only end of element solutions
plot(t, c_num_qb_pl(:,xIndex_qb_0), '-black')

grid on

legend('Numerical Solution (x = 0)', 'Numerical Solution (x = 1)')
ylabel('C')
xlabel('t')

subplot(2,1,2)
plot(BCs.neumann.t(:,1), BCs.neumann.val(:,1), '--black')
hold on
plot(BCs.dirichlet.t(:,end), BCs.dirichlet.val(:,2), '-black')
grid on
ylabel('BC')
xlabel('t')
legend('Neumann BC (x = 0)', 'Dirichlet BC (x = 1)')
