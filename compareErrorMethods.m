clear
%% Define Scheme
intScheme = struct();

intScheme.theta = 1; % theta = 1/2, 0, 1  for crank-nicholson, forward-euler, backward-euler
intScheme.simTime = 1; % How long to run the simulation
intScheme.N = 4400; % Number of time steps (rounded)
intScheme.dt = intScheme.simTime/intScheme.N; % Desired time step

t = (0:1:intScheme.N)*intScheme.dt;
N = 5; % GQ SCHEME N

%% Define elements to analyse
nel = [8 10 12];
L2Norm_val = zeros(length(nel),intScheme.N+1);

%% Material parameters and BCs
D = 1; lambda = 0; f = 0;
BCs = struct(); % Defining BCs
BCs.neumann.pos = [0 0]; BCs.neumann.val = [0 0]; BCs.dirichlet.pos = [1 1]; BCs.dirichlet.val = [0 1];
BCs.neumann.t = 0; BCs.dirichlet.t = 0;

%% Check scheme
mesh = OneDimLinearMeshGen(0,1,nel(end));
assert(D*intScheme.dt/(mesh.nvec(2) - mesh.nvec(1))^2 <= 1/2)

%% Loop through different element numbers
for num = 1:length(nel)
    mesh = OneDimLinearMeshGen(0,1,nel(num));
    mesh.DVec = D*ones(1,mesh.ngn);

    intScheme.theta = 1/2;
    c1 = solveTransientDiffReact_quadraticBasis(mesh, D, lambda, f, BCs, intScheme);
    L2Norm_val1(num,:) = L2Norm_quadratic(N, mesh, t, c1);

    intScheme.theta = 1;
    c2 = solveTransientDiffReact_quadraticBasis(mesh, D, lambda, f, BCs, intScheme);
    L2Norm_val2(num,:) = L2Norm_quadratic(N, mesh, t, c2);

    intScheme.theta = 0;
    c3 = solveTransientDiffReact_quadraticBasis(mesh, D, lambda, f, BCs, intScheme);
    L2Norm_val3(num,:) = L2Norm_quadratic(N, mesh, t, c3);
end

%% Plot
figure()
plot(1./nel, L2Norm_val1(:,end),'black-')
hold on
plot(1./nel, L2Norm_val2(:,end),'black--')
grid on
plot(1./nel, L2Norm_val3(:,end),'black-.')


xlabel('dx')
ylabel('L2 Norm')

legend('Crank Nicholson', 'Backward Euler', 'Forward Euler')