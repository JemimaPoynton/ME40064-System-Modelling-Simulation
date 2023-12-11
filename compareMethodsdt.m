clear
%% Define Scheme
intScheme = struct();

intScheme.theta = 1/2; % theta = 1/2, 0, 1  for crank-nicholson, forward-euler, backward-euler
intScheme.simTime = 0.4; % How long to run the simulation
intScheme.N = 450; % Number of time steps (rounded)
intScheme.dt = intScheme.simTime/intScheme.N; % Desired time step

N = 5; % GQ SCHEME N

%% Define elements to analyse
n = [450 1000 2000 4000 6000 9000];
L2Norm_val = zeros(length(N),intScheme.N+1);

%% Material parameters and BCs
D = 1; lambda = 0; f = 0;
BCs = struct(); % Defining BCs
BCs.neumann.pos = [0 0]; BCs.neumann.val = [0 0]; BCs.dirichlet.pos = [1 1]; BCs.dirichlet.val = [0 1];
BCs.neumann.t = 0; BCs.dirichlet.t = 0;

%% Check scheme
mesh = OneDimLinearMeshGen(0,1,10);
assert(D*intScheme.dt/(mesh.nvec(2) - mesh.nvec(1))^2 <= 1/2)

%% Loop through different element numbers
for num = 1:length(n)
    intScheme.N = n(num);
    t = (0:1:intScheme.N)*intScheme.dt;
    mesh.DVec = D*ones(1,mesh.ngn);

    intScheme.theta = 1/2;
    c1 = solveTransientDiffReact_quadraticBasis(mesh, D, lambda, f, BCs, intScheme);
    
    intScheme.theta = 1;
    c2 = solveTransientDiffReact_quadraticBasis(mesh, D, lambda, f, BCs, intScheme);
    
    L2Norm_val1(num,1:n(num)+1) = L2Norm_quadratic(N, mesh, t, c1);
    L2Norm_val2(num,1:n(num)+1) = L2Norm(N, mesh, t, c2);
end

%% Plot
% figure()
% loglog(1./nel, L2Norm_val)

figure() % To check correct gradient
plot(n, L2Norm_val1(:,end),'black-')
hold on
grid on

plot(n, L2Norm_val2(:,end),'black-')
hold on
grid on

xlabel('log(dx)')
ylabel('log(L2 Norm)')

legend('Quadratic Basis','Gradient = 3','Linear Basis','Gradient = 2')

