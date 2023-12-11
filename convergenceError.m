clear
%% Define Scheme
intScheme = struct();

solver = solverScheme(1/2,1,4400);

t = (0:1:solver.N)*solver.dt;
N = 5; % GQ SCHEME N

%% Define elements to analyse
nel = [4 6 8];
L2Norm_val = zeros(length(nel),solver.N+1);

%% Material parameters and BCs
D = 1; lambda = 0; f = 0;
BCs = struct(); % Defining BCs
BCs.neumann.pos = [0 0]; BCs.neumann.val = [0 0]; BCs.dirichlet.pos = [1 1]; BCs.dirichlet.val = [0 1];
BCs.neumann.t = 0; BCs.dirichlet.t = 0;

%% Check scheme
mesh = OneDimLinearMeshGen(0,1,nel(end));
assert(D*solver.dt/(mesh.nvec(2) - mesh.nvec(1))^2 <= 1/2)

%% Loop through different element numbers
for num = 1:length(nel)
    mesh = OneDimLinearMeshGen(0,1,nel(num));
    mesh.DVec = D*ones(1,mesh.ngn);

    c1 = solveTransientDiffReact_qb(mesh, lambda, f, BCs, solver);
    c2 = solveTransientDiffReact_lb(mesh, lambda, f, BCs, solver);
    L2Norm_val1(num,:) = L2Norm_quadratic(N, mesh, t, c1);
    L2Norm_val2(num,:) = L2Norm(N, mesh, t, c2);
end

%% Plot
% figure()
% loglog(1./nel, L2Norm_val)

figure() % To check correct gradient
plot(log(1./nel), log(L2Norm_val1(:,end)),'black-')
hold on
plot(log(1./nel), 3*log(1./nel) - 3*log(1./nel(2)) + log(L2Norm_val1(2,end)), 'black-.','color',[0.55,0.55,0.55]); 
grid on

plot(log(1./nel), log(L2Norm_val2(:,end)), 'black--')
hold on
plot(log(1./nel), 2*log(1./nel) - 2*log(1./nel(2)) + log(L2Norm_val2(2,end)),'black:', 'color', [0.55,0.55,0.55]); 
grid on

xlabel('log(dx)')
ylabel('log(L2 Norm)')

legend('Quadratic Basis','Gradient = 3','Linear Basis','Gradient = 2')

