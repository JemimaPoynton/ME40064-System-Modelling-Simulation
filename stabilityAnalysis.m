%% Define Material Coefficients
clear
D = 1; lambda = 0; f = 0; el = 20;
neumann = BC([0 0], [0 1], 0, 'neumann');
dirichlet = BC([1 1], [0 1], 0, 'dirichlet');

BCs.dirichlet = dirichlet;
BCs.neumann = neumann;

%% Initialise Mesh
x_length = 2;
mesh = OneDimLinearMeshGen(0,x_length,el);
mesh.DVec = D*ones(1,mesh.ngn);

%% Initialise Time Integration Scheme

%% Iterate through different time steps
bound = [0 0];
non_osc = [0 0];
x = [4 5 6 8 10 12 14 20];
n = [25:10:90 100:50:200 300:400:1000];

for k = 1:length(n)
    N = n(k);
    solver = solverScheme(0, 4, N);
    for i = 1:length(x)
        el = x(i);
        mesh = OneDimLinearMeshGen(0, x_length, el);
        mesh.DVec = D*ones(1,mesh.ngn);

        c = solveTransientDiffReact_quadraticBasis(mesh, D, lambda, f, BCs, solver);
        c = c(:,1:2:end);
        t = (0:1:solver.N)*solver.dt;

        [osc, LB, UB, avg, tavg, avg_cont] = checkOsc(c(:,end-1), 0.0005, t);

        if osc == true
            bound = [bound; solver.dt x_length/el];
        else
            non_osc = [non_osc; solver.dt x_length/el];
        end
    end
end

%% Plot
figure()
scatter(bound(:,2), bound(:,1),'x')
hold on
scatter(non_osc(:,2), non_osc(:,1),'o')
dx = 1./(2:0.1:40);
dt = 0.05*dx.^2;

ylabel('dt')
xlabel('dx')

plot(dx, dt)