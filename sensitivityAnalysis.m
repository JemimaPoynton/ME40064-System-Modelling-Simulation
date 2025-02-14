%% Define Material Coefficients
clear
D = 4e-6; betaVal = 0.01; gammaVal = 0.02; el = 36; %% 18 elements to represent e distance

neumann = BC([0 1], [0 0], 0, 'neumann');
dirichlet = BC([1 0], [90 0], 0, 'dirichlet');

BCs.dirichlet = dirichlet;
BCs.neumann = neumann;

%% Initialise Mesh
E = 0.00166667; % Epidermis UB
De = 0.005; % Dermin UB
B = 0.01; % Sub cutaneous UB

layers = [E, De, B];
elSizes = [10, 10, 10];

mesh = distributedMesh(layers, elSizes);

%% Initialise Time Integration Scheme
solver = solverScheme(1/2, 30, 400);

%% Create a varying field
[mesh, E_idx, De_idx, B_idx] = skinProperties(mesh, layers, [4e-6 5e-6 2e-6], [0 0.01 0.01], [0.02 0.02 0.02]);

%% Choose indices
t = (0:1:solver.N)*solver.dt;
x = mesh.nvec;

analysis_t = 30;
[~, tIndex] = min(abs(t-analysis_t));

%% Vary Coefficients
% vary D, gamma and beta
percentVar = [1.9 1.6, 1.3, 1, 0.7, 0.4 0.1];
len = length(percentVar);
figure()
[~, t30] = min(abs(t-30)); % find index of t = 30

subplot(3,1,1); hold on; grid on

for i = 1:len
    [mesh, E_idx, De_idx, B_idx] = skinProperties(mesh, layers, percentVar(i).*[4e-6 5e-6 2e-6], [0 0.01 0.01], [0.02 0.02 0.02]);

    c =  solveTransientDiffReact_qb_part2(mesh, BCs, solver);
    c = c(:,1:2:end);

    plot(x, c(tIndex,:))
    
    teff = find(c(:,De_idx) > 40, 1);
    K1(i) = trapz(c(teff:t30, De_idx))*solver.dt;
end
legend('160%','130%','100%','70%','40%')

subplot(3,1,2); hold on; grid on

for i = 1:len
    [mesh, E_idx, De_idx, B_idx] = skinProperties(mesh, layers, [4e-6 5e-6 2e-6], percentVar(i).*[0 0.01 0.01], [0.02 0.02 0.02]);

    c =  solveTransientDiffReact_qb_part2(mesh, BCs, solver);
    c = c(:,1:2:end);

    plot(x, c(tIndex,:))
    
    teff = find(c(:,De_idx) > 40, 1);
    K2(i) = trapz(c(teff:t30, De_idx))*solver.dt;
end

subplot(3,1,3); hold on; grid on

for i = 1:len
    [mesh, E_idx, De_idx, B_idx] = skinProperties(mesh, layers, [4e-6 5e-6 2e-6], [0 0.01 0.01], percentVar(i).*[0.02 0.02 0.02]);

    c =  solveTransientDiffReact_qb_part2(mesh, BCs, solver);
    c = c(:,1:2:end);

    plot(x, c(tIndex,:))

    teff = find(c(:,De_idx) > 40, 1);
    K3(i) = trapz(c(teff:t30, De_idx))*solver.dt;
end

%% Plot
% Plot effective delivered dose against variation of each property
figure()
plot((percentVar-1)*100, K1, '-.black')
hold on
plot((percentVar-1)*100, K2, '--black')
plot((percentVar-1)*100, K3, '-black')
grid on

legend('D Variation', 'beta Variation', 'gamma Variation')
xlabel('% Variation')
ylabel('Effective Delivered Dose (K)')

%% Investigate if D has an influence on the gradient of beta and gamma
Dvar = linspace(2,0, 50);
% percentVar = [1.6, 1, 0.4];
% len = length(percentVar);
figure()
grid on
hold on

for n = 1:length(Dvar)    
    for i = 1:len
        [mesh, E_idx, De_idx, B_idx] = skinProperties(mesh, layers, [4e-6 5e-6 2e-6]*Dvar(n), percentVar(i).*[0 0.01 0.01], [0.02 0.02 0.02]);
    
        c =  solveTransientDiffReact_qb_part2(mesh, BCs, solver);
        c = c(:,1:2:end);
        
        teff = find(c(:,De_idx) > 40, 1);
        K_D2(i,n) = trapz(c(teff:t30, De_idx))*solver.dt;
    end

    gradientbeta(n) = (K_D2(1,n) + K_D2(end,n))/(percentVar(1) - percentVar(end));
    
    for i = 1:len
        [mesh, E_idx, De_idx, B_idx] = skinProperties(mesh, layers, [4e-6 5e-6 2e-6]*Dvar(n), [0 0.01 0.01], percentVar(i).*[0.02 0.02 0.02]);
    
        c =  solveTransientDiffReact_qb_part2(mesh, BCs, solver);
        c = c(:,1:2:end);
    
        teff = find(c(:,De_idx) > 40, 1);
        K_D3(i,n) = trapz(c(teff:t30, De_idx))*solver.dt;
    end

    gradientgamma(n) = (K_D3(1,n) + K_D3(end,n))/(percentVar(1) - percentVar(end));
end

%% Plot
plot((percentVar-1)*100, K_D2, '--black')
plot((percentVar-1)*100, K_D3, '-black')

legend('D Variation', 'beta Variation', 'gamma Variation')
xlabel('% Variation')
ylabel('Effective Delivered Dose (K)')

%% Plot
figure()
plot(Dvar, gradientbeta, '-black')
hold on
grid on
plot(Dvar, gradientgamma, '--black')
ylabel('Gradient dK/dF')
xlabel('F_d')
legend('beta Gradient','gamma Gradient')