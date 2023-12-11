function [c] = solveTransientDiffReact_qb_part2(mesh, BCs, solver)
% Function solveTransientDiffReact_qb_part2 solves the transient
% diffusion-reaction equation as defined for drug delivery in Part 2 of the
% brief. Uses quadratic basis functions.
%
% mesh:   mesh of elements through x (depth into the skin)
% BCs:    object of class 'BC' defining the boundary conditions
% solver: object of class solverScheme defining the time-stepping scheme 
%         used for the numerical solver, including the method, step-size and 
%         simulation time
%
% Jemima Poynton 12/23

%% Initialise Timer
overallSolTime = tic;

%% Extract for Readability
dt = solver.dt;
theta = solver.theta;

%% Define solution variable vectors
Ccurrent = zeros(2*mesh.ngn-1,1); % Setup initial conditions
Cnext = zeros(2*mesh.ngn-1,1);
c(1, 1:2*mesh.ngn-1) = Cnext;

%% Formating t = 0 point for BCs
% setup c with inital condition etc. (for plotting only)
if BCs.dirichlet.pos(2) == 1
    c(1, end) = BCs.dirichlet.val(2);
else 
    c(1, end) = 0;
end

if BCs.dirichlet.pos(1) == 1
    c(1, 1) = BCs.dirichlet.val(1);
else 
    c(1, 1) = 0;
end 

%% Solve
for tstep = 2:solver.N+1 % step through starting at 2 to account for t=0 point
    timerVal = tic;
    t = (tstep-1)*dt; % For handling time varying terms

    M = massMat_global_qb(mesh); % Mass matrix
    K = stiffMat_global_qb_part2(mesh); % Stiffness matrix

    globalMat = M + theta*dt*K; % Create global matrix
    globalVec = (M - (1 - theta)*dt*K)*Ccurrent; % Create global vector 
    
    % display time step and time taken for each loop
    disp(['Time step: ' num2str(tstep-1) ', Time Taken: ' num2str(toc(timerVal))]);
    disp(' ')
 
    [globalMat, globalVec] = applyBCs(globalMat, globalVec, BCs, dt, theta, t, mesh);

    % Store solutions
    Cnext = globalMat\globalVec; % Calculate next solution point
    c(tstep, 1:2*mesh.ngn-1) = Cnext;

    Ccurrent = Cnext;

end

disp(['Total Simulation Time: ' num2str(toc(overallSolTime))]);
disp(' ')