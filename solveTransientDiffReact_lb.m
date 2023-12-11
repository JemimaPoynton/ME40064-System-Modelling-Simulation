function [c] = solveTransientDiffReact_lb(mesh, lambda, f, BCs, solver)
% Function solveTransientDiffReact_lb solves the transient
% diffusion-reaction equation for Part 1 of the brief. Uses linear 
% basis functions.
%
% mesh:   mesh of elements through x (depth into the skin)
% lambda: Reaction coefficient
% f:      Source term
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

%% Check for varying parameter
% Check outside loop to reduce if statements, format f if constant

f_size = ones(2, length(f));
f = f_size.*f; % Lines 19 to 18 add a time row set to one if time data is missing (f is assumed constant)

%% Define solution variable vectors
Ccurrent = zeros(mesh.ngn,1); % Setup initial conditions
Cnext = zeros(mesh.ngn,1);
c(1, 1:mesh.ngn) = Cnext;

%% Solve
for tstep = 2:solver.N+1 % step through starting at 2 to account for t=0 point
    timerVal = tic;
    t = (tstep-1)*dt; % For handling time varying terms

    [~, fidx_curr] = min(abs(f(1,:)-t)); % Find closest t to handle the event in which f is not discretised with the same number of points e.g. experimental data
    [~, fidx_next] = min(abs(f(1,:)-(t + dt)));

    M = massMat_global(mesh); % Mass matrix
    K = stiffMat_global(mesh, lambda); % Stiffness matrix

    globalMat = M + theta*dt*K; % Create global matrix
    globalVec = (M - (1 - theta)*dt*K)*Ccurrent; % Create global vector

    disp(['Time step: ' num2str(tstep-1) ', Time Taken: ' num2str(toc(timerVal))]);
    disp(' ')

    for eID = 1:mesh.ne      
        Fcurrent = sourceLEV(f(2,fidx_curr), eID, mesh); % calculate LEV for current f
        Fnext = sourceLEV(f(2,fidx_next), eID, mesh); % calculate LEV for next f

        source = dt*(theta*Fnext + (1 - theta)*Fcurrent); % Source terms matrix
        globalVec(eID:eID+1) = source + globalVec(eID:eID+1); % Add source terms to global vector  
    end
 
    [globalMat, globalVec] = applyBCs(globalMat, globalVec, BCs, dt, theta, t, mesh);
    
    % Store solutions
    Cnext = globalMat\globalVec;
    c(tstep, 1:mesh.ngn) = Cnext;

    Ccurrent = Cnext;
end

disp(['Total Simulation Time: ' num2str(toc(overallSolTime))]);
disp(' ')