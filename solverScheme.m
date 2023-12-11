classdef solverScheme
    % Classdef BC defines the time-stepping scheme used for the numerical
    % solver, including the method, step-size and simulation time 
    %
    % Jemima Poynton 12/23
    properties
        % theta:   in the transient solution equation, defined by the method
        %          used - 1/2, 0, or 1  for crank-nicholson, forward-euler, 
        %          and backward-euler, respectively.
        % simTime: total simulation time
        % N:       number of time steps
        % dt:      time step size

        theta
        simTime 
        N 
        dt 
    end

    methods
        function obj = solverScheme(theta, simTime, N)
            % function solverScheme intialises the solver
            arguments 
                theta(1,1) {mustBeMember(theta, [0 0.5 1])} 
                simTime(1,1) {mustBeReal(simTime), mustBeFinite(simTime)}
                N(1,1) {mustBeReal(N), mustBeInteger(N)}
            end

            obj.theta = theta;
            obj.simTime = simTime;
            obj.N = N;
            obj.dt = simTime/N; % calculate dt based on number of steps
        end
    end
end