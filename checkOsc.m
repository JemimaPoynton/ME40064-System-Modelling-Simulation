function [oscillating, LB, UB, avg, tavg, avg_cont] = checkOsc(c, tol, t)
% function checkOsc evaluates the data defined by c & t for oscillation
% using upper and lower bounds proportional to the timestep, and the
% tolerance, tol. 
%
% Jemima Poynton 12/23

oscillating = false;
% Oscillation is defined by variations of +-tol%

c = c(round(2,0):end); % Ignore first term for averaging
t = t(round(2,0):end);

%% Calculate average between points
for i = 1:length(c)-1
    avg(i) = (c(i) + c(i+1))/2; % Calculate average of midpoints
end

tavg = (t(2) - t(1))/2 + t; % step time into centre point (assumes linear dist)
tavg = tavg(1:end-1); % time associated with the average

%% Interpolate average
avg_cont = interp1(tavg, avg, t(1:end)); % linearly interpolate between the average points

%% Scale tolerance
for i = 2:length(c)-1
    tol_adj(i-1) = 1e3*tol*(i)^-2; % create an adjusted tolerance that scales with timestep
end

%% Set upper and lower bounds and check within range
UB = avg_cont(2:end-1) + tol_adj; % Calculate upper bound before considered oscillating
LB = avg_cont(2:end-1) - tol_adj;

% check if points are within bounds
check = c(2:end-1)' <= LB | c(2:end-1)' >= UB | isnan(c(2:end-1));

if ismember(1,check) % check if any point in c is out of bounds
    oscillating = true;
end