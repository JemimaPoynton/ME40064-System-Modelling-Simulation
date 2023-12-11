function error = errorFun(cDose, mesh, BCs, solver, De_idx, t30)
% function errorFun calculates the objective function defined by error, to
% converge on a value of K of 1000.
%
% cDose: dose applied on skin
% De_idx: index in mesh.nvec of the dermis upper surface
% t30: time index of 30 s

BCs.dirichlet.val = [cDose 0]; % Run with current cDose guess
c =  solveTransientDiffReact_qb_part2(mesh, BCs, solver);
c = c(:,1:2:end);

teff = find(c(:,De_idx) > 40, 1);

K = trapz(c(teff:t30, De_idx))*solver.dt;
error = (1000 - K)^2;