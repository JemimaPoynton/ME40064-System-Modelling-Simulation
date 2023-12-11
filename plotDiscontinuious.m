function plotDiscontinuious(E_idx, De_idx, B_idx, meshVec, x)
% function plotDiscontinuous adds a dashed line between the different
% layers for the property defined by meshVec.

plot(x(1:E_idx), meshVec(1:E_idx), '-black')
hold on

if x(De_idx) ~= x(E_idx) % plot dashed line to show discontinuity
    plot([x(E_idx), x(E_idx)], [min(meshVec(1:E_idx)) max(meshVec(1:De_idx))], '--black')
end

plot(x(E_idx:De_idx), [meshVec(E_idx+1) meshVec(E_idx+1:De_idx) ] , '-black')

if x(De_idx) ~= x(B_idx) % plot dashed line to show discontinuity
    plot([x(De_idx), x(De_idx)], [min(meshVec(E_idx+1:De_idx)) max(meshVec(De_idx+1:B_idx))], '--black')
end

plot(x(De_idx:end), [meshVec(De_idx+1) meshVec(De_idx+1:end) ] , '-black')