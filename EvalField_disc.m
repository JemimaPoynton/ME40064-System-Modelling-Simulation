function [val] = EvalField_disc(meshVec, eID)
% Evaluates the property defined by meshVec (e.g. mesh.betaVec) at eID
% using discontinuous basis functions
val = meshVec(eID);
