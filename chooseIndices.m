function [xIndex, tIndex] = chooseIndices(t, x, tpoint, xpoint)
% function chooseIndices finds index in t and x associated with the values 
% tpoint (which may take vector form) and x point

tIndex = zeros(size(tpoint));
xIndex = zeros(size(xpoint));

for i = 1:length(tpoint)
    [~, tIndex(i)] = min(abs(t-tpoint(i)));
end

for i = 1:length(xpoint)
    [~, xIndex(i)] = min(abs(x-xpoint(i)));
end