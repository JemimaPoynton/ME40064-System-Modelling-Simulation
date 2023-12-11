classdef BC
    % Classdef BC defines a set of boundary conditions at either end of a
    % mesh
    %
    % Jemima Poynton 12/23
    properties
        % pos: boolean for position in the form [start end]
        % val: value of the condition at the [start end] of the mesh 
        %      defined as a gradient for type Neumann.
        % t: time vector associated with the values for time-varying
        %    BCs, defined as 1 for time-invariant
        % type: type label (neumann/dirichlet) for reference

        pos(1,2) {mustBeMember(pos, [0 1]), mustBeInteger(pos)} = [0 0]; 
        val(:,2) {mustBeReal(val), mustBeFinite(val)} = [0 0];
        t(:,1) {mustBeReal(t), mustBeFinite(t)} = 0;
        type = 'undefined';
    end

    methods
        function obj = BC(pos, val, t, type)
            % function BC initialises a set of boundary conditions
            
            obj.pos = pos;
            obj.val = val;
            obj.t = t;
            obj.type = type;
        end

        function obj = setBC(obj, idx, newVal, t)
            % function setBC sets new values of a boundary condition object
            %
            % idx: either zero or one, representing which end of the mesh
            %      [start end] = [0 1] the condition is applied
            % newVal: new BC values
            
            arguments
                obj
                idx(1,1) {mustBeMember(idx, [0 1]), mustBeInteger(idx)}
                newVal(:,1) {mustBeReal(newVal), mustBeFinite(newVal)}
                t(:,1) {mustBeReal(t), mustBeFinite(t)}
            end

            if idx == 0
                obj.val = [newVal obj.val(2).*ones(size(newVal))];
            end

            if idx == 1
               obj.val = [obj.val(1).*ones(size(newVal)) newVal]; 
            end

            obj.t = [t; t(end).*ones(length(obj.val') - length(t'),1)];
        end
    end
end