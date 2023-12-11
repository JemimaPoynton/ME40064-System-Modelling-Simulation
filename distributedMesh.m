function mesh = distributedMesh(layers, elSizes)
% function distributedMesh creates a mesh made of 3 distinct section
% defined by the points in layers, each with different element sizes as
% defined by elSizes
% 
% layers: layers of the skin in the form [E De B]
% elSizes: number of elements in each section in the form [el1 el2 el3]

%% Extract for Readability
E = layers(1); De = layers(2); B = layers(3);
el1 = elSizes(1); el2 = elSizes(2); el3 = elSizes(3);

mesh1 = OneDimLinearMeshGen(0,E,el1);
mesh2 = OneDimLinearMeshGen(E,De,el2);
mesh3 = OneDimLinearMeshGen(De,B,el3);

%% Form one varied density mesh by stitching together 
mesh.ne = mesh1.ne + mesh2.ne + mesh3.ne;
mesh.ngn = mesh.ne + 1;
mesh.nvec = [mesh1.nvec mesh2.nvec(2:end) mesh3.nvec(2:end)];

for i = 1:mesh1.ne
    mesh.elem(i).x = mesh1.elem(i).x;
    mesh.elem(i).n = mesh1.elem(i).n;
end

for i = 1:mesh2.ne
    mesh.elem(i + mesh1.ne).x = mesh2.elem(i).x;
    mesh.elem(i + mesh1.ne).n = mesh2.elem(i).n + mesh1.ne;
end

for i = 1:mesh3.ne
    mesh.elem(i + mesh2.ne + mesh1.ne).x = mesh3.elem(i).x;
    mesh.elem(i + mesh2.ne + mesh1.ne).n = mesh3.elem(i).n + mesh1.ne + mesh2.ne;
end

for i = 1:mesh.ne
    mesh.elem(i).J = 0.5*(mesh.elem(i).x(2) - mesh.elem(i).x(1));
end