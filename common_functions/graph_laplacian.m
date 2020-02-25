function [D,D3D] = graph_laplacian(faces)
Nv  = max(faces(:));
D   = zeros(Nv,Nv);
D3D = zeros(3*Nv,3*Nv);
for node = 1:Nv
    [row,col]               = find(faces == node);
    index                   = faces(row,:);
    index                   = index(:);
    index                   = setdiff(index,node);
    D(node,node)            = size(index,1);
    D(node,index)           = -1;
    D(index,node)           = -1;
    D3D(3*node,3*node)      = size(index,1);
    D3D(3*node-1,3*node-1)  = size(index,1);
    D3D(3*node-2,3*node-2)  = size(index,1);
    D3D(3*node,3*index)     = -1;
    D3D(3*node-1,3*index-1) = -1;
    D3D(3*node-2,3*index-2) = -1;
    D3D(3*index,3*node)     = -1;
    D3D(3*index-1,3*node-1) = -1;
    D3D(3*index-2,3*node-2) = -1;
end
s    = 0.01;

Deg  = diag(diag(D));
A    = abs(D - Deg);
D    = eye(length(A)) - s*A + (s^2)*(Deg - eye(length(A)));
D    = sparse(D); 

Deg  = diag(diag(D3D));
A    = abs(D3D - Deg);
D3D  = eye(length(A)) - s*A + (s^2)*(Deg - eye(length(A)));
D3D  = sparse(D3D); 
end