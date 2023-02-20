%% Load matthew's mesh with 45 rivers
matt = msh('fort.14.45_rivers');

%% Find the river boundary segments (IBTYPE = 22)
nvell = matt.bd.nvell;
ibtype = matt.bd.ibtype;
nbvv = matt.bd.nbvv;

river_segs_num = nvell(ibtype==22);
river_nodes = nbvv(:,ibtype==22);

%% Map the nodes from the old mesh to the new mesh
% Get the x,y of each river node
river_x = river_nodes;
river_y = river_nodes;
[ii,jj,ss] = find(river_nodes);
ii = int64(ii);
jj = int64(jj);

for k = 1:length(ii)
    node = int64(full(river_nodes(ii(k),jj(k))));
    river_x(ii(k),jj(k)) = matt.p(node,1);
    river_y(ii(k),jj(k)) = matt.p(node,2);

    
end