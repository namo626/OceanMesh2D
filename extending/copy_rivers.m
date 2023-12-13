%% Load matthew's mesh with 45 rivers
matt = msh('fort.14.45_rivers');
mv7 = msh('30m_cut_v7.14');
%% Find the river boundary segments (IBTYPE = 22)
nvell = matt.bd.nvell;
ibtype = matt.bd.ibtype;
nbvv = matt.bd.nbvv;

river_segs_num = nvell(ibtype==22);
river_nodes = nbvv(:,ibtype==22);

%% Map the nodes from the old mesh to the new mesh
% Get the x,y of each river node
% river_x = river_nodes;
% river_y = river_nodes;
[ii,jj,ss] = find(river_nodes);
ii = int64(ii);
jj = int64(jj);

river_new = river_nodes;

for k = 1:length(ii)
    node = int64(full(river_nodes(ii(k),jj(k))));
    %river_x(ii(k),jj(k)) = matt.p(node,1);
    %river_y(ii(k),jj(k)) = matt.p(node,2);
    x = matt.p(node,1);
    y = matt.p(node,2);
    river_new(ii(k),jj(k)) = find(mv7.p(:,1) == x & mv7.p(:,2)==y, 1);
    disp(k)
end

%% Put the new rivers into the new mesh
%nbvv2 = [mv7.bd.nbvv river_new];
nbvv2 = sparse(8303,45);
nbvv2(1:size(river_new,1),:) = river_new;
nbvv3 = [mv7.bd.nbvv nbvv2];

% Add this to NBVV
mv7.bd.nbvv = nbvv3;
% Add the count in ibtype and nvell
mv7.bd.ibtype = [mv7.bd.ibtype ibtype(ibtype==22)];
mv7.bd.nbou = mv7.bd.nbou + 45;
mv7.bd.nvell = [mv7.bd.nvell river_segs_num];
mv7.bd.nvel = mv7.bd.nvel + nnz(river_new);

% Pad the rest with zeros...? I don't think so!

%% Write to file
write(mv7, '30m_cut_v8', 'f14')
