%% Load the mesh
%load('Model_30m_Floodplainv2.mat');
%m = msh('fort.14');
shp = shaperead('cut.shp');

%% Carve out from the floodplain

ins = [[shp.X]', [shp.Y]'];
mfp_ins = ExtractSubDomain(mfp, ins);
mfp_ins = Make_Mesh_Boundaries_Traversable(mfp_ins, 0);

mfp_ins = minus(mfp_ins, m);

%% Insert bathymetry for floodplain
mfp_ins = interp(mfp_ins,'../GEBCO_2021.nc');


%% Fix before cat
m2 = fixmeshandcarry(m);
mfp2 = fixmeshandcarry(mfp_ins);
merged2 = cat(m2, mfp2);
merged2 = catBathy(merged2, m2, mfp2);
merged3 = merged2;

%% Test merging
re = [-97, -95.5; 28.3, 30];
plot(merged3, 'proj','none','subdomain',re);
plot(m, 'proj','none','subdomain',re);


%% Internal edge boundary condition
merged3.bd = [];
merged3 = Make_Mesh_Boundaries_Traversable(merged3, 0);
merged3 = makens(merged3,'islands',0);

m_fin = merged3;

%% Plot the boundary
plot(m_fin, 'proj','none','type','bd');

%% Add levees boundary conditions
tic
%mfp_file = 'Model_120m_Floodplainv5.mat';
%ocean_file = 'Model_120m_Oceanv6_Tid.mat';
%centerlines_file = 'centerlines_res200m_snap300m_new2';
mfp_file = '/workspace/OceanMesh2D/namo/Model_30m_Floodplain_Clipped3_Cleaned.mat';
ocean_file = '/workspace/OceanMesh2D/namo/Model_30m_Ocean_Tidv6.mat';
centerlines_file = 'usace_survey_centerline_matlab2';
m_fin = Levees2Islands(m_fin,mfp_file,ocean_file,centerlines_file,70);
toc
m_fin = renum(m_fin);

% !! IMPORTANT !! Save backup
write(m_fin, 'levee_backup','f14');
disp('Finished adding levees.')

%% Add open boundary conditions
m_fin.op = [];
m_fin = makens(m_fin,'outer',0); % Select the init and end point of the bondary and assign boundary option 2 "Elevation BC"
%m = makens(m,'islands',0);

%% Add mainland boundary condition

bnde = extdom_edges2(m_fin.t,m_fin.p);
del = ismember(bnde(:,1),m_fin.bd.nbvv);
bnde(del,:) = []; 
del = ismember(bnde(:,2),m_fin.bd.nbvv);
bnde(del,:) = []; 
del = ismember(bnde(:,1),m_fin.bd.ibconn);
bnde(del,:) = []; 
del = ismember(bnde(:,2),m_fin.bd.ibconn);
bnde(del,:) = []; 
 

cell2 = extdom_polygon(bnde,m_fin.p,-1,0);
lengths = cell2mat(cellfun(@size,cell2,'UniformOutput',false)'); lengths(:,2) = [];
Idx = find(lengths == max(lengths));
ocean_first = unique(cell2{Idx},'rows','stable');
div1 = knnsearch(ocean_first,m_fin.p(m_fin.op.nbdv(1),:),'K',1);
div2 = knnsearch(ocean_first,m_fin.p(m_fin.op.nbdv(end),:),'K',1);
ocean1 = ocean_first(1:min(div1,div2),:);
ocean2 = ocean_first(max(div1,div2):end,:);
nans = find(isnan(ocean1(:,1))); ocean1(nans,:) = [];
nans = find(isnan(ocean2(:,1))); ocean2(nans,:) = [];

if ocean1(1,2)<ocean1(end,2)
    ocean1 = flipud(ocean1);
end
if ocean2(1,2)<ocean2(end,2)
    ocean2 = flipud(ocean2);
end
if mean(ocean1(:,2))>mean(ocean2(:,2))
    ocean_first = [ocean1;ocean2];
else
    ocean_first = [ocean2;ocean1];
end

Idx = knnsearch(m_fin.p,ocean_first,'K',1);
aux = find(~ismember(Idx,m_fin.op.nbdv));
mainland_ids = Idx(aux);


for i = 1:500:length(mainland_ids)

m_fin.bd.nbou = m_fin.bd.nbou + 1;
m_fin.bd.ibtype(m_fin.bd.nbou) = 20;
mini_mainland = mainland_ids(i:min(i+500,length(mainland_ids)));
m_fin.bd.nbvv(1:length(mini_mainland),m_fin.bd.nbou) = mini_mainland;
m_fin.bd.nvell(m_fin.bd.nbou) = length(mini_mainland);
m_fin.bd.barincfsb(end,m_fin.bd.nbou) = 0;
m_fin.bd.barincfsp(end,m_fin.bd.nbou) = 0;
m_fin.bd.barinht(end,m_fin.bd.nbou) = 0;
m_fin.bd.ibconn(end,m_fin.bd.nbou) = 0;
end
m_fin.bd.nvel = length(find(m_fin.bd.nbvv(:)~=0));

%% Add tides and write to file
m_fin = Make_f15(m_fin, '05-Sep-2008 12:00', '14-Sep-2008 06:00', 2, 'const', 'major8','tidal_database','h_tpxo9.v1.nc');


%% Change manning's n in the floodplain to a constant value
% Load the mesh with wrong fort.13
m_fin = msh('fname','30m_cut_v7.14','aux',{'30m_cut_v7.13'});

% Get polygon of floodplain region and assign 0.1 to those points
%pgon = get_boundary_of_mesh(mfp2);
% Assign 0.1 to points in the floodplain
m_fin13 = assignConstManning(m_fin, mfp2.p, 0.1);

% write to fort.13
write(m_new, '30m_cut_v8', 'f13');


