addpath('..')
addpath(genpath('../utilities/'))
addpath(genpath('../datasets/'))
addpath(genpath('../m_map/'))
addpath(genpath('../mesh2d/'))
addpath(genpath('../namo/'))

dom1 = [-100, -95; 23 26];
dom2 = [-85, -82; 25 32];
missipi = [-92 -88.7; 28.5 31.2 ];
bay = [-95.5, -94.5; 29.0, 30.25];
bay2 = [-95.1, -94; 29.25 30];
levees = [-94, -88; 26 32];

%% Extract high-resolution insert (2008 mesh)
% Load High-res mesh % Call it obj1
obj1 = m2008; 
shp = shaperead('~/Documents/cutting2.shp');
%load('Model_120m_Combinedv58_Fixed_Run.mat')
idx = find([shp.id]==2);
ins1 = [[shp(idx).X]',[shp(idx).Y]'];
%[m_ins,ind,~] = extract_subdomain(obj1,ins1);
m_ins = ExtractSubDomain(obj1,ins1);
%m_ins = map_mesh_properties(m_ins,'ind',ind);

%% Clip the outer mesh (120m)
% Load (Coase resolution mesh) % Call it as obj2
obj2 = m30_cut; 
idx = find([shp.id]==1);
ins2 = [[shp(idx).X]',[shp(idx).Y]'];
m_coar = ExtractSubDomain(obj2, ins2);
%[m_coar,ind,~] = extract_subdomain(obj2,ins2);
%m_coar = map_mesh_properties(m_coar,'ind',ind);

%% Extract the nodes along the boundaries
m_ins.op= []; 
m_ins = makens(m_ins,'outer',0);
bd1 = m_ins.p(m_ins.op.nbdv,:);

m_coar.op=[];
m_coar = makens(m_coar,'outer',0);
bd2 = m_coar.p(m_coar.op.nbdv,:);
figure
%plot(bd1(:,1),bd1(:,2),bd2(:,1),bd2(:,2));

%% Define the element 

opts.iter = 200;
opts.kind = 'delfront';
opts.ref1 = 'preserve';
% namo
opts.rho2 = 1;
%bd = [bd1;flipud(bd2);bd1(1,:)];
% THIS IS CRUCIAL !!!
% For a polygon with hole, need to delimit by a row of NaN
bd = [bd1; [nan, nan]; bd2];
[node,edge] = getnan2(bd);
%figure;
%plot(bd(:,1),bd(:,2));

% use local feature size with grade = olfs.dhdx
olfs.dhdx = 0.30;
[vlfs,tlfs, hlfs] = lfshfn2(node,edge,[],olfs) ;

[pc1] = obj1.baryc; [pc2] = obj2.baryc;
tq1 = gettrimeshquan( obj1.p, obj1.t);
tq2 = gettrimeshquan( obj2.p, obj2.t);

bar_length1 = mean(tq1.ds,2);
bar_length2 = mean(tq2.ds,2);
opts.iter = 200;
opts.kind = 'delfront';
opts.ref1 = 'preserve';
% namo
opts.rho2 = 1;
%bd = [bd1;flipud(bd2);bd1(1,:)];
% THIS IS CRUCIAL !!!
% For a polygon with hole, need to delimit by a row of NaN
bd = [bd1; [nan, nan]; bd2];
[node,edge] = getnan2(bd);
%figure;
%plot(bd(:,1),bd(:,2));

% use local feature size with grade = olfs.dhdx
olfs.dhdx = 0.30;
[vlfs,tlfs, hlfs] = lfshfn2(node,edge,[],olfs) ;

[pc1] = obj1.baryc; [pc2] = obj2.baryc;
tq1 = gettrimeshquan( obj1.p, obj1.t);
tq2 = gettrimeshquan( obj2.p, obj2.t);

bar_length1 = mean(tq1.ds,2);
bar_length2 = mean(tq2.ds,2);

[Idx1, Dist1] = knnsearch(pc1,vlfs);
[Idx2, Dist2] = knnsearch(pc2,vlfs);

hlfs = mean([bar_length1(Idx1),bar_length2(Idx2)],2);
[slfs] = idxtri2(vlfs,tlfs) ;
hfun = @trihfn2;
% make the delaunay refinement mesh
[vert,etri,tria,tnum] = refine2(node,edge,[],opts,hfun,vlfs,tlfs,slfs,hlfs);
% smooth mesh to improve quality
[p,etri,t,tnum] = smooth2(vert,etri,tria,tnum);

% put into oceanmesh msh class format and assign projection
mfp = msh(); mfp.p = p; mfp.t = t;
[~,mfp] = setProj(mfp,1,'lam',1);
[Idx1, Dist1] = knnsearch(pc1,vlfs);
[Idx2, Dist2] = knnsearch(pc2,vlfs);

hlfs = mean([bar_length1(Idx1),bar_length2(Idx2)],2);
[slfs] = idxtri2(vlfs,tlfs) ;
hfun = @trihfn2;
% make the delaunay refinement mesh
[vert,etri,tria,tnum] = refine2(node,edge,[],opts,hfun,vlfs,tlfs,slfs,hlfs);
% smooth mesh to improve quality
[p,etri,t,tnum] = smooth2(vert,etri,tria,tnum);

% put into oceanmesh msh class format and assign projection
mfp = msh(); mfp.p = p; mfp.t = t;
[~,mfp] = setProj(mfp,1,'lam',1);

%% Add bathymetry for the buffer zone mesh
dem = 'galveston_13_mhw_2007.nc';
%mfp = interp(mfp,'../GEBCO_2021.nc');
mfp = interp(mfp,dem);
mfp = lim_bathy_slope(mfp,0.1,-1);

%% Join everything
 m_new = cat(m_ins, mfp);
 m_new = catBathy(m_new, m_ins, mfp);
 m_fin = cat(m_new, m_coar);
 m_fin = catBathy(m_fin, m_new, m_coar );

%% Clean the mesh
 m_fin.bd = []; m_fin.op= []; m_fin.f13 = []; m_fin.f15 = []; m_fin.f24=[];
if ~isempty(find(isnan(m_fin.b),1))
    error('Check bathymetry. There are NaN values')
end
m_fin = Make_Mesh_Boundaries_Traversable(m_fin, 0);
 %% Add islands boundary conditions

m_fin = makens(m_fin,'islands',0);
%m = makens(m,'auto',gdat{1},60);

%% Save a backup before adding the levees
write(m_fin,'ship_cut_v3','f14');

%% Add levees boundary conditions
tic
%mfp_file = 'Model_120m_Floodplainv5.mat';
%ocean_file = 'Model_120m_Oceanv6_Tid.mat';
%centerlines_file = 'centerlines_res200m_snap300m_new2';
mfp_file = 'Model_30m_Floodplain_Clipped3_Cleaned.mat';
ocean_file = 'Model_30m_Ocean_Tidv6.mat';
centerlines_file = 'usace_survey_centerline_matlab2';
m_fin = Levees2Islands(m_fin,mfp_file,ocean_file,centerlines_file,70);
toc
m_fin = renum(m_fin);
disp('Finished adding levees.')
%% Save backup
write(m_fin,'ship_cut_v3','f14');
%% Add open boundary conditions

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

%% Add tides (Ike)
name = 'ship_cut_v3';
m_fin = Make_f15(m_fin, '05-Sep-2008 12:00', '14-Sep-2008 06:00', 2, 'const', 'major8','tidal_database','h_tpxo9.v1.nc');
write(m_fin,name,'f14');
write(m_fin,name,'f15');
