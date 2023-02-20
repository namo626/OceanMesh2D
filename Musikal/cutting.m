%% 1. Musikal project. Load the mesh
mglobal = msh('fname','GSTOFSv5.3.14','aux',{'GSTOFSv5.3.13'});


%% 2. Extract the cut from shapefile
% Load High-res mesh % Call it obj1
obj1 = mglobal; 
shp = shaperead('cut_STOFS_mesh.shp');
%idx = find([shp.id]==1);
ins1 = [[shp.X]',[shp.Y]'];
%[m_ins,ind,~] = extract_subdomain(obj1,ins1);
m_fin = ExtractSubDomain(obj1,ins1);
%m_ins = map_mesh_properties(m_ins,'ind',ind);

%% 3. Clean the mesh
m_fin.bd = []; m_fin.op= []; m_fin.f15 = []; 
if ~isempty(find(isnan(m_fin.b),1))
    error('Check bathymetry. There are NaN values')
end

m_fin = Make_Mesh_Boundaries_Traversable(m_fin, 0);
 %% 4. Add islands boundary conditions

m_fin = makens(m_fin,'islands',0);

%% 5. Add levees
mfp_file = '../namo/Model_120m_Floodplainv5.mat';
ocean_file = '../namo/Model_120m_Oceanv6_Tid.mat';
centerlines_file = '../namo/centerlines_res200m_snap300m_new2';
m_fin = Levees2Islands(m_fin,mfp_file,ocean_file,centerlines_file,70);
toc
m_fin = renum(m_fin);
disp('Finished adding levees.')

%% 6. Add open boundary conditions

m_fin = makens(m_fin,'outer',0); % Select the init and end point of the bondary and assign boundary option 2 "Elevation BC"
%m = makens(m,'islands',0);



%% 7. Add mainland boundary condition

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

%% 8. Add tides
m_fin = Make_f15(m_fin, '05-Sep-2008 12:00', '05-Oct-2008 06:00', 2, 'const', 'major8','tidal_database','../namo/h_tpxo9.v1.nc');

%% Save the mesh
write(m_fin, 'musikal_v1', 'f14');
write(m_fin, 'musikal_v1', 'f15');
write(m_fin, 'musikal_v1', 'f13');
