
%% Extract high-resolution insert
% Load High-res mesh % Call it obj1
obj1 = mesh120; 
shp = shaperead('buffer1.shp');
%load('Model_120m_Combinedv58_Fixed_Run.mat')
idx = find([shp.id]==3);
ins1 = [[shp(idx).X]',[shp(idx).Y]'];
[m_ins,ind,~] = extract_subdomain(obj1,ins1);
m_ins = map_mesh_properties(m_ins,'ind',ind);

%% Clipp the low resolution outer mesh
% Load (Coase resolution mesh) % Call it as obj2
obj2 = coarseMesh; 
idx = find([shp.id]==1);
ins2 = [[shp(idx).X]',[shp(idx).Y]'];
[m_coar,ind,~] = extract_subdomain(obj2,ins2);
m_coar = map_mesh_properties(m_coar,'ind',ind);

%% Extract the nodes along the boundaries
m_ins.op= []; 
m_ins = makens(m_ins,'outer',0);
bd1 = m_ins.p(m_ins.op.nbdv,:);

m_coar.op=[];
m_coar = makens(m_coar,'outer',0);
bd2 = m_coar.p(m_coar.op.nbdv,:);
figure
plot(bd1(:,1),bd1(:,2),bd2(:,1),bd2(:,2));

%% Define the element 

opts.iter = 100;
opts.kind = 'delaunay';
opts.ref1 = 'preserve';
bd = [bd1;flipud(bd2);bd1(1,:)];
[node,edge] = getnan2(bd);
figure;
plot(bd(:,1),bd(:,2));

% use local feature size with grade = olfs.dhdx
olfs.dhdx = 0.20;
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

%% Add bathymetry for the buffer zone mesh

mfp = interp(mfp,'../GEBCO_2021.nc');


%% Join everything

 %m_new = plus(m_ins,mfp,'matche');
 %m_new = plus(m_new,m_coar,'match');
 m_new = cat(m_ins, mfp);
 m_new = catBathy(m_new, m_ins, mfp);
 m_fin = cat(m_new, m_coar);
 m_fin = catBathy(m_fin, m_new, m_coar );

