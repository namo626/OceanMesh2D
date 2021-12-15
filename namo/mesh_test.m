% Creating a coarse mesh

addpath('..')
addpath(genpath('..utilities/'))
addpath(genpath('..datasets/'))
addpath(genpath('..m_map/'))

%% STEP 1: set mesh extents and set parameters for mesh.
bbox = [-98 -60;		% lon_min lon_max
        7.9 45.86]; 		% lat_min lat_max
min_el    = 1e3;  		% minimum resolution in meters.
max_el    = 24e3; 		% maximum resolution in meters. 
max_el_ns = 5e3;        % maximum resolution nearshore in meters.
grade     = 0.2; 		% mesh grade in decimal percent.
R         = 2;    		% number of elements to resolve feature width.
slope = 20;
wl = 100;
%% STEP 2: specify geographical datasets and process the geographical data 
%% to be used later with other OceanMesh classes...
coastline = 'GSHHS_f_L1';
dem = 'GEBCO_2021.nc';

gdat = geodata('shp',coastline, ...
                'dem',dem, ...
                'bbox',bbox, ...
                'h0',min_el);
%% STEP 3: create an edge function class
fh = edgefx('geodata',gdat,...
            'fs',R, ...
            'max_el',max_el, ...
            'g',grade, ...
            'wl',wl, ...
            'slp',slope);

%% STEP 4: Pass your edgefx class object along with some meshing options and
% build the mesh...
mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',1,'nscreen',5,'proj','equi');
mshopts = mshopts.build; 

%% STEP 5: Plot it and write a triangulation fort.14 compliant file to disk.
% Get out the msh class and put on nodestrings
coarseMesh = mshopts.grd;
coarseMesh = make_bc(coarseMesh,'auto',gdat); % make the boundary conditions
plot(coarseMesh,'type','bd','proj','equi');


%% STEP 6: Example of plotting a subdomain with bcs

plot(coarseMesh,'type','tri','subdomain',bbox)
