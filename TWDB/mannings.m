%% Load the meshes
m2008 = msh('fname','2008_mesh.14','aux',{'2008_mesh.13'});
m = msh('fname','30m_cut_v2.14','aux',{'30m_cut_v2.13'});
%% Interpolation
% Node numbers of manning-specified nodes
ids = m2008.f13.userval.Atr(4).Val(1,:);
% x and y location of those nodes
xs = m2008.p(ids,1);
ys = m2008.p(ids,2);
vals = m2008.f13.userval.Atr(4).Val(2,:)';
% interpolate the values
F = scatteredInterpolant(xs,ys,vals);

%% Get the polygon of the insert region
% shp = shaperead('~/Documents/cutting2.shp');
% idx = find([shp.id]==2);
shp = shaperead('~/Documents/newOrleans.shp');
idx = find([shp.id]==1);
ins1 = [shp(idx).X]';
ins2 = [shp(idx).Y]';
% Find the indices of points of the 30m mesh inside this polygon
id2 = find(inpolygon(m.p(:,1), m.p(:,2), ins1, ins2));
% interpolate the values
new_values = F(m.p(id2,1), m.p(id2,2));

%% Insert the new values inside the f13 struct
vals1 = m.f13.userval.Atr(3).Val;
vals2 = vals1;

[~,ia,ib] = intersect(vals2(1,:), id2);
vals2(2,ia) = new_values(ib);
[C,ic] = setdiff(id2, vals2(1,:));
vals3 = [vals2 [C'; new_values(ic)']];
[~,ids] = sort(vals3(1,:));
vals3 = vals3(:,ids);

%% Write to the mesh
m_new = m;
m_new.f13.userval.Atr(3).Val = vals3;
m_new.f13.userval.Atr(3).usernumnodes = size(vals3,2);

write(m_new, '30m_cut_v5', 'f13');

%% Second approach; manually set the ship channel from QGIS
shp = shaperead('~/Documents/ship.shp');
idx = find([shp.id]==1);
ins1 = [shp(idx).X]';
ins2 = [shp(idx).Y]';
polg = polyshape(ins1,ins2);