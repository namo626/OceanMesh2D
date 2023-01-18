function m = Levees2Islands(m,mfp_file,ocean_file,centerlines_file,tg_dist)
% clear all 
% close all
% clc


bbox = [-92 -89.3; 29.1 31.5];
m_or = m;
% m = ExtractSubDomain(m,bbox);
% save('-v7.3','test','m') ;

load(mfp_file)
mfp = ExtractSubDomain(mfp,bbox);
mfp = clean(mfp,'passive');

load(ocean_file)
m_tid = ExtractSubDomain(m_tid,bbox);
m_tid = clean(m_tid,'passive');


%% Define levees based on mesh

shp = shaperead(centerlines_file);
levees = [[shp.X]',[shp.Y]'];
[levees2(:,1),levees2(:,2)] = interpm(levees(:,1),levees(:,2),km2deg(10/1000));
levees(isnan(levees(:,1)),:) = [];
levees = unique(levees,'rows');
[Idx,Dist] = ourKNNsearch(m.p',levees',1);
Dist = 1000*deg2km(Dist);
b = find(Dist>tg_dist);
levees(b,:) = [];

[Idx,Dist] = ourKNNsearch(m.p',levees',1);
Dist = 1000*deg2km(Dist);
levees = m.p(Idx,:);
p_idx1 = Idx;

tri  = triangulation(m.t,m.p) ;
edges = tri.edges;
[cn,on] = inpoly(levees2,m.p,edges);
in = find(cn ==1 & on ==0);
TI = unique(pointLocation(tri, levees2(in,:)));
TI(isnan(TI),:)=[];
p_idx = m.t(TI,:);
[Idx,Dist] = ourKNNsearch(m.p(p_idx,:)',levees2(in,:)',1);
p_idx = [p_idx1;p_idx(Idx)];
mdl2 = KDTreeSearcher(m.p(p_idx,:));


levees_ids = cell(length(shp),1);
shp2 = shaperead('usace_survey_centerlines_points_noz_cleaned');
pts = [[shp2.X]',[shp2.Y]',0.3048*[shp2.Elev_good]'];

closed = 0 ;
mdl2_tid = KDTreeSearcher(m_tid.p);
mdl2_mfp = KDTreeSearcher(mfp.p);

for i = 1:length(shp)
    aux1 = [[shp(i).X]',[shp(i).Y]'];
    if aux1(1,1) == aux1(end-1,1) && aux1(1,2) == aux1(end-1,2)
        closed =1;
    end
    [aux(:,1),aux(:,2)] = interpm(aux1(1:end-1,1),aux1(1:end-1,2),km2deg(10/1000));
    [Idx, Dist] = knnsearch(mdl2,aux,'k',1);
    Dist = 1000*deg2km(Dist);
    Idx(Dist>tg_dist) = [];
    Idx = unique(Idx,'stable');
    aux2 = unique(m.p(p_idx(Idx),:),'rows');
    if ~isempty(aux2)
        [Idx, Dist] = ourKNNsearch(aux2', aux1',1);
        extra = find(~ismember(1:length(aux2(:,2)),Idx));
        Dist = 1000*deg2km(Dist);
        Idx(find(isnan(Dist))) = [];
        aux3 = aux2(Idx,:);
        for ii = 1: length(extra)
            tri_p_id1 = ourKNNsearch(tri.Points',aux2(extra(ii),:)',1);
            tri_p_id2 = ourKNNsearch(tri.Points',aux3',1);
            [row,col] = ind2sub(size(edges),find(edges==tri_p_id1));
            ed_all = edges(row,:);
            neig= find(ismember(ed_all,tri_p_id2));
            nods = tri.Points(ed_all(neig),:);
            
            if isempty(nods)
                extra(end+1) = extra(ii);
            elseif length(neig)>1
                nods_ids = min(ourKNNsearch(aux3',nods',1));
                aux3 = [aux3(1:nods_ids,:);aux2(extra(ii),:);aux3(nods_ids+1:end,:)];
            else
                nods_ids = min(ourKNNsearch(aux3',nods',1));
                if (aux2(extra(ii),1)>aux3(nods_ids,1) && aux2(extra(ii),1)<aux3(nods_ids+1,1)) || ...
                        (aux2(extra(ii),2)>aux3(nods_ids,2) && aux2(extra(ii),2)<aux3(nods_ids+1,2))
                    aux3 = [aux3(1:nods_ids,:);aux2(extra(ii),:);aux3(nods_ids+1:end,:)];
                else
                    aux3 = [aux3(1:nods_ids-1,:);aux2(extra(ii),:);aux3(nods_ids:end,:)];
                end
            end
            
        end
        aux4  = ourKNNsearch(pts(:,1:2)',aux3',1);
        levees_ids{i}(:,2:4) = unique([aux3 ,pts(aux4,3)],'rows','stable');
        aux5 = knnsearch(mdl2_tid,levees_ids{i}(:,2:3));
        levees_ids{i}(:,5) = m_tid.b(aux5);
        aux5 = knnsearch(mdl2_mfp,levees_ids{i}(:,2:3));
        levees_ids{i}(:,6) = mfp.b(aux5);
        
        if closed == 1
            levees_ids{i}(end+1,2:4) = levees_ids{i}(1,2:4);
            closed = 0;
        end
        levees_ids{i}(:,1) = [1:length(levees_ids{i}(:,1))]';
    else
        levees_ids{i}= [];
    end
    
    clear aux aux1 Idx
end


%% Open levees BC

centroids = baryc(m);
a = cell2mat(cellfun(@size,levees_ids,'UniformOutput',false));
a = find(a(:,1)<=1);
levees_ids(a) = []; %levees_pos(a)=[];

% Delete levees out of the domain
cell2 = extdom_polygon(extdom_edges2(m.t,m.p),m.p,-1,0);
poly_vec2 = cell2mat(cell2');
[edges2] = Get_poly_edges(poly_vec2);
int_t = []; uni = [];
for i= 1:length(levees_ids)
    del = ismember(levees_ids{i}(:,2:3),poly_vec2,'rows');
    levees_ids{i}(del,:) = [];
end

centroids = baryc(m);
a = cell2mat(cellfun(@size,levees_ids,'UniformOutput',false));
a = find(a(:,1)<=1);
levees_ids(a) = []; %levees_pos(a)=[];
levees_all = [];
for i = 1:length(levees_ids)
    levees_all = [levees_all;levees_ids{i}; nan(1,6)];
end

for i= 1:length(levees_ids)
    disp(['Assigning bc for levees ', num2str(i)])
       
    %If levee is open
    if levees_ids{i}(1,2)~=levees_ids{i}(end,2) && levees_ids{i}(1,3)~=levees_ids{i}(end,3)
        
        intersections = find(levees_all(:,2)==levees_ids{i}(1,2) & levees_all(:,3)==levees_ids{i}(1,3));
        if length(intersections) >= 2
            int_t(end+1,1:2) = levees_ids{i}(1,2:3);
            int_t(end,3) = i;
        elseif length(intersections) ==1
            uni(end+1,1:4) = [levees_ids{i}(1,2:3),i,2];
        end
        
        intersections = find(levees_all(:,2)==levees_ids{i}(end,2) & levees_all(:,3)==levees_ids{i}(end,3));
        if length(intersections) >= 2
            int_t(end+1,1:2) = levees_ids{i}(end,2:3);
            int_t(end,3) = i;
        elseif length(intersections) ==1
            uni(end+1,1:4) = [levees_ids{i}(end,2:3),i,length(levees_ids{i}(:,1))-1];
        end
        
        l1 = length(m.p(:,1))+1;
        l2 = length(m.p(:,1))+length(levees_ids{i}(:,1));
        ids_pos = [l1:1:l2]';
        
        auxA = ismember(m.p(m.t(:,1),:),levees_ids{i}(2:end-1,2:3),'rows');
        auxB = ismember(m.p(m.t(:,2),:),levees_ids{i}(2:end-1,2:3),'rows');
        auxC = ismember(m.p(m.t(:,3),:),levees_ids{i}(2:end-1,2:3),'rows');
        
        aux2 = sum([auxA,auxB,auxC],2);
        t_1c = find(aux2 == 1);
        t_2c = find(aux2 == 2);
        t_3c = find(aux2 == 3);
        
        % Fix elements with one nodes in the levees
        for j =1:length(t_1c)
            
            nodes = m.t(t_1c(j),:)';
            aux4 = ismember(m.p(nodes,:),levees_ids{i}(2:end-1,2:3),'rows');
            out(j,1:2) = nodes(aux4==0);
            x = m.p(out(j,:),1); y = m.p(out(j,:),2);
            
            in = nodes(aux4==1);
            aux4 = ourKNNsearch(levees_ids{i}(2:end-1,2:3)',m.p(in,:)',1);
            aux4 = sort(aux4)+1;
            if aux4(1)==1
                x1 = levees_ids{i}(aux4(1),2); y1 = levees_ids{i}(aux4(1),3);
                x2 = levees_ids{i}(aux4(1)+1,2); y2 = levees_ids{i}(aux4(1)+1,3);
            elseif aux4(end)==length(levees_ids{i}(:,1))
                x1 = levees_ids{i}(aux4(1)-1,2); y1 = levees_ids{i}(aux4(1)-1,3);
                x2 = levees_ids{i}(aux4(1),2); y2 = levees_ids{i}(aux4(1),3);
            else
                x1 = levees_ids{i}(aux4(1)-1,2); y1 = levees_ids{i}(aux4(1)-1,3);
                x2 = levees_ids{i}(aux4(1)+1,2); y2 = levees_ids{i}(aux4(1)+1,3);
            end
            
            cent = centroids (t_1c(j),:);
            d(j,1) = (cent(:,1)-x1)*(y2-y1)-(cent(:,2)-y1)*(x2-x1);
            
            in_id = find(nodes == in);
            id_mod = find(levees_ids{i}(:,2)== m.p(in,1) & levees_ids{i}(:,3)== m.p(in,2));
            
            if d(j,1) >0
                m.t(t_1c(j),in_id) = ids_pos(id_mod);
            end
        end
        clear d x x1 x2 y y1 y2
        
        % Fix elements with 2 nodes in the levees
        for j =1:length(t_2c)
            
            nodes = m.t(t_2c(j),:)';
            aux4 = ismember(m.p(nodes,:),levees_ids{i}(2:end-1,2:3),'rows');
            out(j,1) = nodes(aux4==0);
            x = m.p(out(j),1); y = m.p(out(j),2);
            
            in = nodes(aux4==1);
            aux4 = ourKNNsearch(levees_ids{i}(2:end-1,2:3)',m.p(in,:)',2);
            aux4 = sort(aux4)+1;
            x1 = levees_ids{i}(aux4(1),2); y1 = levees_ids{i}(aux4(1),3);
            x2 = levees_ids{i}(aux4(2),2); y2 = levees_ids{i}(aux4(2),3);
            
            d(j,1) = (x-x1)*(y2-y1)-(y-y1)*(x2-x1)*(x-x1)*(y2-y1)-(y-y1)*(x2-x1);
            
            in_id(1) = find(nodes == in(1)); in_id(2) = find(nodes == in(2));
            id_mod(1) = find(levees_ids{i}(:,2)== m.p(in(1),1) & levees_ids{i}(:,3)== m.p(in(1),2));
            id_mod(2) = find(levees_ids{i}(:,2)== m.p(in(2),1) & levees_ids{i}(:,3)== m.p(in(2),2));
            
            if d(j)>0 && abs(aux4(2)-aux4(1))==1
                m.t(t_2c(j),in_id(1)) = ids_pos(id_mod(1));
                m.t(t_2c(j),in_id(2)) = ids_pos(id_mod(2));
            elseif abs(aux4(2)-aux4(1))~=1
                try
                    x1 = levees_ids{i}(aux4(1),2); y1 = levees_ids{i}(aux4(1),3);
                    x2 = levees_ids{i}(aux4(1)+1,2); y2 = levees_ids{i}(aux4(1)+1,3);
                catch
                    x1 = levees_ids{i}(aux4(1)-1,2); y1 = levees_ids{i}(aux4(1)-1,3);
                    x2 = levees_ids{i}(aux4(1),2); y2 = levees_ids{i}(aux4(1),3);
                end
                d(j,1) = (x-x1)*(y2-y1)-(y-y1)*(x2-x1)*(x-x1)*(y2-y1)-(y-y1)*(x2-x1);
                if d(j,1)>0
                    m.t(t_2c(j),in_id(1)) = ids_pos(id_mod(1));
                    m.t(t_2c(j),in_id(2)) = ids_pos(id_mod(2));
                end
            end
        end
        clear d x x1 x2 y y1 y2
        
        % Fix elements with 3 nodes in the levees
        for j =1:length(t_3c)
            
            nodes = m.t(t_3c(j),:)';
            aux4 = ismember(m.p(nodes,:),levees_ids{i}(2:end-1,2:3),'rows');
            x = centroids(t_3c(j),1); y = centroids(t_3c(j),2);
            
            in = nodes(aux4==1);
            aux4 = ourKNNsearch(levees_ids{i}(2:end-1,2:3)',m.p(in,:)',2);
            aux4 = sort(aux4)+1;
            if abs(aux4(1)-aux4(2))==1
                x1 = levees_ids{i}(aux4(1),2); y1 = levees_ids{i}(aux4(1),3);
                x2 = levees_ids{i}(aux4(2),2); y2 = levees_ids{i}(aux4(2),3);
            elseif abs(aux4(2)-aux4(3))==1
                x1 = levees_ids{i}(aux4(2),2); y1 = levees_ids{i}(aux4(2),3);
                x2 = levees_ids{i}(aux4(3),2); y2 = levees_ids{i}(aux4(3),3);
            else
                x1 = levees_ids{i}(aux4(1),2); y1 = levees_ids{i}(aux4(1),3);
                x2 = levees_ids{i}(aux4(3),2); y2 = levees_ids{i}(aux4(3),3);
            end
            d(j,1) = (x-x1)*(y2-y1)-(y-y1)*(x2-x1)*(x-x1)*(y2-y1)-(y-y1)*(x2-x1);
            
            in_id(1) = find(nodes == in(1)); in_id(2) = find(nodes == in(2)); in_id(3) = find(nodes == in(3));
            id_mod(1) = find(levees_ids{i}(:,2)== m.p(in(1),1) & levees_ids{i}(:,3)== m.p(in(1),2));
            id_mod(2) = find(levees_ids{i}(:,2)== m.p(in(2),1) & levees_ids{i}(:,3)== m.p(in(2),2));
            id_mod(3) = find(levees_ids{i}(:,2)== m.p(in(3),1) & levees_ids{i}(:,3)== m.p(in(3),2));
            
            if d(j)>0
                m.t(t_3c(j),in_id(1)) = ids_pos(id_mod(1));
                m.t(t_3c(j),in_id(2)) = ids_pos(id_mod(2));
                m.t(t_3c(j),in_id(3)) = ids_pos(id_mod(3));
            end
        end
        clear t_3c t_2c t_1c ids_pos aux aux2
        
        m.p = [m.p ; levees_ids{i}(1:end-1,2:3)];
        
    else
        % Fix elements with closed levees
        l1 = length(m.p(:,1))+1;
        l2 = length(m.p(:,1))+length(levees_ids{i}(:,1))-1;
        ids_pos = [l1:1:l2]';
        
        %aux = ismember(m.t,levees_ids{i}(:,1));
        auxA = ismember(m.p(m.t(:,1),:),levees_ids{i}(1:end,2:3),'rows');
        auxB = ismember(m.p(m.t(:,2),:),levees_ids{i}(1:end,2:3),'rows');
        auxC = ismember(m.p(m.t(:,3),:),levees_ids{i}(1:end,2:3),'rows');
        
        aux2 = sum([auxA,auxB,auxC],2);
        t_1c = find(aux2 == 1);
        t_2c = find(aux2 == 2);
        t_3c = find(aux2 == 3);
        
        poly = [levees_ids{i}(:,2:3); nan nan];
        [edges2] = Get_poly_edges(poly);
        
        % Fix elements with one nodes in the levees
        for j =1:length(t_1c)
            
            nodes = m.t(t_1c(j),:)';
            aux4 = ismember(m.p(nodes,:),levees_ids{i}(:,2:3),'rows');
            out(j,1:2) = nodes(aux4==0);
            x = m.p(out(j,:),1); y = m.p(out(j,:),2);
            
            in = nodes(aux4==1);
            
            cent = centroids (t_1c(j),:);
            in2 = inpoly(cent,poly,edges2);
            
            in_id = find(nodes == in);
            id_mod = find(levees_ids{i}(:,2)== m.p(in,1) & levees_ids{i}(:,3)== m.p(in,2));
            
            if in2 == 1
                m.t(t_1c(j),in_id) = ids_pos(id_mod(1));
            end
            
        end
        clear d x x1 x2 y y1 y2
        
        % Fix elements with 2 nodes in the levees
        for j =1:length(t_2c)
            
            nodes = m.t(t_2c(j),:)';
            aux4 = ismember(m.p(nodes,:),levees_ids{i}(:,2:3),'rows');
            out(j,1) = nodes(aux4==0);
            x = m.p(out(j),1); y = m.p(out(j),2);
            
            in = nodes(aux4==1);
            aux4 = ourKNNsearch(levees_ids{i}(:,2:3)',m.p(in,:)',2);
            aux4 = sort(aux4);
            x1 = levees_ids{i}(aux4(1),2); y1 = levees_ids{i}(aux4(1),3);
            x2 = levees_ids{i}(aux4(2),2); y2 = levees_ids{i}(aux4(2),3);
            
            in2 = inpoly([x y],poly,edges2);
            
            in_id(1) = find(nodes == in(1)); in_id(2) = find(nodes == in(2));
            id_mod(1) = find(levees_ids{i}(1:end-1,2)== m.p(in(1),1) & levees_ids{i}(1:end-1,3)== m.p(in(1),2));
            id_mod(2) = find(levees_ids{i}(1:end-1,2)== m.p(in(2),1) & levees_ids{i}(1:end-1,3)== m.p(in(2),2));
            
            if in2 ==1 %d(j)>0 && abs(aux4(2)-aux4(1))==1
                m.t(t_2c(j),in_id(1)) = ids_pos(id_mod(1));
                m.t(t_2c(j),in_id(2)) = ids_pos(id_mod(2));
            end
            
            
        end
        clear d x x1 x2 y y1 y2
        
        % Fix elements with 3 nodes in the levees
        for j =1:length(t_3c)
            
            nodes = m.t(t_3c(j),:)';
            aux4 = ismember(m.p(nodes,:),levees_ids{i}(:,2:3),'rows');
            x = centroids(t_3c(j),1); y = centroids(t_3c(j),2);
            
            in = nodes(aux4==1);
            aux4 = ourKNNsearch(levees_ids{i}(:,2:3)',m.p(in,:)',2);
            aux4 = sort(aux4);
            
            cent = centroids (t_3c(j),:);
            in2 = inpoly(cent,poly,edges2);
            
            in_id(1) = find(nodes == in(1)); in_id(2) = find(nodes == in(2)); in_id(3) = find(nodes == in(3));
            id_mod(1,:) = find(levees_ids{i}(:,2)== m.p(in(1),1) & levees_ids{i}(:,3)== m.p(in(1),2));
            id_mod(2,:) = find(levees_ids{i}(:,2)== m.p(in(2),1) & levees_ids{i}(:,3)== m.p(in(2),2));
            id_mod(3,:) = find(levees_ids{i}(:,2)== m.p(in(3),1) & levees_ids{i}(:,3)== m.p(in(3),2));
            
            if in2==1
                m.t(t_3c(j),in_id(1)) = ids_pos(id_mod(1,1));
                m.t(t_3c(j),in_id(2)) = ids_pos(id_mod(2,1));
                m.t(t_3c(j),in_id(3)) = ids_pos(id_mod(3,1));
            end
            
            
        end
        clear t_3c t_2c t_1c ids_pos aux aux2
        
        m.p = [m.p ; levees_ids{i}(1:end-1,2:3)];
        
    end
end

save('levees_ids.mat','levees_ids')

%% Delete disjoint nodes and reorder the connectivy table

orig = knnsearch(m_or.p,m.p,'K',1);
m.b = []; m.b = m_or.b(orig);

aux = [1:1:length(m.p(:,1))]';
del = find(~ismember(aux,m.t(:,1)) & ~ismember(aux,m.t(:,2)) & ~ismember(aux,m.t(:,3)));

for i =1:length(del)
    
    modi = find(m.t(:,1)>del(i)-i+1);
    m.t(modi,1) = m.t(modi,1)-1;
    modi = find(m.t(:,2)>del(i)-i+1);
    m.t(modi,2) = m.t(modi,2)-1;
    modi = find(m.t(:,3)>del(i)-i+1);
    m.t(modi,3) = m.t(modi,3)-1;
    
end
m.p(del,:) = [];
m.b(del) = [];


vx = m.p(:,1); vy = m.p(:,2);
xt = [vx(m.t(:,1)) vx(m.t(:,2)) vx(m.t(:,3)) vx(m.t(:,1))];
yt = [vy(m.t(:,1)) vy(m.t(:,2)) vy(m.t(:,3)) vy(m.t(:,1))];
dxt = diff(xt,[],2);
dyt = diff(yt,[],2);
xt(:,end) = xt(:,1);
dxt = diff(xt,[],2);
Area = dxt(:,3).*-dyt(:,2) + dxt(:,2).*dyt(:,3);
modi = find(Area<0);
m.t(modi,2:3) = [m.t(modi,3) m.t(modi,2)];



%% Fix t connections
if ~isempty(int_t)
for  i = 1: length(int_t(:,1))
    %      try
    disp(['Triangle number: ',num2str(i)])
    
    % Make sure there are 3 copies of the nodes
    node_ids = find(m.p(:,1) == int_t(i,1) & m.p(:,2) == int_t(i,2));
    if length(node_ids) == 1
        node_ids(end+1) = length(m.p(:,1))+1;
        node_ids(end+1) = length(m.p(:,1))+1;
        m.p(end+1,:) = m.p(node_ids(1),:);
        m.p(end+1,:) = m.p(node_ids(1),:);
    elseif length(node_ids) == 2
        node_ids(end+1) = length(m.p(:,1))+1;
        m.p(end+1,:) = m.p(node_ids(1),:);
    end
    
    % Find all the elements associated to the node
    elem_id1 = find(m.t(:,1)== node_ids(1) | m.t(:,1)== node_ids(2) | m.t(:,1)== node_ids(3));
    elem_id2 = find(m.t(:,2)== node_ids(1) | m.t(:,2)== node_ids(2) | m.t(:,2)== node_ids(3));
    m.t(elem_id2,:) = [m.t(elem_id2,2), m.t(elem_id2,1), m.t(elem_id2,3)];
    elem_id3 = find(m.t(:,3)== node_ids(1) | m.t(:,3)== node_ids(2) | m.t(:,3)== node_ids(3));
    m.t(elem_id3,:) = [m.t(elem_id3,3), m.t(elem_id3,1:2)];
    elem_id = unique([elem_id1;elem_id2;elem_id3]);
    
    
    % Find edges in leeves systems and reorganize conectivity table
    for a = 1:length(levees_ids)
        aux = find(levees_ids{a}(:,2) == m.p(node_ids(1),1) & levees_ids{a}(:,3) == m.p(node_ids(1),2));
        if aux == 1
            seg1 = levees_ids{a}(1:2,2:3);
        elseif aux == length(levees_ids{a}(:,1))
            seg1 = flipud(levees_ids{a}(end-1:end,2:3));
        elseif ~isempty(aux)
            seg2 = flipud(levees_ids{a}(aux-1:aux,2:3));
            seg3 = levees_ids{a}(aux:aux+1,2:3);
        end
    end
    
    seg1_nod = find(seg1(:,1) ~= int_t(i,1) & seg1(:,2)~= int_t(i,2));
    seg1_nod = find(m.p(:,1)==seg1(seg1_nod,1) & m.p(:,2)==seg1(seg1_nod,2));
    [row,col] = ind2sub([length(elem_id) 3],find(ismember(m.t(elem_id,:),seg1_nod)));
    for a = 1:length(col)
        if col(a)~=2
            m.t(elem_id(row(a)),:) = [m.t(elem_id(row(a)),1), fliplr(m.t(elem_id(row(a)),2:3)) ];
        end
    end
    aux = elem_id(row); elem_id(row)=[] ;
    
    seg2_nod = find(seg2(:,1) ~= int_t(i,1) & seg2(:,2)~= int_t(i,2));
    seg2_nod = find(m.p(:,1)==seg2(seg2_nod,1) & m.p(:,2)==seg2(seg2_nod,2));
    [row,col] = ind2sub([length(elem_id) 3],find(ismember(m.t(elem_id,:),seg2_nod)));
    for a = 1:length(col)
        if col(a)~=2
            m.t(elem_id(row(a)),:) = [m.t(elem_id(row(a)),1), fliplr(m.t(elem_id(row(a)),2:3)) ];
        end
    end
    aux = [aux; elem_id(row)]; elem_id(row)=[] ;
    
    seg3_nod = find(seg3(:,1) ~= int_t(i,1) & seg3(:,2)~= int_t(i,2));
    seg3_nod = find(m.p(:,1)==seg3(seg3_nod,1) & m.p(:,2)==seg3(seg3_nod,2));
    [row,col] = ind2sub([length(elem_id) 3],find(ismember(m.t(elem_id,:),seg3_nod)));
    for a = 1:length(col)
        if col(a)~=2
            m.t(elem_id(row(a)),:) = [m.t(elem_id(row(a)),1), fliplr(m.t(elem_id(row(a)),2:3)) ];
        end
    end
    aux = [aux; elem_id(row)]; elem_id(row)=[] ;
    if ~isempty(elem_id)
        elem_id = [aux; elem_id(~ismember(elem_id,aux))];
    else
        elem_id = aux;
    end
    
    
    seg1_nod = find(m.p(:,1)==seg1(2,1) & m.p(:,2)==seg1(2,2));
    [row,col] = ind2sub([length(elem_id) 3],find(ismember(m.t(elem_id,:),seg1_nod(1))));
    if col(1) == 3
        m.t(elem_id(row(1)),2:3) = fliplr(m.t(elem_id(row(1)),2:3));
    end
    
    aux = elem_id(row(1)); elem_id(row(1)) = []; %aux2 = elem_id(row(2));
    u=1;
    while ~isempty(elem_id)
        [row,col] = ind2sub([length(elem_id) 3],find(ismember(m.t(elem_id,:),m.t(aux(end),3))));
        if ~isempty(col)
            if col == 3
                m.t(elem_id(row),2:3) = fliplr(m.t(elem_id(row),2:3));
            end
            aux = [aux;elem_id(row)]; elem_id(row) = [];
        else
            g(u) = length(aux);
            u=u+1;
            aux = [aux;elem_id(1)];
            elem_id(1) = [];
        end
        
    end
    elem_id = aux;
    g1 = elem_id(1:g(1));
    g2 = elem_id(g(1)+1:g(2));
    g3 = elem_id(g(2)+1:end);
    if length(g)~=2
        disp(num2str(i))
    end
    clear g
    
    m.t(g1,1) = node_ids(1);
    m.t(g2,1) = node_ids(2);
    m.t(g3,1) = node_ids(3);
    
    % Error check
    if sum(sum(ismember(m.t(g1,:),m.t(g2,:))))~=0
        warning(['Error T connection ',num2str(i)])
    end
    if sum(sum(ismember(m.t(g1,:),m.t(g3,:))))~=0
        warning(['Error T connection ',num2str(i)])
    end
    if sum(sum(ismember(m.t(g2,:),m.t(g3,:))))~=0
        warning(['Error T connection ',num2str(i)])
    end
    
    %     catch
    %         warning(['fix levees intersection: ',num2str(i)])
    %     end
    
    
    
end
end

%% Fix uni connections
if ~isempty(uni)
for  i = 1: length(uni(:,1))
    
    % try
    disp(['Uni connection number: ',num2str(i)])
    
    % Make sure there are 3 copies of the nodes
    cent = levees_ids{uni(i,3)}(uni(i,4),2:3);
    node_ids1 = find(m.p(:,1) == cent(1,1) & m.p(:,2) == cent(1,2));
    
    % Find all the elements associated to the node
    elem_id1 = find(m.t(:,1)== node_ids1(1) | m.t(:,1)== node_ids1(2));
    elem_id2 = find(m.t(:,2)== node_ids1(1) | m.t(:,2)== node_ids1(2));
    m.t(elem_id2,:) = [m.t(elem_id2,2), m.t(elem_id2,1), m.t(elem_id2,3)];
    elem_id3 = find(m.t(:,3)== node_ids1(1) | m.t(:,3)== node_ids1(2));
    m.t(elem_id3,:) = [m.t(elem_id3,3), m.t(elem_id3,1:2)];
    elem_id = unique([elem_id1;elem_id2;elem_id3]);
    
    if uni(i,4)==2
        seg1 = [uni(i,1:2);cent];
        seg2 = [cent;levees_ids{uni(i,3)}(uni(i,4)+1,2:3)];
    else
        seg1 = [uni(i,1:2);cent];
        seg2 = [cent;levees_ids{uni(i,3)}(uni(i,4)-1,2:3)];
    end
    
    node_ids2 = find(m.p(:,1) == uni(i,1)...
        & m.p(:,2) == uni(i,2));
    
    node_ids3 = find(m.p(:,1) == seg2(end,1)...
        & m.p(:,2) == seg2(end,2));
    
    
    node_ids = [node_ids1;node_ids2;node_ids3];
    for j = 1:length(elem_id)
        aux1 = find(ismember(m.t(elem_id(j),:),node_ids));
        aux2 = find(~ismember(m.t(elem_id(j),:),node_ids));
        m.t(elem_id(j),:) = [m.t(elem_id(j),aux1),m.t(elem_id(j),aux2)];
        if m.p(m.t(elem_id(j),1),1)~=cent(1,1) || m.p(m.t(elem_id(j),1),2)~=cent(1,2)
            m.t(elem_id(j),2:3)= fliplr(m.t(elem_id(j),2:3));
        end
    end

    aux = ismember(m.t(elem_id,:),node_ids);
    aux3 = ismember(m.t(elem_id,:),node_ids3);
    corr_row = find(sum(aux,2) == 2 & sum(aux3,2)==1);
    
    corr2_id = find(m.t(elem_id,1) ~= m.t(elem_id(corr_row(1)),1),1,'first');
    corr2_id = m.t(elem_id(corr2_id),1);
        
    uni_id = find(m.p(:,1) == uni(i,1) & m.p(:,2) == uni(i,2));
    g1 = elem_id(corr_row(1));
    elem_id(corr_row(1)) = [] ;
    while m.t(g1(end),2)~=uni_id && m.t(g1(end),3)~=uni_id
        
        [aux_r, aux_col] = ind2sub(size(m.t(elem_id,2:3)),...
            find(m.t(elem_id,2:3)== m.t(g1(end),3) | m.t(elem_id,2:3)== m.t(g1(end),2)));
        g1 = [g1;elem_id(aux_r)];
        elem_id(aux_r) = [] ;
    end
    
    m.t(g1,1) = m.t(g1(1),1);
    m.t(elem_id,1) = corr2_id;
    
end
end

%% Delete disjoint nodes and reorder the connectivy table

% m_or = load('test.mat');
% m_or = m_or.m;
orig = knnsearch(m_or.p,m.p,'K',1);
m.b = []; m.b = m_or.b(orig);

aux = [1:1:length(m.p(:,1))]';
del = find(~ismember(aux,m.t(:,1)) & ~ismember(aux,m.t(:,2)) & ~ismember(aux,m.t(:,3)));

for i =1:length(del)
    
    modi = find(m.t(:,1)>del(i)-i+1);
    m.t(modi,1) = m.t(modi,1)-1;
    modi = find(m.t(:,2)>del(i)-i+1);
    m.t(modi,2) = m.t(modi,2)-1;
    modi = find(m.t(:,3)>del(i)-i+1);
    m.t(modi,3) = m.t(modi,3)-1;
    
end
m.p(del,:) = [];
m.b(del) = [];


vx = m.p(:,1); vy = m.p(:,2);
xt = [vx(m.t(:,1)) vx(m.t(:,2)) vx(m.t(:,3)) vx(m.t(:,1))];
yt = [vy(m.t(:,1)) vy(m.t(:,2)) vy(m.t(:,3)) vy(m.t(:,1))];
dxt = diff(xt,[],2);
dyt = diff(yt,[],2);
xt(:,end) = xt(:,1);
dxt = diff(xt,[],2);
Area = dxt(:,3).*-dyt(:,2) + dxt(:,2).*dyt(:,3);
modi = find(Area<0);
m.t(modi,2:3) = [m.t(modi,3) m.t(modi,2)];

%% Adjust topo-bathy values

[u,I,J] = unique(m.p, 'rows', 'first');
aux = setdiff(1:size(m.p,1), I);
Dup_Val = m.p(aux,:);
aux = ismember(m.p,Dup_Val);
Dup_ID = find(aux(:,1)==1 & aux(:,2)==1);

centroids = baryc(m);
poly = [bbox(:,1)'; bbox(1,2) bbox(2,1); bbox(:,2)';bbox(1,1) bbox(2,2); bbox(:,1)'; nan nan];
[edges2] = Get_poly_edges(poly);
id_cent = inpoly(centroids,poly,edges2);
id_cent = find(id_cent==1);

mfp_cent = baryc(mfp);

water = knnsearch(m_tid.p,m.p,'K',1);
water = max(0.3,m_tid.b(water(Dup_ID)));

topo = knnsearch(mfp.p,m.p,'K',1);
topo = mfp.b(topo(Dup_ID));


[cell2, cell2_Idx] = extdom_polygon(extdom_edges2(mfp.t,mfp.p),mfp.p,-1,0);
mfp_bd = cellfun(@size,cell2_Idx,'UniformOutput',0);
mfp_bd = cell2mat(cell2');
mfp_bd_id = knnsearch(mfp_bd,m.p,'K',1);


tic
for i =1:length(Dup_ID)
    
    i
    if ~ismember(Dup_ID(i),mfp_bd_id)
        
        neig1 = (m.t(id_cent,1)==Dup_ID(i));
        neig2 = (m.t(id_cent,2)==Dup_ID(i));
        neig3 = (m.t(id_cent,3)==Dup_ID(i));
        neig = max([neig1,neig2,neig3],[],2);
        cent = centroids(id_cent(neig),:);
        
        aux = ~ismember(round(cent,4),round(mfp_cent,4),'rows');
        
        if sum(aux(:))>0
            m.b(Dup_ID(i))= water(i);
        else
            m.b(Dup_ID(i))= min(-0.3,topo(i));
        end
        
    else
        
        neig1 = (m.t(id_cent,1)==Dup_ID(i));
        neig2 = (m.t(id_cent,2)==Dup_ID(i));
        neig3 = (m.t(id_cent,3)==Dup_ID(i));
        neig = max([neig1,neig2,neig3],[],2);
        cent = centroids(id_cent(neig),:);
        
        aux = ~ismember(round(cent,4),round(mfp_cent,4),'rows');
        
        if sum(aux(:))>0
            m.b(Dup_ID(i))= water(i);
        else
            m.b(Dup_ID(i))= topo(i);
        end
        
    end
end
toc

%% Assign boundary conditions
m.bd = []; m.op = [];
% m.bd.barincfsb = zeros(size(m_or.bd.nbvv));
% m.bd.barincfsp = zeros(size(m_or.bd.nbvv));
% m.bd.barinht = zeros(size(m_or.bd.nbvv));
% m.bd.ibconn = zeros(size(m_or.bd.nbvv));

for i = 1:length(levees_ids)
    i
    if isempty(m.bd)
        m.bd.nbou = 1;
    else
        m.bd.nbou = m.bd.nbou+1;
    end
    m.bd.ibtype(m.bd.nbou) = 24;
    nbvv = [];
    ibconn = [];
    barincfsp = [];
    for j = 1:length(levees_ids{i}(:,1))
        
        if sum(sum(ismember(int_t,levees_ids{i}(j,2:3))))>0
            ids  = knnsearch(m.p(Dup_ID,:),levees_ids{i}(j,2:3),'K',3);
            
            [row1,col1]= ind2sub(size(m.t),find(m.t(:)==Dup_ID(ids(1))));
            [row2,col2]= ind2sub(size(m.t),find(m.t(:)==Dup_ID(ids(2))));
            [row3,col3]= ind2sub(size(m.t),find(m.t(:)==Dup_ID(ids(3))));
            rows = [row1;row2;row3]; cols = [col1;col2;col3];
            if j==1
                prev = levees_ids{i}(2,2:3);
            elseif j == length(levees_ids{i}(:,1))
                prev = levees_ids{i}(j-1,2:3);
            else
                prev = [levees_ids{i}(j-1,2:3);levees_ids{i}(j+1,2:3)];
            end
            
            for k = 1:length(prev(:,1))
                prev_ID = find(m.p(:,1) == prev(k,1) & m.p(:,2) == prev(k,2));
                row_prev = [];
                for kk =1:length(prev_ID)
                    [row_prev1,col_prev] = ind2sub(size(m.t(rows,:)),find(m.t(rows,:)==prev_ID(kk)));
                    row_prev = [row_prev;row_prev1];
                end
                
                nbvv(end+1) = m.t(rows(row_prev(1)),cols(row_prev(1)));
                ibconn (end+1) = m.t(rows(row_prev(2)),cols(row_prev(2)));
                barincfsp(end+1) = levees_ids{i}(j,4);
                barincfsp(end) = max(0.2-m.b(nbvv(end)),barincfsp(end));
                barincfsp(end) = max(0.2-m.b(ibconn(end)),barincfsp(end));
                
            end
            
        elseif sum(sum(ismember(uni,levees_ids{i}(j,2:3))))>0
            %             disp('none')
            ids = [];
        else
            ids  = knnsearch(m.p(Dup_ID,:),levees_ids{i}(j,2:3),'K',2);
            nbvv(end+1) = Dup_ID(ids(1));
            ibconn(end+1) = Dup_ID(ids(2));
            barincfsp(end+1) = levees_ids{i}(j,4);
            barincfsp(end) = max(0.2-m.b(nbvv(end)),barincfsp(end));
            barincfsp(end) = max(0.2-m.b(ibconn(end)),barincfsp(end));
        end
        
        
    end
    
    m.bd.nvell(m.bd.nbou) = length(nbvv); % Number of nodes of j levee
    m.bd.barincfsb (1:length(nbvv), m.bd.nbou) = ones(length(nbvv),1); % Coeff 1
    m.bd.barincfsp (1:length(nbvv), m.bd.nbou) = ones(length(nbvv),1); % Coeff 1
    
    
    m.bd.nbvv(1:length(nbvv), m.bd.nbou) = nbvv; % Nodes first side
    m.bd.ibconn(1:length(nbvv), m.bd.nbou) = ibconn; % Nodes second side
    m.bd.barinht(1:length(nbvv), m.bd.nbou) = barincfsp; % Levees elevation
    m.bd.nvel = length(find(m.bd.nbvv~=0));
    
end



for i = find(m.bd.ibtype==24)
    
    nvell = m.bd.nvell(i);
    aux = find(m.bd.barinht(1:nvell,i)<-1*m.b(m.bd.nbvv(1:nvell,i)));
    for j = 1:length(aux)
        if aux(j) == 1
            m.b(m.bd.nbvv(aux(j),i)) = m.b(m.bd.nbvv(aux(j)+1,i));
        elseif aux(j) == nvell
            m.b(m.bd.nbvv(aux(j),i)) = m.b(m.bd.nbvv(aux(j)-1,i));
        else
            m.b(m.bd.nbvv(aux(j),i)) = max(m.b(m.bd.nbvv(aux(j)-1,i)),...
                m.b(m.bd.nbvv(aux(j)+1,i)));
        end
        if m.bd.barinht(aux(j),i)<-1*m.b(m.bd.nbvv(aux(j),i))
            m.b(m.bd.nbvv(aux(j),i)) = -0.8*m.bd.barinht(aux(j),i);
        end
    end
    
    aux = find(m.bd.barinht(1:nvell,i)<-1*m.b(m.bd.ibconn(1:nvell,i)));
    for j = 1:length(aux)
        if aux(j) == 1
            m.b(m.bd.ibconn(aux(j),i)) = m.b(m.bd.ibconn(aux(j)+1,i));
        elseif aux(j) == nvell
            m.b(m.bd.ibconn(aux(j),i)) = m.b(m.bd.ibconn(aux(j)-1,i));
        else
            m.b(m.bd.ibconn(aux(j),i)) = max(m.b(m.bd.ibconn(aux(j)-1,i)),...
                m.b(m.bd.ibconn(aux(j)+1,i)));
        end
        if m.bd.barinht(aux(j),i)<-1*m.b(m.bd.ibconn(aux(j),i))
            m.b(m.bd.ibconn(aux(j),i)) = -0.8*m.bd.barinht(aux(j),i);
        end
    end
    
    b1 = m.b(m.bd.ibconn(1:nvell,i));
    aux = find(m.bd.barinht(1:nvell,i)<=-1*b1 +0.2);
    if ~isempty(aux)
        m.bd.barinht(aux,i)= -1*b1(aux)+0.201;
    end
    
end

tri  = triangulation(m.t,m.p) ;
edges = tri.edges;
for  i = 1:length(uni(:,1))
    id2 = knnsearch(m.p,uni(i,1:2));
    [row,col] = ind2sub(size(edges),find(edges==id2));
    
    aux = find(col==1);
    cand = edges(row(aux),2);
    aux = find(col==2);
    cand = [cand;edges(row(aux),1)];
    [u,I,J] = unique(m.p(cand,:),'rows','stable');
    ixDupRows = setdiff(1:size(cand,1), I);
    id1 = cand(ixDupRows);
    J = find(m.p(cand,1)==m.p(id1,1) & m.p(cand,2)==m.p(id1,2) & cand~=id1);
    id3 = cand(J);
    
    nbou = m.bd.nbou+1;
    m.bd.nbou = nbou;
    m.bd.nvel = m.bd.nvel +3;
    m.bd.ibtype(nbou) = 20;
    m.bd.nbvv(1:3,nbou) = [id1;id2;id3];
    m.bd.nvell(nbou) = 3;
    m.bd.barincfsb(:,nbou) = 0;
    m.bd.barincfsp(:,nbou) = 0;
    m.bd.barinht(:,nbou) = 0;
    m.bd.ibconn(:,nbou) = 0; 
    clear id1 id2 id3
end



%% Update slopes
m.bx = []; m.by = [];
if ~isempty(m_or.bx)
    orig = knnsearch(m_or.p,m.p,'K',1);
    m.bx = []; m.bx = m_or.bx(orig);
    m.by = []; m.by = m_or.by(orig);
end

%% Update boundary conditions

if ~isempty(m_or.bd)
    disp('Updating islands boundary conditions...')
    m.bd.nbou = m.bd.nbou + m_or.bd.nbou;
    m.bd.ibtype = [m.bd.ibtype, 21*ones(size(m_or.bd.ibtype))];
    m.bd.nvell = [m.bd.nvell, m_or.bd.nvell];
    
    id_or = m_or.bd.nbvv;
    id_or = id_or(id_or~=0);
    id_new = ourKNNsearch(m.p',m_or.p(id_or,:)',1);
    
    beg = 1; ends = m_or.bd.nvell(1);
    nbvv_new = zeros(size(m_or.bd.nbvv));
    for i = 1:m_or.bd.nbou
        nvel_or = m_or.bd.nvell(i);
        nbvv_new(1:nvel_or,i) = id_new(beg:ends);
        try
            beg = ends+1; ends = ends + m_or.bd.nvell(i+1);
        catch
            beg = ends+1; ends = length(id_or);
        end
    end
    
    row_max = max(length(m.bd.nbvv(:,1)),length(nbvv_new(:,1)));
    nbvv_new_all = zeros(row_max,m.bd.nbou);
    nbvv_new_all(1:size(m.bd.nbvv,1),1:size(m.bd.nbvv,2))= m.bd.nbvv;
    nbvv_new_all(1:size(nbvv_new,1),size(m.bd.nbvv,2)+1:end) = nbvv_new;
    
    
    m.bd.nbvv = nbvv_new_all;
    m.bd.nvel = length(find(m.bd.nbvv~=0));
end

if length(m.bd.ibconn(1,:))<m.bd.nbou
    m.bd.barincfsb(length(m.bd.nbvv(:,1)),m.bd.nbou) = 0;
    m.bd.barincfsp(length(m.bd.nbvv(:,1)),m.bd.nbou) = 0;
    m.bd.barinht(length(m.bd.nbvv(:,1)),m.bd.nbou) = 0;
    m.bd.ibconn(length(m.bd.nbvv(:,1)),m.bd.nbou) = 0;
end


%% Make sure topo/bathy is constrained
aux = find(m.b>=0 & m.b<0.3);
m.b(aux) = 0.3;
aux = find(m.b<0 & m.b>-0.1); 
m.b(aux) = -0.1; 

% % % %% Update nodal attributes fort.13
% % %
% % % if ~isempty(m.f13)
% % %     aux2 = cell2mat(levees_ids');
% % %     m.f13.NumOfNodes = length(m.p(:,1));
% % %     for i = 1:m.f13.nAttr
% % %
% % %         Val = [m.f13.userval.Atr(i).Val]';
% % %         idx_rep = Val(find(ismember(Val(:,1),aux2(:,1))),:);
% % %
% % %         if ~isempty(idx_rep)
% % %             Val2 = [];
% % %             for j = 1 :length(idx_rep(:,1))
% % %                 aux = find(aux2(:,1)==idx_rep(j,1));
% % %
% % %                 Val2 = [Val2; [aux2(aux,2), idx_rep(j,2)*ones(length(aux),1)]];
% % %
% % %             end
% % %             m.f13.userval.Atr(i).Val = [m.f13.userval.Atr(i).Val, Val2'];
% % %             m.f13.userval.Atr(2).usernumnodes = length(m.f13.userval.Atr(i).Val(1,:));
% % %             disp(['Updating nodal attribute: ', m.f13.userval.Atr(i).AttrName])
% % %
% % %         end
% % %
% % %     end
% % % end

% % % %% Update fort.24
% % %
% % % if ~isempty(m.f24)
% % %     disp('Updating fort.24')
% % %     Val = m.f24.Val;
% % %     for i=1:length(Val(:,1,1))
% % %         Val(i,:,aux2(:,2)) = Val(i,:,aux2(:,1));
% % %     end
% % %     m.f24.Val = Val;
% % % end

end