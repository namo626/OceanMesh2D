function obj = assignConstManning(m, points, val)

%% Find the points inside the polygon
%id = find(inpolygon(m.p(:,1), m.p(:,2), pgon(:,1), pgon(:,2)));
id = dsearchn(m.p, points);
%% Insert the new values inside the f13 struct
vals1 = m.f13.userval.Atr(3).Val;
vals2 = vals1;

[~,ia,ib] = intersect(vals2(1,:), id);
vals2(2,ia) = val;
[C,ic] = setdiff(id, vals2(1,:));
extraVals = val + zeros(1,length(ic));
vals3 = [vals2 [C'; extraVals']];
[~,ids] = sort(vals3(1,:));
vals3 = vals3(:,ids);

obj = m;
obj.f13.userval.Atr(3).Val = vals3;
obj.f13.userval.Atr(3).usernumnodes = size(vals3,2);

end