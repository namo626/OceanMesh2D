mv8 = msh('30m_cut_v8.14');
m45 = msh('fort.14.45_rivers');

%%
nodes8 = nonzeros(mv8.bd.nbvv(:,mv8.bd.ibtype==22));
nodes45 = nonzeros(m45.bd.nbvv(:,m45.bd.ibtype==22));


xy8 = mv8.p(int64(nodes8),:);
xy45 = m45.p(int64(nodes45),:);


if isequal(xy8,xy45)
    disp('The new mesh has the same river coordinates!!');
else
    disp('The new mesh has wrong river coordinates...');
end