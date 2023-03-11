m = msh('30m_cut_v7.14');
plot(m,'type','b','proj','none');

%% Load
m5 = msh('30m_cut_v5.14');

%% Check extension bathy
edge = [-97.9, -93.6; 25.5, 30.6];
re = [-97, -95.5; 28.3, 30];

plot(m,'type','b','proj','none','subdomain',re);
plot(m5,'type','b','proj','none','subdomain',re);

