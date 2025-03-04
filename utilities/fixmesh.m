function [p,t,pix] = fixmesh(p,t,ptol)
%FIXMESH  Remove duplicated/unused nodes and fix element orientation.
%   [P,T]=FIXMESH(P,T)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

if nargin<3, ptol=1024*eps; end
if nargin>=2 && (isempty(p) || isempty(t)), pix=1:size(p,1); return; end

if ptol > 0
    snap      = max(max(p,[],1)-min(p,[],1),[],2)*ptol;
    [~,ix,jx] = unique(round(p/snap)*snap,'rows','stable');
    % namo

    p = p(ix,:);
%     m1 = 7469978;
%     C = p(m1+1:end,:);
%     D = setdiff(C, p(1:m1,:), 'rows', 'stable');
%     p = [p(1:m1,:); D];
else
    ix = [1:length(p)]'; jx = ix;
end

if nargin>=2
    t=reshape(jx(t),size(t));
    
    % namo - comment these
    [pix,~,jx1]=unique(t);
    
    t=reshape(jx1,size(t));
    p=p(pix,:);
    pix=ix(pix);
   
    
    if size(t,2) == size(p,2)+1
        flip          = simpvol(p,t)<0;
        t(flip,[1,2]) = t(flip,[2,1]);
    end
end
