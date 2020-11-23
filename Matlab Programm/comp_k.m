function k = comp_k(x)
%COMP_DX  Compute distance between barycenter and edges of a triangle
%   COMP_DX(X) computes the minimal distance between the barycenter of the 
%   triangle defined by the three vertices X(1,:), X(2,:) and X(3,:) and
%   its three edges.

% Compute barycenter s
s = 1/3*sum(x);

% First edge
vt = x(2,:)-x(1,:);
vn = [-vt(2); vt(1)];
th = [vn -vt']\(x(1,:)-s)';
ps = x(1,:) + th(2)*vt;
k1 = norm(ps-s);

% Second edge
vt = x(3,:)-x(1,:);
vn = [-vt(2); vt(1)];
th = [vn -vt']\(x(1,:)-s)';
ps = x(1,:) + th(2)*vt;
k2 = norm(ps-s);

% Third edge
vt = x(3,:)-x(2,:);
vn = [-vt(2); vt(1)];
th = [vn -vt']\(x(2,:)-s)';
ps = x(2,:) + th(2)*vt;
k3 = norm(ps-s);

% Compute minimal distance
k = min([k1 k2 k3]);