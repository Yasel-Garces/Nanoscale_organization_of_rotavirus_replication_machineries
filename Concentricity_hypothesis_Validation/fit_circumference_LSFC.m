% Adjust a circumference to a set of points in R^2. The algorithm is based in a least square 
% approach. Details about this algorithm is available in this paper:
% "Nanoscale organization of rotavirus replication machineries", eLife
% INPUT:
%     x,y: data points
% OUTPUT:
%     v: implicit form of the circumference (R,cx,cy).
% Author: Yasel Garces (88yasel@gmail.com).

function v = fit_circumference_LSFC(x,y)

% Traslation of the coordinate system
u=x-mean(x);
v=y-mean(y);

% Simplification notation
% Eq. 4
S_uu=sum(u.^2);
S_uuu=sum(u.^3);
S_uv=sum(u.*v);
S_uvv=sum(u.*(v.^2));
% Eq. 5
S_vv=sum(v.^2);
S_vvv=sum(v.^3);
S_vuu=sum(v.*(u.^2));

% Relations matrix (Mx=V)
M=[S_uu S_uv;S_uv S_vv];
V=[0.5*(S_uuu+S_uvv) 0.5*(S_vvv+S_vuu)];

% Value x=inv(M)*V
phi=inv(M)*V';

% Compute the center
cx=mean(x)+phi(1);
cy=mean(y)+phi(2);

% Compute the radius R 
R2=sum(phi.^2)+((S_uu+S_vv)/length(x));
R=sqrt(R2);

% Save parameters in a vector and return function
v=[R cx cy];