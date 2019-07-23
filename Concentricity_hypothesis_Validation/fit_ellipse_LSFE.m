% It adjust a ellipse to points in R2 using the algorithm Direct Least Square Fitting Ellipse (DLSFE).
% The algorithm was described here:
% R. Halir and J. Flusser. Numerically stable direct least-squares fitting of ellipses. In Sixth International Conference of Computer Graphics and Visualization, 1998.
% INPUT:
%     x,y: data points
% OUTPUT:
%     v: implicit form of the ellipse (a,b,cx,cy,theta).

function v = fit_ellipse_LSFE(x,y)

% quadratic part of the design matrix
D1 = [x .^ 2, x .* y, y .^ 2]; 
% linear part of the design matrix
D2 = [x, y, ones(size(x))];    
% quadratic part of the scatter matrix
S1 = D1' * D1;  
 % combined part of the scatter matrix  
S2 = D1' * D2;   
% linear part of the scatter matrix
S3 = D2' * D2;    

if det(S3)==0
    v=[];
    return
end
% for getting a2 from a1
T = - inv(S3) * S2';   
% reduced scatter matrix 
M = S1 + S2 * T;   
% premultiply by inverse matrix of C1
M = [M(3, :) ./ 2; - M(2, :); M(1, :) ./ 2]; 
% solve eigensystem
[evec, eval] = eig(M);    
% evaluate a'Ca (here a' denote the transpose of a)
cond = 4 * evec(1, :) .* evec(3, :) - evec(2, :) .^ 2; 
% eigenvector for min. pos. eigenvalue
a1 = evec(:, cond > 0);    
% ellipse coefficients
coef = [a1; T * a1];    

% Transform the conic for the general form to the implicit form.
% All the coefficients have to be reals.
if ~isreal(coef) || isempty(coef)
    v=[];
    return
end

% Coefficients of the general form of the ellipse
A=coef(1); C=coef(3); B=coef(2);
D=coef(4); E=coef(5); F=coef(6);

% If exist the term xy (rectangle term) we rotate the coordinates axes.
if B~=0
    if A==C
        theta=pi/4;
    else
        % Conic orientation
        theta = atan2(B,A-C)/2 - sign(B)*pi/2;
    end
    % Obtain the scale parameter of the ellipse.
    ct = cos(theta);
    st = sin(theta);
    
    a = A*ct*ct + B*ct*st + C*st*st;
    c = A*st*st - B*ct*st + C*ct*ct;
    d = D*ct + E*st;
    e = -D*st + E*ct;
    f = F;
else
    % If don't exist the term "xy", the rotation not is necessary.
    a=A; c=C; d=D; e=E; f=F;
    theta=pi/2;
    ct = cos(theta);
    st = sin(theta);
end
% Estimate the coefficients of the implicit form of the ellipse.

T=((d^2)/(4*a)) + ((e^2)/(4*c)) - f;

ax = (-d/(2*a))*ct + (e/(2*c))*st;
by = (-d/(2*a))*st - (e/(2*c))*ct;
a=sqrt(T/a);
b=sqrt(T/c);

v=[a,b,ax,by,theta];


if v(2) > v(1),
    v = [b, a, v(3:4), v(end)+(pi/2)];
end
