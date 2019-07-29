% Runs the last step of the algorithm VP-DLSFC. 
% The center of the reference protein is taken 
% as the center of the accompanying protein, 
% and then the radius of the circumference 
% for this second protein is computed.
% INPUT:
%      center: Center of the adjusted circumference to the reference 
%              protein computed through the algorithm DLSFC.
%      x,y: Points of the distribution of the accompanying protein.
% OUTPUT: R: Radius of the adjusted circumference.
% AUTHOR: Yasel Garces(88yasel@gmail.com)

function R=radius_accompanying_protein(center,x,y)
% Data normalization
u=x-mean(x);
v=y-mean(y);

% Simplification (notation)
% Eq. 9 in the Appendix 1
S_uu=sum(u.^2);
S_vv=sum(v.^2);

% Compute the radius 
R2=sum(center.^2)+((S_uu+S_vv)/length(x));
R=sqrt(R2);
