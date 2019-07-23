% It computes the mean and the standard deviation of the distance of 
% a set of points to a circumference.
% INPUT: 
%      circ: circumference representation as [radius center_x center_y]
%    points: set of points as a matrix of n x 2 ([X Y]).
% OUTPUT: [mean_circle,std_circle] --> Vector with the mean distance and
% the standard deviation of the points to the circumference.
% AUTHOR: Yasel Garces (88yasel@gmail.com)
function [mean_circle,std_circle]=circle_std(circ,points)
% Compute the distance
distance=abs(circ(1)^2 - (points(:,1)-circ(2)).^2-(points(:,2)-circ(3)).^2);
% Compute the mean
mean_circle=mean(distance/100);
% Compute the standard deviation
std_circle=std(distance/100);