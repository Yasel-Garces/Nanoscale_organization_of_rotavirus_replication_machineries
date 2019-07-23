% The algorithm DLSFC was compared with other two approaches proposed 
% by Gander et al. (Gander et al., 1994). For a more comprehensive and 
% clear comparison we named these others two methods as 
% "Algebraic Least Square Fitting Circle" (ALSFC) and 
% "Geometric Least Square Fitting Circle" (GLSFC).
% For details see the section "Algorithm Validation" in 
% "Nanoscale organization of rotavirus replication machineries", ELife.
% INPUT: None
% OUTPUT: 
%        SimulationResult.csv: Table with all the information about the 
%                              generated images.
%        Validation_Images: Directory with all the generated images for the
%        validation
% AUTHOR: Yasel Garces (88yasel@gmail.com)
%--------------------------------------------------------------------------
% Create an empty Image (512x512)
image=512*zeros(512);

% Generate 500 random circumferences
R=randi([50,100],1,500); % Radius
CX=randi([100,350],1,500); % x-center
CY=randi([100,350],1,500); % y-center

% Declarate variables
DLSFC=[]; ALSFC=[]; GLSFC=[];
ID=[]; Angle=[]; DistanceStd=[]; DistanceMean=[];

for i=1:length(R)
    % Generation of circle's points and consider partial occlusion
    general_angle=[2*pi 3*pi/2 pi pi/2]; % [360 270 180 90]<<-degrees
    for j=1:length(general_angle)
        angle=general_angle(j); % Generate a specific partial occlusion angle.
        Points=draw_ellipse([R(i) R(i) CX(i) CY(i) 0],[0 angle]);
        % Add gaussian noise to the Points
        general_noise = linspace(0, 10, 20);
        for t=1:length(general_noise)
            noise=general_noise(t);
            n1 = randn(1, numel(Points(:,1))); % noise with mean=0 and std=1;
            Points(:,2) = Points(:,2) + noise*n1'; % Corrupted points
            
            % Computhe the standard deviation of the distance
            [mean_circle(t),std_circle(t)]=circle_std([R(i) CX(i) CY(i)],Points);
            % Adjust circumference
            % DLSFC
            DLSFC_fit(t,:) = fit_circumference_LSFC(Points(:,1),Points(:,2));
            % Find the linear least squares fit
            [ALSFC_zl(t,:), ALSFC_rl(t)] = fitcircle(Points, 'linear');
            % Find the best geometric fit
            [GLSFC_z(t,:), GLSFC_rl(t)] = fitcircle(Points);
            
            % And plot the results
            if  t==5 || t==10 || t==15 || t==20
                x = linspace(0, 2*pi, 100);
                plot(Points(:,1), Points(:,2), '.', ...
                    CX(i) + R(i)* cos(x), CY(i)+R(i)*sin(x),'black',...
                    ALSFC_zl(t,1) + ALSFC_rl(t) * cos(x), ALSFC_zl(t,2) + ALSFC_rl(t) * sin(x), 'b--', ...
                    GLSFC_z(t,1)  + GLSFC_rl(t)  * cos(x), GLSFC_z(t,2)  + GLSFC_rl(t) * sin(x), 'g',...
                    DLSFC_fit(t,2) + DLSFC_fit(t,1)* cos(x), DLSFC_fit(t,3)+DLSFC_fit(t,1)*sin(x),'r',...
                    'markers',20,'LineWidth',2.5)
                axis equal off
                set(gcf,'color','w');
                %                 axis([0 512 0 512])
                %legend('Location','best','Data points','Ground Truth','ALSFC','GLSFC','DLSFC')
                Name=strcat('Validation_Images/ImageFit_','Mean_',num2str(mean_circle(t)),...
                    'Std_',num2str(std_circle(t)),'Angle_',num2str(angle),'.pdf');
                print('-bestfit',Name,'-dpdf','-r1500');
                close all
            end
        end
        % Error -Euclidean distance R^3 -
        DLSFC=[DLSFC; sqrt((sum(([R(i) CX(i) CY(i)]-DLSFC_fit).^2,2)))/100];
        ALSFC=[ALSFC; sqrt((sum(([R(i) CX(i) CY(i)]-[ALSFC_rl', ALSFC_zl]).^2,2)))/100];
        GLSFC=[GLSFC; sqrt((sum(([R(i) CX(i) CY(i)]-[GLSFC_rl', GLSFC_z]).^2,2)))/100];
        
        DistanceMean=[DistanceMean; mean_circle'];
        DistanceStd=[DistanceStd; std_circle'];
        
        % Represent the occlusion angle (for data frame)
        Angle=[Angle; repmat(general_angle(j),length(general_noise),1)];
    end
    ID=[ID; repmat(i,length(general_noise)*length(general_angle),1)];
end

% Concatenate variables
ID=[ID;ID;ID];
Angle=[Angle;Angle;Angle];
DistanceMean=[DistanceMean;DistanceMean;DistanceMean];
DistanceStd=[DistanceStd;DistanceStd;DistanceStd];
Algorithm=[repmat('DLSFC',length(DLSFC),1);repmat('ALSFC',length(ALSFC),1);...
    repmat('GLSFC',length(GLSFC),1)];
Error=[DLSFC; ALSFC; GLSFC];

% Create table of data
T=table(ID,Algorithm,Angle,DistanceMean,DistanceStd,Error);
% Write and save data
writetable(T,'SimulationResult.csv')