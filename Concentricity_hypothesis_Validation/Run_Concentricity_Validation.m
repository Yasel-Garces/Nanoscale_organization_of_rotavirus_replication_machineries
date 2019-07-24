% Program to validate the hypothesis that the central and the accompainy 
% protein are concentric.
% In each experiment we have a couple of proteins (normally NSP2 as
% a central protein). For each combination (for example NSP-NSP4 
% and NSP2-VP6) we adjust the least square circle (algorithm DLSFC) 
% to each protein independently.
% INPUT: 
%      pathname: Cell array that contain the path to read the images for 
%                the validation.
%      otherProteinName: Name of the others proteins.
%      central_protein: Is a vector of lenght = length(pathname). A 
%                component is 1 if NSP2 is red and zero in other case.
%      type: file type of the image
% OUTPUT: The results will be saved in a DistanceCenterProtein.csv file 
%         with the structure:
%   Column 1 (Combination): Name of the protein combined with NSP2
%   Column 2 (Distance): Distance between the center of the adjusted 
%                      circumferences to NSP2 and the accompanying protein.
% AUTHOR: Yasel Garces (88yasel@gmail.com)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% It specifies the directory where the images are
pathname{1} = '/home/yasel/TRABAJO/IBt/Viroplasms/Validation/NSP2rojo-VP4verde/';
pathname{2} = '/home/yasel/TRABAJO/IBt/Viroplasms/Validation/NSP2rojo-VP6verde/';
pathname{3} = '/home/yasel/TRABAJO/IBt/Viroplasms/Validation/NSP2rojo-VP760verde/';
pathname{4} = '/home/yasel/TRABAJO/IBt/Viroplasms/Validation/NSP2rojo-VP7159verde/';
pathname{5} = '/home/yasel/TRABAJO/IBt/Viroplasms/Validation/NSP4rojo-NSP2verde/';
pathname{6} = '/home/yasel/TRABAJO/IBt/Viroplasms/Validation/NSP5rojo-NSP2verde/';
pathname{7} = '/home/yasel/TRABAJO/IBt/Viroplasms/Validation/NSP2rojo-VP1verde/';
pathname{8} = '/home/yasel/TRABAJO/IBt/Viroplasms/Validation/NSP2rojo-VP2verde/';

% Name of the others proteins.
otherProteinName=["VP4","VP6","VP760","VP7159","NSP4","NSP5","VP1","VP2"];
% This variable is 1 of NSP2 is red and zero in other case.
central_protein=[1,1,1,1,0,0,1,1];
% Type of file
type='*.tif';
%--------------------------------------------------------------------------
% Variables Declaration
Distance=[];
Combination=[];

% For each directory (experiment between two different proteins) 
% and each image in the directory, run the algorithms DLSFC and 
% adjust the distribution of both proteins independently.
for j=1:length(pathname)
    % Load files in the directory
    files=get_list_files(pathname{j},type);
    % For each image adjust a circumference to each protein.
    for i=1:length(files)
        % Image name
        nameImg=files(i);     
        % Load the green and red channel.
        image_red=imread(fullfile(pathname{j}, char(nameImg)),1);
        image_green=imread(fullfile(pathname{j}, char(nameImg)),2);
        % Find the pixels with values different to zero
        [x_red, y_red]=find(image_red~=0);
        [x_green, y_green]=find(image_green~=0);
        
        if central_protein(j)==1
            % If the reference protein is labeled in red:
            % Fit a circumference through DLSFC
            CentralP=fit_circumference_LSFC(x_red,y_red);
            OtherP=fit_circumference_LSFC(x_green,y_green);
        else
            % If the reference protein is labeled in green:
            OtherP  =fit_circumference_LSFC(x_red,y_red);
            CentralP=fit_circumference_LSFC(x_green,y_green);
        end
        % Compute the distance between the centers of each circumference
        this_Distance(i,:)=sqrt((CentralP(2)-OtherP(2))^2 + (CentralP(3)-OtherP(3))^2)/100;
    end
    % Save the data
    Distance=[Distance; this_Distance];
    % Save the iformation about the protein combination
    this_Combination=strcat(otherProteinName(j));
    Combination=[Combination; repmat(this_Combination,length(this_Distance(:,1)),1)];  
end
% Create a table with all data
T=table(Combination,Distance);
% Write and save data
writetable(T,'DistanceCenterProtein.csv')

%--------------------------------------------------------------------------
% Nested functions
function out=get_list_files(path,type)
% Get all files
list_dir=dir(fullfile(path,type));
% Extract just the names of the files.
out={list_dir.name};
end