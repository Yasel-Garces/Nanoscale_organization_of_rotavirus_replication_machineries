% Principal scrip to read and process all images in one specific hour
% post infection
% INPUT: 
%      pathname: Cell array that contain the path to read the images for 
%                experimental condition.
%      central_protein: Is a vector of lenght = length(pathname). A 
%                component is 1 if the "nameCentralProtein" is red and 
%                zero in other case.
%      nameCentralProtein: Name of the reference protein.
%      type: file type of the image
% OUTPUT: The results will be saved in a .csv file with the structure:
%   Column 1: Distance between the distribution of both proteins.
%   Column 2: Radius of the circumference that adjust the central protein.
%   Column 3: Radius of the circumference that adjust the other protein.
% AUTHOR: Yasel Garces (88yasel@gmail.com)
% -------------------------------------------------------------------------
% Clear all the variables and clean the workspace 
clear all
clc

% It specifies the directory where the images are
pathname{1} = '/home/yasel/TRABAJO/IBt/Viroplasms/Images Paper/NSP2rojo-VP4verde';
pathname{2} = '/home/yasel/TRABAJO/IBt/Viroplasms/Images Paper/NSP2rojo-VP6verde';
pathname{3} = '/home/yasel/TRABAJO/IBt/Viroplasms/Images Paper/NSP2rojo-VP760verde';
pathname{4} = '/home/yasel/TRABAJO/IBt/Viroplasms/Images Paper/NSP2rojo-VP7159verde';
pathname{5} = '/home/yasel/TRABAJO/IBt/Viroplasms/Images Paper/NSP4rojo-NSP2verde';
pathname{6} = '/home/yasel/TRABAJO/IBt/Viroplasms/Images Paper/NSP5rojo-NSP2verde';
pathname{7} = '/home/yasel/TRABAJO/IBt/Viroplasms/Images Paper/NSP2rojo-VP1verde';
pathname{8} = '/home/yasel/TRABAJO/IBt/Viroplasms/Images Paper/NSP2rojo-VP2verde';
% Central protein declaration vector (=1 if central protein is red and =0 is not)
central_protein=[1,1,1,1,0,0,1,1];

% Central protein (this is the reference protein)
nameCentralProtein='NSP2';
type='*.tif'; % image extension
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% For each directory (experiment between two differents proteins) and each
% image in the directory, run the algorithm VP-DLSFC.
for j=1:length(pathname)
    % Get all the files in the pathname{j} directory
    files=get_list_files(pathname{j},type);
    
    % Variables declaration
    format shortg
    Names=[];
    distancia=[]; ratioCentralP=[]; ratioOther=[];
    
    % Load one by one image and do the processing through the algorithm
    % VP-DLSFC.
    for i=1:length(files)
        nameImg=files(i);
        
        % Load the green and red channel of the image.
        image_red=imread(fullfile(pathname{j}, char(nameImg)),1);
        image_green=imread(fullfile(pathname{j}, char(nameImg)),2);
        % Find the pixels with values different to zero
        [x_red, y_red]=find(image_red~=0);
        [x_green, y_green]=find(image_green~=0);
        
        if central_protein(j)==1
            % If the reference protein is labeled in red:
            % Circumference adjustment
            CentralP=fit_circumference_LSFC(y_red,x_red);
            % Fit the other protein
            centerCentralP=[CentralP(3)-mean(y_red) CentralP(2)-mean(x_red)];
            R_other=closeNSP2(centerCentralP,x_green,y_green);
        else
            % If the reference protein is labeled in green:
            % Circumference adjustment
            CentralP=fit_circunsference_LSFC(y_green,x_green);
            % Fit the other protein
            centerCentralP=[CentralP(2)-mean(y_green) CentralP(3)-mean(x_green)];
            R_other=closeNSP2(centerCentralP,x_red,y_red);
        end
        
        % Distance between the circumferences
        distancia(i)= R_other-CentralP(1);
        ratioCentralP(i)=CentralP(1);
        ratioOther(i)=R_other;
    end
    % Save all the information in a dataframe
    Data=[distancia' ratioCentralP' ratioOther'];
    
    % Take the name of the directory
    NameFile=regexp(pathname{j}, '/', 'split');
    NameFile=char(NameFile(end-1));
    
    % Write results
    pathname='/home/yasel/Dropbox/Paper Viroplasmas/Programs/R/ResultsCSV/NSP2/';
    csvwriteh(strcat(pathname,strcat(NameFile,'.csv')),...
        Data,  {'Distance',strcat('ratio',nameCentralProtein),'ratioOther'})
end

function out=get_list_files(path,type)
% Get all files
list_dir=dir(fullfile(path,type));
% Extract just the names of the files.
out={list_dir.name};
end
