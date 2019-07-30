% To proof the hypothesis that the viral elements of the VPs can be 
% approximated through circumferences, we carried up a series of 
% experiments based on the comparison between the circumference obtained 
% by the algorithm DLSFC, and a least squared ellipse resulting of the 
% ''Direct Least Square Fitting Ellipse'' (DLSFE) algorithm 
% (Fitzgibbon A, Pilu M, Fisher B. Direct Least Squares Fitting of 
%  Ellipses. IEEE Transactions on Pattern Analysis and Machine 
%  Intelligence. 1999; 21(5). doi: 10.1109/ICPR.1996.546029.).
% See the article 
% Garcés et al. Nanoscale organization of rotavirus replication machineries. 
% eLife 2019;8:e42906. https://elifesciences.org/articles/42906, 
% doi: 10.7554/eLife.42906.
% INPUT: 
%      pathname: Cell array that contain the path to read the images for 
%                the validation.
%      reference_protein: Name of the reference protein.
%      otherProteinName: Name of the others proteins.
%      central_protein: Is a vector of lenght = length(pathname). A 
%                component is 1 if the "reference_protein" is red and 
%                zero in other case.
%      type: file type of the image
% OUTPUT: The results will be saved in a Result.csv file with the structure:
%   Column 1 (Combination): Name of the protein combined with NSP2
%   Column 2 (Protein): Adjusted proteins through the algorithms DLSFC and
%                       DLFE
%   Column 3 (MajorSemiAxis): Major semi-axis of the ellipse
%   Column 4 (MinorSemiAxis): Minor semi-axis of the ellipse
%   Column 5 (RadiusCirc): Circumference radius
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

% Name of the reference protein.
reference_protein="NSP2";
% Name of the others proteins.
otherProteinName=["VP4","VP6","VP760","VP7159","NSP4","NSP5","VP1","VP2"];
% Vector of lenght = length(pathname). A component is 1 if the 
% "reference_protein" is red and zero in other case.
central_protein=[1,1,1,1,0,0,1,1];
% Type of file
type='*.tif';
%--------------------------------------------------------------------------
% Variables Declaration
Combination=[]; Protein=[];
MajorSemiAxis=[]; MinorSemiAxis=[];
RadiusCirc=[];

% For each directory (experiment between two differents proteins) and each
% image in the directory, run the algorithms DLSFC and DLSFE.
for j=1:length(pathname)
    % Load files in the directory
    files=get_list_files(pathname{j},type);
    
    % For each image we adjust a ellipse and a circumference 
    % to each protein.
    % Declaration of the storage variables
    Circumference_Central=[];Circumference_Other=[];
    Ellipse_CentralP=[];Ellipse_Other=[];
    for i=1:length(files)
        % Image name
        nameImg=files(i);
        % Load the green and red channel of the image.
        image_red=imread(fullfile(pathname{j}, char(nameImg)),1);
        image_green=imread(fullfile(pathname{j}, char(nameImg)),2);
        % Find the pixels with values different to zero
        [x_red, y_red]=find(image_red~=0);
        [x_green, y_green]=find(image_green~=0);
        
        if central_protein(j)==1
            % If the reference protein is labeled in red:
            % Fit ellipse
            Ellipse_CentralP(i,:)=fit_ellipse_LSFE(x_red,y_red);
            Ellipse_Other(i,:)=fit_ellipse_LSFE(x_green,y_green);
            % Fit circumference
            Circumference_Central(i,:) = fit_circumference_LSFC(x_red,y_red);
            Circumference_Other(i,:) = fit_circumference_LSFC(x_green,y_green);
        else
            % If the reference protein is labeled in green:
            % Fit ellipse
            Ellipse_Other(i,:)  =fit_ellipse_LSFE(x_red,y_red);
            Ellipse_CentralP(i,:)=fit_ellipse_LSFE(x_green,y_green);
            % Fit circumference
            Circumference_Other(i,:) = fit_circumference_LSFC(x_red,y_red);
            Circumference_Central(i,:) = fit_circumference_LSFC(x_green,y_green);
        end    
    end
    
    % Extract the name of the reference & accompanying protein
    NameReferenceProtein=repmat(reference_protein,length(Ellipse_CentralP(:,1)),1);
    NameOtherProtein=repmat(otherProteinName(j),length(Ellipse_CentralP(:,1)),1);
    
    % Concatenate the data
    Protein=[Protein; [string(NameReferenceProtein); string(NameOtherProtein)]];
    MajorSemiAxis=[MajorSemiAxis; [Ellipse_CentralP(:,1); Ellipse_Other(:,1)]];
    MinorSemiAxis=[MinorSemiAxis; [Ellipse_CentralP(:,2); Ellipse_Other(:,2)]];
    % Radius Circumference
    RadiusCirc=[RadiusCirc; [Circumference_Central(:,1); Circumference_Other(:,1)]];
    
    this_Combination=strcat(reference_protein,"_",otherProteinName(j));
    this_Combination=otherProteinName(j);
    Combination=[Combination; repmat(this_Combination,2*length(Ellipse_CentralP(:,1)),1)];
end
% Convert to microns
MajorSemiAxis=MajorSemiAxis/100;
MinorSemiAxis=MinorSemiAxis/100;
RadiusCirc=RadiusCirc/100;

% Create a table with all data
T=table(Combination,Protein,MajorSemiAxis,MinorSemiAxis,RadiusCirc);
% Save data in a csv file
writetable(T,'Result.csv')

%--------------------------------------------------------------------------
% Nested functions
function out=get_list_files(path,type)
% Get all files
list_dir=dir(fullfile(path,type));
% Extract just the names of the files.
out={list_dir.name};
end
