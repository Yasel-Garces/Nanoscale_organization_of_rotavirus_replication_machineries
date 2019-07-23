% Program to validate the hypothesis that the central and accompainy protein are concentric.
% In each experiment we have a couple of proteins (normally NSP2 as
% a central protein). For each combination (for example NSP-NSP4
% and NSP2-VP6) we adjust a least square circle (algorithm DLSFC) and save
% the information relative to the segmentation.
% Author: Yasel Garces

%------------------------------------------------------------
% Load directory
pathname{1} = '/media/yasel/e918b490-ed8b-4ab6-9b85-b300970fb84c/Yasel/TRABAJO/Doctorado/Viroplasmas/Imagenes/PreSeg_Images/SR/Validation/NSP2rojo-dsRNAverde/';
pathname{2} = '/media/yasel/e918b490-ed8b-4ab6-9b85-b300970fb84c/Yasel/TRABAJO/Doctorado/Viroplasmas/Imagenes/PreSeg_Images/SR/Validation/NSP2rojo-PDIverde/';
pathname{3} = '/media/yasel/e918b490-ed8b-4ab6-9b85-b300970fb84c/Yasel/TRABAJO/Doctorado/Viroplasmas/Imagenes/PreSeg_Images/SR/Validation/NSP2rojo-VP4verde/';
pathname{4} = '/media/yasel/e918b490-ed8b-4ab6-9b85-b300970fb84c/Yasel/TRABAJO/Doctorado/Viroplasmas/Imagenes/PreSeg_Images/SR/Validation/NSP2rojo-VP6verde/';
pathname{5} = '/media/yasel/e918b490-ed8b-4ab6-9b85-b300970fb84c/Yasel/TRABAJO/Doctorado/Viroplasmas/Imagenes/PreSeg_Images/SR/Validation/NSP2rojo-VP760verde/';
pathname{6} = '/media/yasel/e918b490-ed8b-4ab6-9b85-b300970fb84c/Yasel/TRABAJO/Doctorado/Viroplasmas/Imagenes/PreSeg_Images/SR/Validation/NSP2rojo-VP7159verde';
pathname{7} = '/media/yasel/e918b490-ed8b-4ab6-9b85-b300970fb84c/Yasel/TRABAJO/Doctorado/Viroplasmas/Imagenes/PreSeg_Images/SR/Validation/NSP4rojo-NSP2verde/';
pathname{8} = '/media/yasel/e918b490-ed8b-4ab6-9b85-b300970fb84c/Yasel/TRABAJO/Doctorado/Viroplasmas/Imagenes/PreSeg_Images/SR/Validation/NSP5rojo-NSP2verde/';
pathname{9} = '/media/yasel/e918b490-ed8b-4ab6-9b85-b300970fb84c/Yasel/TRABAJO/Doctorado/Viroplasmas/Imagenes/PreSeg_Images/SR/Validation/NSP2rojo-VP1verde/';
pathname{10} = '/media/yasel/e918b490-ed8b-4ab6-9b85-b300970fb84c/Yasel/TRABAJO/Doctorado/Viroplasmas/Imagenes/PreSeg_Images/SR/Validation/NSP2rojo-VP2verde/';


% Central protein, 1= red, 0=green
nameCentralProtein="NSP2";
otherProteinName=["dsRNA","PDI","VP4","VP6","VP760","VP7159","NSP4","NSP5","VP1","VP2"];

central_protein=[1,1,1,1,1,1,0,0,1,1];
% Type of file
type='*.tif';

% Declarate loop variables
Distance=[];
Combination=[];

% Run the rutine of fit the ellipse for each directory
for j=1:length(pathname)
    % Load files in the directory
    files=get_list_files(pathname{j},type);
    % For each image we adjust a ellipse to each protein.
    for i=1:length(files)
        % Image name
        nameImg=files(i);
        
        % Load the images in the green and red channel.
        image_red=imread(fullfile(pathname{j}, char(nameImg)),1);
        image_green=imread(fullfile(pathname{j}, char(nameImg)),2);
        
        [x_red, y_red]=find(image_red~=0);
        [x_green, y_green]=find(image_green~=0);
        
        if central_protein(j)==1
            % Ellipses adjustment
            CentralP=fit_circunsference_LSFC(x_red,y_red);
            OtherP=fit_circunsference_LSFC(x_green,y_green);
        else
            % Ellipses adjustment
            OtherP  =fit_circunsference_LSFC(x_red,y_red);
            CentralP=fit_circunsference_LSFC(x_green,y_green);
        end
        this_Distance(i,:)=sqrt((CentralP(2)-OtherP(2))^2 + (CentralP(3)-OtherP(3))^2)/100;
    end
    Distance=[Distance; this_Distance];
    
    this_Combination=strcat(otherProteinName(j));
    Combination=[Combination; repmat(this_Combination,length(this_Distance(:,1)),1)];  
end
% Create a table with all data
T=table(Combination,Distance);
% Write and save data
writetable(T,'DistanceCenterProtein.csv')