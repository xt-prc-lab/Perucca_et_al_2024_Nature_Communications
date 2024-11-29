%%%% AIM %%%%
% Compute the angle toward the edge formed by the fibers that ctFIRE found

%%%% INPUTS %%%%
% CAF image to draw the edge
% ECM images if needed
% The corresponding ctFIRE outputs from these images

%%%% OUTPUTS %%%%
% 'Edgewidth.mat' containing the width of the CAF edge for each analyzed
% folder
% 'MEASURE.mat' containing the fibers angle for each region (IN, EDGE, OUT)
% of each analyzed folder

% Clément Hallopeau 05/2021

%%

% Directory of the codes ("address" of the folder containing them)
codepath =  '/home/challopeau/Documents/PhD/Codes/Paper_MIRO/Matlab_for_ECMV3.2_forAnna_last';

% Directory of the folder containing all the subfolders
path = uigetdir('', 'Select the folder containing the subfolders') ;
path = '/home/challopeau/Downloads/CAF200_COL_1_1/CAF200_COL_1/' ;

folders = dir(path) ; folders = folders(3:end)   ; folders = folders([folders.isdir]) ;% exclude the first 2 folders/files which are always '.' and '..'


%% Draw the Edge of each spheroid 
% Draw a polygon passing along the Fibroblast edge and along the edges of the image inside the spheroid.
% Before closing, move the 4 corners outside of the image

for i = 1:length(folders)
    cd(codepath)

    % Definition of px_to_micro depending on the size of the image
    px_size = 0.3225 ; % pixel size, in microns per pixel, for 2048*2048 field of views
    px_to_micro = 0.3225 * 2048/size(imread([path,filesep,  folders(i).name, filesep, 'CAF.tif']),1) ;   %convert pixel size considering the image binning

    sigma = 50 ; % parameter of the Gaussian filter used to smooth the edge, and avoid sharp angles
    DefineEdge([path, filesep, folders(i).name], sigma, px_to_micro) ;   
end
disp('Edges : DONE')


%% Compute the intensity profiles of each ECM component & CAF

for i = 1:length(folders)
    cd(codepath)
    Window_size = 6 ;  % smoothing parameter, the larger the smoother, see detail in analyze_intensity.mat

    %disp(folders(i).name)
    analyze_intensityECM([path, filesep, folders(i).name], codepath, Window_size) ;   
end
disp('Intensity Profiles : DONE')


%% Compute the position (start/end) of each ECM peak 

for i = 1:length(folders)
    cd(codepath)
    disp(folders(i).name)
    ECMbarrier([path, filesep, folders(i).name]) ;   
end
disp('ECM Peaks : DONE')

%% Compute the Width of the CAF barrier

% The width of the edge is obtained by finding the change points of the  intensity profile. 
% These changepoints usually corresond to the beginning, the top and the end of the 
% intensity peak due to the accumulation of ECM.

Edgewidth = table({}, [], [], [], 'VariableNames', {'Folder', 'Width', 'Start', 'End'}) ; 

for i = 1:length(folders)
    cd([path, filesep, folders(i).name,filesep,'Matlab_data']) ;
    load('PROFILE.mat') ; 
    ipt = findchangepts(PROFILE.CAF_MEAN, 'Statistic', 'linear', 'MaxNumChanges',3) ;;; 
    Edgewidth.Folder{i} = string(folders(i).name) ;
    % since the peak usually increase abruptly, we consider the starting
    % point as the 2nd changepoint
    Edgewidth.Width(i) = PROFILE.R(ipt(end)) - PROFILE.R(ipt(2)) ; ; Edgewidth.Start(i) = ipt(2) ; Edgewidth.End(i) = ipt(end) ;  
end

Edgewidth.Folder = string(Edgewidth.Folder) ;

cd(path) ; save('Edgewidth.mat', 'Edgewidth') ;


%% calculate angle map

path = '/home/challopeau/Downloads/RFN39_COL_1/RFM39_COL_EXP015_Pos12' ;
cd(codepath) ; 

compute_angle_map(path, codepath )

%% Compute the orientation of ECM fibers within each region (in, edge, out)

for i =1:length(folders)
    cd(codepath)
    disp(folders(i).name)
    compute_angle_map = 0 ;  %1 if computation & smoothing of angle map, 0 if the map must be loaded
    Fibers_angles([path, filesep, folders(i).name], codepath, compute_angle_map)  
end
disp('ECM Orientation : DONE')


%% Compute the profile of angles as a function to the distance to the CAF edge
for i =1:length(folders)
    cd(codepath)
    disp(folders(i).name)
    Window_size = 40 ;
    angle_profile([path, filesep, folders(i).name], Window_size) ;
end
disp('Angles profiles : DONE')


%% merge the intensity profiles (PROFILE.mat) and the angle profile (Angle_profile.csv) in one single csv file

for i =1:length(folders)

temp_path = [path, filesep, folders(i).name, filesep, 'Matlab_data'] ; 
load([temp_path, filesep, 'PROFILE.mat']) ; 
angles = readtable([temp_path, filesep, 'Angle_profile.csv']) ; 

col_names = PROFILE.Properties.VariableNames ; 
col_types = cell(1, length(col_names)) ; 
for j = 1:length(col_names) 
 col_types{j} = class(PROFILE.(col_names{j})) ; 
end

new_table = outerjoin(PROFILE(:,find(~contains(col_types, 'cell'))),angles,'Keys','R','MergeKeys', 1) ;
writetable(new_table,[temp_path, filesep, 'Intensity_and_angle_profiles.csv']) ;
end


%% export MEASURE table as a csv

for i = 1:length(folders)
    cd(codepath)
    disp(folders(i).name)
     cd([path, filesep, folders(i).name,filesep,'Matlab_data']) ;
    load('MEASURE') ;
        
pos = find(~cellfun(@isempty,strfind(string(MEASURE.Properties.VariableNames), "THETA")) ==1) ; 

m= 0 ;
for i = 1:length(pos)
    if m < max(cellfun(@length,MEASURE.(MEASURE.Properties.VariableNames{pos(i)})))
        m = max(cellfun(@length,MEASURE.(MEASURE.Properties.VariableNames{pos(i)}))) ;
    end
end

MEASUREtab = ones(m, length(pos)*3)*NaN ;

s=1;; colnames = {} ; pat = ["IN", "EDGE", "OUT"] ; ;;
for i = 1:length(pos)
    for j = 1:3
        colnames = horzcat(colnames, { [MEASURE.Properties.VariableNames{pos(i)}, '_', char(pat(j))] } ) ; 
        MEASUREtab(:,s) = [MEASURE.(MEASURE.Properties.VariableNames{pos(i)}){j}', ones(1,m-size(MEASURE.(MEASURE.Properties.VariableNames{pos(i)}){j},1))*NaN]' ;
        s=s+1 ;
    end
end
   

MEASUREtab = array2table(MEASUREtab, 'VariableNames', colnames) ;

cd(path) ;
writetable(MEASUREtab, 'MEASURE.csv') ;

end

%% For ALL subfolders : merge angle profiles for the different ECM components


% here path is a cell array containing all the "main folders" you wanna use
% (1 "main folder" contains several "subfolders)
path = { '/home/yourusername/Documents/YourCodeFolder/Codes'} %,
path = { '/home/challopeau/Downloads/CAF200_COL_1_1/CAF200_COL_1/' } ; 


MERGE = table() ;
binx = transpose(-1000:2:1000) ;   % discretize the x-axis, in µm. i.e. we divide the x-axis (distance to spheroid center) is very small bins
% if binx = -1000:2:1000, the x-axis will range from -1000 to +1000m with a
% step of 2µm. Visually, on the plot, it's invisible. 
MERGE.('R') = binx ;

% enter each "main folder"
for ii = 1:length(path)
    folders = dir(path{ii}) ; folders = folders([folders.isdir]) ;  % we extract the name of the subfolders
% Then we enter each subfolder
    for i = 3:length(folders)
        % we enter in the "Matlab_data" folder, and load the intensity profiles all the
        % images of the subfolder ('PROFILE.mat')
        cd(strcat(path{ii}, '/',folders(i).name, '/Matlab_data')) ;

        data = readtable('Angle_profile.csv') ; 

        cd([path{ii},'/',folders(i).name]) ; imgs = dir('*.tif') ;

        % Then, we take each image of the subfolder
        for j =1:length(imgs)
            % If the columns of the MERGE dataframe doesn't contain the
            % name of the given image (imgs(j).name), we create the column
            if isempty(cell2mat(  strfind(MERGE.Properties.VariableNames, imgs(j).name(1:end-4))  )) 
                MERGE.(imgs(j).name(1:end-4)) = cell(length(binx),1) ; ; 
            end
        end


                for j =1:length(imgs)
            % Then, we take each tiny bin (binx) 1by1, and try to find the
            % values of PROFILE.R belonging to each bin, and their
            % associated intensities  = find(PROFILE.R-px_to_micro*center(4) < binx(k) & PROFILE.R-px_to_micro*center(4) > binx(k-1))
            for k = 2:length(binx)  % assign values to each tiny bin
               % MERGE.(imgs(j).name(1:end-4)){k-1} = cat(2, MERGE.(imgs(j).name(1:end-4)){k-1}, transpose(PROFILE.([imgs(j).name(1:end-4),'_MEAN'])(find(PROFILE.R < binx(k)+floor(center(4)*px_to_micro)+1 & PROFILE.R > binx(k-1)+floor(center(4)*px_to_micro)+1)))) ;;
               % When we find values of R comprised within the given bin
               % (binx(k) -  binx(k-1)), we append the values of intensity
%                MERGE.(imgs(j).name(1:end-4)){k-1} = cat(2, MERGE.(imgs(j).name(1:end-4)){k-1}, transpose(PROFILE.([imgs(j).name(1:end-4),'_MEAN'])(find(PROFILE.R < binx(k) & PROFILE.R > binx(k-1)))  )) ;;
%                 MERGE.([imgs(j).name(1:end-4), '_', num2str(i-2)]){k-1} = mean(PROFILE.([imgs(j).name(1:end-4),'_MEAN'])(find(PROFILE.R < binx(k) & PROFILE.R > binx(k-1))))  ;
               MERGE.(imgs(j).name(1:end-4)){k-1} = cat(2, MERGE.(imgs(j).name(1:end-4)){k-1}, transpose(data.([imgs(j).name(1:end-4),'_MEAN'])(find(data.R < binx(k) & data.R > binx(k-1)))  )) ;;
                MERGE.([imgs(j).name(1:end-4), '_', num2str(i-2)]){k-1} = mean(data.([imgs(j).name(1:end-4),'_MEAN'])(find(data.R < binx(k) & data.R > binx(k-1))))  ;

            end
            temp =  cell2mat(MERGE.([imgs(j).name(1:end-4), '_', num2str(i-2)])) ;;
            MERGE.([imgs(j).name(1:end-4), '_', num2str(i-2)]) = [temp', nan]' ;;
        end % And we repeat that for all the images
    end % for all the subfolders
end % and for all the main folders
% /!\ /!\ So it's important to name the images in the same way everywhere ! /!\ /!\

    % For each image (COL_4, CAF, FN ...) we concatenate all the values of
% intensity for each value of x. It's now time to compute means and stds
%names = MERGE.Properties.VariableNames
for k = 1:length(imgs)
    MERGE.([imgs(k).name(1:end-4),'_MEAN']) = zeros(length(MERGE.R),1) ;  MERGE.([imgs(k).name(1:end-4),'_STD']) = zeros(length(MERGE.R),1) ; 
    for i = 1:length(MERGE.R)
        MERGE.([imgs(k).name(1:end-4),'_MEAN'])(i) = mean(MERGE.(  imgs(k).name(1:end-4)  ){i}) ; 
        MERGE.([imgs(k).name(1:end-4),'_STD'])(i) = std(MERGE.(  imgs(k).name(1:end-4)  ){i}) ; 
    end
MERGE.([imgs(k).name(1:end-4),'_MEAN'])(isnan(MERGE.([imgs(k).name(1:end-4),'_MEAN']))) = 0 ;
MERGE.([imgs(k).name(1:end-4),'_STD'])(isnan(MERGE.([imgs(k).name(1:end-4),'_STD']))) = 0 ;; 
MERGE.(imgs(k).name(1:end-4)) = [];
end

cd(path{1}) ;
writetable(MERGE, 'PROFILES_ANGLES_MERGED.csv') ;

%% For ALL subfolders : compute mean CAF intensity and orientation at the EDGE (0-50µm) and OUT (200-250µm)

% here path is a cell array containing all the "main folders" you wanna use
% (1 "main folder" contains several "subfolders)
path = { '/home/yourusername/Documents/YourCodeFolder/Codes'} %,
path = { '/home/challopeau/Downloads/CAF200_COL_1_1/CAF200_COL_1/' } ; 


% enter each "main folder"
for ii = 1:length(path)
    folders = dir(path{ii}) ; folders = folders([folders.isdir]) ;  % we extract the name of the subfolders
% Then we enter each subfolder

edge_score_all = table() ; 
    for i = 3:length(folders)
        % we enter in the "Matlab_data" folder, and load the intensity profiles all the
        % images of the subfolder ('PROFILE.mat')
        cd(strcat(path{ii}, '/',folders(i).name, '/Matlab_data')) ;

        data = readtable('Intensity_and_angle_profiles.csv') ; 

        cd([path{ii},'/',folders(i).name]) ; imgs = dir('*.tif') ;

        edge_score = table() ;
        % Then, we take each image of the subfolder
        for j =1:length(imgs) 
            % angle from 0 to 50 µm
            ind = find(contains(data.Properties.VariableNames, [imgs(j).name(1:end-4), '_MEAN_angle']));
            edge_score.([data.Properties.VariableNames{ind}, '_0-50um']) = mean(data.(data.Properties.VariableNames{ind})(find(data.R >0 & data.R < 50))) ;
            % angle from 200 to 250 µm
            edge_score.([data.Properties.VariableNames{ind}, '_200-250um']) = mean(data.(data.Properties.VariableNames{ind})(find(data.R >200 & data.R < 250))) ; 

           %intensity from 0 to 50µm, 200-250µm
            ind2 = find(contains(data.Properties.VariableNames, [imgs(j).name(1:end-4), '_MEAN']));
            ind2 = ind2(ind2 ~= ind)  ;
            edge_score.([data.Properties.VariableNames{ind2}, '_0-50um']) = mean(data.(data.Properties.VariableNames{ind2})(find(data.R >0 & data.R < 50))) ;
            edge_score.([data.Properties.VariableNames{ind2}, '_200-250um']) = mean(data.(data.Properties.VariableNames{ind2})(find(data.R >200 & data.R < 250))) ; 

                    if i==3
         %   edge_score_all = table({}, [], [], [], 'VariableNames', {'Folder', 'Width', 'Start', 'End'}) ; 
                edge_score_all.('folder') = {folders(i).name} ;
                edge_score_all.([data.Properties.VariableNames{ind}, '_0-50um']) = mean(data.(data.Properties.VariableNames{ind})(find(data.R >0 & data.R < 50))) ;
                edge_score_all.([data.Properties.VariableNames{ind}, '_200-250um']) = mean(data.(data.Properties.VariableNames{ind})(find(data.R >200 & data.R < 250))) ; 
                edge_score_all.([data.Properties.VariableNames{ind2}, '_0-50um']) = mean(data.(data.Properties.VariableNames{ind2})(find(data.R >0 & data.R < 50))) ;
                edge_score_all.([data.Properties.VariableNames{ind2}, '_200-250um']) = mean(data.(data.Properties.VariableNames{ind2})(find(data.R >200 & data.R < 250))) ; 
        else
                edge_score_all.('folder'){i-2} = folders(i).name ;
                edge_score_all.([data.Properties.VariableNames{ind}, '_0-50um'])(i-2) = mean(data.(data.Properties.VariableNames{ind})(find(data.R >0 & data.R < 50))) ;
                edge_score_all.([data.Properties.VariableNames{ind}, '_200-250um'])(i-2) = mean(data.(data.Properties.VariableNames{ind})(find(data.R >200 & data.R < 250))) ; 
                edge_score_all.([data.Properties.VariableNames{ind2}, '_0-50um'])(i-2) = mean(data.(data.Properties.VariableNames{ind2})(find(data.R >0 & data.R < 50))) ;
                edge_score_all.([data.Properties.VariableNames{ind2}, '_200-250um'])(i-2) = mean(data.(data.Properties.VariableNames{ind2})(find(data.R >200 & data.R < 250))) ;
        end

        end

     cd(strcat(path{ii}, '/',folders(i).name, '/Matlab_data')) ;
    writetable(edge_score, 'edge_score.csv') ;
    end

    cd(path{ii}) ; 
    writetable(edge_score_all, 'edge_score_all.csv') ;

end





