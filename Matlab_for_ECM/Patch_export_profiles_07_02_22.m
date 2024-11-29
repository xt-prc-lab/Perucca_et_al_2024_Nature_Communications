%% For ALL subfolders : plot intensity profiles for the different ECM components


% here path is a cell array containing all the "main folders" you wanna use
% (1 "main folder" contains several "subfolders)
path = { '/home/yourusername/Documents/YourCodeFolder/Codes'} %,


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
        cd(strcat(path{ii}, '/',folders(i).name, '/Matlab_data')) ;load('PROFILE.mat') ;  ; 
        % We go back to the subfolder and extract the name of the images.
        % The name of the images are used to name (and call) the columns of
        % PROFILE.mat and MERGE
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
               MERGE.(imgs(j).name(1:end-4)){k-1} = cat(2, MERGE.(imgs(j).name(1:end-4)){k-1}, transpose(PROFILE.([imgs(j).name(1:end-4),'_MEAN'])(find(PROFILE.R < binx(k) & PROFILE.R > binx(k-1)))  )) ;;
                MERGE.([imgs(j).name(1:end-4), '_', num2str(i-2)]){k-1} = mean(PROFILE.([imgs(j).name(1:end-4),'_MEAN'])(find(PROFILE.R < binx(k) & PROFILE.R > binx(k-1))))  ;
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
writetable(MERGE, 'PROFILES_MERGED.csv') ;


