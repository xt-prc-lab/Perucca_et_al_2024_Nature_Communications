function [] = morpho_params(path, file)
%%%% AIM %%%%%
% To add morphological parameters to each point of each Trackmate track

%%%% INPUTS %%%%
 %      - Frames for each timepoint of the current file to analyze
 %      - Tracks associated to the analyzed file

 %%%% PROCESSING %%%%
%       - The tracks are formatted as a table 'tracks'. Lines correspond to
%       a given point in space and time. Columns correspond to the timestep
%       between a current point (x,y,t) and the next one (x,y,t+timestep),
%       the X coordinates, the Y coordinates, the ID of the track
%       - For each frame of the analyzed file, the centroid and
%       morphological parameters of each immune cell are extracted in a
%       table 'morpho'
% Then,  each line of 'tracks' is associated to a line of 'morpho' by
% finding the line in 'morpho' having the closest (x,y) coordinates than
% the given line of 'tracks'
% Clément Hallopeau 03/2020

disp(file) 

% import the tracks
cd(path)  ;
tracks = importTrackMateTracks(strcat(file(1:length(file)-4),'_Tracks.xml')) ;   % May not work if the FIJI directory containing diverse scripts is not added to MATLAB path (see MAIN.mat line 5)

% import the frames
cd(strcat(path, '/Frames')) ;
frames = dir('*.tif') ;  % get the name of all the frames inside the folder
frames = {frames(~cellfun('isempty',strfind(cellstr(char(frames.name)),file(1:end-4)))).name} ;  % keep the name of the frames for the current position


%gather the cells coordinates, and their position in the trackmate file (=k)
t = [] ;  %timestep of each timepoint, x, y, track ID
for k = 1:length(tracks)
    c = tracks{k}  ;
    t = cat(1, t, [ [diff(c(:,1))', NaN]', c(:,2), c(:,3), ones(size(c,1),1)*k]);  % timestep between a given timepoint and the next one, position in x, position in y, track ID (= position in the trackmate file)
end


%gather the centroid position, cells morphological parameters, ID of the
%frame, for all the cell of all the frames
morpho = [] ;  
wb = waitbar(0) ;
 for i = 1:length(frames)  % for each frame associated to the analyzed file
    waitbar(i/size(frames,2), wb, 'Frames analysis') ;;
    img = imread(frames{i}) ;  % load the frame whose name is given by frames{i}
    CC = bwconncomp(img) ; % find white components = white connected pixels
    L = labelmatrix(CC) ; %  label the connected components

    stats = regionprops(L,'all', 'table') ; % and analyze their properties

    for k = 1 : size(stats, 1) % we add the parameters extracted from frames{i} to the big dataframe
        % we only extract the following parameters : x, y, timepoint, area,
        % major & minor axis lenth, circularity, orientation
        morpho = cat(1, morpho, [stats(k).Centroid(1),  stats(k).Centroid(2),  i, stats(k).Area ,  stats(k).MajorAxisLength,   stats(k).MinorAxisLength,   stats(k).Circularity,   stats(k).Orientation]) ;
    end
 end
 close(wb) ;

% conversion from array to table
morpho = array2table(morpho, 'VariableNames', {'X', 'Y', 'frame', 'Area', 'Majax', 'Minax', 'Circularity', 'Orientation'}) ;
%morpho = X, Y, frame number, parameters ...

tracks = cat(2, t, zeros(size(t,1),size(morpho,2)-1)) ;  %add empty columns for the extracted parameters
tracks = array2table(tracks, 'VariableNames', {'dt', 'X', 'Y', 'ID','Xm', 'Ym',morpho.Properties.VariableNames{4:end}}) ; 
% tracks = X, Y, cell ID, ...parameters from morpho...


error = [] ;
wb = waitbar(0) ;;
for i =1:size(tracks,1)  %for each position of each cell sent by trackmate
    waitbar(i/size(tracks,1), wb, 'Matching Trackmate & Matlab') ;
    tol = 1 ;%we find the corresponding cell in matlab's output, i.e. the cell having the same/similar position
    % tol  = tolerance -> high tol (=0.1) almost all positions matlab are
    % matched with the ones of trackmate, but high errors. Low tol
    % (=0.0001) low error but positions never match. tol=0.01 looks good
    
    % we first select the positions the most similar to (track.X(i), tracks.Y(i))
    % we decrease the tolerance (tol) while we have more than 10% of the
    % positions selected
    while numel(intersect(find(morpho.X > tracks.X(i)*(1-tol) & morpho.X < tracks.X(i)*(1+tol)), find(morpho.Y > tracks.Y(i)*(1-tol) & morpho.Y < tracks.Y(i)*(1+tol)))) > size(tracks,1)/10
       tol=tol/2 ; 
    end
   tol=2*tol ; 
  
   % if we have found a group of cells sharing similar positions to (track.X(i), tracks.Y(i))
    if numel(intersect(find(morpho.X > tracks.X(i)*(1-tol) & morpho.X < tracks.X(i)*(1+tol)), find(morpho.Y > tracks.Y(i)*(1-tol) & morpho.Y < tracks.Y(i)*(1+tol)))) > 0
        % we extract their index in the matlab output
        pos = intersect(find(morpho.X > tracks.X(i)*(1-tol) & morpho.X < tracks.X(i)*(1+tol)), find(morpho.Y > tracks.Y(i)*(1-tol) & morpho.Y < tracks.Y(i)*(1+tol))) ;
        % allowing to compute the errors in X and Y between the selected
        % positions  and (track.X(i), tracks.Y(i))
        er = abs(tracks.X(i) - morpho.X(pos))/tracks.X(i) + abs(tracks.Y(i) - morpho.Y(pos))/tracks.Y(i) ;
        
        % Among the selected positions, we select the one where the error
        % with (track.X(i), tracks.Y(i)) is the lowest 
        pos = pos(find(er == min(er))) ; pos = pos(1) ;

        % And keep this value of error
        error = cat(1,error, [abs(tracks.X(i) - morpho.X(pos))/tracks.X(i)*100, abs(tracks.Y(i) - morpho.Y(pos))/tracks.Y(i)*100, tol*1000]) ;
        

        tracks(i,7:end) = morpho(pos, 4:end) ;% assign parameters from matlab to the cell "i"
        tracks(i, 5:6) = table(morpho.X(pos), morpho.Y(pos)); % and also export the (x,y) position kept from matlab (to compute the error afterward, to overlay each
        % trackmate track with the reconstitution using matlab's positions
    end
end
close(wb) ;


tracks = standardizeMissing(tracks, 0)  ; %remove zeros due to unassigned positions

% The we plot each trajectory from trackmate, overlaid with the rebuilt trajectory from
%matlab positions'
% And we ask the user if he/she want to keep the track by judging how the
% trackmate (black) and matlab (red) tracks are overlaid, and whether they
% have the same shape or not
disp(sprintf('\n **Track Selection**'))

keep = [] ; % 'keep' contains the ID of the tracks to keep
fig = figure('Position', [2500,500,600,600]) ;  % position to plot on a second screen, if the figure doesn't pop up, replace 2500 by 500
for i  = 1 : length(unique(tracks.ID))  % for each cell
    next = 1 ;
    while next  %plot
    plot(tracks.X(find(tracks.ID == i)), tracks.Y(find(tracks.ID == i)), 'Color', [0 0 0]) ; hold on ;
    plot(tracks.Xm(find(tracks.ID == i)), tracks.Ym(find(tracks.ID == i)), 'Color', 'r') ; hold on ;
    title(['Progress : ',num2str(100*i/length(unique(tracks.ID))), '%']) ; ; 
    ask = input('Keep this track ? [y/n]', 's') ;
    if  isempty(ask) | ask == 'y'  %if 'yes'
        keep = [keep, i] ; % we extract the id of the cell
        clf 
        next = 0 ;
    elseif ask == 'n'
        clf
        next = 0 ;
    end
    if ask ~= 'y' & ask ~= 'n' & ~isempty(ask)
          ask = input('Keep this track ? [y/n]', 's') ;
    end
    end
end
close(fig)

disp(sprintf('\n **Track Selection : DONE**'))


% Here we erase the morphological data of the cells we just discarded
trackscleaned = table2array(tracks) ;% conversion to matrix
exclude = setdiff(trackscleaned(:,4), keep) ; % id of the cells discarded
for k = 1 : length(exclude)
    pos = find(trackscleaned(:,4) == exclude(k)) ; 
    trackscleaned(pos,7:end) = NaN ; 
end
trackscleaned = array2table(trackscleaned, 'VariableNames', tracks.Properties.VariableNames) ;% conversion to table






if ~exist(strcat(path, '/Morpho_params/'), 'dir') 
    mkdir(strcat(path, '/Morpho_params/'))
end

cd([path, '/Morpho_params/'])
save([file(1:end-4),'_Morphoparams.mat'], 'trackscleaned')
    
    clearvars trackscleaned tracks tracksmat keep ask next error pos er morpho
    
end
