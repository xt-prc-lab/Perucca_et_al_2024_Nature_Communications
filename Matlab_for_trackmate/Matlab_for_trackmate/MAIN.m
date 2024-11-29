%% main code to run
% Clément Hallopeau 03/2020

codepath = 'insert.code.path.here' ;                                        %Codes directory
cd(codepath)  ;
 
addpath('/home/Fiji.app/scripts/') ;                                        % add to MATLAB path the FIJI folders containing diverse scripts


% In the same folder, put the first image of each stack (CAF channel)
% Put the corresponding trackmate outputs named as "  nameofyourimage_Tracks.xml  "
path  = uigetdir('/insert.path.to.the.data.to.analyze.here');
tiffiles = struct2table(dir(fullfile(path,'*tif*'))).name  ;                %list of current images in the folder
trackfiles = struct2table(dir(fullfile(path,'*xml*'))).name  ;              %list of the xml/tracks in the folder


%% Draw the edge, measure  the intensity profiles of the given images as respect to the drawn edge
disp(sprintf('\n *****Radius*****'))
%radius in µm


for k =  1:length(tiffiles)                                                 % for image (=position) in the folder
    file = char(tiffiles(k)) ;                                              % get filename
    disp(file)

    sigma=30;                                                               % parameter of the Gaussian filter used to smooth the edge, and avoid sharp angles
    
    % Definition of px_to_micro depending on the image size
    px_size = 0.3225 ;                                                      % pixel size, in microns per pixel, for 2048*2048 field of views
    px_to_micro = px_size*2048/size(imread([path,filesep,file]),1) ;        %convert pixel size considering the image binning

    cd(codepath)
    DefineEdge(path, file,sigma, px_to_micro);  
end

% output in micrometers
disp('*****Radius DONE*****')


%% Measure and save morphological parameters

disp(sprintf('\n *****Morphological Parameters*****'))


for k = 1 : length(tiffiles)
     file = char(tiffiles(k)) ;
    cd(codepath)
    morpho_params(path, file) ;
 end

disp(sprintf('\n *****Morphological Parameters     DONE__',char(datetime('now','Format','HH:mm:ss')),'*****'))

%%  create dataset

% create a new dataframe
% File,dt,X,Y,ID,Area,Majax,Minax,Circularity,Orientation,dist,Max_Area,Bin,Theta,Di,Dm,Vi,Vm,Ri,Rm,AFi,AFm

DATA = table;

leg = {'IN', 'EDGE','OUT'} ; 



for  i=1:length(tiffiles)

    file = char(tiffiles(i)) ;
    disp(file)
    picture=imread(fullfile(path,file));
    img = imadjust(rgb2gray(picture)) ;

    % Definition of px_to_micro depending on the image size
    px_size = 0.3225 ;                                                      % pixel size, in microns per pixel, for 2048*2048 field of views
    px_to_micro = px_size*2048/size(imread([path,filesep,file]),1) ;        %convert pixel size considering the image binning

    % get framerate
    disp([file, sprintf('\n'), 'Input the time between 2 frames (seconds)', sprintf('\n')])
    dt = input('') ;
    disp(sprintf('\n'))
    
    % tracks is the output from the morphological measures
    tracks = load([path,'/Morpho_params/',file(1:end-4),'_Morphoparams.mat']).trackscleaned ;
    
    load([path,'/Radius/',file(1:end-4),'_Radius.mat']) ;                   % load distance matrix
    Mask = zeros(size(Radius)) ;  Mask(find(Radius <0)) = 1 ;
    
    %%%%%%%%% Define the bins %%%%%%%%%%%
    rbins=[min(min(Radius)), 0, 120, 2*120];                                % values in microns
    %%%%%%%%% %%%%%%%%%%%%%%%%%%%


    tracks.Area = px_to_micro^2*tracks.Area ;
    % we create a new table called data, to pick only only the columns from
    % "tracks" that are interesting, and process them. For instance, tracks
    % contains the matlab positions' assigned to each trackmate position
    % (Xm, Ym), we don't use them. 
    data =  tracks(:,[1:4,7:end])  ;
    
    %Distance to the edge for each position
    data.('dist')=ones(size(data,1),1)*NaN;
    for j=1:size(data,1)
        x =round(data.X(j)); y =round(data.Y(j));
        if ~isnan(x)&~isnan(y)
            data.dist(j)=Radius(y,x);
        end
    end

    %Plot each position colorcoded depending on its distance to the edge 
    disp(' ** Plotting Positions distance-colorcoded') ;
    fig = figure('visible', 'off') ; cmap = img ; cmap(:,:,2) = img ; cmap(:,:,3) = img ; imagesc(img, 'CData',cmap); hold on ; cmap = colormap(parula(2+max(abs(data.dist)))) ; for j = 1:size(data,1) ; if ~isnan(data.X(j)) & ~isnan(data.Y(j)) ; scatter(data.X(j), data.Y(j), 4, cmap(round(abs(data.dist(j)))+1, :)) ; hold on ; end ; end ;   contour(Mask)

    if ~exist(strcat(fullfile(path,'Matlab_data')))
        mkdir(strcat(fullfile(path, 'Matlab_data'))) ;
    end
    
    cd([path, filesep, 'Matlab_data']) ; saveas(fig, [file,'_positions_distance2edge.jpg']) ;  close(fig) ; 

    %Compute max cell area 
    IDs = unique(tracks.ID) ; 
    data.('Max_Area') = zeros(size(data,1),1) ; 
    for ii = 1:length(IDs)
        pos = find(tracks.ID == IDs(ii)) ; 
        data.Max_Area(pos) = ones(length(data.Max_Area(pos)),1) * max(tracks.Area(pos)) ;
    end

    % Assign each position to a category
    data.('Bin') = repmat("NaN",size(data,1),1) ; 
    for j = 2:length(rbins)                                                 %take the limits on each bin
        % and find the positions where the cells are within the limits of
        % the considered bin (  [bin(j), bin(j-)]  )
        pos = find(data.dist < rbins(j) & data.dist > rbins(j-1)) ; 
        data.Bin(pos) = repmat(leg{j-1}, length(pos), 1) ; 
    end
   
    
    % Displacement angle toward the edge
    data.('Theta') = zeros(size(data,1),1)*NaN ; 
    IDs = unique(tracks.ID) ; 
    % For smoothing purpose, the toward the edge are computed using the
    % segment formed by the position at time 't' and 't+step'
    %%%%%%%%%% Set the step size to quantify angles %%%%%%%%%%
                                                        step = 4 ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd(codepath) ;
    [Vecfieldx, Vecfieldy] = angle_map(Mask) ;                              % create the maps normal angles to the edge
    
    for j = 1:length(IDs)
        pos = find(tracks.ID == IDs(j)) ; 
        y = round(data.Y(pos)) ; x = round(data.X(pos)) ;
        if length(pos)  > step &  sum(isnan(x)) == 0 &  sum(isnan(y)) == 0  % if the track is longer than "step", and if the positions do not contains NaN (sometimes Trackmate send NaN)
            vedge = [Vecfieldx(sub2ind(size(Vecfieldx), y , x )) , Vecfieldy(sub2ind(size(Vecfieldy), y, x ))] ;   % coordinate of vectors normal to the edge for the given (x,y) coordinates
            vedge = vedge(1:end-step,:) ;                                   % exclude the last ones
            vsegm = [x(1:end-step) - x(1+step), y(1:end-step) - y(1+step) ] ; % segments formed by the track
            ang = rad2deg(acos(dot(vedge, vsegm,2) ./ (vecnorm(vedge,2,2).*vecnorm(vsegm,2,2)) )) ; % angle between each segment and normal vector to the closest edge
            ang = abs(90-ang) ;                                             % convert to have the tangeant angle
            data.Theta(pos(1:end-step)) = ang ; 
        end
    end
    

% assign instant speed & displacements from 1 pos 2 another & mean speed to all pos
    data.('Di') = zeros(size(data,1),1)*NaN ;
    data.('Dm') = zeros(size(data,1),1)*NaN ;
    data.('Vi') = zeros(size(data,1),1)*NaN ;
    data.('Vm') = zeros(size(data,1),1)*NaN ;
    for j = 1:length(IDs)
        pos = find(tracks.ID == IDs(j)) ;
        vsegm=[tracks.X(pos(2:end)) - tracks.X(pos(1:end-1)),    tracks.Y(pos(2:end)) - tracks.Y(pos(1:end-1))];
        data.Di(pos) = px_to_micro .* [vecnorm(vsegm, 2,2)', NaN]' ;
        data.Dm(pos) = px_to_micro .* ones(size(data.Vm(pos),1),1) * mean(vecnorm(vsegm, 2,2)) ;
        data.Vi(pos) = px_to_micro./(dt/60.*tracks.dt(pos)) .* [vecnorm(vsegm, 2,2)', NaN]' ; 
        data.Vm(pos) = ones(size(data.Vm(pos),1),1) *mean(data.Vi(pos), 'omitnan') ; 
    end
    
    
    % Roundess and mean roundess (Xu et al)
    data.('Ri') = zeros(size(data,1),1)*NaN ;
    data.('Rm') = zeros(size(data,1),1)*NaN ;
    
   for j = 1:length(IDs)
        pos = find(tracks.ID == IDs(j)) ;
        circsc = tracks.Area(pos)./(pi*(tracks.Majax(pos)./2).^2) ; 
        data.Ri(pos) = circsc ; 
        data.Rm(pos) = ones(size(data.Rm(pos),1),1) * mean(circsc) ; 
   end
    
    % aspect factor and mean aspect factor (min axis / maj axis)
    data.('AFi') = zeros(size(data,1),1)*NaN ;
    data.('AFm') = zeros(size(data,1),1)*NaN ;
    for j = 1:length(IDs)
        pos = find(tracks.ID == IDs(j)) ;
        circsc = tracks.Minax(pos)./tracks.Majax(pos) ;  
        data.AFi(pos) = circsc ; 
        data.AFm(pos) = ones(size(data.AFm(pos),1),1) * mean(circsc) ; 
    end 
    
    data = [table(repmat(string(file(1:end-4)), size(data,1),1), 'VariableName',{'File'}), data] ; 
    DATA = [DATA; data] ; 
    
end

    % export in csv
cd(path)
name = pwd ;  pos = find(name == filesep ) ;  pos = pos(end) ;name = name(pos+1:end) ;% name of the csv = name of the folder
writetable(DATA, [name,'.csv']) ;

disp(sprintf('\n *****create dataset     DONE__',char(datetime('now','Format','HH:mm:ss')),'*****'))