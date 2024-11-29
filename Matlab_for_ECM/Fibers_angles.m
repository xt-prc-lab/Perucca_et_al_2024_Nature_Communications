function [] = measure_ECM(path, codepath, angmap )
 %%%% AIM %%%%
% Compute the angles of the fibers found by ctFIRE, as respect to the
% Fibroblast Edge drawn previously

%%%% INPUTS %%%%
% - 'Radius.mat' : output from 'DefineEdge.mat'
% - ctFIRE outputs : named as 'ctFIREout_nameoftheimage.mat'

%%%% PROCESSING %%%%
% From the edge, two matrices are derived, where each cell contain the X or Y coordinates 
% of the vector normal to the closest piece of edge. These matrices are smoothed (takes long)
% 
% For each pixel of each fiber, the distance to the edge is derived from
% 'Radius', to further assign a region (IN,EDGE,OUT) for each one of them
% 
% For a given fiber, the local angles that it describes toward the edge are
% computed. 
% For a given region, we extract all the local angles of all the pieces of fibers
% spanning in this region

%%%% OUTPUTS %%%%
% - For each image, a PNG figure with the image, the ctFIRE fibers colorcoded by distance or angle to the edge
% (Overlay-Angle2Edge_nameoftheimage,  Overlay-Distance2Edge_nameoftheimage)
% - 'MEASURE.mat' 

% Clément Hallopeau 05/2021


% for MATLAB versions < R2017, the vecnorm function do not exist, we define
% it here
function p = vecnorm(p,n,d)
    p = sqrt(sum(p.^2,d));
end

cd([path,filesep,'ctFIRE',filesep,'ctFIREout']) ;
matfiles = dir('*.mat') ;


cd([path, filesep, 'Matlab_data']) ; % move into the matlab_data folder
Radius = load('Radius.mat') ;   Radius = Radius.Radius ;%Radius

Mask = zeros(size(Radius)) ; 
Mask(find(Radius <0)) = 1 ; 


%%  For each pixel of the field of view, compute the vector normal to the closest edge

if angmap
    cd(codepath) ;
    [s,s2,Krv,Coord,Normal] = curv(Mask); %compute normal direction
    %here i prepare the variables to make the field of normal
    %direction vectors
    Vecfieldx=zeros(size(Mask));Front=Vecfieldx;Vecfieldy=Vecfieldx;
    Tempx=Vecfieldx;Vecfieldx(:)=NaN;
    Tempy=Vecfieldy;Vecfieldy(:)=NaN;
    %we introduce the values calculated 4 lines above in the matrices
    Tempy(sub2ind(size(Vecfieldx),Coord(2,:),Coord(1,:)))=Normal(:,2);  % Matrices whose size is the size of the image
    % The cells corresponding to the pixels of the edge contains the x or y
    % coordinates of the  normal vectors
    Tempx(sub2ind(size(Vecfieldx),Coord(2,:),Coord(1,:)))=Normal(:,1);
    Vecfieldx=Tempx;Vecfieldy=Tempy;
    
    % we assign to each point the closest vector of the leading edge
    Front(sub2ind(size(Mask),Coord(2,:),Coord(1,:)))=1;
    [DF,LF]=bwdist(Front);  %DF = distance matrix, LF = index of the closest pixel
    [RF,CF]=ind2sub(size(Front),LF);  % RF/CF = Y/X coordinates of the closest edge pixel
    DF=round(DF);   %distance matrix, positive (even inside the Mask)
    for ni=1:max(max(DF));
        Tempx(find(DF==ni))=Tempx(LF(find(DF==ni)));   % Matrices whose size is the size of the image, 
        % each cell of the matrices contains the x or y coordinates of the
        % closest normal vector
        Tempy(find(DF==ni))=Tempy(LF(find(DF==ni)));
    end
   Vecfieldy = Tempy ; Vecfieldx = Tempx ; 
    
    %  mean smoothing of the angle maps
%     wb = waitbar(0) ; 
    for n=1:round(max(max(DF)))  % DF = matrix where each cell contain the distance to the closest edge of the corresponding pixel
        waitbar(n/round(max(max(DF))), wb, 'Smoothing the angle map') ; 
        [im,jm]=ind2sub(size(DF),find(DF==n));  % im, jm = x,y coordinates of the pixels located at a distance "n" from the edge
        for m=1:size(im,1) % for each of these pixel
            in=im(m);
            jn=jm(m);
            temp = Tempx(max(in-DF(in,jn),1):min(in+DF(in,jn),size(DF,1)),...
                max(jn-DF(in,jn),1):min(jn+DF(in,jn),size(DF,2)));
            Vecfieldx(in,jn)=nanmean(nanmean(temp));
            temp = Tempy(max(in-DF(in,jn),1):min(in+DF(in,jn),size(DF,1)),...
                max(jn-DF(in,jn),1):min(jn+DF(in,jn),size(DF,2)));
            Vecfieldy(in,jn)=nanmean(nanmean(temp));
        end
    end
%     close(wb) ;
    
    Vecfield = Vecfieldx ; Vecfield(:,:,2) = Vecfieldy ; 
    cd([path, filesep, 'Matlab_data']) ;
    save('Vecfield.mat', 'Vecfield') ;

end

%%
cd([path, filesep, 'Matlab_data']) ;
load('Vecfield.mat') ;Vecfieldx = Vecfield(:,:,1) ; Vecfieldy = Vecfield(:,:,2) ; 

cd([path, filesep,'ctFIRE',filesep]) ; imgs = dir('*.tif') ;
a = size(imread([path, filesep,'CAF.tif']),1) ;  % a = size of the image used to draw the edge and compute the distance matrix (the CAF image)

%%%%%%%%%%%%%% DEFINE THE REGIONS (in, edge, out) HERE %%%%%%%%%%%%%%%%%
% % Rbins contains the "limits" of each region : 
Rbins = [min(min(Radius)), 0 ;               %IN
    0, 120 ;                                                     %EDGE
    240 , 500] ;                                              %OUT      


MEASURE = table();
for m = 1:length(imgs)
    

b = size(imread([path,filesep,'ctFIRE',filesep,imgs(m).name]),1) ;  % size of the image currently analyzed

fact = b/a ;   % "fact" is a resize factor, because the images to run ctFIRE (contained in the "ctFIRE"folder) may have a different size (e.g. 512*512 px) than the ones in the subfolder
% (e.g. 1048*1048 or 2048*2048) which are use to draw the Fibroblast edge
% and derive the distance matrix 'Radius'
% e.g., fact = 512/1048 = 0.5, or 512/2048=0.25
% If we draw the edge on a 1048 image, the ctFIRE image is
% twice smaller (e.g; 512*512), we need to resize the fibers to size of the CAF image. 
% I.E. divide the fibers by fact

data = load([path,filesep,'ctFIRE',filesep,'ctFIREout',filesep,matfiles(m).name]) ;  %data = fibers from ctFIRE
data = data.data ;;
MEASURE.(['BINS_',imgs(m).name(1:end-4)]) =  num2cell(Rbins,2) ;
MEASURE.(['THETA_',imgs(m).name(1:end-4)]) = cell(size(Rbins,1),1) ;


%% Plot fiber - Distance to edge colorcoded

fig = figure('Position', [200,200,800,800]) ; title('Fibers Distance to Spheroid Edge') ; 

img = imread([path, filesep, imgs(m).name]) ;;  img = imadjust(img) ;;
img = imresize(img, size(imread([path, filesep,'CAF.tif']))) ; 
cmap = img ; cmap(:,:,2) = img ; cmap(:,:,3) = img ; 
imagesc(img, 'CData',cmap) ; hold on ;

cmap = colormap(parula(max(abs(min(min(Radius))), max(max(Radius)) ) )) ; 

 for i=1:length(data.Fa)
    
    fiber = data.Xa(data.Fa(i).v,:);  
    fiber = fiber/fact ;   %rescale the fiber (usually from 512*512px to 1024 or 2048 = size of 'Radius' matrix)
    fiber = fiber( [fiber(:,1) > 0 & fiber(:,1) < size(Radius,2)] & [fiber(:,2) > 0 & fiber(:,2) < size(Radius,1)], :) ; %get rid of fibers out of the image/matrix

    if ~isempty(fiber)
    fiber = [fiber, abs(Radius(sub2ind(size(Radius), round(fiber(:,2)), round(fiber(:,1)))))  ] ; % append the distance to the edge in µm
    
    plot(fiber(:,1), fiber(:,2), 'Color', cmap(1+abs(round(median(fiber(:,end) ))),:), 'LineWidth', 1, 'HandleVisibility', 'off') ; hold on ;
    end
 end
 
 contour(Mask,1,'Color', 'r', 'LineWidth', 1, 'DisplayName', 'CAF Edge')    ; hold on ;%plot the CAF edge
cb = colorbar ;;
caxis([0,max(max(Radius))]) ;  cb.Title.String = sprintf('Distance to \n Edge [µm]') ;

cd([path, filesep, 'Matlab_data']) ; saveas(fig,['Overlay-Distance2Edge-',imgs(m).name(1:end-4),'.png']) ; 
 close(fig) ;



%% Plot fibers - angles colorcoded



% For each pixel of a fiber, the segment formed by
% pixel 'i' and 'i+step' is used to compute the angle at a given location of the fiber, 
% with the vector normal to the edge associated to this location.  
step = 4 ;


fig = figure('Position', [200,200,800,800]) ; title('Fibers Angle to Spheroid Edge') ;  
img = imread([path, filesep, imgs(m).name]) ;;  img = imadjust(img) ;;
img = imresize(img, size(imread([path, filesep,'CAF.tif']))) ; 
cmap = img ; cmap(:,:,2) = img ; cmap(:,:,3) = img ; 
imagesc(img, 'CData',cmap) ; hold on ;

cmap = colormap(parula(91)) ;
thetamax =  0 ;
for i=1:length(data.Fa)
    fiber = data.Xa(data.Fa(i).v,:);  
    fiber = fiber/fact ;
    fiber = fiber( [fiber(:,1) > 0 & fiber(:,1) < size(Radius,1)] & [fiber(:,2) > 0 & fiber(:,2) < size(Radius,1)], :) ; %get rid of fibers out of the image/matrix

    if size(fiber,1) > step+1

        if ~isempty(fiber)                    
            fiber = [fiber, Radius(sub2ind(size(Radius), round(fiber(:,1)), round(fiber(:,2))))] ; ;% append the distance to the edge in µm

            vedge = [Vecfieldx(sub2ind(size(Vecfieldx), fiber(:,2), fiber(:,1))) , Vecfieldy(sub2ind(size(Vecfieldy), fiber(:,2), fiber(:,1)))] ; 
            vedge = vedge(1:end-step,:) ; % exclude the last ones

            vsegm = [fiber(1:end-step,1) - fiber(1+step,1), fiber(1:end-step,2) - fiber(1+step,2) ] ; 

            ang = rad2deg(acos(dot(vedge, vsegm,2) ./ (vecnorm(vedge,2,2).*vecnorm(vsegm,2,2)) )) ; % angle for all the fiber, with the nearest normal vectors
            %ang(isnan(ang)) = [] ; ; %remove nan 
            %ang = median(ang) ; ; % take the median
            ang(ang>90) = 180 - ang(ang>90) ; ; % convert to [0,90] if necessary
            ang = abs(90 - ang)  ; %take the tangent angles


            for j  = 1:size(Rbins,1)  %we take each region 1 by 1
                pos = find(fiber(1:end-step,end) > Rbins(j,1) & fiber(1:end-step,end) < Rbins(j,2)) ; 
                if ~isempty(pos)  % if a piece of the fiber belongs to the given region, append the angles of this piece of fiber
                    MEASURE.(strcat('THETA_',imgs(m).name(1:end-4))){j} = cat(1, MEASURE.(strcat('THETA_',imgs(m).name(1:end-4))){j}, ang(pos)) ;
                end
            end
            ang(isnan(ang)) = [] ; ; %remove nan
             plot(fiber(:,1), fiber(:,2), 'Color', cmap(round(median(ang))+1,:), 'LineWidth', 2, 'HandleVisibility', 'off') ; hold on ;
            thetamax = max(thetamax, round(max(ang))) ; 

        end
    end
end
cb = colorbar ; caxis([0,90]) ; cb.Title.String = sprintf('Angle Toward the Edge [um]') ;

 %contour(Mask,1,'Color', 'w', 'LineWidth', 1, 'DisplayName', 'CAF Edge')
 
   %plot the bins
%   Radius = transpose(Radius) ; 
%  for k = 1:size(Rbins,1)
%      start = zeros(size(Radius))  ;  endd = start ; start(find(Radius < Rbins(k,1))) = 1 ; endd(find(Radius < Rbins(k,2))) = 1 ; 
%      contour(start,1, 'Color', 'g', 'LineStyle', '--', 'LineWidth',2,'DisplayName', 'Region Starts') ; hold on ;contour(endd,1, 'Color', 'r', 'LineStyle', '-',  'LineWidth',2,'DisplayName', 'Region Ends') ;  ; hold on ;
%  end ; legend
 

cd([path, filesep, 'Matlab_data']) ; 
saveas(fig,['Overlay-Angle2Edge-',imgs(m).name(1:end-4),'.png']) ; 
close(fig) ; 
end

save('MEASURE.mat', 'MEASURE') ;


end
