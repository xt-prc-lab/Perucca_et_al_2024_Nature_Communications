function [] = angle_profile(path, Window_size)

% for MATLAB versions < R2017, the vecnorm function do not exist, we define
% it here
function p = vecnorm(p,n,d)
    p = sqrt(sum(p.^2,d));
end

% path = '/home/challopeau/Downloads/CAF200_COL_1_1/CAF200_COL_1/CAFF200_COL_EXP015_Pos0' ;

cd([path,filesep,'ctFIRE',filesep,'ctFIREout']) ;
matfiles = dir('*.mat') ;


cd([path, filesep, 'Matlab_data']) ; % move into the matlab_data folder
Radius = load('Radius.mat') ;   Radius = Radius.Radius ;%Radius
Vecfield = load('Vecfield.mat').Vecfield ;Vecfieldx = Vecfield(:,:,1) ; Vecfieldy = Vecfield(:,:,2) ; 
% load('Vecfield.mat') ;Vecfieldx = Vecfield(:,:,1) ; Vecfieldy = Vecfield(:,:,2) ; 

cd([path, filesep,'ctFIRE',filesep]) ; imgs = dir('*.tif') ;
a = size(imread([path, filesep,'CAF.tif']),1) ;  % a = size of the image used to draw the edge and compute the distance matrix (the CAF image)


% For each pixel of a fiber, the segment formed by
% pixel 'i' and 'i+step' is used to compute the angle at a given location of the fiber, 
% with the vector normal to the edge associated to this location.  
step = 4 ;

R_start = min(min(Radius))  ;                           %distance of the pixel nearer of the centroid (µm)
R_end = max(max(Radius)) ;                            %distance of the pixel the more far away from the centroid (µm)
R_vec = round(R_start):1:round(R_end) ;           % vector of the distance to sphero center, starting at the nearest pixel to the most far away

PROFILE = table() ;
PROFILE.('R') = transpose(R_vec) ;      % Distance to edge in µm


    for m = 1:length(imgs)
        b = size(imread([path,filesep,'ctFIRE',filesep,imgs(m).name]),1) ;  % size of the image currently analyzed

        fact = b/a ;   % "fact" is a resize factor, because the images to run ctFIRE (contained in the "ctFIRE"folder) may have a different size (e.g. 512*512 px) than the ones in the subfolder

        data = load([path,filesep,'ctFIRE',filesep,'ctFIREout',filesep,matfiles(m).name]) ;  %data = fibers from ctFIRE
        data = data.data ;;
    
        angles  = [] ; 
         for i=1:length(data.Fa)
            fiber = data.Xa(data.Fa(i).v,:);  
            fiber = fiber/fact ;   %rescale the fiber (usually from 512*512px to 1024 or 2048 = size of 'Radius' matrix)
            fiber = fiber( [fiber(:,1) > 0 & fiber(:,1) < size(Radius,2)] & [fiber(:,2) > 0 & fiber(:,2) < size(Radius,1)], :) ; %get rid of fibers out of the image/matrix
                
            if size(fiber,1) > step+1
        
                if ~isempty(fiber)
                    fiber = [fiber, Radius(sub2ind(size(Radius), round(fiber(:,2)), round(fiber(:,1))))] ; ;% append the distance to the edge in µm
        
                    vedge = [Vecfieldx(sub2ind(size(Vecfieldx), fiber(:,2), fiber(:,1))) , Vecfieldy(sub2ind(size(Vecfieldy), fiber(:,2), fiber(:,1)))] ; 
                    vedge = vedge(1:end-step,:) ; % exclude the last ones
        
                    vsegm = [fiber(1:end-step,1) - fiber(1+step,1), fiber(1:end-step,2) - fiber(1+step,2) ] ; 
        
                    ang = rad2deg(acos(dot(vedge, vsegm,2) ./ (vecnorm(vedge,2,2).*vecnorm(vsegm,2,2)) )) ; % angle for all the fiber, with the nearest normal vectors
                    %ang(isnan(ang)) = [] ; ; %remove nan 
                    %ang = median(ang) ; ; % take the median
                    ang(ang>90) = 180 - ang(ang>90) ; ; % convert to [0,90] if necessary
                    ang = abs(90 - ang)  ; %take the tangent angles
      
                    % ang = [ang, (fiber(step+1:end,end) + fiber(1:end-step,end))/2] ; % append the mean distance to the edge in µm
                    ang = [ang, fiber(1:end-step,end)] ;
                    
                    % remove nan
                    ang = ang(~isnan(ang(:,1)),:) ; 

                    angles = [angles; ang] ;

                end
            end
         end


         % measure profile of angles
%         R_start = min(angles(:,2))  ;                           %distance of the pixel nearer of the centroid (µm)
%         R_end = max(angles(:,2));                            %distance of the pixel the more far away from the centroid (µm)
%         R_vec = round(R_start):1:round(R_end) ;           % vector of the distance to sphero center, starting at the nearest pixel to the most far away
        num_R=length(R_vec);                                    %number of pixels between the nearest and farer pixels
        data_vec= {} ;%zeros(num_R,1);                           %empty cell with as many cells as pixels between R_start and R_end
        data_mean = zeros(num_R,1);                         % vector as long as the number of pixel, will contain the mean intensity around each pixel
        data_std=zeros(num_R,1);
        data_sem=zeros(num_R,1);

        for i=1:num_R    % for each pixel
            % we find all the pixels around the pixel of interest
            % the "size" of "around" depends on the value of Window_size
            
            % a pixel "i" is at R_vec(i) from the spheroid edge. We wanna find
            % all pixels whose distance is comprised in the range
            % [R_vec(i) - windows_size/2 , R_vec(i) + window_size/2] -> it is like
            % finding all the pixels on a radial slice, centered on R_vec(i), whose
            % width is "windows_size"
            data_index = find(angles(:,2) >=R_vec(i) - Window_size/2 & angles(:,2) <= R_vec(i) + Window_size/2);
            data_vec{i} =  angles(data_index,1) ;%we extract the value of intensity for each pixel of the slice
        end

        % and compute the mean 
        for i = 1:num_R
            data_mean(i) = mean(data_vec{i}, 'omitnan') ;
            data_std(i) = std(data_vec{i}, 'omitnan') ;
            data_sem(i) = std(data_vec{i}, 'omitnan')/sqrt(length(data_vec{i})) ; 
        end


%         figure ; fill([R_vec, fliplr(R_vec)], [data_mean'+ data_std', fliplr(data_mean' - data_std')], [.5,.5,.5], 'linestyle', 'none', 'FaceAlpha', .5) ; hold on ; plot(R_vec,data_mean) ; xlabel('Distance to edge (µm)') ;
% ylabel('Mean fiber angle toward the edge (°)') ; hold on ; plot([R_start, R_end], [0,0], 'linestyle' , '--', 'color', 'k') ; hold on ; 
% plot([0, 0], [min(data_mean-data_std),max(data_mean+data_std)], 'linestyle' , '--', 'color', 'k')


% Gather the normalized radius, the intensity profile along each circles,
% the mean, std and sem of the intensity profile along each circles
PROFILE.(strcat(imgs(m).name(1:end-4),'_MEAN')) = data_mean ;%extract mean
PROFILE.(strcat(imgs(m).name(1:end-4),'_STD')) = data_std ; % standard deviation
PROFILE.(strcat(imgs(m).name(1:end-4),'_SEM')) = data_sem ; %standard error of the mean


    end
    
writetable(PROFILE, strcat(path,filesep,'Matlab_data', filesep, 'Angle_profile.csv')) ;

end

