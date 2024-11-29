function [ ] = analyze_intensity(path, codepath,Window_size)
%%%% AIM %%%%
% Compute the intensity profile of Fibroblast and ECM image, as a function
% of the distance to the Fibroblast edge

%%%% INPUTS %%%%
% The distance matrix 'Radius', each image to analyze

%%%% OUTPUTS %%%%
% 'Radius' : matrix where the value of each cell is the distance to the
% Fibroblast edge ( >0 outside spheroid, <0 inside spheroid, =0 at the
% Fibroblast edge)
% Clément Hallopeau 05/2021


% change directory toward the subfolder
cd(path) ; 

% struct containing the name of all the images within the subfolder
imgs = dir('*.tif') ; 

Radius = load([path, filesep, 'Matlab_data', filesep, 'Radius']) ; Radius = Radius.Radius ; 

PROFILE = table() ; 

%% Intensity profile for each image
for m = 1:length(imgs)   %for each image

%load tiff image
t = Tiff(strcat(path,filesep, imgs(m).name),'r');
imageData = read(t);

%measure mean intensity for all pixels within a certain radius range
data=single(imageData);                                 %convert image data to single


% resize the distance matrix if the size of the image to analyze is
% different than the one used to compute the distance matrix
Radius = imresize(Radius, size(data)) ; 

R_start = min(min((Radius)));                           %distance of the pixel nearer of the centroid (µm)
R_end = max(max((Radius)));                            %distance of the pixel the more far away from the centroid (µm)
R_vec = round(R_start):round(R_end) ;           % vector of the distance to sphero center, starting at the nearest pixel to the most far away
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
    data_index = find(Radius >=R_vec(i) - Window_size/2 & Radius <= R_vec(i) + Window_size/2);
    data_vec{i} =  data(data_index) ;%we extract the value of intensity for each pixel of the slice
end


%compute mean intensity at each position of num_R
for i = 1:num_R
    data_mean(i) = mean(data_vec{i}) ;
end

%normalization of all the measured intensities to have max(mean) =1
for k  = 1 : length(data_vec)  
    data_vec{k} = data_vec{k}/max(data_mean) ; 
end 

% and compute the mean again (this time the mean will be = 1)
for i = 1:num_R
    data_mean(i) = mean(data_vec{i}) ;
    data_std(i) = std(data_vec{i}) ;
    data_sem(i) = std(data_vec{i})/sqrt(length(data_vec{i})) ; 
end

% Gather the normalized radius, the intensity profile along each circles,
% the mean, std and sem of the intensity profile along each circles
PROFILE.('R') = transpose(R_vec) ;      % Distance to edge in µm
PROFILE.(imgs(m).name(1:end-4)) = data_vec'  ; %extract the raw data
PROFILE.(strcat(imgs(m).name(1:end-4),'_MEAN')) = data_mean ;%extract mean normalized to 1
PROFILE.(strcat(imgs(m).name(1:end-4),'_STD')) = data_std ; % standard deviation
PROFILE.(strcat(imgs(m).name(1:end-4),'_SEM')) = data_sem ; %standard error of the mean


end


%SAVE  the profile
cd(strcat(path,filesep,'Matlab_data')) ;
save('PROFILE', 'PROFILE') ;

disp(strcat(path,'        Intensity Profile DONE'))
clearvars PROFILE data_mean data_vec data_std data_sem data_index R_vec R_end  R_start Radius data num_R x_c y_c x y imageData t CD filename xmin xmax img center
end





% m = 0 ;%find max value
% for k = 1: length(data_vec)  ;%find the max value of intensity
%     if max(data_vec{k})> m 
%         m = max(data_vec{k}) ;
%     end
% end


% %find min and max gray values
% [xmax, xmin] = deal(0,0) ;
% for i = 1:num_R
%     if max(data_vec{i}) > xmax ;
%         xmax = max(data_vec{i}) ; 
%     end
%     if min(data_vec{i}) < xmin ;
%         xmin = min(data_vec{i}) ;
%     end
% end
