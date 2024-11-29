function [] = DefineEdge(path,file, sigma, px_to_micro)
%%%% AIM %%%%
% Draw the edge of the Fibroblast

%%%% INPUTS %%%%
% 1 image of the fibrolast of the given file to analyze

%%%% OUTPUTS %%%%
% 'Radius' : matrix where the value of each cell is the distance to the
% Fibroblast edge ( >0 outside spheroid, <0 inside spheroid, =0 at the
% Fibroblast edge)
% Clément Hallopeau 03/2020


picture=imread(fullfile(path,file));
img = imadjust(rgb2gray(picture)) ;
    

%% Draw the CAF edge 
ask = 'n' ; 
next = 0 ; 
while next == 0  % while the user hasn't type 'y' when asking if it's satified : 


% Draw a polygon passing along the Fibroblast edge and along the edges of the image inside the spheroid.
% Before closing, move the 4 corners outside of the image
fig = figure('Position',[200,200,1500,900]) ; 
imshow(img) ;
title(sprintf('Delineate the Edge of Fibroblasts  \n Double click on a point to finish'))
Mask=roipoly;   % to draw, use the roipoly function from MATLAB
close(fig) ; 
 
Mask=double(Mask); % convert the mask as a double matrix
Mask = imbinarize(imgaussfilt(Mask,sigma)) ;   % smooth the mask with gaussian filter

 % plot the image and the edge after smoothing                    
fig = figure('Position', [200, 100, 600,600]) ; 
imshow(img);  ; hold on ;
contour(Mask,1, 'LineWidth', 3, 'Color', 'r')

ask = input('Satisfied by the spheroid edge ? [y/n]', 's') ; 

if  isempty(ask) | ask == 'y'  ; 
    next = 1 ; ;
end
if ask == 'n' ; 
   next = 0 ;  
end
if ask ~= 'y' & ask~='n'  & ~isempty(ask) ; 
    ask = input('Are you satisfied with the current output ? [y/n]', 's')
end

close(fig) ;
end


[Dista1,Linear]=bwdist(~Mask);  %   Distances from  inside the spheroid to the Fibroblast edge
                                        % D = BWDIST(BW) computes the Euclidean distance transform of the
                                       %   binary image BW. For each pixel in BW, the distance transform assigns
                                       %   a number that is the distance between that pixel and the nearest
                                       %   nonzero pixel of BW. ~ is the inverse of the mask. This
                                      
                                     
                            
[Dista2,Linear2]=bwdist(Mask); %   Distances from  outside the spheroid to the Fibroblast edge

Dista1=round(Dista1); % round to get distance in round pixel number.
Dista2=round(Dista2); % round to get distance in round pixel number.
Dista=Dista2-Dista1; % Distances inside the spheroid to the edge will be negative and those of the outside positive.
                     % the border will be 0.


Radius=Dista*px_to_micro;                           %convert Radius from pixel to micrometer

filename=erase(file,".tif");
Radiusname=strcat(filename, "_Radius.mat");

 if ~exist(fullfile(path,'Radius'))
    mkdir(fullfile(path, 'Radius')) ;
end

save(fullfile(path,'Radius',Radiusname),'Radius') ; 


end

