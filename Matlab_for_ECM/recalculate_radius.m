
% 20240627 recalculate radius.mat after having downsized the angle map, or
% the CAF.tif

codepath =  '/home/challopeau/Documents/PhD/Codes/Paper_MIRO/Matlab_for_ECMV3.2_forAnna_last';
cd(codepath) ;

PATH = '/home/challopeau/Downloads/ZS - ALLEXP - SR' ;

folders  = dir(PATH) ; folders = folders(3:end)   ; folders = folders([folders.isdir]) ;% exclude the first 2 folders/files which are always '.' and '..'

   px_size = 0.3225 ; % pixel size, in microns per pixel, for 2048*2048 field of views

for i = 1:length(folders) 
    temp_path = [PATH, '/', folders(i).name] ;

    subfolders =  dir(temp_path) ; subfolders = subfolders(3:end)   ; subfolders = subfolders([subfolders.isdir]) ;% exclude the first 2 folders/files which are always '.' and '..'
    
        for j = 1:length(subfolders)
        path = [temp_path, '/', subfolders(j).name] ; 
        disp(path) 

                     px_to_micro = 0.3225 * 2048/size(imread([path,filesep, 'CAF.tif']),1) ;   %convert pixel size considering the image binning

        cd([path, filesep, 'Matlab_data']) ; % move into the matlab_data folder
        
        Radius = load('Radius.mat') ;   Radius = Radius.Radius ;%Radius
        
        Mask = zeros(size(Radius)) ; 
        Mask(find(Radius <0)) = 1 ; 
        Mask_small = imresize(Mask, [1024,1024]) ;   % 2024/06/27 downsize the mask by 2, otherwise it's too long to smooth. Oriol will have to downsize the CAF image by 2 as well
        Mask = zeros(size(Mask_small)) ; 
        Mask(find(Mask_small ==1)) = 1 ;
                
        
        [Dista1,Linear]=bwdist(~Mask); %   Distances from  inside the spheroid to the Fibroblast edge
                                                % D = BWDIST(BW) computes the Euclidean distance transform of the
                                               %   binary image BW. For each pixel in BW, the distance transform assigns
                                               %   a number that is the distance between that pixel and the nearest
                                               %   nonzero pixel of BW. ~ is the inverse of the mask. 
                                             
                                    
        [Dista2,Linear2]=bwdist(Mask); % distances for the outside of the mask to the border of the mask.
                
        Dista1=round(Dista1); % round to get distance in round pixel number.
        Dista2=round(Dista2); % round to get distance in round pixel number.
        Dista=Dista2-Dista1; % Distances inside the spheroid to the edge will be negative and those of the outside positive.
                             % the border will be 0.
        
        Radius=Dista*px_to_micro;                           %convert Radius from pixel to micrometer
        
        
        cd([path, filesep, 'Matlab_data']) ; save('Radius','Radius') ; 

        end
end
