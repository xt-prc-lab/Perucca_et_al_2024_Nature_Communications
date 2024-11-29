codepath =  '/home/challopeau/Documents/PhD/Codes/Paper_MIRO/Matlab_for_ECMV3.2_forAnna_last';
cd(codepath) ;




PATH = '/home/challopeau/Downloads/ZS - ALLEXP - SR' ;


folders  = dir(PATH) ; folders = folders(3:end)   ; folders = folders([folders.isdir]) ;% exclude the first 2 folders/files which are always '.' and '..'



for i = 1:length(folders) 
    temp_path = [PATH, '/', folders(i).name] ;

    subfolders =  dir(temp_path) ; subfolders = subfolders(3:end)   ; subfolders = subfolders([subfolders.isdir]) ;% exclude the first 2 folders/files which are always '.' and '..'
    
        for j = 1:length(subfolders)
        path = [temp_path, '/', subfolders(j).name] ; 
        disp(path) 
             

        cd([path, filesep, 'Matlab_data']) ; % move into the matlab_data folder
        
        Radius = load('Radius.mat') ;   Radius = Radius.Radius ;%Radius
        
        Mask = zeros(size(Radius)) ; 
        Mask(find(Radius <0)) = 1 ; 
        Mask_small = imresize(Mask, [1024,1024]) ;   % 2024/06/27 downsize the mask by 2, otherwise it's too long to smooth. Oriol will have to downsize the CAF image by 2 as well
        Mask = zeros(size(Mask_small)) ; 
        Mask(find(Mask_small ==1)) = 1 ;
        
        %%  For each pixel of the field of view, compute the vector normal to the closest edge
        
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
%             for n=1:round(max(max(DF)))  % DF = matrix where each cell contain the distance to the closest edge of the corresponding pixel
%                 waitbar(n/round(max(max(DF))), wb, 'Smoothing the angle map') ; 
%                 [im,jm]=ind2sub(size(DF),find(DF==n));  % im, jm = x,y coordinates of the pixels located at a distance "n" from the edge
%                 for m=1:size(im,1) % for each of these pixel
%                     in=im(m);
%                     jn=jm(m);
%                     temp = Tempx(max(in-DF(in,jn),1):min(in+DF(in,jn),size(DF,1)),...
%                         max(jn-DF(in,jn),1):min(jn+DF(in,jn),size(DF,2)));
%                     Vecfieldx(in,jn)=nanmean(nanmean(temp));
%                     temp = Tempy(max(in-DF(in,jn),1):min(in+DF(in,jn),size(DF,1)),...
%                         max(jn-DF(in,jn),1):min(jn+DF(in,jn),size(DF,2)));
%                     Vecfieldy(in,jn)=nanmean(nanmean(temp));
%                 end
%             end
        %     close(wb) ;
            
        %     Vecfield = Vecfieldx ; Vecfield(:,:,2) = Vecfieldy ; 
            cd([path, filesep, 'Matlab_data']) ;
        %     Vecfieldx = Vecfield(:,:,1) ; 
        %     Vecfieldy = Vecfield(:,:,2) ;
            save('Vecfieldx.mat', 'Vecfieldx') ;
            save('Vecfieldy.mat', 'Vecfieldy') ;
        
    end
end

% for MATLAB versions < R2017, the vecnorm function do not exist, we define
% it here
function p = vecnorm(p,n,d)
    p = sqrt(sum(p.^2,d));
end


%% 


PATH = '/home/challopeau/Downloads/ZS - ALLEXP - SR' ;


folders  = dir(PATH) ; folders = folders(3:end)   ; folders = folders([folders.isdir]) ;% exclude the first 2 folders/files which are always '.' and '..'


for i = 1:length(folders) 
    temp_path = [PATH, '/', folders(i).name] ;

    subfolders =  dir(temp_path) ; subfolders = subfolders(3:end)   ; subfolders = subfolders([subfolders.isdir]) ;% exclude the first 2 folders/files which are always '.' and '..'
    
        for j = 1:length(subfolders)
            path = [temp_path, '/', subfolders(j).name] ; 
            disp(path) 

%         temp = load([path, filesep, 'Matlab_data', filesep, 'Vecfield.mat']) ;
%         Vecfield = temp.Vecfieldx_filt ;
%         Vecfield(:,:,2) = temp.Vecfieldy_filt ; 
%         save([path, filesep, 'Matlab_data', filesep, 'Vecfield_4matlab.mat'], 'Vecfield') ;
% 
            clear Vecfield 
            load([path, filesep, 'Matlab_data', filesep, 'Vecfield_4matlab.mat']) ;
            save([path, filesep, 'Matlab_data', filesep, 'Vecfield.mat'], 'Vecfield') ;

        end
end





