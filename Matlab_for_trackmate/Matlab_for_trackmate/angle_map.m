function [Vecfieldx, Vecfieldy] = angle_map(Mask)



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

% we assign to each point the closest vector
%of the leading edge
Front(sub2ind(size(Mask),Coord(2,:),Coord(1,:)))=1;
[DF,LF]=bwdist(Front);  %DF = distance matrix, LF = index of the closest pixel
[RF,CF]=ind2sub(size(Front),LF);  % RF/CF = Y/X coordinates of the closest edge pixel
DF=round(DF);   %distance matrix, positive (even inside the Mask)
%     DF(Dista==0)=0;
wb = waitbar(0) ;
for ni=1:max(max(DF));
    waitbar(ni/max(max(DF)), wb, 'Computing Normal Angles to the Edge')
    Tempx(find(DF==ni))=Tempx(LF(find(DF==ni)));   % Matrices whose size is the size of the image, 
    % each cell of the matrices contains the x or y coordinates of the
    % closest normal vector
    Tempy(find(DF==ni))=Tempy(LF(find(DF==ni)));
end

Vecfieldx  = Tempx ;;
Vecfieldy = Tempy ;;

close(wb) ;
end

