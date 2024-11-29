function [] =   ECM_barrier(path)
 %%%% AIM %%%%
% Find the width of the ECM edges

%%%% INPUTS %%%%
% The intensity profile (PROFILE) of each ECM image

%%%% OUTPUTS %%%%
% PROFILE.mat, where the starting and ending point of the ECM edge are
% added, in microns.
% Clément Hallopeau 05/2021

cd(path) ;   %move into the subfolder
imgs = dir('*.tif') ; %extract the name of all the images (used to call the columns of PROFILE.mat)

cd([path, filesep, 'Matlab_data']) ; % move into the matlab_data folder
load('PROFILE') ;  %load the intensity profile
load('Radius') ;   %Radius

for i = 1:length(imgs)  %for all the images, ie for all the columns of PROFILE.mat
   disp(imgs(i).name(1:3))
    if sum(imgs(i).name(1:3) ~= 'CAF') > 0   % if the image is not CAFs (we don't want to find the edge width of the CAFs)


	MEAN = [imgs(i).name(1:end-4),'_MEAN'] ; 
    STD = [imgs(i).name(1:end-4),'_STD'] ; 
    
    R = PROFILE.R ;  

    % The width of the edge is obtained by finding the change points of the
    % intensity profile. These changepoints usually corresond to the
    % beginning, the top and the end of the intensity peak due to the
    % accumulation of ECM.
    ipt = findchangepts(PROFILE.(MEAN), 'MaxNumChanges', 3, 'Statistic', 'linear') ;  
    ipt(1) = ipt(2) ; ipt(3) = round((ipt(3)+ipt(2))/2) ;  % MATLAB function find changepoint always a bit decayed to the right of the profile, we shift them to the left
    
    mask = zeros(size(Radius)) ; 
    mask(find(Radius > PROFILE.R(ipt(1)) & Radius <     PROFILE.R(ipt(3))   )) = 1 ; 
    maskCAF = zeros(size(Radius)) ; maskCAF(find(Radius > 0 & Radius < 120)) = 1 ; 
    
    % Plot
    fig = figure('Position', [500,400,1800,600]) ; 
    subplot(121) ; cd(path) ; 
    imshow(imresize(imread(imgs(i).name), size(maskCAF))) ; hold on ;
    contour(maskCAF,1,':', 'LineWidth',3) ;  hold on ;
    contour(mask,1, 'LineWidth', 2) ; 
    title([imgs(i).name(1:end-4),' barrier =  ', num2str(PROFILE.R(ipt(1)) - PROFILE.R(ipt(3))), 'µm']) ; 
    
    subplot(122)
    patch([transpose(R), fliplr(transpose(R))], [transpose(PROFILE.(MEAN) - PROFILE.(STD)), fliplr(transpose(PROFILE.(MEAN)+PROFILE.(STD)))], [.9 .9 .9], 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'HandleVisibility', 'off') ; hold on ;
    plot(R, PROFILE.(MEAN), 'DisplayName', 'Mean') ; hold on ;
    plot([PROFILE.R(ipt(1)),PROFILE.R(ipt(1))], [0,1], 'Color', 'k', 'LineStyle', '--') ; hold on ;     plot([PROFILE.R(ipt(3)),PROFILE.R(ipt(3))], [0,1], 'Color', 'k', 'LineStyle', '--') ; hold on ;
    title([imgs(i).name(1:end-4), ' Intensity Profile'])
    legend ;     xlabel('Distance to the Spheroid Edge (µm)') ;    ylabel('Normalized Fluorescence Intensity (a.u.)') ;    grid on

    PROFILE.([imgs(i).name(1:end-4),'_EDGE']) = transpose([PROFILE.R(ipt(1)), PROFILE.R(ipt(3)), zeros(1,length(PROFILE.R)-2)])  ;

    end
end

cd([path, filesep, 'Matlab_data']) ; 
save('PROFILE.mat', 'PROFILE') ; 
saveas(fig,['Edge-',imgs(i).name(1:end-4),'.png']) ; 
close(fig) ;

end 
