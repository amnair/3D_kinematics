function [S,idx] = autotrack_eye(seq_path,seq,imMean,S,idx,frames,cam_view,eye)
% Called by select_eyes.m to autotrack eyes



%% Parameters

invert = 1;


%% Find two blobs closest to selected points

% Define image path
if strcmp(cam_view,'top')
    im_path = [seq_path filesep seq.dir_camT filesep seq.aT(idx).name];
else
    im_path = [seq_path filesep seq.dir_camS filesep seq.aS(idx).name];
end

% load image with background subtracton
im = give_im(im_path,invert,imMean);

% roi binary
imROI   = roipoly(im,S.x_roi(idx,:),S.y_roi(idx,:));

% Threshold image
imBW = im2bw(imcomplement(im),S.tVal(idx));

% Crop image
imBW    = ~imBW & imROI;

% Get boundary data
[B,L] = bwboundaries(imBW,'noholes');
 
if isempty(B)
    warning('Please adjust the threshold')
    return
end

% Calc dist of each blob to each eye
for i = 1:length(B)
    cntrs(i,:) = [mean(B{i}(:,2)) mean(B{i}(:,1))];
    dist_L(i) =  sqrt((cntrs(i,1)-S.leftEye(idx,1)).^2 + ...
                      (cntrs(i,2)-S.leftEye(idx,2)).^2);
    dist_R(i) =  sqrt((cntrs(i,1)-S.rightEye(idx,1)).^2 + ...
                      (cntrs(i,2)-S.rightEye(idx,2)).^2); 
    perim(i) = length(B{i});              
end

% Select closest centers
idx_R = find(dist_R==min(dist_R),1,'first');
idx_L = find(dist_L==min(dist_L),1,'first');

% Make sure centers aren't the same
if idx_R == idx_L
    
    warning(['Left and right eyes at same point: ' ...
             'Reselect eyes and/or adjust threshold']);
    return

% Store results
else
    
    S.leftEye(idx,:)   = cntrs(idx_L,:);
    S.rightEye(idx,:)  = cntrs(idx_R,:);
    left_perim  = perim(idx_L);
    right_perim = perim(idx_R);
    
end

clear cntrs perim idx_R idx_L B L


%% Loop trhough next frames

% Define prior coordinates
L_next = [S.leftEye(idx,1) S.leftEye(idx,2)];
R_next = [S.rightEye(idx,1) S.rightEye(idx,2)];

% Make figure
f2 = figure('DoubleBuffer','on','CurrentCharacter','1');

% Advance frame index
idx = idx + 1;

while true
    
    % Center roi around eyes from last frame
    S.x_roi(idx,:) = S.x_roi(idx-1,:)-mean(S.x_roi(idx-1,:)) + ...
                     mean([L_next(1) R_next(1)]);
    S.y_roi(idx,:) = S.y_roi(idx-1,:)-mean(S.y_roi(idx-1,:)) + ...
                     mean([L_next(2) R_next(2)]);

    % Define image path
    if strcmp(cam_view,'top')
        im_path = [seq_path filesep seq.dir_camT filesep seq.aT(idx).name];
    else
        im_path = [seq_path filesep seq.dir_camS filesep seq.aS(idx).name];
    end
    
    % load image with background subtracton
    im = give_im(im_path,invert,imMean);
    
    % roi binary
    imROI   = roipoly(im,S.x_roi(idx,:),S.y_roi(idx,:));
    
    % Threshold image
    imBW = im2bw(imcomplement(im),S.tVal(idx));
    
    % Crop image
    imBW    = ~imBW & imROI;
    
    % Get boundary data
    [B,L] = bwboundaries(imBW,'noholes');
    
    % Check for at leat 2 blobs
    if isempty(B) || length(B)<2
        disp(' ')
        disp('autotracking failed')
        break
    end
    
    clear imBW imROI
    
    % List blobs of most similar periphery
    for i = 1:length(B)
        diff_perim(i) = abs(length(B{i})-left_perim);
    end
    
    % Sort blobs by similarity to past periphery
    [tmp,i_sort] = sort(diff_perim);
    
    % Keep only 3 blobs of similar periphery
    B_tmp{1} = B{i_sort(1)};
    B_tmp{2} = B{i_sort(2)};
    if length(B)>2
        B_tmp{3} = B{i_sort(3)};
    end
    B = B_tmp;
    
    clear tmp i_sort diff_perim B_tmp
    
    % Define centers of the three similar blobs
    cntrs(1,:) = [mean(B{1}(:,2)) mean(B{1}(:,1))];
    cntrs(2,:) = [mean(B{2}(:,2)) mean(B{2}(:,1))];
    if length(B)>2
        cntrs(3,:) = [mean(B{3}(:,2)) mean(B{3}(:,1))];
    end
    
    % Difference in left position btwn prior coordinates and current centers
    L_dist = sqrt( (cntrs(:,1)-L_next(1)).^2 + (cntrs(:,2)-L_next(2)).^2 );
    
    % Assign blob closest to left eye
    L_idx = find(L_dist==min(L_dist),1,'first');
    
    % Difference in right eye position btwn prior coordinates and current centers
    R_dist = sqrt( (cntrs(:,1)-R_next(1)).^2 + (cntrs(:,2)-R_next(2)).^2 );
    
    % Assign blob closest to left eye
    R_idx = find(R_dist==min(R_dist),1,'first');
     
    % Store results
    S.leftEye(idx,:)   = cntrs(L_idx,:);
    S.rightEye(idx,:)  = cntrs(R_idx,:);
    
    clear cntrs R_dist L_dist
    
    
    % Perimeter values for next loop
    left_perim  = length(B{L_idx});
    right_perim = length(B{R_idx});
    
    % Display frame
    figure(f2)
    warning off
    imshow(im)
    title(['Frame ' num2str(frames(idx)) '/' num2str(frames(end)) ...
           '  Hit key to stop'])
    warning on
    hold on
    
    for i = 1:length(B)
        
        plot(B{R_idx}(:,2),B{R_idx}(:,1),'r')
        plot(S.rightEye(idx,1),S.rightEye(idx,2),'r+')
        
        plot(B{L_idx}(:,2),B{L_idx}(:,1),'g')
        plot(S.leftEye(idx,1),S.leftEye(idx,2),'g+')
        
        axis([min(S.x_roi(idx,:)) max(S.x_roi(idx,:)) ...
            min(S.y_roi(idx,:)) max(S.y_roi(idx,:))]);
        pause(.5)
    end
    
    % Define predicted next coordinates
    L_diff = [S.leftEye(idx,1)-S.leftEye(idx-1,1) ...
              S.leftEye(idx,2)-S.leftEye(idx-1,2)];
    R_diff = [S.rightEye(idx,1)-S.rightEye(idx-1,1) ...
              S.rightEye(idx,2)-S.rightEye(idx-1,2)];      
    L_next = S.leftEye(idx,:) + L_diff;
    R_next = S.rightEye(idx,:) + R_diff;
    
    clear L_diff R_diff
    
    % Check if at the end of the movie
    if idx == length(frames)
        break
    else
        idx = idx + 1;
    end
    
    clear B L cntrs L_diff L_idx R_idx im
    
    % Stop if button pressed
    if ~strcmp(get(f2,'CurrentCharacter'),'1') 
        disp(' ')
        disp('Autotracking stopped')
        break
        
    elseif ~isnan(S.leftEye(idx,1))
        disp(' ')
        disp('You have reached coordinates already defined')
        break
        
    end
    
end

close(f2)
 