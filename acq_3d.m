function acq_3d(seq_path)
% Code used to acquire the 3-d kinematics of larval fish in 3d from 2
% camera views.  This code is designed with impulse chamber experiments in
% mind, but could be adapted to a variety of experimental situations.
% It assumes a right-handed coordinate system, where the z-axis points away
% from gravity and the x-y plane is horizontal.  It also assumes that the
% cameras are directed with perpendicular orientations.
% 
% Each sequence directory should contain the following directories:
%   'top view'   - video recording from the cam1 view as individual tiffs
%   'side view'  - video recording from the cam2 view, as individual tiffs
%
% Each sequence directory should contain the following calibration images:
%   'scale top.tif' 
%   'scale side.tif'


%% Parameters

% Suffix for video frame filenames (e.g. 'tif' or 'TIFF')
im_suf = 'tif';

% Top cam directory name
seq.dir_camT = 'top view';

% Side cam prefix
seq.dir_camS = 'side view';

% Top cam calibration filename
seq.fname_camT_cal = 'scale top.tif';

% Side cam calibration filename
seq.fname_camS_cal = 'scale side.tif';

% Max number of frames for creating the mean image
maxFrames = 1000;


%% Specifying sequence information 

% Prompt for directory
if nargin < 1
    seq_path = uigetdir(pwd,'Choose directory for a particular sequence');
end

if isempty(dir([seq_path filesep 'seq_info.mat']))
    
    % Check for sequence directories
    if isempty(dir([seq_path filesep seq.dir_camT]))
        error(['The sequence directory must have a directory called "'...
            dir_camT '"']);
        
    elseif isempty(dir([seq_path filesep seq.dir_camS]))
        error(['The sequence directory must have a directory called "'...
            dir_camS '"']);
    end
    
    % Look for images in video directories
    seq.aT = dir([seq_path filesep seq.dir_camT filesep '*.' im_suf]);
    seq.aS = dir([seq_path filesep seq.dir_camS filesep '*.' im_suf]);
    
    % Verify images
    if isempty(seq.aT)
        error(['No images identified in ' seq_path filesep seq.dir_camT])
    elseif isempty(seq.aS)
        error(['No images identified in ' seq_path filesep seq.dir_camS])
    end
    
    % Check frame number
    if length(seq.aT) > length(seq.aS)
        warning(['more frames from top view than side view: ' ...
            'trimming duration to shorter sequence'])
        
    elseif length(seq.aT) > length(seq.aS)
        warning(['more frames from side view than top view: '...
            'trimming duration to shorter sequence'])
    end
    
    % Specify number of frames
    seq.num_frames = min([length(seq.aS) length(seq.aT)]);
    
    % Check for calibration files
    if isempty(dir([seq_path filesep seq.fname_camT_cal]))
        [seq.fname_camT_cal,path,fIdx] = uigetdir(seq_path,...
            'Choose top view calibration file');
    end
    
    if isempty(dir([seq_path filesep seq.fname_camS_cal]))
        [seq.fname_camS_cal,path,fIdx] = uigetdir(seq_path,...
            'Choose side view calibration file');
    end
    
    % Prompt for sequence information
    answer = inputdlg({'Fish number','Sequence number'},...
        'Sequence information',1,{'1','1'});
    
    % Store sequence info
    seq.fish_num = str2num(answer{1});
    seq.seq_num  = str2num(answer{2});
    
    % Clean up
    clear answer im_suf
    
    save([seq_path filesep 'seq_info.mat'],'seq')
    
else
    disp(' ')
    disp(' Loading sequence info')
    load([seq_path filesep 'seq_info.mat'])
end


%% Create mean images

% Calculate mean image does not exist
if isempty(dir([seq_path filesep 'meanImage_top.tif'])) || ...
   isempty(dir([seq_path filesep 'meanImage_side.tif']))  

    pathT = [seq_path filesep seq.dir_camT];
    pathS = [seq_path filesep seq.dir_camS];
    
    % Define list of frame numbers, depending on max number of frames
    % requested
    if seq.num_frames > maxFrames
        dframe = floor(seq.num_frames/maxFrames);
        frIdx = 1:dframe:seq.num_frames;
        clear dframe
    else
        frIdx = 1:seq.num_frames;
    end
    
    % Create waitbar
    h = waitbar(0,...
            ['Mean image: ' num2str(1)],...
             'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    % Create sum image based on first frame
    [imCurr,tmp] = imread([pathT filesep seq.aT(frIdx(1)).name]);
    imSumT = double(imCurr);
    
    [imCurr,tmp] = imread([pathS filesep seq.aS(frIdx(1)).name]);
    imSumS = double(imCurr);
    
    clear imCurr tmp
    
    
    % Loop through frames 
    for i = 1:length(frIdx)
        
        % Add current frame to sum image
        [imCurr,tmp] = imread([pathT filesep seq.aT(frIdx(i)).name]);
        imSumT       = imSumT + double(imCurr);
        
        [imCurr,tmp] = imread([pathS filesep seq.aS(frIdx(i)).name]);
        imSumS       = imSumS + double(imCurr);
        
        clear tmp imCurr
        
        % Update status bar
        h = waitbar(i/length(frIdx),h,...
            ['Mean image: ' num2str(i) ' of ' num2str(length(frIdx)) ' frames']);
        
        % Quit m-file, if cancel button pushed
        if getappdata(h,'canceling')
            close force
            return
        end
        
    end
    
    % Calculate mean from sum images
    imMeanT = uint8(round(imSumT./(length(frIdx)+1)));
    imMeanT = imMeanT(:,:,1);
    
    imMeanS = uint8(round(imSumS./(length(frIdx)+1)));
    imMeanS = imMeanS(:,:,1);
    
    % Write image to movie dir
    imwrite(imMeanT,[seq_path filesep 'meanImage_top.tif'],'tif',...
            'Compression','none');
    imwrite(imMeanS,[seq_path filesep 'meanImage_side.tif'],'tif',...
            'Compression','none');
    
    close force
    clear frIdx h i imSumT imSumS
        
    %imMean = rgb2gray(imMean);
      
    
% Load mean image, if present
else
    
    disp(' ')
    disp('Loading mean images . . .');
    imMeanT = imread([seq_path filesep 'meanImage_top.tif']);
    imMeanS = imread([seq_path filesep 'meanImage_side.tif']);
    
end


%% Acquire calibration data

if isempty(dir([seq_path filesep 'cal_data.mat']))
    
    f = figure;
    
    disp(' ')
    disp('This code assumes a right-handed coordinate system.')
    disp(' The cameras should be arranged as follows:')
    disp(' ')
    disp(' Top cam: y-axis points down on screen')
    disp('          x-axis points left on screen')
    disp( ' ')
    disp(' Side cam: x-axis points left on screen')
    disp('           z-axis points up on screen')
    disp( ' ')
    
    % Calibrate magnification in top cam -------------------------------
    im_cal_camT = imread([seq_path filesep seq.fname_camT_cal]); 
    
    % Adjust contrast
    im_cal_camT = imadjust(im_cal_camT);
    
    cal.camT.const = calibrate(im_cal_camT,...
            'Top cam: choose 2 points over known disance, press return');
    
    clear im_cal_camT
      
     % Calibrate magnificantion in cam2 -------------------------------  
    im_cal_camS = imread([seq_path filesep seq.fname_camS_cal]); 
    
    % Adjust contrast
    im_cal_camS = imadjust(im_cal_camS);
    
    cal.camS.const = calibrate(im_cal_camS,...
              'Side cam: choose 2 points over known disance, press return');
    
     clear im_cal_camS
       
    % Select same point in the two views ------------------------------
    
    % Point selection
    [xT,yT,xS,yS] = select_common(seq_path,seq,imMeanT,imMeanS);
    cal.camT.pnt = [xT yT];
    cal.camS.pnt = [xS yS];
    
    % Save 'cal'
    save([seq_path filesep 'cal_data.mat'],'cal')
    
    clear im_cam1 im_cam2
    close
       
else
    disp(' ')
    disp('Loading calibration data . . .')
    disp(' ')
    
    % Load 'cal' structure
    load([seq_path filesep 'cal_data.mat'])
end


%% Establish directory structure for current sequence

% Prompt for sequence
if ~isempty(seqs)
   seq_path = uigetdir(seq_path,'Choose sequence to analyze');
   disp(' '); disp('Choose a sequence to analyze');
else
    error('No sequence directories present')
end

% Get video file info
a_camT = dir([seq_path filesep pre_camT '*.' im_suf]);
a_camS = dir([seq_path filesep pre_camS '*.' im_suf]);

% Check for images
if isempty(a_camT)
    error('No top image files present')
elseif isempty(a_camS)
    error('No side image files present')
elseif length(a_camT) ~= length(a_camS)
    error('There is an unequal number of side and top frames');
end


%% Prep data structure, set zoom

if isempty(dir([seq_path filesep 'coord_data.mat']))
    
    % Prompt for start and end frames
    prompt = {'Start frame','End frame','Frame rate'};
    defaults = {'1',num2str(length(a_camT)),'500'};
    
    answer = inputdlg(prompt,' ',1,defaults);
    
    d.start_frame = str2num(answer{1});
    d.end_frame   = str2num(answer{2}); 
    d.frame_rate  = str2num(answer{3});
    
    clear answer prompt defaults
    
    % Define number of frames
    d.frame_skip = 1; 
    d.frames     = d.start_frame:d.frame_skip:d.end_frame;
    d.num_frames = length(d.frames);
    
    % Create empty fields in 'd'
    d.top.xLim  = nan(d.num_frames,2);
    d.top.yLim  = nan(d.num_frames,2);
    d.side.xLim = nan(d.num_frames,2);
    d.side.yLim = nan(d.num_frames,2);
    
    d.top.leftEye   = nan(d.num_frames,2);
    d.top.rightEye  = nan(d.num_frames,2);
    d.side.leftEye  = nan(d.num_frames,2);
    d.side.rightEye = nan(d.num_frames,2);
    
    d.x_coord       = nan(d.num_frames,1);
    d.y_coord       = nan(d.num_frames,1);
    d.z_coord       = nan(d.num_frames,1);
    
    % Create figure window
    f = figure;
    %set(f,'WindowStyle','modal')
    set(f,'DoubleBuffer','on')
    set(f,'Name','Point selection')
    set(f,'Menubar','figure')
    
    % Read first frame
    imT = imread2([seq_path filesep a_camT(d.start_frame).name]);
    imS = imread2([seq_path filesep a_camS(d.start_frame).name]);
    
    % Prompt for zoom
    warning off
    hT(1) = imshow(imT);
    title('Top view: select zoom level, press return')
    add_labels(cal.camT)
    zoom on
    pause
    
    % Store top axes in 'd'
    d.top.xLim(1,:) = xlim;
    d.top.yLim(1,:) = ylim;
    
    % Prompt for zoom
    hS(1) = imshow(imS);
    title('Side view: select zoom level, press return')
    add_labels(cal.camS)
    zoom on
    pause
    
    % Store side axes in 'd'
    d.side.xLim(1,:) = xlim;
    d.side.yLim(1,:) = ylim;
    
    
    warning on
    clear hS hT   
    close
    
    save([seq_path filesep 'coord_data.mat'],'d')

else
    disp(' ')
    disp('Loading coordinate data . . .')
    disp(' ')
    
    % Load 'd' structure
    load([seq_path filesep 'coord_data.mat'])

end


%% Acquire coordinates

% Commands
disp('Commands ')
disp('   Left click   - collect coordinate')
disp('   Right click    - delete coordinate')
disp('   Return        - finish point collecting')
disp('   ESC           - stop acq_3d')
disp(' ')

% Define index number
idx = find(~isnan(d.top.leftEye(:,1)),1,'last');
if isempty(idx)
    idx = 1;
else
    idx = idx + d.frame_skip;
end

% Define frame limits
xlim_range_s = range(d.side.xLim(idx,:));
ylim_range_s = range(d.side.yLim(idx,:));
xlim_range_t = range(d.top.xLim(idx,:));
ylim_range_t = range(d.top.yLim(idx,:));

% Loop through frames from the two perspectives -------------------------
while true
    
    % Read image files for current time 
    cFrame = d.frames(idx);
    imT = imread2([seq_path filesep a_camT(idx).name]);
    imS = imread2([seq_path filesep a_camS(idx).name]);
     
    % Define ralative frames
    next_idx = min([idx+d.frame_skip d.num_frames]);
    last_idx = max([idx-d.frame_skip 1]);
    
    % Define handles
    hS = [];
    hT = [];
    
    skip_side = 0;
    
    %% ACQUIRE HEAD CENTER, TOP VIEW 
    
    % Acquire, only if no top point
    %if isnan(d.top.leftEye(idx,1))
        
        % Display side view
        warning off
        subplot(1,2,2)
        hS(1) = imshow(imS);
        
        % Zoom by setting xlim, ylim
        xlim(d.side.xLim(idx,:));
        ylim(d.side.yLim(idx,:));
        
        % Title and labels
        htxt(2) = title(['Side (frame ' num2str(idx) ...
               ' of ' num2str(d.num_frames) ')']);
        add_labels(cal.camS)
        
        % Display top view
        subplot(1,2,1)
        hT(1) = imshow(imT);
        xlim(d.top.xLim(idx,:));
        ylim(d.top.yLim(idx,:));
        
        % Draw boarder
        hold on
        bd(1) = plot([d.top.xLim(idx,1) d.top.xLim(idx,2) ...
                      d.top.xLim(idx,2) d.top.xLim(idx,1) ...
                      d.top.xLim(idx,1)],[d.top.yLim(idx,1) ...
                      d.top.yLim(idx,1) d.top.yLim(idx,2) ...
                      d.top.yLim(idx,2) d.top.yLim(idx,1)],'r-');
        set(bd(1),'LineWidth',3);
        hold off
               
%         hold on
%         xVals = [d.top.xLim(idx,1) d.top.xLim(idx,2) ...
%                  d.top.xLim(idx,2) d.top.xLim(idx,1) ...
%                  d.top.xLim(idx,1)];
%         yVals = [d.top.yLim(idx,1) d.top.yLim(idx,1) ...
%                  d.top.yLim(idx,2) d.top.yLim(idx,2) ...
%                  d.top.yLim(idx,1)];  
%              
%         plot(xVals,yVals,'w-')   
%         hold off
        
      
        %title(['Top (frame ' num2str(idx) ')'])
        htxt(1) = title(['Top: click on the left eye (frame ' num2str(idx) ...
               ' of ' num2str(d.num_frames) ')']);
        set(htxt(1),'Color','r')
        add_labels(cal.camT)
        
        warning on
        
        % Prompt to input a point
        [x_tmp,y_tmp,b_tmp] = ginput(1);
        
        % Left click
        if b_tmp ==1
            d.top.leftEye(idx,:) = [x_tmp y_tmp];
            
            % Adjust image limits for next frame
            y_cntr = y_tmp;
            y_min  = max([1 y_cntr-ylim_range_t/2]);
            y_max  = min([size(imT,1) y_cntr+ylim_range_t/2]);
            x_min  = max([1 x_tmp-xlim_range_t/2]);
            x_max  = min([size(imT,2) x_tmp+xlim_range_t/2]);
            
            d.top.xLim(next_idx,:) = [x_min x_max];
            d.top.yLim(next_idx,:) = [y_min y_max];
                         
            clear y_ctnr y_min y_max x_min x_max
            
        % Right click
        elseif b_tmp ==3
            d.top.leftEye(idx,:) = nan(1,2);
            d.side.leftEye(idx,:) = nan(1,2);
            
            skip_side = 1;
            
            % Jump back current frame number
            idx = last_idx;
            
        % 'Return'
        elseif isempty(b_tmp)
            close
            % Clear coord values
            clear x_tmp y_tmp b_tmp
            break
            
        % esc
        elseif b_tmp == 27
            return
        end
        
        % Clear coord values
        clear x_tmp y_tmp b_tmp
        delete(bd(1))
        set(htxt(1),'Color','k')
    
    
    %% ACQUIRE HEAD CENTER, SIDE VIEW 
     if ~skip_side   
        % Overlay selected point
        subplot(1,2,1)
        hold on

        hT(2) = plot(d.top.leftEye(idx,1),d.top.leftEye(idx,2),'r+');
        htxt(1) = title(['Top (frame ' num2str(idx) ...
                  ' of ' num2str(d.num_frames) ')']);
        set(htxt(1),'Color','k')
        hold off
        
        % Plot line on side view
        subplot(1,2,2)
        
        % Draw boarder
        hold on
        bd(2) = plot([d.side.xLim(idx,1) d.side.xLim(idx,2) ...
                      d.side.xLim(idx,2) d.side.xLim(idx,1) ...
                      d.side.xLim(idx,1)],[d.side.yLim(idx,1) ...
                      d.side.yLim(idx,1) d.side.yLim(idx,2) ...
                      d.side.yLim(idx,2) d.side.yLim(idx,1)],'r-');
        set(bd(2),'LineWidth',3);
        hold off
        
        
%         hold on
%         
%         xVals = [d.side.xLim(idx,1) d.side.xLim(idx,2) ...
%                  d.side.xLim(idx,2) d.side.xLim(idx,1) ...
%                  d.side.xLim(idx,1)];
%         yVals = [d.side.yLim(idx,1) d.side.yLim(idx,1) ...
%                  d.side.yLim(idx,2) d.side.yLim(idx,2) ...
%                  d.side.yLim(idx,1)];  
%              
%         plot(xVals,yVals,'w-')   
%         
%         %hS(2) = plot([x_side_tmp x_side_tmp],ylim,'r--');
%         hold off
        htxt(2) = title(['Side: click on the left eye (frame ' num2str(idx) ...
               ' of ' num2str(d.num_frames) ')']);
        set(htxt(2),'Color','r') 
        pause(.01)
        
        % Prompt for z-coordinate of point
        [x_tmp,y_tmp,b_tmp] = ginput(1);
        
        % If left click
        if b_tmp ==1
            d.side.leftEye(idx,:) = [x_tmp y_tmp];
            
            % Plot on side view
            subplot(1,2,2)
            hold on
            hS(3) = plot(d.side.leftEye(idx,1),...
                         d.side.leftEye(idx,2),'r+');
            hold off
            htxt(2) = title(['Side: (frame ' num2str(idx) ...
                   ' of ' num2str(d.num_frames) ')']);
            set(htxt(2),'Color','k')   
            
            % Adjust image limits for next frame
            y_cntr = y_tmp;
            y_min  = max([1 y_cntr-ylim_range_s/2]);
            y_max  = min([size(imS,1) y_cntr+ylim_range_s/2]);
            x_min  = max([1 x_tmp-xlim_range_s/2]);
            x_max  = min([size(imS,2) x_tmp+xlim_range_s/2]);
            
            d.side.xLim(next_idx,:) = [x_min x_max];
            d.side.yLim(next_idx,:) = [y_min y_max];
                         
            clear y_ctnr y_min y_max x_min x_max
            
            % Advance current frame number
            idx = next_idx;
            
        % If right click
        elseif b_tmp ==3
            d.side.leftEye(idx,:) = nan(1,2);
            d.top.leftEye(idx,:)  = nan(1,2);
            
            % Delete guideline for side view
            delete(hT(2))
            
        % If 'Return'
        elseif isempty(b_tmp)
            clear x_tmp y_tmp b_tmp
            close
            break
            
        %If esc
        elseif b_tmp == 27
            return
        end
        
        clear x_tmp y_tmp b_tmp
        delete(bd(2))
        set(htxt(2),'Color','k')
        
        % Calculate coordinates in global system
        %d.x_coord = cal.camT.const .* (d.top.leftEye(idx,1)-cal.camT.pnt(1));
        %d.y_coord = cal.camT.const .* (d.top.leftEye(idx,2)-cal.camT.pnt(2));
        %d.z_coord = cal.camS.const .* (d.side.leftEye(idx,2)-cal.camS.pnt(2));
     end
   % end
    
    % Clear variables for next loop
    clear imT imS hT hS next_idx last_idx
    
    % Save data
    save([seq_path filesep 'coord_data.mat'],'d')
    
    pause(.1)
    
    
end

clear xlim_range_s ylim_range_s xlim_range_t ylim_range_t



function im = give_im(imPath,invert,imMean,x_roi,y_roi)

% Load images
im = imread(imPath);

% Adjust grayscale values
im   = (imadjust(im));

%img = adapthisteq(img,'clipLimit',0.02,'Distribution','rayleigh');

% Subtract background
warning off
im = imsubtract(imadjust(imMean),im);
warning on

%im(find(im>255))  = 255;

if ~invert
    im_T = imcomplement(im_T);
    im_S = imcomplement(im_S);
end
% 
% % Use roi to crop image
% if (nargin > 3) && ~isempty(x_roi)
%     roiI = roipoly(im,x_roi,y_roi);
%     img = uint8(255.*ones(size(im,1),size(im,2)));
%     img(roiI(:)) = im(roiI(:));
%     im = img;
% end
%     im_T = uint8(255.*ones(size(im_T,1),size(im_T,2)));
%     im_S = uint8(255.*ones(size(im_S,1),size(im_S,2)));
% 
%     


function frame_to_global(cal,pts,cam_num)
% Transforms video frame coordinates into global coordinate

% Find common axis dimension
if isfield(cal.cam1,'x') && isfield(cal.cam2,'x')
    com_dim = 1;
elseif isfield(cal.cam1,'y') && isfield(cal.cam2,'y')
    com_dim = 2;
elseif isfield(cal.cam1,'z') && isfield(cal.cam2,'z')
    com_dim = 3;
end

function global_to_frame(cal,pts,cam_num)
% Transforms global coordinates into video frame coordinates


function S = localSystem(P1,P2,P3)
% Defines a transformation vector for a local coordinate system in an
% inertial frame of reference.  Uses P1 as the xaxis and P2 as the origin, and 
% P3 as the z-axis. Coordinates must be (1x3) vectors. Note: if theses axes 
% are not orthogonal, the z-axis direction is assumed to be more accurate
% than the x-axis and the x-axis direction is adjusted to make the coordinates 
% orthoganal.
 
% Check dimensions of inputs
if size(P1,1)~=1 || size(P1,2)~=3 ||...
   size(P2,1)~=1 || size(P2,2)~=3 ||...
   size(P3,1)~=1 || size(P3,2)~=3
    error('Coordinates must be 1x3 vectors');
end
 
% Define units vectors for x and y axes
xAxis   = (P1-P2)./norm(P1-P2);
zAxis   = (P3-P2)./norm(P3-P2);
 
% Define yaxis from the cross product
yAxis   = cross(zAxis,xAxis);
yAxis   = yAxis./norm(yAxis);
 
% Redefine the xaxis, so all axes are orthoganal
xAxis   = cross(yAxis,zAxis);
 
% Define transformation matrix
S       = [xAxis' yAxis' zAxis'];
 
function [xn,yn,zn] = localToGlobal(x,y,z,origin,S)
% Transforms coordinates from the local coordinate system to the global
% system. Coordinates may be given as nxm matricies of equal dimensions.
 
% Check dimensions of inputs
if ~( (size(origin,1)==1 && size(origin,2)==3) ||...
      (size(origin,1)==3 && size(origin,2)==1) )   
    error('Origin must be a 1x3 or 3x1 vector');
    
elseif size(S,1)~=3 || size(S,2)~=3
    error('S must be a 3x3 matrix');
    
elseif ~min(size(x)==size(y)) || ~min(size(x)==size(z)) || ~min(size(y)==size(z))
    error('x, y, & z must have the same dimensions')
    
end
 
% Loop through column to complete transformation
for i = 1:size(x,2)
    pts     = [x(:,i) y(:,i) z(:,i)];
    pts     = [inv(S)'*pts']';
 
    xn(:,i) = pts(:,1) + origin(1);
    yn(:,i) = pts(:,2) + origin(2);
    zn(:,i) = pts(:,3) + origin(3);
    
    clear pts 
end
 
function [xn,yn,zn] = globalToLocal(x,y,z,origin,S)
% Transforms coordinates from the global coordinate system to the local
% system. Coordinates may be given as nxm matricies of equal dimensions.
 
% Check dimensions of inputs
if ~( (size(origin,1)==1 && size(origin,2)==3) ||...
      (size(origin,1)==3 && size(origin,2)==1) )       
    error('Origin must be a 1x3 or 3x1 vector');
       
elseif size(S,1)~=3 || size(S,2)~=3
    error('S must be a 3x3 matrix');
    
elseif ~min(size(x)==size(y)) || ~min(size(x)==size(z)) || ~min(size(y)==size(z))
    error('x, y, & z must have the same dimensions')
    
end
 
% Loop through column to complete transformation
for i = 1:size(x,2)
    pts         = [x(:,i) y(:,i) z(:,i)];    
    pts(:,1)    = x(:,i)-origin(1);
    pts(:,2)    = y(:,i)-origin(2);
    pts(:,3)    = z(:,i)-origin(3);
    pts         = [S'*pts']';
    
    xn(:,i)     = pts(:,1);
    yn(:,i)     = pts(:,2);
    zn(:,i)     = pts(:,3);
    
    clear pts
end

function add_labels(cam)
% Adds axes to current view, according to calibration data

if isfield(cam,'xaxis')
    if cam.xaxis == [1 0]
        xlabel('x ->')
    elseif cam.xaxis == [-1 0]
        xlabel('<- x')
    elseif cam.xaxis == [0 -1]
        ylabel('x ->')
    elseif cam.xaxis == [0 1]
        ylabel('<- x')
    end
end

if isfield(cam,'yaxis')
    if cam.yaxis == [1 0]
        xlabel('y ->')
    elseif cam.yaxis == [-1 0]
        xlabel('<- y')
    elseif cam.yaxis == [0 -1]
        ylabel('y ->')
    elseif cam.yaxis == [0 1]
        ylabel('<- y')
    end
end

if isfield(cam,'zaxis')
    if cam.zaxis == [1 0]
        xlabel('z ->')
    elseif cam.zaxis == [-1 0]
        xlabel('<- z')
    elseif cam.zaxis == [0 -1]
        ylabel('z ->')
    elseif cam.zaxis == [0 1]
        ylabel('<- z')
    end
end

