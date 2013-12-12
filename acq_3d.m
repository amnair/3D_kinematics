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
%   'scale T.tif' 
%   'scale S.tif'


%set(0,'DefaultFigurePointer','arrow')

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
            seq.dir_camT '"']);
        
    elseif isempty(dir([seq_path filesep seq.dir_camS]))
        error(['The sequence directory must have a directory called "'...
            seq.dir_camS '"']);
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
        [seq.fname_camT_cal,path,fIdx] = uigetfile([seq_path filesep '*.*'],...
            'Choose top view calibration file');
    end
    
    if isempty(dir([seq_path filesep seq.fname_camS_cal]))
        [seq.fname_camS_cal,path,fIdx] = uigetfile([seq_path filesep '*.*'],...
            'Choose side view calibration file');
    end
    
    % Prompt for sequence information
    answer = inputdlg({'Fish number','Sequence number','Frame rate (fps)'},...
        'Sequence information',1,{'1','1','1000'});
    
    % Store sequence info
    seq.fish_num   = str2num(answer{1});
    seq.seq_num    = str2num(answer{2});
    seq.frame_rate = str2num(answer{3});
    
    % Clean up
    clear answer im_suf
    
    save([seq_path filesep 'seq_info.mat'],'seq')
    
else
    disp(' ')
    disp(' Loading sequence info')
    load([seq_path filesep 'seq_info.mat'])
end

clear im_suf 


%% Create mean images

% Calculate mean image does not exist
if isempty(dir([seq_path filesep 'meanImage_T.tif'])) || ...
   isempty(dir([seq_path filesep 'meanImage_S.tif']))  

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
    imwrite(imMeanT,[seq_path filesep 'meanImage_T.tif'],'tif',...
            'Compression','none');
    imwrite(imMeanS,[seq_path filesep 'meanImage_S.tif'],'tif',...
            'Compression','none');
    
    close force
    clear frIdx h i imSumT imSumS
        
    %imMean = rgb2gray(imMean);
      
    
% Load mean image, if present
else
    
    disp(' ')
    disp('Loading mean images . . .');
    imMeanT = imread([seq_path filesep 'meanImage_T.tif']);
    imMeanS = imread([seq_path filesep 'meanImage_S.tif']);
    
end

clear maxFrames


%% Acquire/load calibration data

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


%% Prep data structure, set zoom 

if isempty(dir([seq_path filesep 'coord_data.mat']))
    
    % Define start and end frames  
    d.start_frame = 1;
    d.end_frame   = seq.num_frames; 
    d.frame_rate  = seq.frame_rate;
    
    clear answer prompt defaults
    
    % Define frames
    d.frame_skip = 1; 
    d.frames     = d.start_frame:d.frame_skip:d.end_frame;
    d.num_frames = length(d.frames);
    
    % Create empty fields in 'd'
    d.T.x_roi  = nan(d.num_frames,5);
    d.T.y_roi  = nan(d.num_frames,5);
    d.S.x_roi  = nan(d.num_frames,5);
    d.S.y_roi  = nan(d.num_frames,5);
    
    d.T.mod_roi = zeros(d.num_frames,1);
    d.S.mod_roi = zeros(d.num_frames,1);
    
    d.T.leftEye   = nan(d.num_frames,2);
    d.T.rightEye  = nan(d.num_frames,2);
    d.T.SB        = nan(d.num_frames,2);
    
    d.S.leftEye   = nan(d.num_frames,2);
    d.S.rightEye  = nan(d.num_frames,2);
    d.S.SB        = nan(d.num_frames,2); 
    
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
    imT = imread([seq_path filesep seq.dir_camT filesep seq.aT(d.start_frame).name]);
    imS = imread([seq_path filesep seq.dir_camS filesep seq.aS(d.start_frame).name]);
    
     % Prompt for zoom
     warning off
%     hT(1) = imshow(imT);
%     title('Top view: select zoom level, press return')
%     add_labels(cal.camT)
%     zoom on
%     pause

    hT(1) = imshow(imT);
    title('Top view: select zoom level, press return')
    add_labels(cal.camT)
    [x_tmp,y_tmp,h_im] = select_roi;
    
    d.T.x_roi = repmat(x_tmp,size(d.T.x_roi,1),1);
    d.T.y_roi = repmat(y_tmp,size(d.T.y_roi,1),1);
    
    % Determine default threshold value
    x_tmp = floor(min(d.T.x_roi(1,:))):ceil(max(d.T.x_roi(1,:)));
    y_tmp = floor(min(d.T.y_roi(1,:))):ceil(max(d.T.y_roi(1,:)));
    im_tmp = imT(y_tmp,x_tmp);
            
    d.T.tVal = ones(size(d.T.y_roi,1),1).*graythresh(im_tmp)/2;
    d.T.mod_tVal = zeros(size(d.T.y_roi,1),1);
    
    clear h_im x_tmp y_tmp im_tmp
    
    % Prompt for zoom
    hS(1) = imshow(imS);
    title('Side view: select zoom level, press return')
    add_labels(cal.camS)
    [x_tmp,y_tmp,h_im] = select_roi;
    
    d.S.x_roi = repmat(x_tmp,size(d.S.x_roi,1),1);
    d.S.y_roi = repmat(y_tmp,size(d.S.y_roi,1),1);
    
    % Determine default threshold value
    x_tmp = floor(min(d.S.x_roi(1,:))):ceil(max(d.S.x_roi(1,:)));
    y_tmp = floor(min(d.S.y_roi(1,:))):ceil(max(d.S.y_roi(1,:)));
    im_tmp = imS(y_tmp,x_tmp);

    d.S.tVal = ones(size(d.S.y_roi,1),1).*graythresh(im_tmp)/2;
    d.S.mod_tVal = zeros(size(d.S.y_roi,1),1);
    
    clear x_roi y_roi h_im hS hT h_im x_tmp y_tmp
   
    warning on  
    close
    
    % Set current frame index
    d.curr_frame_idx = 1;
    
    save([seq_path filesep 'coord_data.mat'],'d')

else
    disp(' ')
    disp('Loading coordinate data . . .')
    disp(' ')
    
    % Load 'd' structure
    load([seq_path filesep 'coord_data.mat'])

end

if isempty(dir([seq_path filesep 'tail_data.mat']))
    % Set up structure to hold tail coordinate data 
    t.S.dorsalEdge = cell(seq.num_frames,2);
    t.S.ventralEdge = cell(seq.num_frames,2);
    t.T.dorsalEdge = cell(seq.num_frames,2);
    t.T.ventralEdge = cell(seq.num_frames,2);
    
    % Save tail coordinate data
    save([seq_path filesep 'tail_data.mat'],'t')
else
    disp(' ')
    disp('Loading tail coordinate data . . .')
    disp(' ')
    
    % Load 'd' structure
    load([seq_path filesep 'tail_data.mat'])
    
end


%% Acquire coordinates


select_eyes(seq_path,d,t,seq,cal,imMeanT,imMeanS)


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

function [x_roi,y_roi,h_im] = select_roi
hold on

% Prompt for first point
disp(' ');
disp('Select first point')
[x1,y1,but1] = ginput(1);
h1 = plot(x1,y1,'m+');

% Prompt for second
disp(' ');
disp('Select second point')
[x2,y2,but2] = ginput(1);

% Define roi
delete(h1)
x_roi = [x1 x2 x2 x1 x1];
y_roi = [y1 y1 y2 y2 y1];

% Display briefly
h2 = plot(x_roi,y_roi,'m-');
pause(.3)
delete(h2)

% Specify current figure handle
h_im = gca;

hold off


        


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

