function acq_3d(root_path)
% Code used to acquire the 3-d kinematics of larval fish in 3d from 2
% camera views.  This code is designed with impulse chamber experiments in
% mind, but could be adapted to a variety of experimental situations.
% It assumes a right-handed coordinate system, where the z-axis points away
% from gravity and the x-y plane is horizontal.  It also assumes that the
% cameras are directed with perpendicular orientations.
% 
% Each sequence directory should contain the following directories:
%   cal    - should have 2 video frames, each with an image of a
%            ruler/micrometer from the two camera views ('cam1_ruler.tif' 
%            & 'cam2_ruler.tif') and 2 video frames, each including a sharp
%            point that can be viewed from both cams
%            ('cam1_pnt.tif','cam2_pnt.tif')
%   cam1   - video recording from the cam1 view as individual tiffs
%   cam2    - video recording from the cam2 view, as individual tiffs
%


%% General parameters

% Suffix for video frame filenames (e.g. 'tif' or 'TIFF')
im_suf = 'tif';

% Top cam  prefix
pre_camT = 'cam_top_frame_';

% Side cam prefix
pre_camS = 'cam_side_frame_';

% Top cam calibration filename
camT_cal_name = 'cam_top_ruler_001.tif';

% Side cam calibration filename
camS_cal_name = 'cam_side_ruler_001.tif';

% Top cam view of image of same point
camT_pnt_name = 'cam_top_pnt_001.tif';

% Side cam view of image of same point
camS_pnt_name = 'cam_side_pnt_001.tif';


%% Establish directory for the group of sequences

% Prompt for directory
if nargin < 1
   root_path = uigetdir(pwd,'Choose directory for a particular fish');
end

% Initialize variables
a = dir(root_path);
j = 1;
calT_path = []; calS_path = []; pntT_path = []; pntS_path = []; seqs = [];

% Step through each item in directory
for i = 3:length(a)
    
    % Sequence directory
    if a(i).isdir && a(i).name(1)=='0'
        seqs(j).name = a(i).name;
        j = j + 1;
        
    % Top calibration    
    elseif strcmp(a(i).name,camT_cal_name)
        calT_path = [root_path filesep a(i).name];
        
    % Side calibration    
    elseif strcmp(a(i).name,camS_cal_name)
        calS_path = [root_path filesep a(i).name];
    
    % Top point    
    elseif strcmp(a(i).name,camT_pnt_name)
        pntT_path = [root_path filesep a(i).name];
        
    % Side point    
    elseif strcmp(a(i).name,camS_pnt_name)
        pntS_path = [root_path filesep a(i).name];
    
    end
end

% Check for all calibration files
if isempty(calT_path)
    error('No top calibration file present')
elseif isempty(calS_path)
    error('No side calibration file present')
elseif isempty(pntT_path)
    error('No top point file present')
elseif isempty(pntS_path)
    error('No side point file present')
end

% Check for sequences
if isempty(seqs)
    error('No sequence directories present')
end


%% Acquire calibration data

if isempty(dir([root_path filesep 'cal_data.mat']))
    
    f = figure;
    
    disp(' ')
    disp('This code assumes a right-handed coordinate system.')
    disp(' The cameras should be arranged as follows:')
    disp(' ')
    disp(' Top cam: y-axis points down on screen')
    disp('          x-axis points left on screen')
    disp( ' ')
    disp('Side cam: x-axis points left on screen')
    disp('          z-axis points up on screen')
    disp( ' ')
    
    % Calibrate magnificantion in top cam -------------------------------
    im_cal_camT = imread2([root_path filesep camT_cal_name]); 
    
    cal.camT.const = calibrate(im_cal_camT,...
            'Top cam: choose 2 points over known disance, press return');
    
    clear im_cal_camT
      
     % Calibrate magnificantion in cam2 -------------------------------  
    im_cal_camS = imread2([root_path filesep camS_cal_name]); 
    
    cal.camS.const = calibrate(im_cal_camS,...
              'Side cam: choose 2 points over known disance, press return');
    
    clear im_cal_camS
       
    % Select same point in the two views ------------------------------
    % Load images
    im_pnt_camT = imread2([root_path filesep camT_pnt_name]); 
    im_pnt_camS = imread2([root_path filesep camS_pnt_name]); 
    
    % Point selection
    [xT,yT] = point_select(im_pnt_camT, 'Top cam: Click on point');
    cal.camT.pnt = [xT yT];
    
    [xT,yT] = point_select(im_pnt_camS, 'Side cam: Click on point');
    cal.camS.pnt = [xT,yT];
    
    % Save 'cal'
    save([root_path filesep 'cal_data.mat'],'cal')
    
    clear im_cam1 im_cam2
    close
       
else
    disp(' ')
    disp('Loading calibration data . . .')
    disp(' ')
    
    % Load 'cal' structure
    load([root_path filesep 'cal_data.mat'])
end


%% Establish directory structure for current sequence

% Prompt for sequence
if ~isempty(seqs)
   seq_path = uigetdir(root_path,'Choose sequence to analyze');
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

    

function [x,y] = point_select(im,txt)

figure;
warning off
imshow(im);
warning on
title(txt)

set(gcf,'DoubleBuffer','on');
disp(' '); disp(' ');
disp('Left mouse button picks point.');disp(' ');
disp('z - zooms');
disp('Press return when done.')

% Loop for interactive input
but = 1; h = []; x = []; y = [];
while 1 == 1
    [xi,yi,but] = ginput(1);
    
    % Return pressed (quit out)
    if isempty(but)
        if isempty(x)
            warning('Select point before pressing return')
        else
            break
        end
        
    % Right click (add point)
    elseif but==1            
        x = xi;
        y = yi;
        
    % Zoom
    elseif but==122 
        zoom on
        pause
        
    end

    hold on
    delete(h)
    h = plot(x,y,'ro-.'); 
    hold off
end
close


    


function im = imread2(im_path)
% modifeid version of imread

im = imread(im_path); 
% Adjust contrast
im = imadjust(im);

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

