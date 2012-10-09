function [xT,yT,xS,yS] = select_common(seq_path,seq,imMeanT,imMeanS)
% Used by acq_3d to have a user select a common point from two views


% Invert image?
invert = 1;

% Allows the selection of a single point in a movie
figure;

% Give instructions
set(gcf,'DoubleBuffer','on');
disp(' '); disp(' ');
disp('Left mouse button picks point.');disp(' ');
disp('z - zoom');
disp('u - unzoom');
disp('j - jump to frame')
disp('Press return when done.')

% Initialize variables
but = 1; 
subplot(1,2,1)
hT = plot([],[]); 
hold on
hS = plot([],[]); 
xT = []; yT = []; xS = []; yS = [];
x_roiT = []; y_roiT = [];
x_roiS = []; y_roiS = [];

% Starting image number
imNum = 1;

while 1 == 1
  
    % Define image paths
    imPath_T = [seq_path filesep seq.dir_camT filesep seq.aT(imNum).name];
    imPath_S = [seq_path filesep seq.dir_camS filesep seq.aS(imNum).name];
    
    % Read images
    imT = give_im(imPath_T,invert,imMeanT,x_roiT,y_roiT);
    imS = give_im(imPath_S,invert,imMeanS,x_roiS,y_roiS);
     
    % Turn off warning during imshow
    warning off
    
    % Plot the top view
    h_T = subplot(1,2,1);
      imshow(imT);
      cMap = colormap('gray');
      hold on
      title(['Frame ' num2str(imNum) '/' num2str(seq.num_frames) ' Top'])
      
      % Set zoom
      if ~isempty(x_roiT)
          axis([min(x_roiT) max(x_roiT) min(y_roiT) max(y_roiT)])
      end
      
      % Plot point
      delete(hT)
      hT = plot(xT,yT,'r+');
      xlabel('X ->')
      ylabel('Y ->')
      
    % Plot the side view
    h_S = subplot(1,2,2);
      imshow(imS);
      cMap = colormap('gray');
      hold on
      title(['Frame ' num2str(imNum) '/' num2str(seq.num_frames) ' Side'])
      
      % Set zoom
      if ~isempty(x_roiS)
        axis([min(x_roiS) max(x_roiS) min(y_roiS) max(y_roiS)])
      end
      
      % Plot point
      delete(hS)
      hS = plot(xS,yS,'r+');
      xlabel('Y ->')
      ylabel('Z ->')
      
    % Turn warning back on
    warning on
    
    % G-INPUT!
    [xi,yi,but] = ginput(1);
    
    % Return pressed (quit out)
    if isempty(but)
        if isempty(xS) || isempty(xT)
            warning('Select points before pressing return')
        else
            break
        end
        
    % Left click (add point)
    elseif but==1            
        x = xi;
        y = yi;
        
        % Figure out which window
        if gca==h_T
            xT = xi;
            yT = yi;
            
        elseif gca==h_S
            xS = xi;
            yS = yi;
        end
    
    % Jump to frame 
    elseif but==106
        answer = inputdlg('Frame number','Jump to frame',1, ...
            {num2str(round(seq.num_frames/2))});
        
        % Define image number
        imNum = min([str2num(answer{1}) seq.num_frames]);  
        imNum = max([1 imNum]);
        
        clear answer
        
    % Zoom (select roi)
    elseif but==122 
        
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
        
        % Figure out which window
        if gca==h_T
            x_roiT = x_roi;
            y_roiT = y_roi;
            
        elseif gca==h_S
            x_roiS = x_roi;
            y_roiS = y_roi;
        else
            error('not sure of image handle')
        end
        
        clear x_roi y_roi h1 h2 but1 but2
        
        % Left click (add point)
        if but==1
            x = xi;
            y = yi;
            
            %pause
        end
        
    % unzoom ("u")
    elseif but == 117
        if gca==h_T
            x_roiT = [];
            y_roiT = [];
            
        else
            x_roiS = [];
            y_roiS = [];
        end
        
    % Left arrow
    elseif but == 28
        
        % Back up frame
        imNum = max([1 imNum-1]);
        
    % Right arrow
    elseif but == 29
        
        % Advance frame
        imNum = min([imNum+1 seq.num_frames]);
        
    end

end
close


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
    im = imcomplement(im);
end
