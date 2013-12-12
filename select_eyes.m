
function select_eyes(seq_path,d,t,seq,cal,imMeanT,imMeanS)
% Function called by acq_3d to select eyes of larvae in 3d




%% Parameters

% ??
skip_side = 0;

% Invert the image
invert = 1;

% Acquisiton mode
mode = 0;

sub_back = 1;
    

%% Display commands
disp('Commands ')
disp('   Left click   - collect coordinate')
disp('   Right click  - delete coordinate')
disp(' ')
disp('   m    - toggle mode (left eye/right eye/swim bladder/dorsal edge/ventral edge)')
disp('   i    - invert image')
disp('   s    - subtract background (default)')
disp('   a    - autotracking')
disp('   t    - adjust threshold (for autot tracking)')
disp('   j    - jump to specific frame')
disp('   1-9  - jump to position in video')
disp('   z    - adjust zoom')
disp('   c    - clear points on frame')
disp(' ')
disp('   Return  - finish point collecting')
disp('   ESC     - stop acq_3d')
disp(' ')



%% Loop through frames from the two perspectives 

% Define frame index number to begin
idx = d.curr_frame_idx;

% Make figure window
f = figure;

% Image handle
h_im = [];

while true
    
    
    %% Initialize
    
    % Define image paths
    pathT = [seq_path filesep seq.dir_camT filesep seq.aT(idx).name];
    pathS = [seq_path filesep seq.dir_camS filesep seq.aS(idx).name];
    
    % Read image files for current time
    cFrame = d.frames(idx);

    % load image with (or without) background subtracton
    if sub_back
        imT = give_im(pathT,invert,imMeanT);
        imS = give_im(pathS,invert,imMeanS);
    else
        imT = give_im(pathT,invert);
        imS = give_im(pathS,invert);
    end
    
    % Store current frame (index)
    d.curr_frame_idx = idx;
    
    % Define relative frames
    next_idx = min([idx+d.frame_skip d.num_frames]);
    last_idx = max([idx-d.frame_skip 1]);
    
    % Define handles
    hS = [];
    hT = [];
    
    
    %% Display current frame
    
    % Display side view
    figure(f)
    warning off
    hS_sm = subplot(3,2,1);
    imshow(imS);
    hold on
    plot(d.S.x_roi(idx,:),d.S.y_roi(idx,:),'y-')
    xlabel('y ->')
    ylabel('z ->')
    
    title(['Side (' num2str(cFrame) '/' num2str(d.end_frame) ')'])
    
    % Zoomed-in version
    hS_lg = subplot(3,2,[3 5]);
    imshow(imS);
    hold on
    axis([min(d.S.x_roi(idx,:)) max(d.S.x_roi(idx,:)) ...
        min(d.S.y_roi(idx,:)) max(d.S.y_roi(idx,:))]);
    
    % Plot data
    hold on
    hL(1) = plot(d.S.leftEye(idx,1),d.S.leftEye(idx,2),'g+');
    hR(1) = plot(d.S.rightEye(idx,1),d.S.rightEye(idx,2),'r+');
    hS(1) = plot(d.S.SB(idx,1),d.S.SB(idx,2),'m+');
    if ~isempty(t.S.dorsalEdge{idx,1}) && ~isempty(t.S.dorsalEdge{idx,2})
        hD(1) = plot(t.S.dorsalEdge{idx,1},t.S.dorsalEdge{idx,2},'bs-','MarkerFaceColor','b','MarkerSize',3);
    end
    if ~isempty(t.S.ventralEdge{idx,1}) && ~isempty(t.S.ventralEdge{idx,2})
        hD(1) = plot(t.S.ventralEdge{idx,1},t.S.ventralEdge{idx,2},'ys-','MarkerFaceColor','y','MarkerSize',3);
    end
    %Plot previous dorsal and ventral edge from previous frame
    if idx > 1
        if ~isempty(t.S.dorsalEdge{idx-1,1}) && ~isempty(t.S.dorsalEdge{idx-1,2})
            hDLprev(1) = plot(t.S.dorsalEdge{idx-1,1},t.S.dorsalEdge{idx-1,2},'Color',[0 0 0.25]);
            hVprev(1) = plot(t.S.ventralEdge{idx-1,1},t.S.ventralEdge{idx-1,2},'Color',[0.25 0.25 0]);
        end
        if ~isempty(t.S.ventralEdge{idx-1,1}) && ~isempty(t.S.ventralEdge{idx-1,2})
            hVprev(1) = plot(t.S.ventralEdge{idx-1,1},t.S.ventralEdge{idx-1,2},'Color',[0.25 0.25 0]);
        end
    end
    hold off
    
    % Title
    if mode==0
        title('Select left eye')
    elseif mode==1
        title('Select right eye')
    elseif mode==2
        title('Select swim bladder')
    elseif mode==3
        title('Select dorsal edge')
    elseif mode==4
        title('Select ventral edge')
    end
    
    % Display top view
    warning off
    hT_sm = subplot(3,2,2);
    imshow(imT);
    hold on
    plot(d.T.x_roi(idx,:),d.T.y_roi(idx,:),'y-')
    xlabel('x ->')
    ylabel('y ->')
    
    title(['Top (' num2str(cFrame) '/' num2str(d.end_frame) ')'])
    
    % Zoomed-in version
    hT_lg = subplot(3,2,[4 6]);
    imshow(imT);
    axis([min(d.T.x_roi(idx,:)) max(d.T.x_roi(idx,:)) ...
        min(d.T.y_roi(idx,:)) max(d.T.y_roi(idx,:))]);
    
    % Plot data
    hold on
    hL(2) = plot(d.T.leftEye(idx,1),d.T.leftEye(idx,2),'g+');
    hR(2) = plot(d.T.rightEye(idx,1),d.T.rightEye(idx,2),'r+');
    hS(2) = plot(d.T.SB(idx,1),d.T.SB(idx,2),'m+');
    if ~isempty(t.T.dorsalEdge{idx,1}) && ~isempty(t.T.dorsalEdge{idx,2})
        hD(2) = plot(t.T.dorsalEdge{idx,1},t.T.dorsalEdge{idx,2},'bs-','MarkerFaceColor','b','MarkerSize',3);
        
    end
    if ~isempty(t.T.ventralEdge{idx,1}) && ~isempty(t.T.ventralEdge{idx,2})
        hV(2) = plot(t.T.ventralEdge{idx,1},t.T.ventralEdge{idx,2},'ys-','MarkerFaceColor','y','MarkerSize',3);
        
    end
    %Plot previous dorsal and ventral edge from previous frame
    if idx > 1
        if~isempty(t.T.dorsalEdge{idx-1,1}) && ~isempty(t.T.dorsalEdge{idx-1,2})
            hDLprev(2) = plot(t.T.dorsalEdge{idx-1,1},t.T.dorsalEdge{idx-1,2},'Color',[0 0 0.25]);
        end
        if~isempty(t.T.ventralEdge{idx-1,1}) && ~isempty(t.T.ventralEdge{idx-1,2})
            hVprev(2) = plot(t.T.ventralEdge{idx-1,1},t.T.ventralEdge{idx-1,2},'Color',[0.25 0.25 0]);
        end
    end
    hold off
    
    % Title
    if mode==0
        title('Select left eye')
    elseif mode==1
        title('Select right eye')
    elseif mode==2
        title('Select swim bladder')
    elseif mode==3
        title('Select dorsal edge')
    elseif mode==4
        title('Select ventral edge')
    end
    
    % Turn warnings back on
    warning on
    
    % Prompt to input a point
    figure(f)
    [x_tmp,y_tmp,b_tmp] = ginput(1);
    
    
    %% Left click - Store value or adjust roi
    if b_tmp ==1
        
        % Zoomed top view
        if gca==hT_lg  
            
            if mode==0
                d.T.leftEye(idx,:) = [x_tmp y_tmp];
            elseif mode==1
                d.T.rightEye(idx,:) = [x_tmp y_tmp];
            elseif mode==2
                d.T.SB(idx,:) = [x_tmp y_tmp];
            elseif mode==3
                t.T.dorsalEdge{idx,1} = [t.T.dorsalEdge{idx,1} x_tmp];
                t.T.dorsalEdge{idx,2} = [t.T.dorsalEdge{idx,2} y_tmp];
            else
                t.T.ventralEdge{idx,1} = [t.T.ventralEdge{idx,1} x_tmp];
                t.T.ventralEdge{idx,2} = [t.T.ventralEdge{idx,2} y_tmp];
            end
            
        % Zoomed side view
        elseif gca==hS_lg
            
            if mode==0
                d.S.leftEye(idx,:) = [x_tmp y_tmp];
            elseif mode==1
                d.S.rightEye(idx,:) = [x_tmp y_tmp];
            elseif mode==2
                d.S.SB(idx,:) = [x_tmp y_tmp];
            elseif mode==3
                t.S.dorsalEdge{idx,1} = [t.S.dorsalEdge{idx,1} x_tmp];
                t.S.dorsalEdge{idx,2} = [t.S.dorsalEdge{idx,2} y_tmp];
            else
                t.S.ventralEdge{idx,1} = [t.S.ventralEdge{idx,1} x_tmp];
                t.S.ventralEdge{idx,2} = [t.S.ventralEdge{idx,2} y_tmp];
            end
            
        % Unzoomed top view
        elseif gca==hT_sm
            
            new_x = d.T.x_roi(idx,:)-mean(d.T.x_roi(idx,:))+x_tmp;
            new_y = d.T.y_roi(idx,:)-mean(d.T.y_roi(idx,:))+y_tmp;
            
            mod_remain = [[1:length(d.frames)] > idx]';  
            mod_next = min([find(mod_remain & d.T.mod_roi,1,'first') ...
                            length(d.frames)]);
            
            d.T.x_roi(idx:mod_next,:) = repmat(new_x,length(idx:mod_next),1);
            d.T.y_roi(idx:mod_next,:) = repmat(new_y,length(idx:mod_next),1);
            
            d.T.mod_roi(idx) = 1;
            
            clear new_x new_y mod_remain mod_next
            
        % Unzoomed side view
        elseif gca==hS_sm  
            
            new_x = d.S.x_roi(idx,:)-mean(d.S.x_roi(idx,:))+x_tmp;
            new_y = d.S.y_roi(idx,:)-mean(d.S.y_roi(idx,:))+y_tmp;
            
            mod_remain = [[1:length(d.frames)] > idx]';  
            mod_next = min([find(mod_remain & d.S.mod_roi,1,'first') ...
                            length(d.frames)]);
            
            d.S.x_roi(idx:mod_next,:) = repmat(new_x,length(idx:mod_next),1);
            d.S.y_roi(idx:mod_next,:) = repmat(new_y,length(idx:mod_next),1);
            
            d.S.mod_roi(idx) = 1;
            
            clear new_x new_y mod_remain mod_next
            
        end
        
        h_im = gca;
        
        clear y_ctnr y_min y_max x_min x_max
        
    %% Right click
    
    elseif b_tmp ==3
        
        % Overwrite data with nans
        if gca==hT_lg     
            if mode==0
                d.T.leftEye(idx,:) = nan(1,2);
                
            elseif mode==1
                d.T.rightEye(idx,:) = nan(1,2);
                
            elseif mode==2
                d.T.SB(idx,:) = nan(1,2);
                
            elseif mode==3
                t.T.dorsalEdge{idx,1} = t.T.dorsalEdge{idx,1}(1:end-1);
                t.T.dorsalEdge{idx,2} = t.T.dorsalEdge{idx,2}(1:end-1);
                
            else
                t.T.ventralEdge{idx,1} = t.T.ventralEdge{idx,1}(1:end-1);
                t.T.ventralEdge{idx,2} = t.T.ventralEdge{idx,2}(1:end-1);
                
            end
               
        elseif gca==hS_lg         
            if mode==0
                d.S.leftEye(idx,:) = nan(1,2);
                
            elseif mode==1
                d.S.rightEye(idx,:) = nan(1,2);
                
            elseif mode==2
                d.S.SB(idx,:) = nan(1,2);
                
            elseif mode==3
                t.S.dorsalEdge{idx,1} = t.S.dorsalEdge{idx,1}(1:end-1);
                t.S.dorsalEdge{idx,2} = t.S.dorsalEdge{idx,2}(1:end-1);
                
            else
                t.S.ventralEdge{idx,1} = t.S.ventralEdge{idx,1}(1:end-1);
                t.S.ventralEdge{idx,2} = t.S.ventralEdge{idx,2}(1:end-1);
                
            end
            
        else
            warning('You can only erase coordinates from the zoomed view')
        end
            
        %skip_side = 1;
        
        % Jump back current frame number
        %idx = last_idx;
        
        
   %% Autotrack ('a')
    elseif b_tmp == 97    
        
         % Top view
        if isempty(h_im)
            warning('you first need to click on a window')
            
        elseif mode > 1
            warning('Change mode -- you can only track the eyes');
            
        % Track from top view
        elseif (h_im==hT_lg) || (h_im==hT_sm)  
            
            if isnan(d.T.leftEye(idx,1)) || isnan(d.T.rightEye(idx,1))
                warning(['You first need to select points for both eyes'])
                
            elseif mode==0
                [d.T,idx] = autotrack_eye(seq_path,seq,imMeanT,d.T,idx,d.frames,'top','left');
                d.curr_frame_idx = idx;
                
            elseif mode==1
                [d.T,idx] = autotrack_eye(seq_path,seq,imMeanT,d.T,idx,d.frames,'top','right');
                
            end
           
            d.curr_frame_idx = idx;
            
        % Side view
        elseif (h_im==hS_lg) || (h_im==hS_sm)
            
            if isnan(d.S.leftEye(idx,1)) || isnan(d.S.rightEye(idx,1))
                warning(['You first need to select points for both eyes'])
            
            elseif mode==0
                [d.S,idx] = autotrack_eye(seq_path,seq,imMeanS,d.S,idx,d.frames,'side','left');
                d.curr_frame_idx = idx;
                
            elseif mode==1
                [d.S,idx] = autotrack_eye(seq_path,seq,imMeanS,d.S,idx,d.frames,'side','right');
                
            end
            
            d.curr_frame_idx = idx;
            
        end
        
        
       
        
        
    %% Adjust threshold ('t')
    elseif b_tmp == 116
        
        if ~invert
            warning('You first need to invert the image')
        end
           
        % Top view
        if isempty(h_im)
            warning('you first need to click on a window')
        
        elseif (h_im==hT_lg) || (h_im==hT_sm)
            
            x_tmp = floor(min(d.T.x_roi(idx,:))):ceil(max(d.T.x_roi(idx,:)));
            y_tmp = floor(min(d.T.y_roi(idx,:))):ceil(max(d.T.y_roi(idx,:)));
            
            im_tmp = imcomplement(imT(y_tmp,x_tmp));
            
            % Matlab guesses a threshold value
            waitfor(threshFinder(im_tmp,d,'top',seq_path))
            
            clear x_tmp y_tmp im_tmp
            
        % Side view
        elseif (h_im==hS_lg) || (h_im==hS_sm)
            
            x_tmp = floor(min(d.S.x_roi(idx,:))):ceil(max(d.S.x_roi(idx,:)));
            y_tmp = floor(min(d.S.y_roi(idx,:))):ceil(max(d.S.y_roi(idx,:)));
            
            im_tmp = imcomplement(imS(y_tmp,x_tmp));
            
            % Matlab guesses a threshold value
            waitfor(threshFinder(im_tmp,d,'side',seq_path))
            
            clear x_tmp y_tmp im_tmp
        end
        
        % Load updated 'd' structure
        load([seq_path filesep 'coord_data.mat'])
        
        
    %% Other commands
    
    % Left arrow
    elseif b_tmp == 28
        
        % Back up frame
        idx = last_idx;
        
    % Right arrow
    elseif b_tmp == 29
        
        % Advance frame
        idx = next_idx;
        
    % Jump to frame 
    elseif b_tmp==106
        
        answer = inputdlg('Frame number','Jump to frame',1, ...
            {num2str(round(seq.num_frames/2))});
        
        % Define frame number
        idx = find(str2num(answer{1})==d.frames);
        
        clear answer  
    
    % Change mode ("m")
    elseif b_tmp == 109
        
        if mode==4
            mode = 0;
        else 
            mode = mode+1;
        end
        
    % Invert image ("i")
    elseif b_tmp == 105
        
        invert = abs(invert - 1);
        
    % Subtract background ("s")
    elseif b_tmp == 115
        
        sub_back = abs(sub_back - 1);
        
    % Clear points
    elseif b_tmp == 99
        
         if isempty(h_im)
            warning('you first need to click on a window')
        
        elseif (h_im==hT_lg) || (h_im==hT_sm)
            
            d.T.leftEye(idx,:) = [nan nan];
            d.T.rightEye(idx,:) = [nan nan];
            
        elseif (h_im==hS_lg) || (h_im==hS_sm)
            
            d.S.leftEye(idx,:) = [nan nan];
            d.S.rightEye(idx,:) = [nan nan]; 
            
         end
    
    % Jump to a different location in the video with numbers
    elseif and((b_tmp>=49),(b_tmp<=57))
        
        % Define frame number
        idx = max([1 round(length(d.frames)*(b_tmp-49)/8)]);
    
    % Zoom (select roi)
    elseif b_tmp==122
        
        [new_x,new_y,h_im] = select_roi;

        mod_remain = [[1:length(d.frames)] > idx]';
       
        if h_im==hT_sm     
            
            mod_next = min([find(mod_remain & d.T.mod_roi,1,'first') ...
                length(d.frames)]);
            
            d.T.x_roi(idx:mod_next,:) = repmat(new_x,length(idx:mod_next),1);
            d.T.y_roi(idx:mod_next,:) = repmat(new_y,length(idx:mod_next),1);
            
            d.T.mod_roi(idx) = 1;
            
            clear mod_remain mod_next
            
        elseif h_im==hS_sm
            
            mod_next = min([find(mod_remain & d.S.mod_roi,1,'first') ...
                length(d.frames)]);
            
            d.S.x_roi(idx:mod_next,:) = repmat(new_x,length(idx:mod_next),1);
            d.S.y_roi(idx:mod_next,:) = repmat(new_y,length(idx:mod_next),1);
            
            d.S.mod_roi(idx) = 1;
            
            clear mod_remain mod_next
            
        else
            warning('You can only adjust zoom in full-frame images')
        end
        
        clear x_roi y_roi h_im num_remain
        

       
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
    %delete(bd(1))
    %set(htxt(1),'Color','k')
    
    % Save data
    save([seq_path filesep 'coord_data.mat'],'d')
    save([seq_path filesep 'tail_data.mat'],'t')
    
    %pause(.05)
end    

close



function [x_roi,y_roi,h_im] = select_roi
hold on

% Prompt for first point
disp(' ');
disp('Select first point')
[x1,y1,but1] = ginput(1);
h1 = plot(x1,y1,'y+');

% Prompt for second
disp(' ');
disp('Select second point')
[x2,y2,but2] = ginput(1);

% Define roi
delete(h1)
x_roi = [x1 x2 x2 x1 x1];
y_roi = [y1 y1 y2 y2 y1];

% Display briefly
h2 = plot(x_roi,y_roi,'y-');
pause(.3)
delete(h2)

% Specify current figure handle
h_im = gca;

hold off

