function varargout = densitometry(varargin)
% DENSITOMETRY MATLAB code for densitometry.fig
%      Richard Klein (2014)
%      densitometry@rklein.me
%
%      DENSITOMETRY, by itself, creates a new DENSITOMETRY or raises the existing
%      singleton*.
%
%      H = DENSITOMETRY returns the handle to a new DENSITOMETRY or the handle to
%      the existing singleton*.
%
%      DENSITOMETRY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DENSITOMETRY.M with the given input arguments.
%
%      DENSITOMETRY('Property','Value',...) creates a new DENSITOMETRY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before densitometry_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to densitometry_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help densitometry

% Last Modified by GUIDE v2.5 04-Mar-2015 16:41:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @densitometry_OpeningFcn, ...
                   'gui_OutputFcn',  @densitometry_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before densitometry is made visible.
function densitometry_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to densitometry (see VARARGIN)
    handles.output = hObject;
    handles.std = 1;
    handles.stdp = [];
    handles.density = figure('Visible','off');
    handles.segments = 0;
    handles.mode = 0;
    handles.mergeregions = {};
    % Update handles structure
    guidata(hObject, handles);
    
    [filename, path] = uigetfile('*.jpg');
    filename = fullfile([path,filename]);
    handles.filename = filename;

    image = imread(filename);
    image = rgb2gray(image);
    
    h = figure;
    imshow(image);title('Crop the image.');
    image1 = imcrop(h);
    if(size(image1, 1) ~= 0)
        image = image1;
        close(h);
    end

    handles.image = image;
    intensity = sum(sum(255 - handles.image));
    pixels    = size(handles.image, 1) * size(handles.image, 2);
    average = intensity/pixels;
    
    newim = 255 - handles.image;
    newim = newim - average;
    handles.fg = imadjust(newim);   
    
    guidata(hObject, handles);
    
    getSegmentation1(hObject, eventdata, handles);
    handles = guidata(hObject);
    % drawDensity(hObject, eventdata, handles);
    redrawAll(hObject, handles);
    % UIWAIT makes densitometry wait for user response (see UIRESUME)
    % uiwait(handles.figure1);

function getSegmentation1(hObject, ~, ~)
    % Calculate segments
    handles = guidata(hObject);
    handles.segments = getSegments(handles.image);
    guidata(hObject, handles);


function segmentation = getSegments(image)
    % Calculate segmentations
    I = image;
    [~, threshold] = edge(I, 'sobel');
    fudgeFactor = .5;
    BWs = edge(I,'sobel', threshold * fudgeFactor);
    %figure, imshow(BWs), title('binary gradient mask');

    %se90 = strel('line', 3, 90);
    %se0 = strel('line', 3, 0);

    BWsdil = imclose(BWs, strel('disk',2));
    %figure, imshow(BWsdil), title('dilated gradient mask');

    BWdfill= imfill(BWsdil,'holes');
    %figure, imshow(BWdfill);
    %title('binary image with filled holes');

    seD = strel('disk',2);
    segmentation = imopen(BWdfill,seD);   

function segout=setOutlines(hObject, handles)
    % Calculate outline image
    BWoutline = bwperim(handles.segments);
    handles.outline = handles.image;
    handles.outline(BWoutline) = 255;
    guidata(hObject, handles);
    
    axes(handles.axes2), imshow(handles.outline);
    segout = handles.outline;

function calcDensityDarkness(hObject, ~, handles)
    


    
function drawDensity(hObject, eventdata, handles)
% hObject    handle to mergeregions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %handles =  guidata(hObject);

    cc = bwconncomp(handles.segments, 4);
    %lbl = bwlabel(handles.segments, 4);
    [Lrgb, lbl] = getLabels(handles);
    

    data = regionprops(cc,'Area');

    STD = handles.std;
    if(size(handles.stdp,2) ~= 2)
        STD_value = 1;
    else
        STD_value = lbl(round(handles.stdp(2)), round(handles.stdp(1)));
    end
    density = zeros(1, cc.NumObjects);
    grey = zeros(1, cc.NumObjects);

    if(cc.NumObjects > 0)
        curr = 1;
        for i = 1:cc.NumObjects
            %segment = false(size(handles.segments));
            %segment(cc.PixelIdxList{i}) = true;
            local = sum(sum(lbl==i));

            if(local < 5)
                if(i == STD)
                    STD = STD + 1;
                end
                continue;
            else
                density(curr) = local; %data(i).Area; %/data(STD).Area*100;
                if(handles.mode == 1)
                    grey(curr) = sum(sum(handles.image(lbl == i)));
                end
                curr = curr + 1;
            end
            if(STD_value == i)
                STD = curr-1;
            end
        end
        density = density(1:curr-1);
        grey = grey(1:curr-1);
        
        
        
        if(handles.mode == 0)
            STD_density = density(STD);
            for i = 1:curr-1
                density(i) = density(i)/STD_density * 100;
            end
        else
            STD_density = grey(STD);
            %STD_density = grey(STD)/density(STD);
            for i = 1:curr-1
                density(i) = (grey(i))/STD_density * 100;
                %density(i) = (grey(i)/density(i))/STD_density * 100;
            end
        end
        density
        figure(handles.density); bar(density);
        labels = cell(size(density,2));
        for i = 1:curr-1
            labels{i} = i; 
        end
        labels{STD} = 'STD';
        set(gca,'XTickLabel',labels);
        grid;
    end

% --- Outputs from this function are returned to the command line.
function varargout = densitometry_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
    varargout{1} = handles.output;


% --- Executes on button press in button_unmask_segments.
function button_unmask_segments_Callback(hObject, eventdata, handles)
% hObject    handle to button_unmask_segments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    maxRowSeg = size(handles.segments,1);
    maxColSeg = size(handles.segments,2);
    
    rect = getrect(handles.axes1);
    rowS = round(max(rect(2),1));
    rowE = round(min(rect(2)+rect(4), maxRowSeg));
    colS = round(max(rect(1),1));
    colE = round(min(rect(1)+rect(3), maxColSeg));
    
    handles.segments(rowS:rowE, colS:colE) = 0;
    Lrgb = label2rgb(bwlabel(handles.segments), 'jet', 'k');
     Lrgb = getLabels(handles);
    axes(handles.axes1), imshow(Lrgb);
    guidata(hObject, handles);
    redrawAll(hObject, handles);


% --- Executes on button press in button_select_std.
function button_select_std_Callback(hObject, eventdata, handles)
% hObject    handle to button_select_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    cc = bwconncomp(handles.segments, 4);

    data2 = regionprops(cc,'PixelIdxList','PixelList');

    axes(handles.axes1);
    STD = -1;
    while STD == -1
        [x,y] = ginput(1);
        handles.stdp = [x,y];
    
        for i = 1:cc.NumObjects
            if(any(abs(data2(i).PixelList(:,1) - x) < 1) && any(abs(data2(i).PixelList(:,2) - y) < 1))
                STD = i;
            end
        end
    end

    handles.std = STD;
    guidata(hObject, handles);
    redrawAll(hObject, handles);

% --- Executes on button press in button_unmask_outline.
function button_unmask_outline_Callback(hObject, eventdata, handles)
% hObject    handle to button_unmask_outline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    maxRowSeg = size(handles.segments,1);
    maxColSeg = size(handles.segments,2);
    
    rect = getrect(handles.axes2);
    rowS = round(max(rect(2),1));
    rowE = round(min(rect(2)+rect(4), maxRowSeg));
    colS = round(max(rect(1),1));
    colE = round(min(rect(1)+rect(3), maxColSeg));

    handles.segments(rowS:rowE, colS:colE) = 0;
    guidata(hObject, handles);    
    redrawAll(hObject, handles);

% --- Executes on button press in button_mask_outline.
function button_mask_outline_Callback(hObject, eventdata, handles)
% hObject    handle to button_mask_outline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    maxRowSeg = size(handles.outline,1);
    maxColSeg = size(handles.outline,2);
    
    rect = getrect(handles.axes2);
    rowS = round(max(rect(2),1));
    rowE = round(min(rect(2)+rect(4), maxRowSeg));
    colS = round(max(rect(1),1));
    colE = round(min(rect(1)+rect(3), maxColSeg));
    
    handles.segments(rowS:rowE, colS:colE) = 1;
    guidata(hObject, handles);    
    redrawAll(hObject, handles);

% --- Executes on button press in button_autosegment.
function button_autosegment_Callback(hObject, eventdata, handles)
% hObject    handle to button_autosegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    maxRowSeg = size(handles.segments,1);
    maxColSeg = size(handles.segments,2);
    
    rect = getrect(handles.axes2);
    rowS = round(max(rect(2),1));
    rowE = round(min(rect(2)+rect(4), maxRowSeg));
    colS = round(max(rect(1),1));
    colE = round(min(rect(1)+rect(3), maxColSeg));
    
    handles.segments(rowS:rowE, colS:colE) = getSegments(handles.image(rowS:rowE, colS:colE));
    guidata(hObject, handles);    
    redrawAll(hObject, handles);


% --- Executes on button press in button_save_images.
function button_save_images_Callback(~, ~, handles)
% hObject    handle to button_save_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [Path, Name, ext] = fileparts(handles.filename);
    outline_name = fullfile([Path,'/', Name, '_outline.jpg']);
    imwrite(handles.outline,outline_name);
    segments_name = fullfile([Path, '/', Name, '_segments.jpg']);
    imwrite(handles.segments,segments_name);
    cropped_name = fullfile([Path, '/', Name, '_cropped.jpg']);
    imwrite(handles.image,cropped_name);
    density_name = fullfile([Path, '/', Name, '_density.jpg']);
    saveas(handles.density, density_name);


% --- Executes on button press in button_densities.
function button_densities_Callback(hObject, eventdata, handles)
% hObject    handle to button_densities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    drawDensity(hObject, eventdata, handles);


% --- Executes on button press in mergeregions.
function mergeregions_Callback(hObject, eventdata, handles)
% hObject    handle to mergeregions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    maxRowSeg = size(handles.segments,1);
    maxColSeg = size(handles.segments,2);
    
    rect = getrect(handles.axes1);
    rowS = round(max(rect(2),1));
    rowE = round(min(rect(2)+rect(4), maxRowSeg));
    colS = round(max(rect(1),1));
    colE = round(min(rect(1)+rect(3), maxColSeg));
    lbl = bwlabel(handles.segments);
    
    handles.mergeregions{end+1} = [rowS, rowE, colS, colE];
    guidata(hObject, handles);
    redrawAll(hObject,handles);
    
function [Lrgb, lbl] = getLabels(handles)
    lbl = bwlabel(handles.segments);

    regions = size(handles.mergeregions, 2);
    for i = 1:regions
        rect = handles.mergeregions{i};
        rowS = rect(1);
        rowE = rect(2);
        colS = rect(3);
        colE = rect(4);
        
        lblS = lbl(rowS:rowE, colS:colE);
        newLbl = max(max(lblS));
        lblS(lblS > 0) = newLbl;
        lbl(rowS:rowE, colS:colE) = lblS;
    end

    Lrgb = label2rgb(lbl, 'jet', 'k');
    
    



% --- Executes on button press in clearmerge.
function clearmerge_Callback(hObject, eventdata, handles)
% hObject    handle to clearmerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.mergeregions = {};
    guidata(hObject, handles);
    redrawAll(hObject, handles);


function drawMergedRegions(handles)
     n = size(handles.mergeregions, 2);
     for i = 1:n
        rect = handles.mergeregions{i};
        rowS = rect(1);
        rowE = rect(2);
        colS = rect(3);
        colE = rect(4);
        
        rMin = min(rowS, rowE);
        rMax = max(rowS, rowE);
        cMin = min(colS, colE);
        cMax = max(colS, colE);
        
        axes(handles.axes1),  rectangle('Position', [cMin, rMin, cMax-cMin, rMax-rMin], 'LineWidth',2, 'EdgeColor', 'r');
     end                
 
 
function redrawAll(hObject, handles)
    Lrgb = getLabels(handles);
    axes(handles.axes1), imshow(Lrgb);
    drawMergedRegions(handles);
    guidata(hObject, handles);
    setOutlines(hObject, handles);
    drawDensity(hObject, 0, handles);


% --- Executes on button press in unmaskp.
function unmaskp_Callback(hObject, eventdata, handles)
% hObject    handle to unmaskp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    while(1==1)
        [col, row, button] = ginput(1);
        if (button ~= 1)
            break;
        end
        col = round(col);
        row = round(row);
                    
        handles.segments(row-2:row+2, col-2:col+2) = 0;
        Lrgb = getLabels(handles);
        axes(handles.axes1), imshow(Lrgb);
        %drawMergedRegions(handles);
        %guidata(hObject, handles);
        setOutlines(hObject, handles);
        
        
    end

    Lrgb = getLabels(handles);
    axes(handles.axes1), imshow(Lrgb);
    guidata(hObject, handles);
    redrawAll(hObject, handles);
        
    


% --- Executes on button press in maskpaint.
function maskpaint_Callback(hObject, eventdata, handles)
% hObject    handle to maskpaint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    while(1==1)
        [col, row, button] = ginput(1);
        if (button ~= 1)
            break;
        end
        col = round(col);
        row = round(row);
                    
        handles.segments(row-2:row+2, col-2:col+2) = 1;
        Lrgb = getLabels(handles);
        axes(handles.axes1), imshow(Lrgb);
        %drawMergedRegions(handles);
        %guidata(hObject, handles);
        setOutlines(hObject, handles);
        
        
    end

    Lrgb = getLabels(handles);
    axes(handles.axes1), imshow(Lrgb);
    guidata(hObject, handles);
    redrawAll(hObject, handles);
        


% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)
% hObject    handle to back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    tmp = handles.image;
    handles.image = handles.fg;
    handles.fg = tmp;
    handles.mode = 1 - handles.mode;
    guidata(hObject, handles);
    redrawAll(hObject, handles);
    
    
