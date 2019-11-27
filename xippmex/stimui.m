function varargout = stimui(varargin)
% STIMUI MATLAB code for stimui.fig
%      STIMUI, by itself, creates a new STIMUI or raises the existing
%      singleton*.
%
%      H = STIMUI returns the handle to a new STIMUI or the handle to
%      the existing singleton*.
%
%      STIMUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STIMUI.M with the given input arguments.
%
%      STIMUI('Property','Value',...) creates a new STIMUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stimui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stimui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stimui

% Last Modified by GUIDE v2.5 03-Jan-2014 15:38:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stimui_OpeningFcn, ...
                   'gui_OutputFcn',  @stimui_OutputFcn, ...
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

% --- Executes just before stimui is made visible.
function stimui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stimui (see VARARGIN)

% Choose default command line output for stimui
handles.output = hObject;

handles.amp_scale = [];
handles.amp_scale = 3.75;
if ~isempty(varargin)
    handles.amp_scale = varargin{1};
end

% set the maximum amplitude to reflect the current headstage scaling.
set(handles.amp_slider, 'max', 127 * handles.amp_scale);

% DEFAULT_AMP = 10.0;
DEFAULT_AMP = 10 * handles.amp_scale;
DEFAULT_FREQ = 30.0;
DEFAULT_DUR = 0.2;
DEFAULT_TL = 1000.0;
DEFAULT_TD = 0.0;
DEFAULT_FS = 1.0;

set(handles.amp_slider, 'Value', DEFAULT_AMP);
set(handles.amp_edit, 'String', sprintf('%4.2f', DEFAULT_AMP));
set(handles.freq_slider, 'Value', DEFAULT_FREQ);
set(handles.freq_edit, 'String', sprintf('%4.2f', DEFAULT_FREQ));
set(handles.dur_slider, 'Value', DEFAULT_DUR);
set(handles.dur_edit, 'String', sprintf('%4.2f', DEFAULT_DUR));
set(handles.tl_slider, 'Value', DEFAULT_TL);
set(handles.tl_edit, 'String', sprintf('%4.2f', DEFAULT_TL));
set(handles.td_slider, 'Value', DEFAULT_TD);
set(handles.td_edit, 'String', sprintf('%4.2f', DEFAULT_TD));
set(handles.fs_slider, 'Value', DEFAULT_FS);
set(handles.fs_edit, 'String', sprintf('%4.2f', DEFAULT_FS));

set(handles.enable_checkbox, 'Value', 1);
set(handles.pol_checkbox, 'Value', 1);

handles.stim_rounding = [];
handles.stim_rounding.dur = 0.005;
handles.stim_rounding.amp = handles.amp_scale;
handles.stim_rounding.tl = 0.03333333;
handles.stim_rounding.td = 0.03333333;
handles.stim_rounding.freq = 1;

handles.stim_params = [];

handles.fs_elec = [];

if xippmex()
    % find stim channels and population list
    stim_elec_list = xippmex('elec', 'stim');
    labels = {};
    for i=1:length(stim_elec_list)
        labels = [labels,  sprintf('%03d', stim_elec_list(i))];
    end
    set(handles.found_list, 'String', cellstr(labels));
else
    % quit...
end
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = stimui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on slider movement.
function amp_slider_Callback(hObject, eventdata, handles)
% hObject    handle to amp_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.amp_slider, 'Value');
% ensuring we have integer divisions of stim applitude to match
% what is seen on the stim chip
% set(handles.amp_slider, 'Value', round(value));
% switching this to integer divisions.  See above
% put the slider value in increments of the uA divisions of the stim chip
set(handles.amp_slider, 'Value', ...
    round(value / handles.stim_rounding.amp)*handles.stim_rounding.amp);

set(handles.amp_edit, 'String', sprintf('%4.2f', get(handles.amp_slider, 'Value')));

% --- Executes during object creation, after setting all properties.
function amp_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amp_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on selection change in found_list.
function found_list_Callback(hObject, eventdata, handles)
% hObject    handle to found_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% val = get(handles.found_list, 'Value');
% found_labels = get(handles.found_list, 'String');
% elec = str2double(found_labels(val));
% 
% DEFAULT_AMP = 560.0;
% DEFAULT_FREQ = 200.0;
% DEFAULT_DUR = 0.20;
% DEFAULT_TL = 200.0;
% DEFAULT_TD = 0.0;
% 
% set(handles.amp_slider, 'Value', DEFAULT_AMP);
% set(handles.amp_edit, 'String', sprintf('%4.2f', DEFAULT_AMP));
% set(handles.freq_slider, 'Value', DEFAULT_FREQ);
% set(handles.freq_edit, 'String', sprintf('%4.2f', DEFAULT_FREQ));
% set(handles.dur_slider, 'Value', DEFAULT_DUR);
% set(handles.dur_edit, 'String', sprintf('%4.2f', DEFAULT_DUR));
% set(handles.tl_slider, 'Value', DEFAULT_TL);
% set(handles.tl_edit, 'String', sprintf('%4.2f', DEFAULT_TL));
% set(handles.td_slider, 'Value', DEFAULT_TD);
% set(handles.td_edit, 'String', sprintf('%4.2f', DEFAULT_TD));
% set(handles.enable_checkbox, 'Value', 1);

% --- Executes during object creation, after setting all properties.
function found_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to found_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on selection change in stim_list.
function stim_list_Callback(hObject, eventdata, handles)
% hObject    handle to stim_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value_list = get(handles.stim_list, 'Value');
stim_labels = get(handles.stim_list, 'String');
found_labels = get(handles.found_list, 'String');

selected_elec_labels = stim_labels(value_list);
% values_list
% find this entry in the found_list
% for iValue=1:length(value_list)
%     elec = {stim_labels(value_list(iValue))}
% end
found_values = zeros(1, length(value_list));

% go through all the saved stim parameters and find which entry they
% are found in in the found labels list
for i=1:length(selected_elec_labels)
    % strcmp(selected_elec_labels{i}, found_la
    search_labels = find(strcmp(selected_elec_labels(i), found_labels));
    if ~isempty(search_labels)
        found_values(i) = search_labels;
    end
end

set(handles.found_list, 'Value', found_values);
if length(value_list) == 1
   param_index = ...
       [handles.stim_params(:).elec]==str2double(selected_elec_labels);
   set_stim_params_display(hObject, eventdata, handles, ...
       handles.stim_params(param_index));
end

function set_stim_params_display(hObject, eventdata, handles, params)
set(handles.amp_slider, 'Value', params.amp);
set(handles.amp_edit, 'String', sprintf('%4.2f', params.amp));
set(handles.dur_slider, 'Value', params.dur);
set(handles.dur_edit, 'String', sprintf('%4.2f', params.dur));
set(handles.freq_slider, 'Value', params.freq);
set(handles.freq_edit, 'String', sprintf('%4.2f', params.freq));
set(handles.tl_slider, 'Value', params.tl);
set(handles.tl_edit, 'String', sprintf('%4.2f', params.tl));
set(handles.td_slider, 'Value', params.td);
set(handles.td_edit, 'String', sprintf('%4.2f', params.td));
set(handles.enable_checkbox, 'Value', params.enabled);
set(handles.pol_checkbox, 'Value', params.pol);

% --- Executes during object creation, after setting all properties.
function stim_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function amp_edit_Callback(hObject, eventdata, handles)
% hObject    handle to amp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(handles.amp_edit, 'string'));
% Check that this was a valid numeric string
if value ~= value
    error('invalid duration: %s', get(handles.dur_edit, 'string'));
end
value = round(value / handles.stim_rounding.amp)*handles.stim_rounding.amp;
if value > get(handles.amp_slider, 'max') ...
        || value < get(handles.amp_slider, 'min')
    error('duration value: %d is out of bounds', value);
end
set(handles.amp_slider, 'Value', value);
set(handles.amp_edit, 'String', ...
    sprintf('%4.2f', get(handles.amp_slider, 'Value')));

% --- Executes during object creation, after setting all properties.
function amp_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on slider movement.
function dur_slider_Callback(hObject, eventdata, handles)
% hObject    handle to dur_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.dur_slider, 'Value');
set(handles.dur_slider, 'Value', round(value / 0.05)*0.05);
set(handles.dur_edit, 'String', sprintf('%4.2f', get(handles.dur_slider, 'Value')));

% --- Executes during object creation, after setting all properties.
function dur_slider_CreateFcn(hObject, ~, handles)
% hObject    handle to dur_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function dur_edit_Callback(hObject, eventdata, handles)
% hObject    handle to dur_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(handles.dur_edit, 'string'));
% Check that this was a valid numeric string
if value ~= value
    error('invalid duration: %s', get(handles.dur_edit, 'string'));
end
value = round(value / handles.stim_rounding.dur)*handles.stim_rounding.dur;
if value > get(handles.dur_slider, 'max') ...
        || value < get(handles.dur_slider, 'min')
    error('duration value: %d is out of bounds', value);
end
set(handles.dur_slider, 'Value', value);
set(handles.dur_edit, 'String', sprintf('%4.3f', get(handles.dur_slider, 'Value')));

% --- Executes during object creation, after setting all properties.
function dur_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dur_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on slider movement.
function freq_slider_Callback(hObject, eventdata, handles)
% hObject    handle to freq_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.freq_slider, 'Value');
value = round(value/handles.stim_rounding.freq)*handles.stim_rounding.freq;
set(handles.freq_slider, 'Value', value);

set(handles.freq_edit, ...
    'String', sprintf('%4.2f', value));
if get(handles.single_pulse_checkbox, 'Value')
    period = 1000/value;
    set(handles.tl_slider, 'Value', period);
    set(handles.tl_edit, 'string', sprintf('%4.2f', period));
end

% --- Executes during object creation, after setting all properties.
function freq_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function freq_edit_Callback(hObject, eventdata, handles)
% hObject    handle to freq_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(handles.freq_edit, 'string'));
% Check that this was a valid numeric string
if value ~= value
    error('invalid duration: %s', get(handles.freq_edit, 'string'));
end
value = round(value / handles.stim_rounding.freq)*handles.stim_rounding.freq;
if value > get(handles.freq_slider, 'max') ...
        || value < get(handles.freq_slider, 'min')
    error('duration value: %d is out of bounds', value);
end
set(handles.freq_slider, 'Value', value);
set(handles.freq_edit, 'String', sprintf('%4.2f', value));

if get(handles.single_pulse_checkbox, 'Value')
    period = ceil(1e5/value)/100;
    set(handles.tl_slider, 'Value', period);
    set(handles.tl_edit, 'string', sprintf('%4.2f', period));
end
% --- Executes during object creation, after setting all properties.
function freq_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on slider movement.
function tl_slider_Callback(hObject, eventdata, handles)
% hObject    handle to tl_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.tl_slider, 'Value');
% force the train length to be an integer multiple of the stim period.  It
% will be any way, why not make it explicit on the ui.
period = ceil(1.0e5/get(handles.freq_slider, 'Value'))/100;
value = floor(value / period)*period;
% value = round(value / period)*period;
set(handles.tl_slider, 'Value', value);
set(handles.tl_edit, 'String', sprintf('%4.2f', get(handles.tl_slider, 'Value')));

% --- Executes during object creation, after setting all properties.
function tl_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tl_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function tl_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tl_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(handles.tl_edit, 'string'));
% Check that this was a valid numeric string
if value ~= value
    error('invalid duration: %s', get(handles.tl_edit, 'string'));
end
if value > get(handles.tl_slider, 'max') ...
        || value < get(handles.tl_slider, 'min')
    error('duration value: %d is out of bounds', value);
end
period = ceil(1e5/get(handles.freq_slider, 'Value'))/100;
% value = round(value / period)*period;
value = floor(value / period)*period;
set(handles.tl_slider, 'Value', value); 
set(handles.tl_edit, 'String', sprintf('%4.2f', value));

% -----------------------------------------------------------------------------
function tl_edit_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.

% -----------------------------------------------------------------------------
function process_stim(hObject, eventdata, handles, stim_params)



amplitudes = [stim_params(:).amp]/handles.stim_rounding.amp;
fs_value = get(handles.fs_slider, 'Value');
fs_elec = [];
if get(handles.fs_checkbox, 'Value')
    fs_elec = handles.fs_elec; 
end

if get(handles.fs_checkbox, 'Value')
    xippmex('fastsettle','stim',[stim_params(:).elec], 2, ones(length(stim_params(:).elec),1)*fs_value);
else    
    xippmex('fastsettle','stim',[stim_params(:).elec], 1);
end

stim_str = stim_param_to_string([stim_params(:).elec], ...
    [stim_params(:).tl], [stim_params(:).freq], [stim_params(:).dur], ...
    amplitudes, [stim_params.td], [stim_params.pol]);
% If check box to print string is selected, then display stimulation string
if get(handles.print_checkbox, 'Value')
    fprintf('%s\n', stim_str);
end

handles = guidata(hObject);

% clear any existing stim data for the desired electrodes
% xippmex('spike', [stim_params(:).elec], 1);

% send out stim
xippmex('stim', stim_str);

% removing the drawing feature of stimulation for now.  Will bring 
% it back in a new way that can handle asycronous drawing in a new GUI
% wait until the longest stim train length has finished
% pause(max([stim_params(:).tl])/1000);
% bail if we are not interested in drawing the stim.
% if get(handles.display_checkbox, 'Value')
%     stim_display('stimui', handles.main_stim_figure, ...
%     'stim_params', stim_params);
% end

guidata(hObject, handles);

% --- Executes on button press in stim_saved_button.
function stim_saved_button_Callback(hObject, eventdata, handles)
% hObject    handle to stim_saved_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% goes through the list of saved stim channels and produces the 
% stim string for the saved stim channels.
enabled_stim_params = ...
    handles.stim_params([handles.stim_params(:).enabled] == 1);

process_stim(hObject, eventdata, handles, enabled_stim_params)
guidata(hObject, handles);

% -----------------------------------------------------------------------------
function stim_single_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in 'stim single' button.
% Finds the stim parameters for the currently selected stim channel
% and produces the stim string to be sent to the stim interface
% codes
val = get(handles.found_list, 'Value');
found_labels = get(handles.found_list, 'String');
elec = str2double(found_labels(val));

stim_params = [];
stim_params.elec = elec;
stim_params.tl = get(handles.tl_slider, 'Value');
stim_params.dur = get(handles.dur_slider, 'Value');
stim_params.td = get(handles.td_slider, 'Value');
stim_params.freq = get(handles.freq_slider, 'Value');
stim_params.amp = get(handles.amp_slider, 'Value');
stim_params.pol = get(handles.pol_checkbox, 'Value');
%handles.stim_params = [handles.stim_params, stim_params];
process_stim(hObject, eventdata, handles, stim_params);
    
guidata(hObject, handles);

% -----------------------------------------------------------------------------
function td_slider_Callback(hObject, eventdata, handles)
% Executes on slider movement for train delay
% hObject    handle to td_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.td_slider, 'Value');
max_value = get(handles.td_slider, 'Max');
rounded_value = ...
    round(value / handles.stim_rounding.td)*handles.stim_rounding.td;
if rounded_value > max_value
    rounded_value = max_value;
end
set(handles.td_slider, 'Value', rounded_value);
set(handles.td_edit, 'String', sprintf('%4.2f', get(handles.td_slider, 'Value')));

% -----------------------------------------------------------------------------
function td_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to td_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% -----------------------------------------------------------------------------
function td_edit_Callback(hObject, eventdata, handles)
% --- Executes on editting of the train delay edit box
value = str2double(get(handles.td_edit, 'string'));
% Check that this was a valid numeric string
if value ~= value
    error('invalid duration: %s', get(handles.td_edit, 'string'));
end
value = round(value / handles.stim_rounding.td)*handles.stim_rounding.td;
if value > get(handles.td_slider, 'max') ...
        || value < get(handles.td_slider, 'min')
    error('duration value: %d is out of bounds', value);
end
set(handles.td_slider, 'Value', value);
set(handles.td_edit, 'String', sprintf('%4.2f', value));

% -----------------------------------------------------------------------------
function td_edit_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.

% -----------------------------------------------------------------------------
function enable_checkbox_Callback(hObject, eventdata, handles)
% --- Executes on button press in enable_checkbox.

% -----------------------------------------------------------------------------
function add_button_Callback(hObject, eventdata, handles)
% --- Executes on button press of the 'save' button
value_list = get(handles.found_list, 'Value');
for iValue=1:length(value_list)
    val = value_list(iValue);
    found_labels = get(handles.found_list, 'String');
    elec = str2double(found_labels(val));
    % find if this electrode is already listed and replace if so replace it
    index = 0;
    if ~isempty(handles.stim_params)
        index = find([handles.stim_params(:).elec] == elec);
    end
    if isempty(index) || ~index
        index = length(handles.stim_params) + 1;
    end
    handles.stim_params(index).elec = elec;
    handles.stim_params(index).amp = get(handles.amp_slider, 'Value');
    handles.stim_params(index).dur = get(handles.dur_slider, 'Value');
    handles.stim_params(index).freq = get(handles.freq_slider, 'Value');
    handles.stim_params(index).tl = get(handles.tl_slider, 'Value');
    handles.stim_params(index).td = get(handles.td_slider, 'Value');
    handles.stim_params(index).enabled = get(handles.enable_checkbox, 'Value');
    handles.stim_params(index).pol = get(handles.pol_checkbox, 'Value');
    labels = {};

    for i=1:length(handles.stim_params)
        labels = [labels sprintf('%03d', handles.stim_params(i).elec)];
    end
    set(handles.stim_list, 'String', labels);
    % handles.emg_norm = [];
end
% Update handles structure
guidata(hObject, handles);

% -----------------------------------------------------------------------------
function remove_button_Callback(hObject, eventdata, handles)
% --- Executes on button press of the 'remove' button
value_list = get(handles.stim_list, 'Value');
found_labels = get(handles.stim_list, 'String');

for value_index=1:length(value_list)
    elec = str2double(found_labels(value_index));
    handles.stim_params = ...
        handles.stim_params([handles.stim_params(:).elec] ~= elec);
end
labels = {};

for i=1:length(handles.stim_params)
    labels = [labels sprintf('%03d', handles.stim_params(i).elec)];
end
set(handles.stim_list, 'String', labels);
set(handles.stim_list, 'Value', 1);

handles.emg_norm = [];
% Update handles structure
guidata(hObject, handles);

% -----------------------------------------------------------------------------
function single_pulse_checkbox_Callback(hObject, eventdata, handles)
% --- Executes on button press in 'single pulse' checkbox.

value = get(handles.single_pulse_checkbox, 'Value');
if value
    period = ceil(1e5/get(handles.freq_slider, 'Value'))/100;

    set(handles.tl_slider, 'value', period);
    set(handles.tl_edit, 'string', sprintf('%4.2f', period));
    set(handles.tl_slider, 'enable', 'off');
    set(handles.tl_edit, 'enable', 'off');
else
    set(handles.tl_slider, 'enable', 'on');
    set(handles.tl_edit', 'enable', 'on');
end

% -----------------------------------------------------------------------------
function fs_checkbox_Callback(hObject, eventdata, handles)
% --- Executes on button press in 'fast settle' checkbox.

% -----------------------------------------------------------------------------
function fs_slider_Callback(hObject, eventdata, handles)
% --- Executes on move of 'fast settle' slider
% rounds fast settle vaues and updates the fast settle edit box
fs_str = sprintf('%4.2f', get(handles.fs_slider, 'Value'));
set(handles.fs_edit, 'String', fs_str);

% --- Executes during object creation, after setting all properties.
function fs_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function fs_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fs_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function fs_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in display_checkbox.
function display_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to display_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in fast_settle_button.
function fast_settle_button_Callback(hObject, eventdata, handles)
% hObject    handle to fast_settle_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fast_settle_settings('stimui', hObject);

% --- Executes on button press in recruit_button.
function recruit_button_Callback(hObject, eventdata, handles)
% hObject    handle to recruit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
recruit('stimui', hObject);


% --- Executes on button press in pol_checkbox.
function pol_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to pol_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pol_checkbox


% --- Executes on button press in bipolar_button.
function bipolar_button_Callback(hObject, eventdata, handles)
% hObject    handle to bipolar_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value_list = get(handles.found_list, 'Value');
if length(value_list) ~= 2
    error('bipolar stim requires two electrodes');
end

for value_index=1:length(value_list)
    val = value_list(value_index);
    found_labels = get(handles.found_list, 'String');
    elec = str2double(found_labels(val));
    % find if this electrode is already listed and replace if so replace it
    index = 0;
    if ~isempty(handles.stim_params)
        index = find([handles.stim_params(:).elec] == elec);
    end
    if isempty(index) || ~index
        index = length(handles.stim_params) + 1;
    end
    handles.stim_params(index).elec = elec;
    handles.stim_params(index).amp = get(handles.amp_slider, 'Value');
    handles.stim_params(index).dur = get(handles.dur_slider, 'Value');
    handles.stim_params(index).freq = get(handles.freq_slider, 'Value');
    handles.stim_params(index).tl = get(handles.tl_slider, 'Value');
    handles.stim_params(index).td = get(handles.td_slider, 'Value');
    handles.stim_params(index).enabled = get(handles.enable_checkbox, 'Value');
    pol_select = get(handles.pol_checkbox, 'Value');
    if mod(value_index, 2)
        handles.stim_params(index).pol = pol_select;
    else
        handles.stim_params(index).pol = ~pol_select;
    end
    
    labels = {};

    for i=1:length(handles.stim_params)
        labels = [labels sprintf('%03d', handles.stim_params(i).elec)];
    end
    set(handles.stim_list, 'String', labels);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in print_checkbox.
function print_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to print_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of print_checkbox
