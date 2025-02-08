function varargout = WT(varargin)
% WT MATLAB code for WT.fig
%      WT, by itself, creates a new WT or raises the existing
%      singleton*.
%
%      H = WT returns the handle to a new WT or the handle to
%      the existing singleton*.
%
%      WT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WT.M with the given input arguments.
%
%      WT('Property','Value',...) creates a new WT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WT

% Last Modified by GUIDE v2.5 14-Jul-2020 14:48:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WT_OpeningFcn, ...
                   'gui_OutputFcn',  @WT_OutputFcn, ...
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


% --- Executes just before WT is made visible.
function WT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WT (see VARARGIN)

% Choose default command line output for WT
handles.output = hObject;

%initialize all variables
handles.sr = 500;
handles.samples_to_save = 20;
handles.num_sent = 0;
handles.da = 0;
handles.dv = 0;
handles.control_airspeed = 1;
handles.act_airspeed = 10*rand(1);
handles.act_voltage = rand(1)/10;
handles.sim = 0;
handles.ref1_1 = 'T1';
handles.ref1_2 = 'T1';
handles.ref2_1 = 'T1';
handles.ref2_2 = 'T1';
handles.ref3_1 = 'T1';
handles.ref3_2 = 'T1';
handles.ref4_1 = 'T1';
handles.ref4_2 = 'T1';
handles.p_atm = 82737.7; %[pa]
handles.rho = .95; %air density [kg/m^3]
handles.temp = 300; %atmospheric temperature [K]
handles.max_num_save = 100; %maximum number of samples to save per button press
handles.max_sr = 1000; %maximum sample rate allowable [Hz]


%get the file the user wishes to save data to.  Write the header info
[file, path] = uiputfile('*.txt');
handles.fname = strcat(path, file);
fid = fopen(handles.fname, 'w');
fprintf(fid, 'Time [sec], Pressure 1 [psi], Pressure 2 [psi], Airspeed [m/s], Voltage [V]\n');
fclose(fid);

% Update handles structure
guidata(hObject, handles);

%This looked bad so removing it.  Keeping code here for reference.
% Add the image of the WT
% axes(handles.wt_image);
% matlabImage = imread('wt_photoshop.png');
% image(matlabImage);
% axis off
% axis image


% UIWAIT makes WT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = WT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in send_samples.
function send_samples_Callback(hObject, eventdata, handles)
% hObject    handle to send_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function sample_rate_Callback(hObject, eventdata, handles)
% hObject    handle to sample_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sample_rate as text
%        str2double(get(hObject,'String')) returns contents of sample_rate as a double
handles.sr = str2double(get(hObject,'String'));

%check that the sample rate doesn't exceed the limit of 1kHz, less than 1,
%or other bad format
handles.sr = check_sample_format(handles.sr, handles.max_sr);
set(hObject, 'String', num2str(handles.sr));

% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function sample_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function samples_per_set_Callback(hObject, eventdata, handles)
% hObject    handle to samples_per_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of samples_per_set as text
%        str2double(get(hObject,'String')) returns contents of samples_per_set as a double
handles.samples_to_save = str2double(get(hObject,'String'));

%check that the number of samples to save is also well formatted, similarly
%to the sample rate
handles.samples_to_save = check_sample_format(handles.samples_to_save, ...
    handles.max_num_save);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function samples_per_set_CreateFcn(hObject, eventdata, handles)
% hObject    handle to samples_per_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function samples_sent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to samples_sent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




function desired_airspeed_Callback(hObject, eventdata, handles)
% hObject    handle to desired_airspeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of desired_airspeed as text
%        str2double(get(hObject,'String')) returns contents of desired_airspeed as a double


% --- Executes during object creation, after setting all properties.
function desired_airspeed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to desired_airspeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function desired_voltage_Callback(hObject, eventdata, handles)
% hObject    handle to desired_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of desired_voltage as text
%        str2double(get(hObject,'String')) returns contents of desired_voltage as a double
%handles.desired_voltage = str2double(get(hObject,'String'));

% Update handles structure
%guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function desired_voltage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to desired_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in airspeed_voltage_toggle.
function airspeed_voltage_toggle_Callback(hObject, eventdata, handles)
% hObject    handle to airspeed_voltage_toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of airspeed_voltage_toggle


% --- Executes on button press in zero_airspeed.
function zero_airspeed_Callback(hObject, eventdata, handles)
% hObject    handle to zero_airspeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in start_sim.
function start_sim_Callback(hObject, eventdata, handles)
% hObject    handle to start_sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%grab all of the relative text field entries.  This is added here instead
%of the regular way since this is more of a 'real-time' application so it
%has to be set up a little differently.
handles.samples_to_save = check_sample_format( ... 
    str2double(get(handles.samples_per_set, 'String')), handles.max_num_save);
handles.sr = check_sample_format( ... 
    str2double(get(handles.sample_rate, 'String')), handles.max_sr);
set(handles.send_samples, 'value', 0); %set to 0 in case it was pressed before the sim began
handles.p_atm = str2double(get(handles.atm_pressure, 'String'));
handles.rho = str2double(get(handles.atm_density, 'String'));
handles.temp = str2double(get(handles.atm_temp, 'String'));


%grab all of the sensor references
contents = cellstr(get(handles.ref1_first,'String'));
handles.ref1_1 = contents{get(handles.ref1_first,'Value')};
contents = cellstr(get(handles.ref1_second,'String'));
handles.ref1_2 = contents{get(handles.ref1_second,'Value')};

contents = cellstr(get(handles.ref2_first,'String'));
handles.ref2_1 = contents{get(handles.ref2_first,'Value')};
contents = cellstr(get(handles.ref2_second,'String'));
handles.ref2_2 = contents{get(handles.ref2_second,'Value')};

contents = cellstr(get(handles.ref3_first,'String'));
handles.ref3_1 = contents{get(handles.ref3_first,'Value')};
contents = cellstr(get(handles.ref3_second,'String'));
handles.ref3_2 = contents{get(handles.ref3_second,'Value')};

contents = cellstr(get(handles.ref4_first,'String'));
handles.ref4_1 = contents{get(handles.ref4_first,'Value')};
contents = cellstr(get(handles.ref4_second,'String'));
handles.ref4_2 = contents{get(handles.ref4_second,'Value')};

%create a bunch of relevant variables for the simulation
stop_sim = 0;
pause_time = 0.2;
display_time = 5; %show 5 seconds of data at all times
offset1 = randn(1) * 50;
offset2 = randn(1) * 75;
pressure_data1 = zeros(1,display_time*handles.sr);
pressure_data2 = zeros(1,display_time*handles.sr);
time = linspace(1/handles.sr,display_time,handles.sr*display_time);
num_points = floor(handles.sr*pause_time);
axes(handles.pressure_plot);
fid = fopen(handles.fname, 'a+');
sts = handles.samples_to_save-1; %samples to save
as_offset = randn(1)*10; %random airspeed offset displayed to the measured airspeed field


while(~stop_sim)
    %plot the current values for pressure
    pressure_data = [pressure_data1; pressure_data2];
    set(gca,'color', 'k')
    set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'w','w','w'})
    plot(handles.pressure_plot, time, pressure_data);
    xlabel(handles.pressure_plot,'Time (s)');
    ylabel(handles.pressure_plot,'Pressure (Pa)');
    legend('Reference 1', 'Reference 2', 'location', 'NorthWest', 'color', 'w');
    pause(pause_time);
    
    %get the updated values if they have changed
    handles.dv = str2double(get(handles.desired_voltage, 'String'));
    handles.da = str2double(get(handles.desired_airspeed, 'String'));
    stop_sim = get(handles.stop_simulation, 'Value');
    
    
    %check if we are controlling the airspeed or the voltage.  Update the
    %desired airspeed accordingly
    handles.control_airspeed = ~get(handles.airspeed_voltage_toggle, 'value');
    if(~handles.control_airspeed)
        desired_airspeed = handles.dv * 6.09 - 1.78;
    else
        desired_airspeed = handles.da;
    end
    
    %check that the desired airspeed is not too high or negative
    max_as = 50; %maximum airspeed is 50 m/s
    if(desired_airspeed > max_as) 
        desired_airspeed = 50;
        set(handles.desired_airspeed, 'String', num2str(desired_airspeed));
        set(handles.desired_voltage, 'String', num2str(8.5));
    end
    
    if(desired_airspeed < 0)
        desired_airspeed = 0;
        set(handles.desired_airspeed, 'String', '0');
        set(handles.desired_voltage, 'String', '0');
    end
    
    %generate the next set of pressure points
    time = time+pause_time;
    pressure_data1(1:num_points) = [];
    pressure_data2(1:num_points) = [];
    new_pressure1 = gen_pressure(pressure_data1, desired_airspeed, ...
        handles.ref1_1, handles.ref1_2, pause_time, num_points, offset1, ...
        handles.p_atm, handles.rho);
    new_pressure2 = gen_pressure(pressure_data2, desired_airspeed, ...
        handles.ref2_1, handles.ref2_2, pause_time, num_points, offset2, ... 
        handles.p_atm, handles.rho);
    
    %and update the vectors
    pressure_data1 = [pressure_data1 new_pressure1];
    pressure_data2 = [pressure_data2 new_pressure2];
    
    %update the actual airspeed and actual voltage
    handles.act_airspeed = desired_airspeed + randn(1)/50;
    handles.act_voltage = (handles.act_airspeed+1.78)/6.09 + randn(1)/1000;
    
    set(handles.actual_airspeed, 'String', handles.act_airspeed+as_offset);
    set(handles.actual_voltage, 'String', handles.act_voltage);
    
    %update some other values
    set(handles.atm_temp, 'String', num2str(handles.temp+randn(1)/10));
    set(handles.atm_pressure, 'String', num2str(handles.p_atm+randn(1)*10));
    set(handles.atm_density, 'String', num2str(handles.rho+randn(1)/100));
    
    %check if the save samples button was pressed
    save_samps = get(handles.send_samples, 'value');
    if(save_samps)
        handles.num_sent = handles.num_sent + 1;
        set(handles.samples_sent, 'String', num2str(handles.num_sent));
        set(handles.send_samples, 'value', 0);
        data = [time(end-sts:end)' pressure_data1(end-sts:end)' pressure_data2(end-sts:end)' ... 
            handles.act_airspeed+randn(sts+1,1) handles.act_voltage+randn(sts+1,1)]';
        fprintf(fid, '%f, %f, %f, %f, %f\n', data);
    end
    
    %check if the zero airspeed was pressed
    zero_stuff = get(handles.zero_airspeed, 'value');
    if(zero_stuff)
        pause(1)
        offset1 = randn(1);
        offset2 = randn(1);
        as_offset = 0;
        set(handles.zero_airspeed, 'value', 0);
    end
    
    %deal with the utube manometer plots
    p3 = compute_goal(handles.ref3_1, handles.ref3_2, desired_airspeed, ...
        handles.p_atm, handles.rho);
    p4 = compute_goal(handles.ref4_1, handles.ref4_2, desired_airspeed, ...
        handles.p_atm, handles.rho);
    [x1, y1, diff1] = gen_manometer(p3, 1000); %the second number is rho of water [kg/m^3]
    [x2, y2, diff2] = gen_manometer(p4, 1000);
    u1 = plot(handles.utube1, x1, y1, '.r');
    grid(handles.utube1, 'on');
    ylabel('Height [mm]');
    u2 = plot(handles.utube2, x2, y2, '.g');
    grid(handles.utube2, 'on');
    ylabel(handles.utube1, 'Height [mm]');
    
    %update the pressure and difference values in the text boxes
    set(handles.pressure1_output, 'String', num2str(pressure_data1(end)));
    set(handles.pressure2_output, 'String', num2str(pressure_data2(end)));
    set(handles.utube1_diff_box, 'String', num2str(diff1*1000));
    set(handles.utube2_diff_box, 'String', num2str(diff2*1000));
    
end

set(handles.stop_simulation, 'value', 0);
fclose(fid);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in stop_simulation.
function stop_simulation_Callback(hObject, eventdata, handles)
% hObject    handle to stop_simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stop_simulation


% --- Executes on selection change in ref1_first.
function ref1_first_Callback(hObject, eventdata, handles)
% hObject    handle to ref1_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ref1_first contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ref1_first
%contents = cellstr(get(hObject,'String'));
%handles.ref1_1 = contents{get(hObject,'Value')};

% Update handles structure
%guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ref1_first_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref1_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ref1_second.
function ref1_second_Callback(hObject, eventdata, handles)
% hObject    handle to ref1_second (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ref1_second contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ref1_second
contents = cellstr(get(hObject,'String'));
handles.ref1_2 = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ref1_second_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref1_second (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ref2_first.
function ref2_first_Callback(hObject, eventdata, handles)
% hObject    handle to ref2_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ref2_first contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ref2_first
contents = cellstr(get(hObject,'String'));
handles.ref2_1 = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ref2_first_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref2_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ref2_second.
function ref2_second_Callback(hObject, eventdata, handles)
% hObject    handle to ref2_second (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ref2_second contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ref2_second
contents = cellstr(get(hObject,'String'));
handles.ref2_2 = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ref2_second_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref2_second (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ref3_first.
function ref3_first_Callback(hObject, eventdata, handles)
% hObject    handle to ref3_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ref3_first contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ref3_first
contents = cellstr(get(hObject,'String'));
handles.ref3_1 = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ref3_first_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref3_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ref3_second.
function ref3_second_Callback(hObject, eventdata, handles)
% hObject    handle to ref3_second (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ref3_second contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ref3_second
contents = cellstr(get(hObject,'String'));
handles.ref3_2 = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ref3_second_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref3_second (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ref4_first.
function ref4_first_Callback(hObject, eventdata, handles)
% hObject    handle to ref4_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ref4_first contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ref4_first
contents = cellstr(get(hObject,'String'));
handles.ref4_1 = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ref4_first_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref4_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ref4_second.
function ref4_second_Callback(hObject, eventdata, handles)
% hObject    handle to ref4_second (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ref4_second contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ref4_second
contents = cellstr(get(hObject,'String'));
handles.ref4_2 = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ref4_second_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref4_second (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pressure1_output_Callback(hObject, eventdata, handles)
% hObject    handle to pressure1_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pressure1_output as text
%        str2double(get(hObject,'String')) returns contents of pressure1_output as a double


% --- Executes during object creation, after setting all properties.
function pressure1_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pressure1_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pressure2_output_Callback(hObject, eventdata, handles)
% hObject    handle to pressure2_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pressure2_output as text
%        str2double(get(hObject,'String')) returns contents of pressure2_output as a double


% --- Executes during object creation, after setting all properties.
function pressure2_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pressure2_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function utube1_diff_box_Callback(hObject, eventdata, handles)
% hObject    handle to utube1_diff_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of utube1_diff_box as text
%        str2double(get(hObject,'String')) returns contents of utube1_diff_box as a double


% --- Executes during object creation, after setting all properties.
function utube1_diff_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to utube1_diff_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function utube2_diff_box_Callback(hObject, eventdata, handles)
% hObject    handle to utube2_diff_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of utube2_diff_box as text
%        str2double(get(hObject,'String')) returns contents of utube2_diff_box as a double


% --- Executes during object creation, after setting all properties.
function utube2_diff_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to utube2_diff_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
