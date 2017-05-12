clear; clc

% CRC=27051 --- contains CRC checksum
% Version=1.6 --- version number of export file
% Data_layout=amplitude,phase,timestamp,validflag,ttlin,mark,ttlout
% Patient_ID=000000
% Date=08-30-2007 --- (month-day-year)
% Total_study_time=247.582 --- actual recording length in sec to 0.001s
% Samples_per_second=30 --- NTSC 30 or 15, CCIR 25 or 12.5
% Scale_factor=10.0 --- scale factor for signal data in mm
%  For ex: if value_resp_wave (amplitude) in cm, the scale factor is 10 for mm rep.

% file1 = 'TT4DSIN_6322_2.vxp';
file1 = uigetfile('.VXP');
[amplitude,phase,timestamp,validflag,ttlin,ttlout,mark,header] = vxp_reader(file1);

%Related to time...
sample_time  = header{1,6}; %total sample time (s)
sample_freq  = header{1,7}; %sampling frequency
sample_inte  = 1/sample_freq;
sample_total = sample_time * sample_freq;
times        = linspace(0,sample_time, sample_total);
timestamp    = cell2mat(timestamp) - timestamp{1,1};

%Phase, amplitude, and mark marticies .
phase     = cell2mat(phase);
amplitude = cell2mat(amplitude);
mark      = cell2mat(mark);

%Recentering ampltude around mean value...
% amplitude = amplitude - mean(amplitude);

%Determination of the indicies when the phase value is closest to phase 0
%or phase pi.
phase_zero = find(mark == 1);
phase_pi   = find(mark == 2);

%Amplitude plot.
fig00 = figure;
figure(fig00);
plot(timestamp, amplitude)
hold on
plot(timestamp(phase_zero), amplitude(phase_zero),'-ro')
plot(timestamp(phase_pi), amplitude(phase_pi),'-go')

%Phase plot.
fig0 = figure;
figure(fig0);
plot(timestamp,phase);
hold on
plot(timestamp(phase_zero), phase(phase_zero),'-ro');
plot(timestamp(phase_pi), phase(phase_pi),'-go');


%Other figures.
fig1 = figure;
fig2 = figure;

figure(fig1);
plot(timestamp/1000,phase);
xlim([0,timestamp(end)/1000]);
ylim([0,2*pi]);
title('Phase')

figure(fig2);
plot(timestamp/1000,amplitude);
xlim([0,timestamp(end)/1000]);
title('Amplitude')



% fid = fopen(file1,'r');
% 
% %Header container.
% header = {};
% 
% %Data sructures to be populated.
% amplitude = {};
% phase = {};
% timestamp = {};
% validflag = {};
% ttlin = {};
% mark = {};
% ttlout = {};
% 
% for i = 1:10
%     if i == 1 || i == 10
%         %Don't read [Header] or [Data] lines.
%         tchar = fgetl(fid);
%         continue
%     else
%         tchar = fgetl(fid);
%         ts = regexp(tchar,'=','split');        
%         if any(strcmp(ts,'Date')) || any(strcmp(ts,'Data_layout'))
%             %String-type data fields.
%             header = [header, ts{1,2}];
%         else
%             %Numeric-type data fields.
%             header = [header, str2double(ts{1,2})];
%         end
%     end
% end 
% 
% %Loop for the data structures.
% while ischar(tchar)
%     tchar = fgetl(fid);
%     if ischar(tchar) == 1
%         ts = textscan(tchar, '%s', 'Delimiter', ',');
% 
%         st1 = str2double(ts{1}{1});
%         st2 = str2double(ts{1}{2});
%         st3 = str2double(ts{1}{3});
%         st4 = str2double(ts{1}{4});
%         st5 = str2double(ts{1}{5});
%         st7 = str2double(ts{1}{7});
%         
%         if any(strcmp(ts{1}{6},''))
%             st6 = 0;
%         else
%             st6 = ts{1}{6};
%         end
% 
%         amplitude = [amplitude ,st1];
%         phase = [phase, st2];
%         timestamp = [timestamp, st3];
%         validflag = [validflag, st4];
%         ttlin = [ttlin, st5];
%         mark = [mark, st6];
%         ttlout = [ttlout, st7];
%     end
% end
% 
% fclose(fid);