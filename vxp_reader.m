function [amplitude,phase,timestamp,validflag,ttlin,ttlout,mark,header] = vxp_reader(file1)
%Function to call when reading export of "motion data files" for
%postprocessing of respiration motion data (v1.6 and newer). Files are
%recorded in a test format after CT session. File is written in text
%format.

%File consists of two main sections: [Header], [Data].

%The structure of the [Header] is as follows:
% [Header] --- indicates the beginning of the header section
% CRC=27051 --- contains CRC checksum
% Version=1.6 --- version number of export file
% Data_layout=amplitude,phase,timestamp,validflag,ttlin,mark,ttlout
% Patient_ID=000000
% Date=08-30-2007 --- (month-day-year)
% Total_study_time=247.582 --- actual recording length in sec to 0.001s
% Samples_per_second=30 --- NTSC 30 or 15, CCIR 25 or 12.5
% Scale_factor=10.0 --- scale factor for signal data in mm

% [Data] %indicates the beginning of the header section
% Data structure is specified in the Data_layout field.
% amplitude,phase,timestamp,validflag,ttlin,mark,ttlout

fid = fopen(file1,'r');

%Header container.
header = {};

%Data sructures to be populated.
amplitude = {};
phase     = {};
timestamp = {};
validflag = {};
ttlin     = {};
mark      = {};
ttlout    = {};

for i = 1:10
    if i == 1 || i == 10
        %Don't read [Header] or [Data] lines.
        tchar = fgetl(fid);
        continue
    else
        tchar = fgetl(fid);
        ts = regexp(tchar,'=','split');        
        if any(strcmp(ts,'Date')) || any(strcmp(ts,'Data_layout')) || any(strcmp(ts,'Patient_ID'))
            %String-type data fields.
            header = [header, ts{1,2}];
        else
            %Numeric-type data fields.
            header = [header, str2double(ts{1,2})];
        end
    end
end 

%Loop for the data structure.
while ischar(tchar)
    tchar = fgetl(fid);
    if ischar(tchar) == 1
        ts = textscan(tchar, '%s', 'Delimiter', ',');

        st1 = str2double(ts{1}{1});
        st2 = str2double(ts{1}{2});
        st3 = str2double(ts{1}{3});
        st4 = str2double(ts{1}{4});
        st5 = str2double(ts{1}{5});
        st7 = str2double(ts{1}{7});
        
        if any(strcmp(ts{1}{6},''))
            st6 = 0;
        elseif any(strcmp(ts{1}{6},'Z'))
            %If Z; zero phase.
            st6 = 1;
        else
            %If P; pi phase.
            st6 = 2;
        end

        amplitude = [amplitude ,st1];
        phase     = [phase, st2];
        timestamp = [timestamp, st3];
        validflag = [validflag, st4];
        ttlin     = [ttlin, st5];
        mark      = [mark, st6];
        ttlout    = [ttlout, st7];
    end
end

fclose(fid);