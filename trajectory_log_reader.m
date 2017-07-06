%  Copyright (C) 2017 Paul Ionele - All Rights Reserved
%  You may NOT use, distribute, or modify this code unless express
%  written and signed consent has been given by the author Paul Ionele.

%%
% This is the primary script for reading TLF binary files, performing scale
% transformations, processing TLF data (phase sorting, then arc sorting
% machine axis data), and plotting.
%%

clear;clc

%Integers and floats are stored in 'little endian' format.
%ens = 'ieee-le'; %order for numerical interpretation of byte sequence

directory = 'data_in'; %change this if required
pathname  = fullfile(pwd,directory,filesep);

try
    [filename,pathname] = uigetfile([pathname,'*.bin'],'Select the TLF Binary File');
    filepath = fullfile(pathname, filename);
    fid1 = fopen(filepath,'rb');
catch exception
    %Can add exceptions here later.
    disp(exception.identifier)
    break
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READING TLF

%Header.
header.signiture   = fread(fid1, 16, 'char=>char')'; %*char is shorthand
header.version     = fread(fid1, 16, '*char')';
header.hsize       = fread(fid1, 1, '*int')';
header.sampling    = fread(fid1, 1, '*int')';
header.num_axes    = fread(fid1, 1, '*int')';

%Machine parameters/axis added to header.
header.axes_enumeration = fread(fid1, header.num_axes, '*int');
header.axes_sampling    = fread(fid1, header.num_axes, '*int');
header.axis_scale       = fread(fid1, 1, '*int');
header.num_sbeams       = fread(fid1, 1, '*int');
header.truncation       = fread(fid1, 1, '*int');
header.num_snaps        = fread(fid1, 1, '*int');
header.mlc_model        = fread(fid1, 1, '*int');
header.reserved1        = fread(fid1, 1024 - (64+header.num_axes*8));

%%%Subbeams.
for i = 1:header.num_sbeams
    subbeam(i).cp         = fread(fid1, 1, '*int');
    subbeam(i).mu         = fread(fid1, 1, '*float');
    subbeam(i).rad_time   = fread(fid1, 1, '*float');
    subbeam(i).sbeam_seq  = fread(fid1, 1, '*int');
    subbeam(i).sbeam_name = fread(fid1, 512, '*char')';
    subbeam(i).reserved2  = fread(fid1, 32, '*char')';
end

total_samples = sum(header.axes_sampling); %parameters per snapshot

%axis_data is where the machine information extracted from the TLF is
%stored. It is a 3D array, with the first dimension (rows) as the number of
%snapshots (recordings, 1 recording every 20 ms), the second dimension
%(cols) as the number of machine parameters or axes that are
%monitored/recorded (this is a fixed number), and the third dimension
%having a length of 2. In the third dimension, the first values indexed for
%any given m or n, represent the EXPECTED machine values. The second values
%indexed represent the ACTUAL machine values recorded. Overall, we can
%think of the the axis_data array as two overlapping planes where the first
%place consists of expected values and the second plane of actual values.
axis_data = zeros(header.num_snaps, total_samples, 2); %2 -> expected and actual

%Data is read from TLF and stored immediately in axis_data array.
for i = 1:header.num_snaps  %snapshot every 20 ms
    for j = 1:total_samples %parameters per snapshot
        for k = 1:2         %expected (1) and actual (2) values
            axis_data(i,j,k) = fread(fid1, 1, '*float');
        end
    end
end
fclose(fid1);

%Removing recordings from TLF; based on initial beam hold. While the beam
%is initially held, those recordings (all cooresponding axis/parameters)
%are removed. When the beam counter is not held, the MU for the
%cooresponding recording is reset to zero.
rc = 0; %to count number of lines removed
while 1
    if axis_data(1, 14, 2) == 2
        axis_data(1,:,:) = [];
        rc = rc + 1;
    else
        break
    end
end
axis_data(1,13,:) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SCALE TRANFORMATIONS.

%First snap is at 0 ms (?), and each preceeding snap is at 20*n ms; where n
%is the snap number. So an array of times:
tlf_times = single(0:20:(20*header.num_snaps - 1)); %do we lose a sample time off the end???
tlf_times = tlf_times(1:end-rc); %tlf truncated

%Collimator Rotation
collrot_e = axis_data(:,1,1)';
collrot_a = axis_data(:,1,2)';

collrot_a_iec121 = arrayfun( @(x) mod(180 - x, 360), collrot_a);

%Gantry Rotation
gantrot_e = axis_data(:,2,1)';
gantrot_a = axis_data(:,2,2)';

gantrot_a_iec121 = arrayfun( @(x) mod(180 - x, 360), gantrot_a);

%Y1
y1_e = axis_data(:,3,1)';
y1_a = axis_data(:,3,2)';

y1_a_iec121 = arrayfun( @(x) - x, y1_a);

%Y2
y2_e = axis_data(:,4,1)';
y2_a = axis_data(:,4,2)';

y2_a_iec121 = y2_a;

%X1
x1_e = axis_data(:,5,1)';
x1_a = axis_data(:,5,2)';

x1_a_iec121 = arrayfun( @(x) - x, x1_a);

%X2

x2_e = axis_data(:,6,1)';
x2_a = axis_data(:,6,2)';

x2_a_iec121 = x2_a;

%Couch Vrt
couchvrt_e = axis_data(:,7,1)';
couchvrt_a = axis_data(:,7,2)';

couchvrt_a_iec121 = arrayfun( @(x) -(x - 100), couchvrt_a);

%Couch Lng
couchlng_e = axis_data(:,8,1)';
couchlng_a = axis_data(:,8,2)';

couchlng_a_iec121 = couchlng_a;

%Couch Lat
couchlat_e = axis_data(:,9,1)';
couchlat_a = axis_data(:,9,2)';

couch_lat_iec121 = arrayfun( @(x) x - 100, couchlat_a);

%Couch Rtn
couchrot_e = axis_data(:,10,1)';
couchrot_a = axis_data(:,10,2)';

couchrot_a_iec121 = arrayfun( @(x) mod(180 - x, 360), couchrot_a);

%Couch Pit
%Unused; couch does not support this function.
%Value stored is largest FP number represented in single precsn 32-bits.
couchpit_e = axis_data(:,11,1)';
couchpit_a = axis_data(:,11,2)';

%Couch Rol
%Unused; couch does not support this function.
couchrol_e = axis_data(:,12,1)';
couchrol_a = axis_data(:,12,2)';

%Monitor Units
mu_e = axis_data(:,13,1)';
mu_a = axis_data(:,13,2)';

%Beam Hold
%When integer value is 2, the beam is held.
%Bits are flipped b/c want beam on as integer value 1; not() op. gives
%logical array.
beamh_e = not( axis_data(:,14,1)' / 2 );
beamh_a = not( axis_data(:,14,2)' / 2 );

%Control Pt
cp_e = axis_data(:,15,1)';
cp_a = axis_data(:,15,2)';

%%% MLC Things
%Carriage A Positions
cara_e = axis_data(:,16,1)';
cara_a = axis_data(:,16,2)';
%Carriage B Positions
carb_e = axis_data(:,17,1)';
carb_a = axis_data(:,17,2)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER PROMPT FOR MW OR VXP SELECTION.

[file1, PATHNAME] = uigetfile(fullfile(pwd,'data_in','*.VXP;*MW*'),'select file');
file1 = fullfile(PATHNAME, file1);

if strfind(file1, 'MW') > 0
    %File selected is MW file. Check isn't robust for all MW variations.
    [amplitude,phase,rpm_times,beamenable] = mw_reader(file1);
    
    %Removing recordings from MW; based on initial beam hold. While the
    %beam is initially held, those recordings (all cooresponding
    %axis/parameters) are removed.
    
    try
        while 1
            if beamenable(1) == 0
                amplitude(1)  = [];
                phase(1)      = [];
                rpm_times(1)  = [];
                beamenable(1) = [];
            else
                break
            end
        end
    catch exception
        if strcmp(exception.identifier, 'MATLAB:badsubscript')
            disp('No "beam on" instances present in MW recording.')
            disp('Exiting program...')
            break
        else
            disp(exception.identifier)
            error('An unexpected error has occurred.')
        end
    end
else
    %VXP file.
    [amplitude,phase,timestamp,validflag,ttlin,ttlout,mark,headerv] = vxp_reader(file1);
    rpm_times = cell2mat(timestamp) - timestamp{1,1};
    clearvars timestamp validflag ttlin ttlout
    phase     = cell2mat(phase);
    % amplitude = cell2mat(amplitude);
    % mark      = cell2mat(mark);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SYNCHRONIZE MW AND TLF.

%Resampling MW at increasing sampling frequency (TLF; 20 ms). tlf_times are
%the query points for 'beamenable'. If tlf_times(end) > rpm_times(end), the
%interp1 function will return NaN values for beamenable_q beyond
%rpm_times(end). These NaN values are removed from beamenable_q, AND
%tlf_times is shorted an amount cooresponding to the NaN values removed.
%This is probably unecessary as it seems the TLF recording is always
%shorter than the RPM recording.


% beamenable_q = interp1(rpm_times, beamenable, tlf_times); %queried points
% tlf_times(find(isnan(beamenable_q)))    = 0;
% beamenable_q(find(isnan(beamenable_q))) = 0;
% 
% %Any points that do NOT have bitvals 0 or 1, need to be reassigned. This is
% %accomplished using a function based on that by Steven Thomas and
% %implemented in the trajectory_log_phase_sort function.
% 
% % Now beamenable_q is in tlf_time. We need the match such that the last
% % beam on point is at tlf_time = 1.91*10^4 ms.

%%% Implemented process for synchronization. The TLF 'leads' the MW in
%%% time, so a simple shift applied to to the TLF can sync fairly well.

tlf_times = tlf_times + 4700; %yup, that's it.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION CALLS FOR PHASE SORTING FOLLOWED BY ARC SORTING.

%Sorting the TLF file information into 10 different phases.
[sorted_phase, phase_tlf2] = trajectory_log_phase_sort(tlf_times, rpm_times, phase);

%Sorting the 10 different phases into n (1, 2, or 3) arcs.
[sorted_phase_arc, intra_arc] = arc_separator(cp_a, subbeam, sorted_phase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTTING.
figure;

plot(rpm_times,beamenable,'b')
hold on
plot(tlf_times,beamh_a,'r')
ylim([-0.5,1.5])
legend('RPM/MW Recording', 'TLF Recording')

% t_tlf = 1.91*10^4;
% t_mw  = 2.366*10^4;
%
% delta = t_mw - t_tlf;










%%%Plotting
% sz = 25;
% scatter(tlf_times(sorted_phase_arc{1,1}), mu_a(sorted_phase_arc{1,1}),sz,'r','filled')
% hold on
% scatter(tlf_times(sorted_phase_arc{2,1}), mu_a(sorted_phase_arc{2,1}),sz,'b','filled')
% scatter(tlf_times(sorted_phase_arc{3,1}), mu_a(sorted_phase_arc{3,1}),sz,'c','filled')
% scatter(tlf_times(sorted_phase_arc{5,1}), mu_a(sorted_phase_arc{5,1}),sz,'g','filled')
% scatter(tlf_times(sorted_phase_arc{10,1}), mu_a(sorted_phase_arc{10,1}),sz,'m','filled')
%
% scatter(tlf_times(sorted_phase_arc{1,2}), mu_a(sorted_phase_arc{1,2}),sz,'r','filled')
% scatter(tlf_times(sorted_phase_arc{2,2}), mu_a(sorted_phase_arc{2,2}),sz,'b','filled')
% scatter(tlf_times(sorted_phase_arc{3,2}), mu_a(sorted_phase_arc{3,2}),sz,'c','filled')
% scatter(tlf_times(sorted_phase_arc{5,2}), mu_a(sorted_phase_arc{5,2}),sz,'g','filled')
% scatter(tlf_times(sorted_phase_arc{10,2}), mu_a(sorted_phase_arc{10,2}),sz,'m','filled')


% plot(tlf_times(subbeam(1).arc(1):subbeam(1).arc(2)), cp_a(subbeam(1).arc(1):subbeam(1).arc(2)))
% hold on
% plot(tlf_times(subbeam(2).arc(1):subbeam(2).arc(2)), cp_a(subbeam(2).arc(1):subbeam(2).arc(2)),'-g')
% plot(tlf_times(intra_arc{1}(1):intra_arc{1}(2)), cp_a(intra_arc{1}(1):intra_arc{1}(2)),'-r')