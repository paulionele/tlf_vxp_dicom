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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READING TLF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Processing the TLF binary file...')

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
axis_data = zeros(header.num_snaps, total_samples, 2);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We will proceed with using the expected values "_e", since that is what is
%used in Tony's code.

%First snap is at 0 ms (?), and each preceeding snap is at 20*n ms; where n
%is the snap number. So an array of times:
tlf_times = single(0:20:(20*header.num_snaps - 1)); %do we lose a sample time off the end???
tlf_times = tlf_times(1:end - rc); %tlf truncated

%Collimator Rotation
collrot_e = axis_data(:,1,1)';
collrot_a = axis_data(:,1,2)';

collrot_a_iec121 = arrayfun( @(x) mod(180 - x, 360), collrot_a);
collrot_e_iec121 = arrayfun( @(x) mod(180 - x, 360), collrot_e);

%Gantry Rotation
gantrot_e = axis_data(:,2,1)';
gantrot_a = axis_data(:,2,2)';

gantrot_a_iec121 = arrayfun( @(x) mod(180 - x, 360), gantrot_a);
gantrot_e_iec121 = arrayfun( @(x) mod(180 - x, 360), gantrot_e);

%Y1
y1_e = axis_data(:,3,1)';
y1_a = axis_data(:,3,2)';

y1_a_iec121 = arrayfun( @(x) - x, y1_a);
y1_e_iec121 = arrayfun( @(x) - x, y1_e);

%Y2
y2_e = axis_data(:,4,1)';
y2_a = axis_data(:,4,2)';

y2_a_iec121 = y2_a;
y2_e_iec121 = y2_e;

%X1
x1_e = axis_data(:,5,1)';
x1_a = axis_data(:,5,2)';

x1_a_iec121 = arrayfun( @(x) - x, x1_a);
x1_e_iec121 = arrayfun( @(x) - x, x1_e);

%X2

x2_e = axis_data(:,6,1)';
x2_a = axis_data(:,6,2)';

x2_a_iec121 = x2_a;
x2_e_iec121 = x2_e;

%Couch Vrt
couchvrt_e = axis_data(:,7,1)';
couchvrt_a = axis_data(:,7,2)';

couchvrt_a_iec121 = arrayfun( @(x) -(x - 100), couchvrt_a);
couchvrt_e_iec121 = arrayfun( @(x) -(x - 100), couchvrt_e);

%Couch Lng
couchlng_e = axis_data(:,8,1)';
couchlng_a = axis_data(:,8,2)';

couchlng_a_iec121 = couchlng_a;
couchlng_e_iec121 = couchlng_e;

%Couch Lat
couchlat_e = axis_data(:,9,1)';
couchlat_a = axis_data(:,9,2)';

couchlat_a_iec121 = arrayfun( @(x) x - 100, couchlat_a);
couchlat_e_iec121 = arrayfun( @(x) x - 100, couchlat_e);

%Couch Rtn
couchrot_e = axis_data(:,10,1)';
couchrot_a = axis_data(:,10,2)';

couchrot_a_iec121 = arrayfun( @(x) mod(180 - x, 360), couchrot_a);
couchrot_e_iec121 = arrayfun( @(x) mod(180 - x, 360), couchrot_e);

%Couch Pit
%Unused! Couch does not support this function.
%Value stored is largest FP number represented in single precsn 32-bits.
couchpit_e = axis_data(:,11,1)';
couchpit_a = axis_data(:,11,2)';

%Couch Rol
%Unused! Couch does not support this function.
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
%USER PROMPT FOR PREPROGRAMMED WAVEFORM OR MW-VXP SELECTION.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
user_query_waveform = input('Press 1 for pre-programmed waveform option, ENTER key for MW-VXP: ');

if isempty(user_query_waveform)
    %Proceed with processing selected MW-VXP file.
    disp('User selected MW-VXP.')
    
    [file1, PATHNAME] = uigetfile(fullfile(pwd,'data_in','*.VXP;*MW*'),'select file');
    FULLPATH = fullfile(PATHNAME, file1);
    
    if strfind(file1, 'MW') > 0
        %File selected is MW file. Check isn't robust for all MW variations.
        [amplitude,phase,rpm_times,beamenable] = mw_reader(FULLPATH);
        
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
            else
                disp(exception.identifier)
                error('An unexpected error has occurred.')
            end
        end
    else
        %VXP file.
        [amplitude,phase,timestamp,validflag,ttlin,ttlout,mark,headerv] = vxp_reader(FULLPATH);
        rpm_times = cell2mat(timestamp) - timestamp{1,1};
        clearvars timestamp validflag ttlin ttlout
        phase     = cell2mat(phase);
        % amplitude = cell2mat(amplitude);
        % mark      = cell2mat(mark);
        
    end
else
    %Call for waveform generator function.
    period = input('Enter a period in seconds or hit ENTER key for default (4s): ');
    [rpm_times, phase] = waveform_generator(period);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SYNCHRONIZE MW AND TLF.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implemented process for synchronization. The TLF 'leads' the MW in
% time, so a simple shift applied to to the TLF can sync fairly well. The
% quality of the synchronization may be variable. In the future, a more
% sophisticated matching algorithm may be required.

if isempty(user_query_waveform)
    %The MW-VXP option has been selected.
    tlf_times_shifted = tlf_times + 4700; %yup, that's it.
else
    %Copy and rename so sorting functions work.
    %tlf_times will be used for plotting purposes from now on.
    tlf_times_shifted = tlf_times;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION CALLS FOR PHASE SORTING FOLLOWED BY ARC SORTING.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sorting the TLF file information into 10 different phases.
[sorted_phase, ~] = phase_sorter(tlf_times_shifted, rpm_times, phase);

%Sorting the 10 different phases into n (1, 2, or 3) arcs.
[sorted_phase_arc, intra_arc, arc_tlf_indicies] = arc_sorter(cp_e, subbeam, sorted_phase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PROCESSING FOR RP-PLAN CONSTRUCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The data retrieved from the TLF and sorted by the selected waveform needs
% to be further processed. Processed data needs to be used to modify 10
% pre-existing DICOM RP files so that 10 or 20 phased RP files are
% produced. These "phased-plans" can be imported into Eclipse and a
% sum-dose computed.

% Basically, the overall process is:
% - record TLF and waveform, sort TLF data using the waveform, and process
%   into a format following DICOM standards and accepted by Eclipse.
% - expos
% For each phased plan to be constructed, the MU at each CP must be
% specified in the form of "CumulativeMetersetWeight" (CMW) as found under
% "BeamSequence" -> "Item_x" (beam number) -> "ControlPointSequence" -> 
% "Item_xxx" (hundred+ CPs) -> "CumulativeMetersetWeight". CMU is equal to
% the cumulative MU delivered up to and including that CP in that field,
% divided by the total MU to be delivered in that field.

% The data needs to be reprocessed so that it has a suitable form for
% writing to the DICOM RP file. Determining the cumulative MU delivered for
% each field is the first step.



cumulative_mu = zeros(1, size(sorted_phase_arc, 2));

for i = 1:length(cumulative_mu)
    for j = 1:size(sorted_phase_arc, 1)
        mu_e_phase_arc = mu_e(sorted_phase_arc{j,i});
        cumulative_mu(i) = cumulative_mu(i) + mu_e_phase_arc(end);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WRITING OUT DICOM RT PLANS (10 or 20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Will this be a separate function? If we want it as a seperate function
% and many inputs are required, then we will probably have to pass the
% arguments as a structure or cell array.
%TBD TBD TBD


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTTING.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tlf_times = tlf_times/1000; %converting ms to s

% figure;
% plot(rpm_times,beamenable,'b')
% hold on
% plot(tlf_times,beamh_e,'r')
% ylim([-0.5,1.5])
% legend('RPM/MW Recording', 'TLF Recording')


f1 = figure();
figure(f1);
sz = 25;
scatter(tlf_times(sorted_phase_arc{1,1}), mu_e(sorted_phase_arc{1,1}),sz,'r','filled')
title('TLF Times vs. Total MU Delivered During Treatment')
xlabel('TLF Times [s]')
ylabel('Total MU Delivered')

hold on
scatter(tlf_times(sorted_phase_arc{2,1}), mu_e(sorted_phase_arc{2,1}),sz,'b','filled')
scatter(tlf_times(sorted_phase_arc{3,1}), mu_e(sorted_phase_arc{3,1}),sz,'c','filled')
scatter(tlf_times(sorted_phase_arc{5,1}), mu_e(sorted_phase_arc{5,1}),sz,'g','filled')
scatter(tlf_times(sorted_phase_arc{10,1}), mu_e(sorted_phase_arc{10,1}),sz,'m','filled')

scatter(tlf_times(sorted_phase_arc{1,2}), mu_e(sorted_phase_arc{1,2}),sz,'r','filled')
scatter(tlf_times(sorted_phase_arc{2,2}), mu_e(sorted_phase_arc{2,2}),sz,'b','filled')
scatter(tlf_times(sorted_phase_arc{3,2}), mu_e(sorted_phase_arc{3,2}),sz,'c','filled')
scatter(tlf_times(sorted_phase_arc{5,2}), mu_e(sorted_phase_arc{5,2}),sz,'g','filled')
scatter(tlf_times(sorted_phase_arc{10,2}), mu_e(sorted_phase_arc{10,2}),sz,'m','filled')


% plot(tlf_times(subbeam(1).arc(1):subbeam(1).arc(2)), cp_e(subbeam(1).arc(1):subbeam(1).arc(2)))
% hold on
% plot(tlf_times(subbeam(2).arc(1):subbeam(2).arc(2)), cp_e(subbeam(2).arc(1):subbeam(2).arc(2)),'-g')
% plot(tlf_times(intra_arc{1}(1):intra_arc{1}(2)), cp_a(intra_arc{1}(1):intra_arc{1}(2)),'-r')