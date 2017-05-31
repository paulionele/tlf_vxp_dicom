%  Copyright (C) 2017 Paul Ionele - All Rights Reserved
%  You may NOT use, distribute and modify this code unless express
%  written and signed consent has been given by the author Paul Ionele.

clear;clc
%TLF reader.
%Integers and floats are stored in 'little endian' format.

file1 = uigetfile('.bin');
fid1 = fopen(file1,'rb');
%ens = 'ieee-le'; %order for numerical interpretation of byte sequence

%%%Header.
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
axis_data = zeros(header.num_snaps, total_samples, 2); %2 -> actual and measured

for i = 1:header.num_snaps
    for j = 1:total_samples
        for k = 1:2
            axis_data(i,j,k) = fread(fid1, 1, '*float');
        end
    end
end

fclose(fid1);

%First snap is at 0 ms (?), and each preceeding snap is at 20*n ms; where n
%is the snap number. So an array of times:
tlf_times = single(0:20:(20*header.num_snaps - 1)); %do we lose a sample time off the end???

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
beamh_e = axis_data(:,14,1)';
beamh_a = axis_data(:,14,2)';

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
%%% (Option 1) Getting VXP file information.
file1 = uigetfile('.VXP');
[amplitude,phase,timestamp,validflag,ttlin,ttlout,mark,headerv] = vxp_reader(file1);
vxp_times = cell2mat(timestamp) - timestamp{1,1}; 
clearvars timestamp validflag ttlin ttlout
phase     = cell2mat(phase);
% amplitude = cell2mat(amplitude);
% mark      = cell2mat(mark);

%%%(Option 2) Getting MW file information.
% file1 = uigetfile('*MW*');
% [aa] = trajectory_log_reader(file1);

%%%Sorting the TLF file information.
[sorted_phase, phase_tlf2] = trajectory_log_phase_sort(tlf_times, vxp_times, phase);

%%%Excising index ranges for seperate arcs.
[sorted_phase_arc, intra_arc] = arc_separator(cp_a, subbeam, sorted_phase);

%%%Plotting
sz = 25;
scatter(tlf_times(sorted_phase_arc{1,1}), mu_a(sorted_phase_arc{1,1}),sz,'r','filled')
hold on
scatter(tlf_times(sorted_phase_arc{2,1}), mu_a(sorted_phase_arc{2,1}),sz,'b','filled')
scatter(tlf_times(sorted_phase_arc{3,1}), mu_a(sorted_phase_arc{3,1}),sz,'c','filled')
scatter(tlf_times(sorted_phase_arc{5,1}), mu_a(sorted_phase_arc{5,1}),sz,'g','filled')
scatter(tlf_times(sorted_phase_arc{10,1}), mu_a(sorted_phase_arc{10,1}),sz,'m','filled')

scatter(tlf_times(sorted_phase_arc{1,2}), mu_a(sorted_phase_arc{1,2}),sz,'r','filled')
scatter(tlf_times(sorted_phase_arc{2,2}), mu_a(sorted_phase_arc{2,2}),sz,'b','filled')
scatter(tlf_times(sorted_phase_arc{3,2}), mu_a(sorted_phase_arc{3,2}),sz,'c','filled')
scatter(tlf_times(sorted_phase_arc{5,2}), mu_a(sorted_phase_arc{5,2}),sz,'g','filled')
scatter(tlf_times(sorted_phase_arc{10,2}), mu_a(sorted_phase_arc{10,2}),sz,'m','filled')


% plot(tlf_times(subbeam(1).arc(1):subbeam(1).arc(2)), cp_a(subbeam(1).arc(1):subbeam(1).arc(2)))
% hold on
% plot(tlf_times(subbeam(2).arc(1):subbeam(2).arc(2)), cp_a(subbeam(2).arc(1):subbeam(2).arc(2)),'-g')
% plot(tlf_times(intra_arc{1}(1):intra_arc{1}(2)), cp_a(intra_arc{1}(1):intra_arc{1}(2)),'-r')