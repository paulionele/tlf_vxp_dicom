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
times = 0:20:(20*header.num_snaps - 1);

%Collimator Rotation
collrot_e = axis_data(:,1,1)';
collrot_a = axis_data(:,1,2)';
%Gantry Rotation
gantrot_e = axis_data(:,2,1)';
gantrot_a = axis_data(:,2,2)';
%Y1
y1_e = axis_data(:,3,1)';
y1_a = axis_data(:,3,2)';
%Y2
y2_e = axis_data(:,4,1)';
y2_a = axis_data(:,4,2)';
%X1
x1_e = axis_data(:,5,1)';
x1_a = axis_data(:,5,2)';
%X2
x2_e = axis_data(:,6,1)';
x2_a = axis_data(:,6,2)';
%Couch Vrt
couchvrt_e = axis_data(:,7,1)';
couchvrt_a = axis_data(:,7,2)';
%Couch Lng
couchlng_e = axis_data(:,8,1)';
couchlng_a = axis_data(:,8,2)';
%Couch Lat
couchlat_e = axis_data(:,9,1)';
couchlat_a = axis_data(:,9,2)';
%Couch Rtn
couchrot_e = axis_data(:,10,1)';
couchrot_a = axis_data(:,10,2)';
%Couch Pit
couchpit_e = axis_data(:,11,1)';
couchpit_a = axis_data(:,11,2)';
%Couch Rol
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

%Getting all leaf positions...
% for i = 18:137
%     mlc_leaf_e = axis_data(:,i,1);
%     mlc_leaf_a = axis_data(:,i,2);
% end

%Getting VXP file information.
file1 = uigetfile('.VXP');
[amplitude,phase,timestamp,validflag,ttlin,ttlout,mark,header] = vxp_reader(file1);
timestamp = cell2mat(timestamp) - timestamp{1,1};
phase     = cell2mat(phase);
% amplitude = cell2mat(amplitude);
% mark      = cell2mat(mark);

%Calling function to sort timestamps based on phases from VXP...
sorted = sinusoidal_trace(timestamp, trace);
