function [ output_args ] = sinusoidal_trace(phase)
%Function for sorting TLF, based on simple pre-programmed respiratory
%trace. This function does not use any VXP or MW information.

%Period and offset can be determined from the VXP file phase information.
%period = 3; %period in seconds
offset = 0; %[0,2*pi)

p = linspace(0, 2*pi, 11); %setting evenly distributed phases

for i = 1:length(phase)
    [~,t_index] = min(abs(p - phase(i)));
    

end