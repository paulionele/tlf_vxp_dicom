function [ output_args ] = sinusoidal_trace( input_args )
%Function for sorting TLF, based on simple pre-programmed respiratory
%trace. This function does not use any VXP or MW information.

%Period and offset can be determined from the VXP file phase information.
%period = 3; %period in seconds
offset = 0; %[0,2*pi)

phases_u = linspace(0, 2*pi, 11); %phases unshifted
shift = (phases_u(2) - phases_u(1))/2; %half phase shift
phases_s = phases_u - shift;

phases = phases_s - offset;


end

