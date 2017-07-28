function [rpm_times, phase] = waveform_generator(period)
%waveform_generator is a function that replaces the mw_reader and
%vxp_reader functions, for testing purposes.
%   A baseline performance assesment of the 4D VMAT calculation can be
%   performed using a preprogrammed repiratory waveform. No phase shift
%   applied; phase begins at zero phase.

period_predef = 4; %4 second respiratory period

if length(period) == 0
   period = period_predef;
   message = sprintf('No period specified. Using predefined period: %1.2f s.', period);
   disp(message)
else
   message = sprintf('Period: %1.2f seconds.', period);
   disp(message)
end

sampling_interval = 30; %in ms
sampling_duration = 200000; %1000 ms = 1 s
rpm_times = 0:sampling_interval:sampling_duration;

theta = 2*pi*rpm_times/(period*1000); %working in ms
phase = mod(theta - 2*pi, 2*pi); %returned phases on [0,2*pi)
end

