function [rpm_times, phase] = waveform_generator(period)
%This function can be called instead of the mw_reader and vxp_reader
%functions, for testing purposes. The user specifies a
%period in seconds, and the function returns an array of phase values
%between 0 and 2*pi as a function of rpm_times.

%INPUT: user specified period in seconds.
%OUTPUT: rpm_times with pre-defined time interval and duration. phase array
%with discrete values between 0 and 2*pi, with a 1:1 mapping to rpm_times.

%   A baseline performance assesment of the 4D VMAT calculation can be
%   performed using a preprogrammed repiratory waveform. No phase shift
%   applied; phase begins at zero phase.

period_predef = 4; %4 second respiratory period
sampling_interval = 30; %in ms
sampling_duration = 200000; %1000 ms = 1 s
rpm_times = 0:sampling_interval:sampling_duration;

if length(period) == 0
   period = period_predef;
   message = sprintf('No period specified. Using predefined period: %1.2f s.', period);
   disp(message)
else
   message = sprintf('Period: %1.2f seconds.', period);
   disp(message)
end

theta = 2*pi*rpm_times/(period*1000); %working in ms
phase = mod(theta - 2*pi, 2*pi); %returned phases on [0,2*pi)
end
