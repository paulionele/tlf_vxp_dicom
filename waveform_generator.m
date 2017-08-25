function [rpm_times, phase] = waveform_generator(period, NUMPHASES)
%This function can be called instead of the mw_reader and vxp_reader
%functions, for testing purposes. The user specifies a period in seconds,
%and the function returns an array of phase values between 0 and 2*pi as a
%function of rpm_times. A baseline performance assesment of the 4D VMAT
%calculation can be performed using a preprogrammed repiratory waveform. No
%phase shift applied; phase begins at zero phase.

%INPUT: user specified period in seconds.
%OUTPUT: rpm_times with pre-defined time interval and duration. phase array
%with discrete values between 0 and 2*pi, with a 1:1 mapping to rpm_times.

period_predef     = 4.0; %4 second respiratory period
sampling_interval = 20.0; %in ms

%Period determination.
if length(period) == 0
   period = period_predef;
   message = sprintf('No period specified. Using predefined period: %1.2f s.', period);
   disp(message)
else
   period = double(period);
   message = sprintf('Period: %1.2f seconds.', period);
   disp(message)
end

%Preparation.
period_repetition = 200; %number of respiratory cycles/periods
period_samples    = period/(sampling_interval/1000); %samples per period (respiratory cycle)
phases_samples    = period_samples/NUMPHASES; %samples per phase
rpm_times         = 0:sampling_interval:(period_repetition*period_samples*sampling_interval);

%Produces discrete wrapped-phase.
i = 1; %scratch var for preallocation
pp = 1;
phase = zeros(1,period_repetition*period_samples);
while pp < period_repetition+1
    %This outer while loop repeats a single built phase.
    
    %These inner two loops build the phase values for a single period.
    for j = 0:(NUMPHASES-1)
        %For each phase value...
        for k = 1:phases_samples
            %...that phase value is repeated.
            phase(i) = j;
            i = i+1;
        end
    end
    pp = pp+1;
end


%Produces continuous wrapped-phase.
% theta = 2*pi*rpm_times/(period*1000); %working in ms
% phase = mod(theta - 2*pi, 2*pi); %returned phases on [0,2*pi)
end