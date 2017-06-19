clear;clc

function [sorted_phase, phase_tlf2] = trajectory_log_phase_sort(tlf_times, rpm_times, phase)
%Function for phase sorting TLF data into 10 phases, based on phase
%information from the VXP or MW file.

% INPUTS: tlf_times (20 ms sampling interval), rpm_times and phase
% (sampling interval and phase specific to MW or VXP, depending on user
% selection).

% OUTPUTS: sorted_phase (phase-sorted indicies where indicies can be used
% to ref. axis), phase_tlf2 : 1xM array of phase values (1:1 coorelation
% with tlf_time), not currently used!

% Complete description of PART 1 and PART 2.

%This function can be improved by preallocating matrix then removing empty
%entries.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1.

%Interpolated phase values, based on TLF times.
phase_tlf = interp1(rpm_times, phase, tlf_times);
phase_tlf(isnan(phase_tlf)) = []; %removing NaN entries
length_tlf = length(phase_tlf); %truncating the time array

%Defining 10 equall spaced phases.
aa = linspace(0, 2*pi, 11);

%Array storing phase values for each index (tlf_time).
sorted_phase = [];

%Individual arrays for storing indices belonging to each phase.
p0  = [];
p10 = [];
p20 = [];
p30 = [];
p40 = [];
p50 = [];
p60 = [];
p70 = [];
p80 = [];
p90 = [];

%Preliminary sorting based on continuous phase. Membership testing is
%performed for each phase value in phase_tlf; i.e. if phase_tlf(i) has a
%value 0.1, then it is assigned a phase value "0" and this value is added
%to the sorted_phase array.
for i = 1:length_tlf
    if (phase_tlf(i) >= aa(1)) && (phase_tlf(i) < aa(2))
        sorted_phase = [sorted_phase, 0];
    elseif (phase_tlf(i) >= aa(2)) && (phase_tlf(i) < aa(3))
        sorted_phase = [sorted_phase, 1];
    elseif (phase_tlf(i) >= aa(3)) && (phase_tlf(i) < aa(4))
        sorted_phase = [sorted_phase, 2];
    elseif (phase_tlf(i) >= aa(4)) && (phase_tlf(i) < aa(5))
        sorted_phase = [sorted_phase, 3];
    elseif (phase_tlf(i) >= aa(5)) && (phase_tlf(i) < aa(6))
        sorted_phase = [sorted_phase, 4];
    elseif (phase_tlf(i) >= aa(6)) && (phase_tlf(i) < aa(7))
        sorted_phase = [sorted_phase, 5];
    elseif (phase_tlf(i) >= aa(7)) && (phase_tlf(i) < aa(8))
        sorted_phase = [sorted_phase, 6];
    elseif (phase_tlf(i) >= aa(8)) && (phase_tlf(i) < aa(9))
        sorted_phase = [sorted_phase, 7];
    elseif (phase_tlf(i) >= aa(9)) && (phase_tlf(i) < aa(10))
        sorted_phase = [sorted_phase, 8];
    else
        sorted_phase = [sorted_phase, 9];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2.

%Due to different sampling rates between the TLF and MW/VXP file, the
%interp1 function returns some phase values that "wrap around" from 2*pi to
%0; that is, these values exist between 2*pi and 0, when they cannot. A
%function was devised to reassign such phase values to either the previous
%phase value (2*pi) or next phase value (0).

%%% Function by Steven Thomas (BCCA).
% Reassigns phase points in the discontinuities.
mode_window = 21;
phase_tlf2 = sorted_phase;
for jmode = ceil(mode_window/2):length_tlf-ceil(mode_window/2)
    phase_tlf2(jmode) = mode(sorted_phase(jmode-floor(mode_window/2):jmode+floor(mode_window/2)));
end
%%% End function by Steven Thomas.

%Cell array structure for containing 10 phases and indicies that reference
%TLF data within a specific phase.
for i = 1:length_tlf
    if (phase_tlf2(i) == 0)
        p0 = [p0,i];
    elseif (phase_tlf2(i) == 1)
        p10 = [p10,i];
    elseif (phase_tlf2(i) == 2)
        p20 = [p20,i];
    elseif (phase_tlf2(i) == 3)
        p30 = [p30,i];
    elseif (phase_tlf2(i) == 4)
        p40 = [p40,i];
    elseif (phase_tlf2(i) == 5)
        p50 = [p50,i];
    elseif (phase_tlf2(i) == 6)
        p60 = [p60,i];
    elseif (phase_tlf2(i) == 7)
        p70 = [p70,i];
    elseif (phase_tlf2(i) == 8)
        p80 = [p80,i];
    else
        p90 = [p90,i];
    end
end

%Returned is indicies belonging to each phase.
sorted_phase = {p0; p10; p20; p30; p40; p50; p60; p70; p80; p90};

%Original function by S. Thomas.
%
% mode_window=21
% Phase_tlf2=Phase_tlf;
% for jmode=ceil(mode_window/2):length_tlf-ceil(mode_window/2);
%     Phase_tlf2(jmode)=mode(Phase_tlf(jmode-floor(mode_window/2):jmode+floor(mode_window/2)));
% end