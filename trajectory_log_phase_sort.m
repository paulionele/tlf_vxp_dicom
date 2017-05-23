function [sorted_phase, phase_tlf2] = trajectory_log_phase_sort(tlf_times, vxp_times, phase)
%Function for sorting TLF, based on phase information from the VXP file.
%sorted_phase : phase-sorted indicies (indicies can be used to ref. axis)
%phase_tlf2 : 1xM array of phase values (1:1 coorelation with tlf_time)

%This function can be improved by preallocating matrix then removing empty
%entries.

%Interpolated phase values, based on TLF times.
phase_tlf = interp1(vxp_times, phase, tlf_times);
phase_tlf(isnan(phase_tlf)) = []; %removing NaN entries
length_tlf = length(phase_tlf); %truncating the time array

aa = linspace(0, 2*pi, 11); %10 phases equally spaced

sorted_phase = [];

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

%Prelim sorting based on continuous phase.
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

%Function by Stephen. Reassigns phase points in the discontinuities.
mode_window = 21;
phase_tlf2 = sorted_phase;
for jmode = ceil(mode_window/2):length_tlf-ceil(mode_window/2)
    phase_tlf2(jmode) = mode(sorted_phase(jmode-floor(mode_window/2):jmode+floor(mode_window/2)));
end

%Cell array structure for containing 10 phases and indicies that
%reference TLF points within a specific phase.
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
        p90 = [p90,tlf_times(i)];
    end
end

%Returned is indicies belonging to each phase.
sorted_phase = {p0; p10; p20; p30; p40; p50; p60; p70; p80; p90};

%Original function.
%
% mode_window=21
% Phase_tlf2=Phase_tlf;
% for jmode=ceil(mode_window/2):length_tlf-ceil(mode_window/2);
%     Phase_tlf2(jmode)=mode(Phase_tlf(jmode-floor(mode_window/2):jmode+floor(mode_window/2)));
% end