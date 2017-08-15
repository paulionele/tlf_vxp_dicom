function [sorted_phase] = phase_sorter(tlf_times, rpm_times, phase, num_phases)
%Function for phase sorting TLF data into 10 (or 20 phases, code would need
%to be re-written for that though) phases, based on phase information from
%the VXP or MW file.

% INPUTS: tlf_times (20 ms sampling interval), rpm_times (~30 ms sampling
% interval), phase (sampling interval and phase specific to MW, VXP, or
% preprogrammed waveform, depending on user selection), and the number of
% phases.

% OUTPUTS: sorted_phase (phase-sorted indicies where indicies can be used
% to reference TLF recordings).

% PART 1 - The phase values arrive with a sampling interval different than
% the sampling interval of the TLF and since a phase value must be mapped
% to each TLF recording, it is necessary to calculate phase values for TLF
% recordings by interpolation. Once a phase value on [0, 2*pi) has been
% associated with each TLF recording, the floating point phase values on
% [0,2*pi] must be converted to an interger value on [0,9] (or [0,n-1],
% where n is the number of phases). Because of the different sampling
% rates, some TLF recordings have a sample that 'wraps' from 2*pi to 0.
% These cannot exist and in PART 2, are reassigned to either 2*pi or 0.

% PART 2 - The values that 'wrap' from 2*pi to 0 are removed by a function
% by Steven Thomas (BCCA). After that, a for-loop and if-else construct
% sort the 1D array of phase values [0,n-1] into n, 1D arrays for each
% phase value. In addition, if the next i+1 TLF recording is associated
% with a different phase, then that next TLF recording is added to the
% current phase and also added to the next cooresponding phase. Basically,
% there is a double addition of a recording every time there is a phase
% change.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1.

%Interpolated phase values, based on TLF times.
phase_tlf = interp1(rpm_times, phase, tlf_times);
phase_tlf(isnan(phase_tlf)) = []; %removing NaN entries
length_tlf = length(phase_tlf); %number of snapshots/recordings

%Defining n equally spaced phases on 0 to 2*pi.
aa = linspace(0, 2*pi, num_phases + 1);

%A preallocated array storing phase values for each index (tlf_time).
phase_array = zeros(1,length_tlf);

%The for-loop below assigns each continuous phase value on [0,2*pi) an
%integer on [0,n-1]. Logical indexing is used below. Testing is performed
%by relational operators, then the logical arrays returned are compared by
%logical AND. Note, the short AND (&&) does not work for whatever reason
%because "&& operator must be convertible to scalar values". In any case,
%logical indexing is used to reference the entire phase_array and set the
%referenced values to the appropriate values.
for i = 1:num_phases
    phase_array( (phase_tlf >= aa(i)) & (phase_tlf < aa(i+1)) ) = i-1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2.

%Due to different sampling rates between the TLF and MW/VXP file, the
%interp1 function returns some phase values that 'wrap' from 2*pi to 0;
%that is, these values exist between 2*pi and 0, when they cannot. A
%function was devised by Steven Thomas to reassign such phase values to
%either the previous phase value (2*pi) or next phase value (0).

%Function by Steven Thomas (BCCA).
%Reassigns phase points in the discontinuities.
mode_window = 21;
phase_tlf2 = phase_array;
for jmode = ceil(mode_window/2):length_tlf-ceil(mode_window/2)
    phase_tlf2(jmode) = mode(phase_array(jmode-floor(mode_window/2):jmode+floor(mode_window/2)));
end
% End function by Steven Thomas.

%Below, the phase_tlf2 array (a processed version of phase_array) is sorted
%into n different arrays, cooresponding to n different phases.
x = ones(1,num_phases); %used to maintain indicies to ref sorted_phase subarrays
sorted_phase = cell(num_phases,1);
for i = 1:num_phases
    sorted_phase{i} = NaN(1,length(phase_array));
end

i = 1; %increment for while loop
entries_added = 0;
while i < length_tlf
    %Loop through entire TLF, expect for last (ref outside array dim).
    for j = 1:num_phases
        %Loop through numerical representation of phase (ex: 0 to 9).
        if phase_tlf2(i) == (j-1)
            %If match to current recorded phase.
            
            if (phase_tlf2(i+1) == (j-2) || phase_tlf2(i+1) == (j))
                %Special case here. Add current recording to current phase
                %AND add next recording to current and next phase. New rec
                %could be either a 'higher' or 'lower' phase hence we have
                %j-2 for one less than j-1 and j as one more.
                
                %Add the current recording.
                sorted_phase{j}(x(j)) = i;
                x(j) = x(j) + 1;
                
                %Add the next recording to the current phase.
                sorted_phase{j}(x(j)) = i+1;
                x(j) = x(j) + 1;
                
                %Add the next recording to the next phase.
                %phase_tlf2 returns 0 to 9 but need to increment
                %for index use.
                q = phase_tlf2(i+1) + 1; %easy to see index
                sorted_phase{q}(x(q)) = i+1;
                x(q) = x(q) + 1;
                
                %Mark to skip next recording. This increment in combination
                %with the increment at the end causes a skip.
                i = i+1;
                entries_added = entries_added + 1;
                break
            else
                %If test returns true, then add the index
                %cooresponding to the phase entry that tested true,
                %to the array cooresponding to that phase.
                sorted_phase{j}(x(j)) = i;
                x(j) = x(j) + 1;
                break
            end %special test
        end %phase test
    end %phase loop
    i = i+1;
end %TLF loop

for i = 1:num_phases
    %Removing NaN array entries.
    sorted_phase{i}(isnan(sorted_phase{i})) = [];
end

%Original function by S. Thomas for reference.
%
% mode_window=21
% Phase_tlf2=Phase_tlf;
% for jmode=ceil(mode_window/2):length_tlf-ceil(mode_window/2);
%     Phase_tlf2(jmode)=mode(Phase_tlf(jmode-floor(mode_window/2):jmode+floor(mode_window/2)));
% end