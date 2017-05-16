clear;clc
load('testvars')

%Returns interpolated Phase_tlf2-values for each TLF time.
Phase_tlf = interp1(vxp_times, phase, tlf_times);
Phase_tlf(isnan(Phase_tlf)) = []; %removing NaN entries

length_tlf = length(Phase_tlf); %truncating the time array
tlf_times = tlf_times(1:length_tlf);

%Function by Stephen. It does something.
mode_window = 21;
Phase_tlf2 = Phase_tlf;
for jmode = ceil(mode_window/2):(length_tlf - ceil(mode_window/2))
    Phase_tlf2(jmode) = mode( Phase_tlf( (jmode - floor(mode_window/2)) : (jmode + floor(mode_window/2)) ) );
end

%Can introduce phase offset.
%offset = pi/10; %[0,2*pi)

aa = linspace(0, 2*pi, 11); %10 phases equally spaced

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

for i = 1:length_tlf
    %[~,t_index] = min(abs(p-Phase_tlf2(i)));
    
    if (Phase_tlf2(i) >= aa(1)) && (Phase_tlf2(i) < aa(2))
        p0 = [p0,i];
        
    elseif (Phase_tlf2(i) >= aa(2)) && (Phase_tlf2(i) < aa(3))
        p10 = [p10,i];
        
    elseif (Phase_tlf2(i) >= aa(3)) && (Phase_tlf2(i) < aa(4)) 
        p20 = [p20,i];
        
    elseif (Phase_tlf2(i) >= aa(4)) && (Phase_tlf2(i) < aa(5))
        p30 = [p30,i];
        
    elseif (Phase_tlf2(i) >= aa(5)) && (Phase_tlf2(i) < aa(6))
        p40 = [p40,i];
        
    elseif (Phase_tlf2(i) >= aa(6)) && (Phase_tlf2(i) < aa(7))
        p50 = [p50,i];
        
    elseif (Phase_tlf2(i) >= aa(8)) && (Phase_tlf2(i) < aa(9))
        p60 = [p60,i];
        
    elseif (Phase_tlf2(i) >= aa(9)) && (Phase_tlf2(i) < aa(10))
        p70 = [p70,i];
        
    elseif (Phase_tlf2(i) >= aa(10)) && (Phase_tlf2(i) < aa(11))
        p80 = [p80,i];
        
    else
        p90 = [p90,tlf_times(i)];
    end
end

sorted = {p0; p10; p20; p30; p40; p50; p60; p70; p80; p90};

