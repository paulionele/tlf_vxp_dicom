function [sorted] = sinusoidal_trace(timestamp, phase)
%Function for sorting TLF, based on phase information from the VXP file.
%Note that the cooresponding times are also sorted with the phase; it is
%the timestamps that are returned.

%Can introduce phase offset.
%offset = pi/10; %[0,2*pi)

aa = linspace(0, 2*pi, 11); %10 phases equally spaced
lp = length(phase);

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

for i = 1:lp
    %[~,t_index] = min(abs(p-phase(i)));
    
    if (phase(i) >= aa(1)) && (phase(i) < aa(2))
        p0 = [p0,timestamp(i)];
        
    elseif (phase(i) >= aa(2)) && (phase(i) < aa(3))
        p10 = [p10,timestamp(i)];
        
    elseif (phase(i) >= aa(3)) && (phase(i) < aa(4)) 
        p20 = [p20,timestamp(i)];
        
    elseif (phase(i) >= aa(4)) && (phase(i) < aa(5))
        p30 = [p30,timestamp(i)];
        
    elseif (phase(i) >= aa(5)) && (phase(i) < aa(6))
        p40 = [p40,timestamp(i)];
        
    elseif (phase(i) >= aa(6)) && (phase(i) < aa(7))
        p50 = [p50,timestamp(i)];
        
    elseif (phase(i) >= aa(8)) && (phase(i) < aa(9))
        p60 = [p60,timestamp(i)];
        
    elseif (phase(i) >= aa(9)) && (phase(i) < aa(10))
        p70 = [p70,timestamp(i)];
        
    elseif (phase(i) >= aa(10)) && (phase(i) < aa(11))
        p80 = [p80,timestamp(i)];
        
    else
        p90 = [p90,timestamp(i)];
    end
end

sorted = {p0; p10; p20; p30; p40; p50; p60; p70; p80; p90};
% 
% mode_window=21
% Phase_tlf2=Phase_tlf;
% for jmode=ceil(mode_window/2):length_tlf-ceil(mode_window/2);
%     Phase_tlf2(jmode)=mode(Phase_tlf(jmode-floor(mode_window/2):jmode+floor(mode_window/2)));
% end
