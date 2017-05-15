function [sorted] = sinusoidal_trace(phases)
%Function for sorting TLF, based on simple pre-programmed respiratory
%trace. This function does not use any VXP or MW information.

%Can introduce phase offset.
%offset = pi/10; %[0,2*pi)

aa = linspace(0, 2*pi, 11); %10 phases equally spaced
lp = length(phases);

p0 = [];
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
    %[~,t_index] = min(abs(p-phases(i)));
    
    if (phases(i) >= aa(1)) && (phases(i) < aa(2))
        p0 = [p0,phases(i)];
    elseif (phases(i) >= aa(2)) && (phases(i) < aa(3))
        p10 = [p10,phases(i)];
    elseif (phases(i) >= aa(3)) && (phases(i) < aa(4)) 
        p20 = [p20,phases(i)];
    elseif (phases(i) >= aa(4)) && (phases(i) < aa(5))
        p30 = [p30,phases(i)];
    elseif (phases(i) >= aa(5)) && (phases(i) < aa(6))
        p40 = [p40,phases(i)];
    elseif (phases(i) >= aa(6)) && (phases(i) < aa(7))
        p50 = [p50,phases(i)];
    elseif (phases(i) >= aa(8)) && (phases(i) < aa(9))
        p60 = [p60,phases(i)];
    elseif (phases(i) >= aa(9)) && (phases(i) < aa(10))
        p70 = [p70,phases(i)];
    elseif (phases(i) >= aa(10)) && (phases(i) < aa(11))
        p80 = [p80,phases(i)];
    else
        p90 = [p90,phases(i)];
    end
end

sorted = {p0; p10; p20; p30; p40; p50; p60; p70; p80; p90};