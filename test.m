clear;clc

phases = linspace(0,2*pi,5);
ph = [0,0.5,0.7,1.0,3.14,pi,6,6.14,2*pi];

lp = length(ph);

p0 = [];
p10 = [];
p20 = [];
p30 = [];
% p40 = [];
% p50 = [];
% p60 = [];
% p70 = [];
% p80 = [];
% p90 = [];

for i = 1:lp
    %[~,t_index] = min(abs(p-ph(i)));
    
    if (ph(i) >= phases(1)) && (ph(i) < phases(2))
        p0 = [p0,ph(i)];
    elseif (ph(i) >= phases(2)) && (ph(i) < phases(3))
        p10 = [p10,ph(i)];
    elseif (ph(i) >= phases(3)) && (ph(i) < phases(4))
        p20 = [p20,ph(i)];
    else
        p30 = [p30,ph(i)];
    end
end

sorted = {p0; p10; p20; p30};