function [arc_1, arc_2, intra_arc] = arc_separator(beamh_a, subbeam)
%Function for seperating arcs. Based on beam hold.
%   Detailed explanation goes here

aa = cell(length(subbeams), 1);
lt = length(beamh_a);

r = 1; %beam counter
q = 0; %

for i = 1:lt
    if beamh_a(i) == 0
        %Add index to beam(r).
        aa{r} = [aa{r}, i];
        q = 0;
    elseif beamh_a(i) == 2 && (i < 10 || i > (lt - 10))
        %Add to beam(r).
        %Beam hold activated, but at start or end of treatment.
        aa{r} = [aa{r}, i];
        q = 0;
    else
        %Beam hold activated, do not add to any beam.
        %Increment beam counter.
        if q == 0
            r = r + 1;
            q = 1;
        else
            continue
        end
        
        %Recording intrarc indicies.
        intra_arc = [intra_arc, i];
    end
end

subbeam(1).arc = [];

for i = 1:length(subbeam)
    subbeam(i).arc = aa{i};
end

%Need to modify below and the objects returned by the function if plan
%consists of more than one arc.
arc_1 = subbeam(1).arc;
arc_2 = subbeam(2).arc;