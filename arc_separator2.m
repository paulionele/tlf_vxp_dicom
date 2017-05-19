function [subbeam, intra_arc] = arc_separator2(cp_a, subbeam)
%Function for seperating arcs. Based on control points.

number_subbeams  = length(subbeam);
aa = cell(number_subbeams, 1); %store indicies for each subbeam

number_intra_arc = length(subbeam) - 1;
intra_arc = cell(number_intra_arc, 1);

for i = 1:number_subbeams
    %Get time/index cooresp. to final CP in subbeam sequence.
    if i ~= length(subbeam)
        %Not final subbeam; get initial CP.
        if i == 1
            %First beam of beam sequence. ex: [0,0,0,0.0002,...] -> 3 -> 60 ms
            ai = find(cp_a == 0);
            ai = ai(end);
        else
            %Second or later (except very last) subbeam.
            %The starting subbeam CP is always the one listed for the beam
            %plus 1. Refer to OneNote comments.
            ai = find(cp_a == subbeam(i).cp + 1);
            ai = ai(end);
        end
        %Final CP of beam i; get from beam i+1
        af = find(cp_a == subbeam(i+1).cp);
        af = af(end); %index of time of final CP in sb. seq.
    else
        %Final subbeam.
        ai = find(cp_a == subbeam(i).cp + 1);
        ai = ai(end);
        
        af = max(cp_a); %max # CP in for entire delivery
        af = find(cp_a == af); %indicies for max # CP
        af = af(1); %index of final CP in final sb. seq. Note the MU tapers so best to take first instance.
    end
    aa{i} = [ai, af];        
end

for i = 1:number_intra_arc
    intra_arc{i} = [aa{i}(2) + 1, aa{i+1}(1) - 1];
end

%Storing subbeam index information into pre-existing subbeam structure.
%New field (arc) added to each sub-structure.
%I wonder if this may cause issues...
subbeam(1).arc = [];
for i = 1:length(subbeam)
    subbeam(i).arc = aa{i};
end

%Need to modify below and the objects returned by the function if plan
%consists of more than one arc.
% arc_1 = subbeam(1).arc;
% arc_2 = subbeam(2).arc;