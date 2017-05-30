function [sorted_phase_arc, intra_arc, arc1_tlf_indicies, arc2_tlf_indicies] = arc_separator(cp_a, subbeam, sorted_phase)
%Function for seperating arcs. Based on control points.
%Limited to 2 arc/subbeam treatments; in various areas.

number_subbeams  = length(subbeam);
aa = cell(number_subbeams, 1); %store indicies for each subbeam (2 subbeams for now)

number_intra_arc = length(subbeam) - 1;
intra_arc = cell(number_intra_arc, 1);

for i = 1:number_subbeams
    %Get time/index cooresp. to final CP in subbeam sequence.
    if i ~= number_subbeams
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
        if number_subbeams == 1
            %Only 1 subbeam.
            ai = find(cp_a == 0);
            ai = ai(end);
            
            af = max(cp_a); %max # CP in for entire delivery
            af = find(cp_a == af); %indicies for max # CP
            af = af(1); %index of final CP in final sb. seq. Note the MU tapers so best to take first instance.
        else
            %2 or more subbeams.
            ai = find(cp_a == subbeam(i).cp + 1);
            ai = ai(end);

            af = max(cp_a); %max # CP in for entire delivery
            af = find(cp_a == af); %indicies for max # CP
            af = af(1); %index of final CP in final sb. seq. Note the MU tapers so best to take first instance.
        end
    end
    aa{i} = [ai, af];
end

for i = 1:number_intra_arc
    %Will run once for 2 subbeams.
    intra_arc{i} = [aa{i}(2) + 1, aa{i+1}(1) - 1];
end

%Storing subbeam index information into pre-existing subbeam structure.
%New field (arc) added to each sub-structure.
% subbeam(1).arc = [];
% for i = 1:number_subbeams
%     subbeam(i).arc = aa{i};
% end

%%%Below this, code only works for 2 subbeams.
%Set of indicies for each arc.
arc1_tlf_indicies = aa{1}(1):aa{1}(2);
arc2_tlf_indicies = aa{2}(1):aa{2}(2);
sorted_phase_arc = cell(10,2);

%Preallocating. Col.1 is ARC 1, Col. 2 is ARC 2.
for i = 1:size(sorted_phase_arc, 1)
    sorted_phase_arc{i,1} = NaN(1, length(sorted_phase{i}));
    sorted_phase_arc{i,2} = NaN(1, length(sorted_phase{i}));
end

%Sorting each phase into 2 arcs.
for i = 1:length(sorted_phase) %(1:10)
    for j = 1:length(sorted_phase{i}) %(1:500 ish)
        if ismember(sorted_phase{i}(j), arc1_tlf_indicies)
            %Arc 1.
            sorted_phase_arc{i,1}(j) = sorted_phase{i}(j);
            
        elseif ismember(sorted_phase{i}(j), arc2_tlf_indicies)
            %Arc 2.
            sorted_phase_arc{i,2}(j) = sorted_phase{i}(j);
            
        else
            %Intra-arc.
        end
    end
end

%Removing NaN entries.
for i = 1:size(sorted_phase_arc, 1)
    sorted_phase_arc{i,1}(isnan(sorted_phase_arc{i,1})) = [];
    sorted_phase_arc{i,2}(isnan(sorted_phase_arc{i,2})) = [];
end

% x = {'xy','yz',nan}
% x(cellfun(@(x) any(isnan(x)),x)) = []

%Need to modify below and the objects returned by the function if plan
%consists of more than one arc.
% arc_1 = subbeam(1).arc;
% arc_2 = subbeam(2).arc;