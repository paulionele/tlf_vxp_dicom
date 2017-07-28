function [sorted_phase_arc, intra_arc] = arc_separator(cp_a, subbeam, sorted_phase)
%Function for seperating arcs. The function works by identifying the number
%of subbeams and using the control point (CP) information stored with each
%subbeam, the starting control point. Additionally, the total number of CPs
%are known. 

% INPUTS: cp_a and subbeam structure from TLF, sorted_phase from
% trajectory_log_phase_sort function (phase sorted TLF 'data'; a list of
% indicies sorted by phase).

% OUTPUTS: same structure as sorted_phase, except additional columns may be
% present if multiple subbeams exist. 

% PART 1 - Identifying the index range for each subbeam.
%The information is unusual in the sense that the starting CP for any beam
%that is NOT the first subbeam in the sequence, is given by the CP listed
%for that subbeam, plus one. EX: if subbeam 2 of 3 lists a CP of 113, then
%its actual starting CP is 113 + 1 = 114. Since no final CP is listed for
%any subbeam other than the final one (where final CP = total CP), then to
%find the final CP for any (except for last) subbeam, we look at the NEXT
%subbeam and use the CP listed for that one. EX: if subbeam 3 of 3 lists a
%CP of 113, then this is the final CP of subbeam 2 of 3.

% PART 2 - Sorting the phase_sorted indicies into n different arcs.
%The goal is now to seperate each of the 10 (or 20!) phases in n different
%arcs. I.e. if there is 2 arcs, then phase 0 will be seperated into arc 1
%and arc 2. The cell array aa contains the absolute index range for each
%subbeam. From a previous function, this current function receives a
%phase-sorted array of indicies (tlf indicies, number of indicies = number
%of 20ms snapshots); an array for each phase. We now pass through each of
%the 10 phases, and for a given phase, test membership of an index to the
%index range for each arc as determined in PART 1. Returned is a cell array
%sorted_phase_arc that maintains the sorted phase_phase structure, but now
%has n cols for the n arcs delivered.

number_subbeams  = length(subbeam); %1, 2, or 3 subbeams usually
aa = cell(number_subbeams, 1); %array to store index data for subbeams

number_intra_arc = length(subbeam) - 1;
intra_arc = cell(number_intra_arc, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART 1.

for i = 1:number_subbeams
    %Get time/index cooresp. to final CP in subbeam sequence.
    if i ~= number_subbeams
        %Not final subbeam; get initial CP of current subbeam.
        if i == 1
            %First beam of beam sequence. ex: [0,0,0,0.0002,...] -> 3 -> 60 ms
            %The beam starts with the final zero CP in the CP sequence.
            ai = find(cp_a == 0);
            ai = ai(end);
        else
            %Second or later (except very last) subbeam.
            %The starting subbeam CP is always the one listed for the beam
            %plus 1. EX: 113 listed but actual starting CP is 114.
            ai = find(cp_a == subbeam(i).cp + 1);
            ai = ai(end);
        end
        %Final CP of subbeam i; get from the next subbeam i+1.
        af = find(cp_a == subbeam(i+1).cp);
        af = af(1); %index of time of final CP in sb. seq. [<113>, 113, 113,...]
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART 2.

%Storing subbeam index information into pre-existing subbeam structure.
%New field (arc) added to each sub-structure.
%DO NOT NEED THIS CURRENTLY.
% subbeam(1).arc = [];
% for i = 1:number_subbeams
%     subbeam(i).arc = aa{i};
% end

%Set of indicies for each arc.
%Used for membership testing in the triple loop below.
arc_tlf_indicies  = cell(number_subbeams, 1);
for i = 1:number_subbeams  
    %Range of indicies belonging to each arc.
    arc_tlf_indicies{i} = aa{i}(1):aa{i}(2);
end

%Preallocating with NaN (easy to remove). Col.1 is ARC 1, Col. 2 is ARC 2.
sorted_phase_arc = cell(10, number_subbeams);
for i = 1:number_subbeams    
    for j = 1:size(sorted_phase_arc, 1)
        sorted_phase_arc{j,i} = NaN(1, length(sorted_phase{j}));
    end
end

%Sorting each phase into n arcs.
%This can be reduced from 3 loops (probably).
for n = 1:number_subbeams %1, 2, or 3
    for i = 1:length(sorted_phase) %(1:10)
        for j = 1:length(sorted_phase{i}) %(1:500 ish)
            
            %Pass thru sorted_phase and test membership in an arc.
            if ismember(sorted_phase{i}(j), arc_tlf_indicies{n})
                sorted_phase_arc{i,n}(j) = sorted_phase{i}(j);
            end
            
        end
    end
end

%Removing NaN entries.
for i = 1:number_subbeams
    for j = 1:size(sorted_phase_arc, 1)
        sorted_phase_arc{j,i}(isnan(sorted_phase_arc{j,i})) = [];
    end
end