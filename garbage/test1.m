clear;clc

load('www')

%Function for seperating arcs. Based on control points.

number_subbeams  = length(subbeam);
aa = cell(number_subbeams, 1);

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
        %sorted_phase_arc{i,2}(isnan(sorted_phase_arc{i,2})) = [];
    end
end

%%%Plotting
sz = 25;
scatter(tlf_times(sorted_phase_arc{1,1}), mu_a(sorted_phase_arc{1,1}),sz,'r','filled')
hold on
scatter(tlf_times(sorted_phase_arc{2,1}), mu_a(sorted_phase_arc{2,1}),sz,'b','filled')
scatter(tlf_times(sorted_phase_arc{3,1}), mu_a(sorted_phase_arc{3,1}),sz,'c','filled')
scatter(tlf_times(sorted_phase_arc{5,1}), mu_a(sorted_phase_arc{5,1}),sz,'g','filled')
scatter(tlf_times(sorted_phase_arc{10,1}), mu_a(sorted_phase_arc{10,1}),sz,'m','filled')

scatter(tlf_times(sorted_phase_arc{1,2}), mu_a(sorted_phase_arc{1,2}),sz,'r','filled')
scatter(tlf_times(sorted_phase_arc{2,2}), mu_a(sorted_phase_arc{2,2}),sz,'b','filled')
scatter(tlf_times(sorted_phase_arc{3,2}), mu_a(sorted_phase_arc{3,2}),sz,'c','filled')
scatter(tlf_times(sorted_phase_arc{5,2}), mu_a(sorted_phase_arc{5,2}),sz,'g','filled')
scatter(tlf_times(sorted_phase_arc{10,2}), mu_a(sorted_phase_arc{10,2}),sz,'m','filled')


% plot(tlf_times(subbeam(1).arc(1):subbeam(1).arc(2)), cp_a(subbeam(1).arc(1):subbeam(1).arc(2)))
% hold on
% plot(tlf_times(subbeam(2).arc(1):subbeam(2).arc(2)), cp_a(subbeam(2).arc(1):subbeam(2).arc(2)),'-g')
% plot(tlf_times(intra_arc{1}(1):intra_arc{1}(2)), cp_a(intra_arc{1}(1):intra_arc{1}(2)),'-r')