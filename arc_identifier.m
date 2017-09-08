function [intra_arc] = arc_identifier(cp_e, subbeam)
%Function for identifying arcs, based on the arc_sorter function. PART 1
%from that function has been stripped and placed into its own separate
%function.

% INPUTS: cp_e and subbeam structure from TLF

% OUTPUTS: intra_arc, index range BETWEEN the end of the first beam and
% start of second beam. Add or remove to these indicies if the actual end
% or start of the subbeams is desired.

% PART 1 - Identifying the index range for each subbeam.
%The information is unusual in the sense that the starting CP for any beam
%that is NOT the first subbeam in the sequence, is given by the CP listed
%for that subbeam, plus one. EX: if subbeam 2 of 3 lists a CP of 113, then
%its actual starting CP is 113 + 1 = 114. Since no final CP is listed for
%any subbeam other than the final one (where final CP = total CP), then to
%find the final CP for any (except for last) subbeam, we look at the NEXT
%subbeam and use the CP listed for that one. EX: if subbeam 3 of 3 lists a
%CP of 113, then this is the final CP of subbeam 2 of 3.

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
            ai = find(cp_e == 0);
            ai = ai(end);
        else
            %Second or later (except very last) subbeam.
            %The starting subbeam CP is always the one listed for the beam
            %plus 1. EX: 113 listed but actual starting CP is 114.
            ai = find(cp_e == subbeam(i).cp + 1);
            ai = ai(end);
        end
        %Final CP of subbeam i; get from the next subbeam i+1.
        af = find(cp_e == subbeam(i+1).cp);
        af = af(1); %index of time of final CP in sb. seq. [<113>, 113, 113,...]
        
    else
        %Final subbeam.
        if number_subbeams == 1
            %Only 1 subbeam.
            ai = find(cp_e == 0);
            ai = ai(end);
            
            af = max(cp_e); %max # CP in for entire delivery
            af = find(cp_e == af); %indicies for max # CP
            af = af(1); %index of final CP in final sb. seq. Note the MU tapers so best to take first instance.
        else
            %2 or more subbeams.
            ai = find(cp_e == subbeam(i).cp + 1);
            ai = ai(end);
            
            af = max(cp_e); %max # CP in for entire delivery
            af = find(cp_e == af); %indicies for max # CP
            af = af(1); %index of final CP in final sb. seq. Note the MU tapers so best to take first instance.
        end
        
    end %subbeam number test
    aa{i} = [ai, af];
end %subbeam loop

for i = 1:number_intra_arc
    %Will run once for 2 subbeams.
    intra_arc{i} = [aa{i}(2) + 1, aa{i+1}(1) - 1];
end

end %function