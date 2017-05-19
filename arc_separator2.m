function [arc_1, arc_2, intra_arc] = arc_separator2(cp_a, subbeam)
%Function for seperating arcs. Based on control points and mu???
%   Detailed explanation goes here


length_time     = length(cp_a);
number_subbeams = length(subbeam);

cps = zeros(1,number_subbeams); %starting CP for each beam
cpt = max(cp_a); %maximum CP of entire delivery (ex:227)

aa = cell(number_subbeams, 1); %store info for each subbeam

for i = 1:number_subbeams
    %Determine starting subbeam control points.
    cps(i) = subbeam(i).cp; %(ex:0, 113)
    
    %Get time cooresp. to final CP in subbeam sequence.
    if i ~= length(subbeam)
        %Not final subbeam.
        a1 = find(cp_a == subbeam(i+1).cp);
        a1 = a1(end); %index of time of final CP in sb. seq.
    else
        %Final subbeam.
        a1 = find(cp_a == subbeam(i).cp + 1);
        a1 = a1(end); %index of final CP in final sb. seq.
    end
    
end



for i = 1:length_time

end


%Storing subbeam information.
subbeam(1).arc = [];
for i = 1:length(subbeam)
    subbeam(i).arc = aa{i};
end

%Need to modify below and the objects returned by the function if plan
%consists of more than one arc.
arc_1 = subbeam(1).arc;
arc_2 = subbeam(2).arc;