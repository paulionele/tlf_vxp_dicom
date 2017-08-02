%Data is directly from TLF !!!!!
%Data is 2 D matrix of all the TLF data (20K lines and ~270 col)


%Looks like this section constructs the Phase structure by referencing the
%data read from the tlf (Data).

%%%%This part sets this time point as both the end of the last phase and
%%%%the beginning of the next phase.
    Phase{phase_number}.ARC_num(counter(phase_number))=ARC_num(i); % Steven Added
    Phase{phase_number}.Col(counter(phase_number))=Data(i,1); % Steven Added
    Phase{phase_number}.Control_Point(counter(phase_number))=Data(i,29);
    Phase{phase_number}.Gantry(counter(phase_number))=Data(i,3);
    Phase{phase_number}.MU(counter(phase_number))=Data(i,25);
    %MU_shift=Data(i,21);
    Phase{phase_number}.MU_shift(counter(phase_number))=MU_shift;
    %need to re-organize the MLC data to match DICOM format
    for j=1:60  %  loop through 60 leaf
        Phase{phase_number}.MLC{counter(phase_number)}.MLC_A(1,j) = Data(i,33+2*j)*10.;
        Phase{phase_number}.MLC{counter(phase_number)}.MLC_B(1,j) = -Data(i,33+120+2*j)*10.;
    end




%I guess this has to do with preprocessing of data.


%***********************************
%   Steven added This section divides one arc into two arcs
%***********************************

%phase2 is a structure that contains all data for all 10 phases

for i=1:num_phases
    
 
 
    ARC1_index=find(Phase2{i}.ARC_num==1)
    
    PhaseARC1{i}.ARC_num=Phase2{i}.ARC_num(ARC1_index);
    PhaseARC1{i}.Control_Point=Phase2{i}.Control_Point(ARC1_index);
    PhaseARC1{i}.Gantry_IEC=Phase2{i}.Gantry_IEC(ARC1_index);
    PhaseARC1{i}.MLC=Phase2{i}.MLC(ARC1_index);
    PhaseARC1{i}.MU_final=Phase2{i}.MU_final(ARC1_index);
    
    ARC2_index=find(Phase2{i}.ARC_num==2)
    PhaseARC2{i}.ARC_num=Phase2{i}.ARC_num(ARC2_index);
    PhaseARC2{i}.Control_Point=Phase2{i}.Control_Point(ARC2_index);
    PhaseARC2{i}.Gantry_IEC=Phase2{i}.Gantry_IEC(ARC2_index);
    PhaseARC2{i}.MLC=Phase2{i}.MLC(ARC2_index);
    PhaseARC2{i}.MU_=Phase2{i}.MU(ARC2_index);
    PhaseARC2{i}.MU_final=Phase2{i}.MU_final(ARC2_index)-Phase2{i}.MU_final(ARC2_index(1));
end



%Now a different section. 



    



