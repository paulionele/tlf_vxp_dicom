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









%***********************************
%   Steven added This section divides one arc into two arcs
%***********************************

%phase2 is a structure that contains all data for all 10 phases

for i=1:num_phases
    
    ARC1_index = find(Phase2{i}.ARC_num==1)
    
    PhaseARC1{i}.ARC_num=Phase2{i}.ARC_num(ARC1_index);
    PhaseARC1{i}.Control_Point=Phase2{i}.Control_Point(ARC1_index);
    PhaseARC1{i}.Gantry_IEC=Phase2{i}.Gantry_IEC(ARC1_index);
    PhaseARC1{i}.MLC=Phase2{i}.MLC(ARC1_index);
    PhaseARC1{i}.MU_final=Phase2{i}.MU_final(ARC1_index);
    
    ARC2_index = find(Phase2{i}.ARC_num==2)
    PhaseARC2{i}.ARC_num=Phase2{i}.ARC_num(ARC2_index);
    PhaseARC2{i}.Control_Point=Phase2{i}.Control_Point(ARC2_index);
    PhaseARC2{i}.Gantry_IEC=Phase2{i}.Gantry_IEC(ARC2_index);
    PhaseARC2{i}.MLC=Phase2{i}.MLC(ARC2_index);
    PhaseARC2{i}.MU_=Phase2{i}.MU(ARC2_index);
    PhaseARC2{i}.MU_final=Phase2{i}.MU_final(ARC2_index)-Phase2{i}.MU_final(ARC2_index(1));
end




















% ******************************************************
%  Sorting done now Creating the DICOM planning files
%            To import into Eclipse TPS  2 ARC version
% *******************************************************

%List of items to be modified for field 1 (Item_1), same thing for field 2:
%RP.FractionGroupSequence(1).Item_1.ReferencedBeamSequence.Item_1.BeamMeterset
%RP.BeamSequence.Item_1.NumberOfControlPoints
%RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).ControlPointIndex
%RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).GantryAngle
%RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).BeamLimitingDevicePositionSequence.Item_3.LeafJawPositions ---> IF CP = 1
%RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).BeamLimitingDevicePositionSequence.Item_1.LeafJawPositions ---> IF CP > 1
%RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).BeamLimitingDevicePositionSequence.Item_1.RTBeamLimitingDeviceType
%RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).GantryRotationDirection = 'CW';
%RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).CumulativeMetersetWeight
%RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).ReferencedDoseReferenceSequence.Item_1.CumulativeDoseReferenceCoefficient ---> same as CumulativeMetersetWeight?
%RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).ReferencedDoseReferenceSequence.Item_1.ReferencedDoseReferenceNumber

%Select directory containing DICOM RT plans that were previously exported.
Directory_name = uigetdir('','Select Folder containing DICOM RT Plan');
dir_struct = dir(Directory_name);
cd(Directory_name);
[sorted_names,sorted_index] = sortrows({dir_struct.name}')

counter = 1; %this is an index for enumeration of the DICOM_RTplan_names

for i = 1:length(sorted_names)
    %Loop through files within the selected directory.
    [pathstr,name,ext] = fileparts(sorted_names{i});
    if strcmp(ext,'.dcm')
        %If file is a DICOM file...
        DICOM_RTplan_names{counter} = strcat(name,ext);
        counter = counter + 1;
    end
end



% [FileName_dicom,PathName_dicom] = uigetfile('*.dcm','Select folder containing DICOM RT Plan files for each phase');
%
% % before creating DICOM we will create the MLC as a matrix
% %  TT=horzcat(Phase{1}.MLC{1}.MLC_A,Phase{1}.MLC{1}.MLC_B);

for j = 1:num_phases
    %Loop through each of the phases, 10 or 20. Within this loop, the two
    %(or more) arcs are processed. Basically, a new 'phased plan' is
    %constructed for each phase using information from the phase and arc
    %sorted data. Once a phased plan has been created, the next phased plan
    %is processed until all 10 or 20 plans are generated.
    
    %Load current plan information into structure named RP. Names are
    %obtained from DICOM_RTplan_names cell array.
    RP = dicominfo(DICOM_RTplan_names{j});
    %         RP_init=dicominfo('RP.TTQUASAR.2 arcsST4.dcm')
    
    %I suspect that PhaseARC1_2 is the phase and arc sorted data from the
    %TLF-MW/VXP/PW. Simple variable rename here then...
    Phase = PhaseARC1_2;
    
    %Testing on phase 0 only for now.
    
    %Below, we begin to modify selected parts of the DICOM RP file. These
    %files were previously exported from the planning system. Some things
    %we don't want to change (i.e. isocenter coordinates already determined
    %for that phase, etc.???) and other things we need to. Essentially, we
    %replace whatever is necessary to make a particular plan a 'phased
    %plan' entirely... which included filling it with information that
    %from the cooresponding respiratory phase.
    
    %Max value of array(?) MU_final for cooresponding phase and for the 1st
    %arc. Not sure what MU_final would be an array?
    RP.FractionGroupSequence(1).Item_1.ReferencedBeamSequence.Item_1.BeamMeterset = max(Phase{j}.MU_final);
    
    %The number of control points (typically around 180 per field) is set
    %equal to the length of this array (which I don't know what is) Hmm
    %why? Also, this will lead to an issue when we merge (bring all 10/20
    %plans back together in the TPS), since we could easily end up with
    %more than 500 CP per field, which is not allowed. At some point we
    %will need to downsample, reduce the CP.
    num_control_points = length(Phase{j}.MU_final);
    RP.BeamSequence.Item_1.NumberOfControlPoints = num_control_points;
    
    for i = 1:num_control_points
        %Loop through all the CPs. Items within the ControlPointSequence are
        %enumerated starting with 1; the names are variable with 'Item_'.
        
        %Specifying some preliminary things. First, the CP index within each
        %item counts from 0. Second thing is the gantry rotation angle in
        %IEC scale.
        RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)])...
            .ControlPointIndex = i - 1;
        RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)])...
            .GantryAngle = Phase{j}.Gantry_IEC(i);
        
        if i == 1
            %If first CP in ControlPointSequence, we will do *this*.
            %horzcat concatenates two arrays horizontally, ex:
            %horzcat([1,2,3], [5,6,7]) -> [1,2,3.5,6,7]. Item_3 is listed
            %below, Item_1/2 refer to jaw positions.
            
            %Individual MLC positions represented as a m × 1 array where m
            %is the sum of leaves in leaf banks A and B. Leaf bank A
            %positions are indexed 0 to m/2-1. Leaf bank B positions are
            %indexed m/2 to m-1. The positions are specified in
            %millimeters. NOTE: it seems that MLC_B are specified first?
            RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)])...
                .BeamLimitingDevicePositionSequence.Item_3.LeafJawPositions =...
                (horzcat(Phase{j}.MLC{i}.MLC_B, Phase{j}.MLC{i}.MLC_A))';
            %Phase j, MLCs for CP i, and bank A or B.
        else
            %If not first CP...
            RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)])...
                .BeamLimitingDevicePositionSequence.Item_1.RTBeamLimitingDeviceType = 'MLCX';
            RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)])...
                .BeamLimitingDevicePositionSequence.Item_1.LeafJawPositions =...
                (horzcat(Phase{j}.MLC{i}.MLC_B, Phase{j}.MLC{i}.MLC_A))';
            
            %To recap, normally Item_3 contains the MLC positions but in
            %the above clause, they are listed in Item_1? Unless the 1st
            %CP, then they are in Item_3. In addition, I only see MLCX
            %specified as the device type, but MLCX and MLC positions
            %always seem to be specified together.
        end
        
        %Below is always done, regardless of the CP.
        RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)])...
            .GantryRotationDirection = 'CW';
        
        RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)])...
            .CumulativeMetersetWeight = Phase{j}.MU_final(i)/max(Phase{j}.MU_final);
        
        RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)])...
            .ReferencedDoseReferenceSequence.Item_1.CumulativeDoseReferenceCoefficient =...
            Phase{j}.MU_final(i)/max(Phase{j}.MU_final);
        
        RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)])...
            .ReferencedDoseReferenceSequence.Item_1.ReferencedDoseReferenceNumber = 1;
        
    end
    
    
    Phase = PhaseARC2_2;
    %Testing on phase 0 only for now
    RP.FractionGroupSequence(1).Item_1.ReferencedBeamSequence.Item_2.BeamMeterset=max(Phase{j}.MU_final);
    num_control_points=length(Phase{j}.MU_final);
    RP.BeamSequence.Item_2.NumberOfControlPoints=num_control_points;
    for i=1:num_control_points
        RP.BeamSequence.Item_2.ControlPointSequence.(['Item_' num2str(i)]).ControlPointIndex=i-1;
        RP.BeamSequence.Item_2.ControlPointSequence.(['Item_' num2str(i)]).GantryAngle=Phase{j}.Gantry_IEC(i);
        if i == 1
            RP.BeamSequence.Item_2.ControlPointSequence.(['Item_' num2str(i)]).BeamLimitingDevicePositionSequence.Item_3.LeafJawPositions=(horzcat(Phase{j}.MLC{i}.MLC_B,Phase{j}.MLC{i}.MLC_A))';
        else
            RP.BeamSequence.Item_2.ControlPointSequence.(['Item_' num2str(i)]).BeamLimitingDevicePositionSequence.Item_1.RTBeamLimitingDeviceType='MLCX';
            RP.BeamSequence.Item_2.ControlPointSequence.(['Item_' num2str(i)]).BeamLimitingDevicePositionSequence.Item_1.LeafJawPositions=(horzcat(Phase{j}.MLC{i}.MLC_B,Phase{j}.MLC{i}.MLC_A))';
        end
        RP.BeamSequence.Item_2.ControlPointSequence.(['Item_' num2str(i)]).GantryRotationDirection='CC';
        RP.BeamSequence.Item_2.ControlPointSequence.(['Item_' num2str(i)]).CumulativeMetersetWeight=Phase{j}.MU_final(i)/max(Phase{j}.MU_final);
        RP.BeamSequence.Item_2.ControlPointSequence.(['Item_' num2str(i)]).ReferencedDoseReferenceSequence.Item_1.CumulativeDoseReferenceCoefficient=...
            Phase{j}.MU_final(i)/max(Phase{j}.MU_final);
        RP.BeamSequence.Item_2.ControlPointSequence.(['Item_' num2str(i)]).ReferencedDoseReferenceSequence.Item_1.ReferencedDoseReferenceNumber=1;
        
    end
    %     dicom_output=['STAR5.dcm'];
    %     dicomwrite(1,dicom_output,RP,'CreateMode', 'copy');
    dicom_output=['TT_' DICOM_RTplan_names{j}];
    dicomwrite(1,dicom_output,RP,'CreateMode', 'copy');
    
end