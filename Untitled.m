Phase{phase_number}.ARC_num(counter(phase_number))=ARC_num(i); % Steven Added
Phase{phase_number}.Col(counter(phase_number))=Data(i,1); % Steven Added
Phase{phase_number}.Control_Point(counter(phase_number))=Data(i,29);
Phase{phase_number}.Gantry(counter(phase_number))=Data(i,3);
Phase{phase_number}.MU(counter(phase_number))=Data(i,25);
%MU_shift=Data(i,21);
Phase{phase_number}.MU_shift(counter(phase_number))=MU_shift;
%need to re-organize the MLC data to match DICOM format
for j=1:60
    Phase{phase_number}.MLC{counter(phase_number)}.MLC_A(1,j)=Data(i,33+2*j)*10.;
    Phase{phase_number}.MLC{counter(phase_number)}.MLC_B(1,j)=-Data(i,33+120+2*j)*10.;
end




for j=1:num_phases
    RP=dicominfo(DICOM_RTplan_names{j});
    %         RP_init=dicominfo('RP.TTQUASAR.2 arcsST4.dcm')
    
    Phase=PhaseARC1_2;
    %Testing on phase 0 only for now
    RP.FractionGroupSequence(1).Item_1.ReferencedBeamSequence.Item_1.BeamMeterset=max(Phase{j}.MU_final);
    num_control_points=length(Phase{j}.MU_final);
    RP.BeamSequence.Item_1.NumberOfControlPoints=num_control_points;
    for i=1:num_control_points
        RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).ControlPointIndex=i-1;
        RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).GantryAngle=Phase{j}.Gantry_IEC(i);
        if i==1
            RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).BeamLimitingDevicePositionSequence.Item_3.LeafJawPositions=(horzcat(Phase{j}.MLC{i}.MLC_B,Phase{j}.MLC{i}.MLC_A))';
        else
            RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).BeamLimitingDevicePositionSequence.Item_1.RTBeamLimitingDeviceType='MLCX';
            RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).BeamLimitingDevicePositionSequence.Item_1.LeafJawPositions=(horzcat(Phase{j}.MLC{i}.MLC_B,Phase{j}.MLC{i}.MLC_A))';
        end
        RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).GantryRotationDirection='CW';
        RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).CumulativeMetersetWeight=Phase{j}.MU_final(i)/max(Phase{j}.MU_final);
        RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).ReferencedDoseReferenceSequence.Item_1.CumulativeDoseReferenceCoefficient=...
            Phase{j}.MU_final(i)/max(Phase{j}.MU_final);
        RP.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(i)]).ReferencedDoseReferenceSequence.Item_1.ReferencedDoseReferenceNumber=1;
        
    end
    
    
    Phase=PhaseARC2_2;
    %Testing on phase 0 only for now
    RP.FractionGroupSequence(1).Item_1.ReferencedBeamSequence.Item_2.BeamMeterset=max(Phase{j}.MU_final);
    num_control_points=length(Phase{j}.MU_final);
    RP.BeamSequence.Item_2.NumberOfControlPoints=num_control_points;
    for i=1:num_control_points
        RP.BeamSequence.Item_2.ControlPointSequence.(['Item_' num2str(i)]).ControlPointIndex=i-1;
        RP.BeamSequence.Item_2.ControlPointSequence.(['Item_' num2str(i)]).GantryAngle=Phase{j}.Gantry_IEC(i);
        if i==1
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


