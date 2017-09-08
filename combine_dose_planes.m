clear; clc 

Directory_name=uigetdir('','Select Folder containing DICOM RT Dose Planes');

dir_struct=dir(Directory_name);
cd(Directory_name);
[sorted_names,sorted_index] = sortrows({dir_struct.name}')

counter=1;

for i=1:length(sorted_names)
    [pathstr,name,ext] = fileparts(sorted_names{i});
    if ( strcmp(ext,'.dcm')  )
        DICOM_RTplan_names{counter}=strcat(name,ext);
        counter=counter+1;
    end
end

SumDose=zeros(301,301); %export dimensions +1
for i1=1:length(DICOM_RTplan_names)
    info_all{i1}=dicominfo(DICOM_RTplan_names{i1});
    phase_dose=double(dicomread(DICOM_RTplan_names{i1}));
    SumDose=SumDose+phase_dose;
end

write_info=info_all{1};
SumDoseWrite=uint32(SumDose);
%
plan_name='SUMDose.dcm';
%
dicomwrite(SumDoseWrite,plan_name,write_info,'CreateMode', 'copy');