

function util_remake_nifti_vols()

    pt     = 'P0001';
    new_pt = 'V0001';
    
    in_nifti_dir  = fullfile('C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN', pt,     'Segmentations');
    
    out_nifti_dir = fullfile('C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN', new_pt, 'Segmentations');

    if ~exist(out_nifti_dir, 'dir'), mkdir(out_nifti_dir); end
    
    switch pt
        case 'P0000'
            prone.nifti_dis  = fullfile(in_nifti_dir, '7094_2050_ADC_10_6_mm_s_Prone_VP.nii.gz');
            prone.nifti_ref  = fullfile(in_nifti_dir, '7094_2151_Water_DCE_Prone_FirstPhase_VP.nii.gz');
            supine.nifti_dis = fullfile(in_nifti_dir, '7094_950_ADC_10_6_mm_s_Supine_VP.nii.gz');
            supine.nifti_ref = fullfile(in_nifti_dir, '7094_1751_Water_DCE_Supine_LastPhase_VP.nii.gz');
        case 'P0001'
            prone.nifti_dis  = fullfile(in_nifti_dir, '7094_2050_ADC_10_6_mm_s_Prone.nii.gz');
            prone.nifti_ref  = fullfile(in_nifti_dir, '7094_2151_Water_DCE_Prone_FirstPhase.nii.gz');
            supine.nifti_dis = fullfile(in_nifti_dir, '7094_950_ADC_10_6_mm_s_Supine.nii.gz');
            supine.nifti_ref = fullfile(in_nifti_dir, '7094_1751_Water_DCE_Supine_LastPhase.nii.gz');
    end
    
    [V, h] = read_nifti(prone.nifti_dis);
    
    data_type = h.Datatype;
    data_fn_h = str2func(data_type);
    
    V1 = (V == 1);
    niftiwrite(data_fn_h(V1), fullfile(out_nifti_dir, 'Prone_ADC_Prostate.nii'), h, 'Compressed',true);
 
    [V, h] = read_nifti(supine.nifti_dis);
     
    data_type = h.Datatype;
    data_fn_h = str2func(data_type);
     
    V1 = (V == 1);
    niftiwrite(data_fn_h(V1), fullfile(out_nifti_dir, 'Supine_ADC_Prostate.nii'), h, 'Compressed',true);
 
    [V, h] = read_nifti(prone.nifti_ref);
    
    data_type = h.Datatype;
    data_fn_h = str2func(data_type);
    
    V1 = (V == 1);
    niftiwrite(data_fn_h(V1), fullfile(out_nifti_dir, 'Prone_T1W_Prostate.nii'), h, 'Compressed',true);
 
    V2 = (V == 2);
    niftiwrite(data_fn_h(V2), fullfile(out_nifti_dir, 'Prone_T1W_Air.nii'), h, 'Compressed',true);
 
    V3 = (V == 2);
    niftiwrite(data_fn_h(V3), fullfile(out_nifti_dir, 'Prone_T1W_Rectum.nii'), h, 'Compressed',true);
    
    [V, h] = read_nifti(supine.nifti_ref);
    
    data_type = h.Datatype;
    data_fn_h = str2func(data_type);
      
    V1 = (V == 1);
    niftiwrite(data_fn_h(V1), fullfile(out_nifti_dir, 'Supine_T1W_Prostate.nii'), h, 'Compressed',true);
 
    V2 = (V == 2);
    niftiwrite(data_fn_h(V2), fullfile(out_nifti_dir, 'Supine_T1W_Air.nii'), h, 'Compressed',true);
 
    V3 = (V == 2);
    niftiwrite(data_fn_h(V3), fullfile(out_nifti_dir, 'Supine_T1W_Rectum.nii'), h, 'Compressed',true);
    
    return
    
    
    function [V, h] = read_nifti(nii_file_in)
        V = niftiread(nii_file_in);
        h = niftiinfo(nii_file_in);
    end
    
end