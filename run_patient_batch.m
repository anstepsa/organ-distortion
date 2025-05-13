



function run_patient_batch()

    %{ 
    ---------------------------------------------------------------------------
    
    [1] Arrange DICOMs and segmentations as follows:-

          <BASE-DIR>\<PT_ID>\<DICOM-ADC-SUPINE>
                            \<DICOM-DCE-SUPINE>
                            \<DICOM-ADC-PRONE>
                            \<DICOM-DCE-PRONE>
                            \Segmentations\<*.nii.gz>   x 8
                    \<PT_ID>\ etc...
                              etc...

    [2] Check (by viewing) ROI coding numbers for NIFTI outlines is:-

            PG     = 1
            Rectum = 2
            Air    = 3

        [ If not amend 'supp_dist_control_script.m' accordingly at top... ]

    [3] [ May need to edit lines 573 + 575 in supp_dist_control_script.m
          accordingly if ROI file names aren't standard... ]

    [3] Input <PT_ID>s into config.PT_ID_LIST below....

    ---------------------------------------------------------------------------
    
    For QC:-

    [i]   Set BASE_DIR to '.\qc-processing'
    
    [ii]  Change input PT_ID_LIST
    
    [iii] Turn on VERBOSE
    
    [iv]  In 'dist_par_compute_radial_distortion_by_slice.m', comment out error
          on VERBOSE and remove 'parfor'
    
    [v]   In 'supp_dist_control_script_v4.m', set DILATE_BY_MM to just '10'...
          ... and in set_up_table_out(), replace stats.REF.AIR_VOLUME(2) & (3)
          with NaNs.

    [vi]  Optionally, un-comment figure saveas() lines in 'dist_par_compute_radial_distortion_by_slice.m' 
          and 'supp_show_mask_and_dicom.m'.

    %} 


    % ----------------------------------------------------------------------------
       
    config.DO_CALC_DISTORTION = true;
    config.DO_CALC_AIR_VOLUME = true;
   
    % ---- Set-up ----------------------------------------------------------------
   
    config.PT_ID_LIST         = {'V0001', ...
                                 'V0002', ...
                                 'V0003', ...
                                 'V0004', ...
                                 'V0005', ...
                                 'V0006', ...
                                 'V0007', ...
                                 'V0008', ...
                                 'V0009', ...
                                 'V0010', ...
                                 'V0011', ...
                                 'V0012', ...
                                 'V0013', ...
                                 'V0000' ...   % virtual phantom
                                 };
     
    config.PT_ID_LIST         = {'V0014', ...
                                 'V0015', ...
                                 'V0016', ...
                                 'V0017', ...
                                 'V0018' ...
                                 };
           
    config.PT_ID_LIST         = {'V0018'};
    
    config.PT_ID_LIST         = {'V0019', ...
                                 'V0020', ...
                                 'V0021', ...
                                 'V0022' ...
                                 };
                             
    config.PT_ID_LIST         = {'V0023', ...
                                 'V0024', ...
                                 'V0025', ...
                                 'V0026' ...
                                 };
       
    config.PT_ID_LIST         = {'V0027', ...
                                 'V0028', ...
                                 'V0029', ...
                                 'V0030', ...
                                 'V0031', ...
                                 'V0032' ...
                                 };
                              
    config.PT_ID_LIST         = {'V0033', ...
                                 'V0034', ...
                                 'V0035', ...
                                 'V0036' ...
                                };

    config.PT_ID_LIST         = {'V0036', ...
                                 'V0037', ...
                                 'V0038', ...
                                 'V0039', ...
                                 'V0040', ...
                                 'V0041', ...
                                 'V0042', ...
                                 'V0043' ...
                                };

    config.PT_ID_LIST         = {'V0044', ...
                                 'V0045' ...
                                };
                            
    config.PT_ID_LIST         = {'V0044'};
    config.PT_ID_LIST         = {'V0045'};

    config.PT_ID_LIST         = {'V0046', ...
                                 'V0047'};
                             
    config.PT_ID_LIST         = {'V0048', ...
                                 'V0049', ...
                                 'V0050', ...
                                 'V0051', ...                                 
                                 'V0052'...                                 
                                 };

    
    config.PT_ID_LIST         = {'V0053'};
   
    config.BASE_DIR           = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\anon-processing';
    
    config.OUT_DIR            = fullfile(config.BASE_DIR, 'results');

    config.WORK_DIR           = fullfile(config.BASE_DIR, 'work');
    
    config.TIME_STAMP         = datestr(datetime('now'), 'yyyy-mm-dd-HH-MM-SS'); %#ok<DATST>
  
    config.VERBOSE            = false;
    
    config.OUT_FILE           = [config.TIME_STAMP '-Distortion-Results.csv'];
   
    % ----------------------------------------------------------------------------
       
    if ~exist(config.OUT_DIR,  'dir'), mkdir(config.OUT_DIR);  end
    if ~exist(config.WORK_DIR, 'dir'), mkdir(config.WORK_DIR); end
     
    T_RESULTS_ALL = [];
    
    np = length(config.PT_ID_LIST);
    
    for p = 1:np

        config.PT_ID = config.PT_ID_LIST{p};

        disp( ' ');
        disp( '-------------------------------------');
        disp(['       Processing ' config.PT_ID]);
        disp( '-------------------------------------');
        disp( ' ');
        
        T_RESULTS_BOTH = supp_dist_control_script_v4(config);

        temp_out_file = fullfile(config.WORK_DIR, [config.PT_ID '-Distortion-Results.csv']);
        
        write_results_for_all(temp_out_file, T_RESULTS_BOTH);
        
        T_RESULTS_ALL = add_lines_to_table_out(T_RESULTS_ALL, T_RESULTS_BOTH);

        close('all');
        
    end
    
    out_file = fullfile(config.OUT_DIR, config.OUT_FILE);
    
    write_results_for_all(out_file, T_RESULTS_ALL)
    
    return
    
    
    function write_results_for_all(results_out_file, T_OUT)
       
        try
            writetable(T_OUT, results_out_file);
        catch ME
            warning('Error on writing results!');
            disp(ME.message);
        end
        
    end

    function T_OUT = add_lines_to_table_out(T_OUT, T_new)
       
        % Start the output table if blank or add the new line to it...
        
        if isempty(T_OUT)
            T_OUT = T_new;
        else
            T_OUT = cat(1, T_OUT, T_new);
        end
        
   end

end
