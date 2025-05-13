



function util_sort_dist_files_for_analysis()

%{

    Sort all DICOM and ROI files for supine/prone prostate distortion study
    into arrangement suitable for analysis...

    Does some more anonymisation too...

%}

    config.START_PT_ID_AT = 52;  % 0  % 2  % 6  % The last 'Vxxxx' number allocated so far...

    config.RAW_DATA_PATH = 'C:\Users\...\ANON_DICOM_DATA\CURRENT\ANON_PROS\raw-data\8790';

    config.BASE_WRITE_PATH = 'C:\Users\...\ANON_DICOM_DATA\CURRENT\ANON_PROS\temp';
    
    config.PT_KEY_FILE = fullfile(config.BASE_WRITE_PATH, 'PatientKey.txt');
    
    config.BASE_SEG_PATH = fullfile(config.RAW_DATA_PATH, 'Segmentation');
    
    sd = dir(fullfile(config.RAW_DATA_PATH, '*.'));
    
    count = 0;
    
    for i = 1:length(sd)
        
        if contains(sd(i).name, '.') || contains(sd(i).name, 'seg', 'IgnoreCase',true), continue; end
        
        count = count + 1;
        
        pt_ref = sd(i).name;
        
        config.PT_DICOM_DIR = fullfile(config.RAW_DATA_PATH, pt_ref);
        config.PT_SEG_DIR   = fullfile(config.BASE_SEG_PATH, pt_ref);        
        
        pt_id_num = config.START_PT_ID_AT + count;
        
        config.PT_ID = ['V' num2str(pt_id_num, '%04d')];
        
        config.PT_WRITE_DIR  = fullfile(config.BASE_WRITE_PATH, config.PT_ID);
        config.SEG_WRITE_DIR = fullfile(config.PT_WRITE_DIR, 'Segmentations');
       
        if ~exist(config.PT_WRITE_DIR, 'dir'),  mkdir(config.PT_WRITE_DIR); end
        if ~exist(config.SEG_WRITE_DIR, 'dir'), mkdir(config.SEG_WRITE_DIR); end
        
        copy_dicom_series(config.PT_DICOM_DIR, config.PT_ID, config.PT_WRITE_DIR);
        copy_segmentations(config.PT_SEG_DIR, pt_ref, config.SEG_WRITE_DIR);
        
        note_key(pt_ref, config.PT_ID, config.PT_KEY_FILE);
        
    end

    return
    
    
    function note_key(pt_ref, pt_id, out_file)
       
        disp([pt_id ' == ' pt_ref]);
        
        fid = fopen(out_file, 'at');
        if fid < 0
            error('File failed to open!');
        end
        
        fprintf(fid, '%s == %s\n', pt_id, pt_ref);
        
        fclose(fid);
        
    end
    
    function copy_dicom_series(dicom_dir, pt_id, out_dir)
        
        ssd = dir(fullfile(dicom_dir, '*.'));
    
        nsd = length(ssd);
        c = 1;
        for d = 1:nsd
            if ~contains(ssd(d).name, '.') 
                these_dicom_dirs{c} = ssd(d).name;   %#ok<AGROW>
                c = c + 1;
            end
        end
        
        SUPINE_DICOM_DIRS = get_match(these_dicom_dirs, 'SUPINE', {{'ADC'},{'DCE'}});
        PRONE_DICOM_DIRS  = get_match(these_dicom_dirs, 'PRONE',  {{'ADC'},{'DCE'}});
        
        for s = 1:numel(SUPINE_DICOM_DIRS)
            my_copyfile(pt_id, fullfile(dicom_dir, SUPINE_DICOM_DIRS{s}), fullfile(out_dir, SUPINE_DICOM_DIRS{s}));
        end
        for s = 1:numel(PRONE_DICOM_DIRS)
            my_copyfile(pt_id, fullfile(dicom_dir, PRONE_DICOM_DIRS{s}),  fullfile(out_dir, PRONE_DICOM_DIRS{s}));
        end
        
    end

    function my_copyfile(pt_id, in_dir, out_dir)
        
        if ~exist(out_dir, 'dir'), mkdir(out_dir); end
        
        df = dir(fullfile(in_dir, '*.dcm'));
        
        for dd = 1:length(df)
            
            H = dicominfo(fullfile(in_dir, df(dd).name));
            I = dicomread(H);
            
            H.PatientName.FamilyName = pt_id;
            H.PatientName.GivenName  = '';
        
            dicomwrite(I, fullfile(out_dir, df(dd).name), H, 'CreateMode','Create', 'WritePrivate',true);
            
        end
        
    end


    function copy_segmentations(seg_dir, pt_ref, out_dir)
        
        rf = dir(fullfile(seg_dir, '*.nii.gz'));
        
        nrf = length(rf);
        
        assert(nrf == 8);
        
        pr_sz = strlength(pt_ref);

        for r = 1:nrf
            these_roi_files{r} = rf(r).name; %#ok<AGROW>
        end
        
        % A {string} or {set of strings} that defines : i.e.: ADC-PROSTATE,  WATER-PROSTATE,  WATER-AIR,  WATER-RECTUM
        SUPINE_ROI_FILES = get_match(these_roi_files, 'SUPINE', ...
                                          {{'ADC','PG'}, {'DCE','PG'}, {'DCE','AIR'}, {'DCE','RECTUM'}});
        PRONE_ROI_FILES  = get_match(these_roi_files, 'PRONE',  ...
                                          {{'ADC','PG'}, {'DCE','PG'}, {'DCE','AIR'}, {'DCE','RECTUM'}});
       
        for rr = 1:numel(SUPINE_ROI_FILES)
            new_name = SUPINE_ROI_FILES{rr};
            if contains(new_name, pt_ref)
                new_name = new_name(pr_sz+2:end);  % Assume pt_ref is [n] digits and at front of name...
            end
            copyfile(fullfile(seg_dir, SUPINE_ROI_FILES{rr}), fullfile(out_dir, new_name));
        end
        for rr = 1:numel(PRONE_ROI_FILES)
            new_name = PRONE_ROI_FILES{rr};
            if contains(new_name, pt_ref)
                new_name = new_name(pr_sz+2:end);  % Assume pt_ref is [n] digits and at front of name...
            end
            copyfile(fullfile(seg_dir, PRONE_ROI_FILES{rr}),  fullfile(out_dir, new_name));
        end
        
    end
      
        
    function str_out = get_match(C, pos, im_type)
        str_out = cell([numel(im_type),1]);
        for ii = 1:length(im_type)
            num_matches = 0;
            for cc = 1:length(C)
                if contains(C{cc}, pos, 'IgnoreCase',true)
                    found_str = zeros(numel(im_type{ii}),1);
                    for jj = 1:numel(im_type{ii})
                        if ~contains(C{cc}, im_type{ii}{jj}, 'IgnoreCase',true)
                            break
                        else
                            found_str(jj) = 1;
                        end
                    end
                    if all(found_str)
                        match_vals{ii}{num_matches+1} = C{cc}; %#ok<AGROW>
                        num_matches = num_matches + 1;
                    end
                end
            end
            if num_matches < 1
                warning('No matches in search for files or directories!');
                str_out{ii} = '<no matches>';               
            elseif num_matches > 1
                warning('More than one match in search for files or directories!');
                str_out{ii} = '<not resolved unambiguously>';
            else
                str_out{ii} = match_vals{ii}{1};
            end
        end
    end
             

end