

function stats_on_scan_times_ss()
% Produce summary statistics on scan times for prostate supine vs. prone study...

    OUT_DIR           = 'C:\Users\...\my_data\UNI_RESEARCH_PROJECTS\ANON_PROS';

    OUT_FILE_NAME     = 'Prostate-Scan-Times.xlsx';

    T_IN = readtable(fullfile(OUT_DIR, OUT_FILE_NAME));

    T_CHECK = readtable("C:\Users\...\ANON_DICOM_DATA\CURRENT\ANON_PROS\Distortion-Results-ALL-With-Graphs-With-Exam-Num.xlsx");

    % -------------------------------------------------------------------------------------------------

    VALID_EXAMS = T_CHECK(1:105,1);

    VALID_EXAMS = VALID_EXAMS(~isnan(VALID_EXAMS{:,1}),:);

    VALID_EXAMS = sortrows(VALID_EXAMS, "Var1");

    valid_exams = string(VALID_EXAMS{:,1});

    valid_exam_str = strcat("Exam", valid_exams);
    
    exam = T_IN{1,1};
    c = 1;
    count = 1;

    for r = 2:height(T_IN)
    
        if r == height(T_IN)
            c = c + 1;
            disp(char(T_IN{r,2}));
            disp([char(exam) ' : ' num2str(c)]);
        elseif strcmpi(T_IN{r,1}, exam)
            c = c + 1;
        else
            disp(char(T_IN{r-1,2}));
            disp([char(exam) ' : ' num2str(c)]);
            exam = T_IN{r,1};
            c = 1;
            count = count + 1;
        end
    end

    disp(['Count = ' num2str(count)])

    idx = ismember(T_IN.exam_name, valid_exam_str);

    T_IN_GOOD = T_IN(idx,:);

    T_IN_GOOD(559,:) = [];  % A duplicate series with contrast # Exam 8924

    good_exams = sortrows(unique(T_IN_GOOD{:,1}));

    disp(['Number of good exams = ' num2str(length(good_exams))]);


    idx = contains(T_IN_GOOD.series_name,'_test_Prone') | contains(T_IN_GOOD.series_name, 'ARDL_L_PRONE') | contains(T_IN_GOOD.series_name, 'T2W_axial_prone');

    T_TOTAL = T_IN_GOOD(idx, :);
  
    good_tot = sortrows(unique(T_TOTAL{:,1})); %#ok<NASGU>

  
    tot_times = T_TOTAL.t_elap_s + T_TOTAL.t_int_s;

    disp([ 'mean total exam time = (n = ' num2str(length(tot_times)) ') : ' num2str(floor(mean(tot_times) / 60)) ' minutes ' num2str(round(rem(mean(tot_times), 60))) ' seconds'...
           ' +/-' num2str(round(std(tot_times))) ' seconds']);

    disp([ 'range lo = ' num2str(floor(min(tot_times) /60)) ' minutes ' num2str(round(rem(min(tot_times), 60))) ' seconds' ]);
    disp([ 'range hi = ' num2str(floor(max(tot_times) /60)) ' minutes ' num2str(round(rem(max(tot_times), 60))) ' seconds' ]);

    idx_1 = contains(T_IN_GOOD.series_name,'_WATER_DCE_Lava_Flex_C');
    idx_2 = contains(T_IN_GOOD.series_name,'_WATER_DCE_MRI');

    T_SUPINE = T_IN_GOOD((idx_1 | idx_2), :);

    supine_times = T_SUPINE.t_elap_s + T_SUPINE.t_ser_s;

    disp([ 'mean supine exam time = (n = ' num2str(length(supine_times)) ') : ' num2str(floor(mean(supine_times) / 60)) ' minutes ' num2str(round(rem(mean(supine_times), 60))) ' seconds' ...
           ' +/-' num2str(round(std(supine_times))) ' seconds']);

    disp([ 'range lo = ' num2str(floor(min(supine_times) /60)) ' minutes ' num2str(round(rem(min(supine_times), 60))) ' seconds' ]);
    disp([ 'range hi = ' num2str(floor(max(supine_times) /60)) ' minutes ' num2str(round(rem(max(supine_times), 60))) ' seconds' ]);


    change_times = T_SUPINE.t_int_s - T_SUPINE.t_ser_s;

    disp([ 'mean change-over time = (n = ' num2str(length(change_times)) ') : ' num2str(floor(mean(change_times) / 60)) ' minutes ' num2str(round(rem(mean(change_times), 60))) ' seconds' ...
           ' +/-' num2str(round(std(change_times))) ' seconds']);

    disp([ 'range lo = ' num2str(floor(min(change_times) /60)) ' minutes ' num2str(round(rem(min(change_times), 60))) ' seconds' ]);
    disp([ 'range hi = ' num2str(floor(max(change_times) /60)) ' minutes ' num2str(round(rem(max(change_times), 60))) ' seconds' ]);

    idx_loc = contains(T_IN_GOOD.series_name,'_3Plane_Localizer_H_Res_Prone');
    idx_dwi = contains(T_IN_GOOD.series_name,'_SoC_Prone') | contains(T_IN_GOOD.series_name, 'ARDL_L_PRONE');

    idx = idx_loc | idx_dwi;

    T_DWI_PRONE = T_IN_GOOD(idx, :);

    count = 0;
    for r = 1:height(T_DWI_PRONE)

        if contains(T_DWI_PRONE.series_name{r}, 'DWI')
            count = count + 1;
            dwi_times(count, 1) = T_DWI_PRONE{r,"t_int_s"} + T_DWI_PRONE{r-1,"t_int_s"}; %#ok<AGROW>
        end

    end

    good_dwi = sortrows(unique(T_DWI_PRONE{:,1})); %#ok<NASGU>

    disp([ 'mean loc + DWI time = (n = ' num2str(length(dwi_times)) ') : ' num2str(floor(mean(dwi_times) / 60)) ' minutes ' num2str(round(rem(mean(dwi_times), 60))) ' seconds' ...
           ' +/-' num2str(round(std(dwi_times))) ' seconds']);
  
    disp([ 'range lo = ' num2str(floor(min(dwi_times) /60)) ' minutes ' num2str(round(rem(min(dwi_times), 60))) ' seconds' ]);
    disp([ 'range hi = ' num2str(floor(max(dwi_times) /60)) ' minutes ' num2str(round(rem(max(dwi_times), 60))) ' seconds' ]);

    added_times = change_times + dwi_times;

    disp([ 'mean added time = (n = ' num2str(length(added_times)) ') : ' num2str(floor(mean(added_times) /60)) ' minutes ' num2str(round(rem(mean(added_times), 60))) ' seconds' ...
          ' +/-' num2str(round(std(added_times))) ' seconds']);
  
    disp([ 'range lo = ' num2str(floor(min(added_times) /60)) ' minutes ' num2str(round(rem(min(added_times), 60))) ' seconds' ]);
    disp([ 'range hi = ' num2str(floor(max(added_times) /60)) ' minutes ' num2str(round(rem(max(added_times), 60))) ' seconds' ]);

end