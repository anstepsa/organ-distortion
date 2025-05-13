


function stats_log_regression_air()

    config.RESULTS_DIR = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN';
    
    config.DIST_RESULTS_XLS = 'Distortion-Results-ALL-With-Graphs-With-Exam-Num.xlsx';

    config.LOADING_RESULTS_XLS = 'RECTAL LOADING_Prone Imaging patients_TB_anon.xlsx';
    
    config.PLOT_OUT_DIR = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\anon-processing\results\stats';
   
    config.THRESHOLD_SCORES = 4; %3:4;

    config.AIR_DISTANCES = 20; %[10, 15, 20];

    config.AIR_INDEX     = 3;  %[1,2,3]

    config.LEAVE_LAST_OUT = true;

    % ------------------------------------------------------------------------------------------------

    if ~exist(config.PLOT_OUT_DIR, 'dir'), mkdir(config.PLOT_OUT_DIR); end
    
    T_DIST_RESULTS = readtable(fullfile(config.RESULTS_DIR, config.DIST_RESULTS_XLS), 'Sheet','ALL-Corrected-d-post');
       
    T_LOADING_RESULTS = readtable(fullfile(config.RESULTS_DIR, config.LOADING_RESULTS_XLS));

    T_DIST_PRONE  = table();
    T_DIST_SUPINE = table();

    nr = height(T_DIST_RESULTS);

    if config.LEAVE_LAST_OUT
        nr = nr - 2;
    end

    for r = 1:nr
           
        if     strcmpi(T_DIST_RESULTS{r, 'scan_position'}, 'PRONE')
            T_DIST_PRONE  = cat(1, T_DIST_PRONE,  T_DIST_RESULTS(r,:));
        elseif strcmpi(T_DIST_RESULTS{r, 'scan_position'}, 'SUPINE')
            T_DIST_SUPINE = cat(1, T_DIST_SUPINE, T_DIST_RESULTS(r,:));
        else
            error('Unknown scan position...');
        end

    end

    x = NaN(height(T_DIST_SUPINE), 1);
    y = NaN(height(T_DIST_SUPINE), 1);  
    n = NaN(height(T_DIST_SUPINE), 1);  

    p = {};

    this_supine_score = NaN(height(T_DIST_SUPINE), 1);  

    for a = 1:length(config.AIR_DISTANCES)

        a_val = config.AIR_DISTANCES(a);
    
        for thr = config.THRESHOLD_SCORES
    
            for r = 1:height(T_DIST_SUPINE)
        
                this_pt_exam = T_DIST_SUPINE{r, 'pt'};

                p{r} = T_DIST_SUPINE{r, 'pt_id'}; %#ok<AGROW>

                x(r) = T_DIST_SUPINE{r, ['ref_prox_air_volume_' num2str(config.AIR_INDEX(a))]};
        
                found_this_pt = false;
        
                for i = 1:height(T_LOADING_RESULTS)
        
                    this_loading_pt = T_LOADING_RESULTS{i, 'ExamID'};
        
                    if this_loading_pt == this_pt_exam
                        
                        found_this_pt = true;

                        this_supine_score(r) = T_LOADING_RESULTS{i, 'SupineRectalScore'};
        
                        n(r) = T_LOADING_RESULTS{i, 'No'};

                        break
        
                    end
        
                end
        
                if found_this_pt
                    y(r) = (this_supine_score(r) >= thr);
                else
                    error('Patient exam not found!');
                end
           
            end
        
            do_logistic_reg(x, y, a_val, thr);

            % Compile a spreadsheet for export...

            T_OUT = [];

            for r = 1:height(T_DIST_SUPINE)
                C = {n(r), p{r}, x(r), this_supine_score(r)};
                if isempty(T_OUT)
                    T_OUT = cell2table(C, 'VariableNames',{'No', 'PtID', 'SupineAirVolume', 'SupineAirScore'});
                else
                    T_OUT = cat(1, T_OUT, cell2table(C, 'VariableNames',{'No', 'PtID', 'SupineAirVolume', 'SupineAirScore'}));
                end
            end
 
            writetable(T_OUT, fullfile(config.PLOT_OUT_DIR, 'Supine-Air-Results.xlsx'));

            
            % Quick QC check...

            x_n = T_OUT{:,'SupineAirVolume'};
            y_n = T_OUT{:,'SupineAirScore'};

            for r = 1:length(y_n)
                y_n(r) = (y_n(r) >= 4);
            end
    
            do_logistic_reg(x_n, y_n, a_val, thr);
           
        end
        
    end


    function fit_y = do_logistic_reg(x, y, a_val, thr)
    
        b = glmfit(x, y, 'binomial', 'link','logit');
       
        x_new = linspace(-10,20,100);
        z = b(1) + (x_new * b(2));
        fit_y = 1 ./ (1 + exp(-z));
    
        figure();
        plot(x, y, 'or', x_new, fit_y,':k');
    
        xlabel(['Air within ' num2str(a_val) ' mm of prostate [cm^3]']);
        ylabel(['Qualitative ''significant air'' from Likert >= ' num2str(thr)]);
    
        yticks([0,1]);
        yticklabels([0,1]);
    
        air_threshold = -b(1) / b(2);
    
        hold on
        plot([air_threshold, air_threshold], [0,1], 'b-');
    
        text(0.55, 0.1, ['Threshold air volume = ' num2str(air_threshold, '%.2f') ' cm^3'], 'Units','normalized');

        disp(['Threshold = ' num2str(air_threshold)]);
        
    end

end
