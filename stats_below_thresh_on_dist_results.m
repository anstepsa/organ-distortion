

function stats_below_thresh_on_dist_results()

    config.SAVE_FIGS = true;%false;

    addpath('C:\Users\abg28\my_data\UNI_SOURCE\Statistics\swtest');

    config.RESULTS_DIR = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\anon-processing\results';
    
    config.RESULTS_XLS = 'Distortion-Results-ALL-With-Graphs-8.xlsx';
    
    config.PLOT_OUT_BASE = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\anon-processing\results\stats-n-33';
    
    config.AIR_THRESHOLD_LEVELS = [NaN, NaN, 4.0]; %[   2.0,     3.0,     5.0]; %4.0];   % [cm^3] within {air_distance}
    
    config.AIR_DISTANCES        = {'10 mm', '15 mm', '20 mm'};

    config.DIST_MEASURES = {'POSTERIOR'};

    % ------------------------------------------------------------------------------------
        
    for a = 3:length(config.AIR_DISTANCES)

        for d = 1:length(config.DIST_MEASURES)

            air_measure = ['ref_prox_air_volume_' num2str(a)];

            air_distance = config.AIR_DISTANCES{a};
            
            config.AIR_THRESHOLD = config.AIR_THRESHOLD_LEVELS(a);

            config.DIST_MEASURE_TO_USE = config.DIST_MEASURES{d};

            if strcmpi(config.DIST_MEASURE_TO_USE, 'POSTERIOR')
                distortion_measure = 'd_post_res';
            else
                distortion_measure = 'd_res';
            end

            config.PLOT_OUT_DIR = fullfile(config.PLOT_OUT_BASE, [config.DIST_MEASURE_TO_USE '-Air-Within-' strrep(air_distance,' ','-')]);

            do_these_plots(config, air_measure, air_distance, distortion_measure);

            %pause();

        end

    end

    return


    function do_these_plots(config, air_measure, air_distance, dist_measure)

        SAVE_FIGS = config.SAVE_FIGS;

        if ~exist(config.PLOT_OUT_DIR, 'dir'), mkdir(config.PLOT_OUT_DIR); end
    
        T_RESULTS = readtable(fullfile(config.RESULTS_DIR, config.RESULTS_XLS), 'Sheet','ALL-Corrected-d-post');
       
        T_PRONE  = table();
        T_SUPINE = table();
        
        T_PRONE_BELOW_SUPINE_AIR_THRESH  = table();
        T_SUPINE_BELOW_SUPINE_AIR_THRESH = table();
    
        for r = 1:height(T_RESULTS)
           
            if     strcmpi(T_RESULTS{r, 'scan_position'}, 'PRONE')
                T_PRONE  = cat(1, T_PRONE,  T_RESULTS(r,:));
            elseif strcmpi(T_RESULTS{r, 'scan_position'}, 'SUPINE')
                T_SUPINE = cat(1, T_SUPINE, T_RESULTS(r,:));
            else
                error('Unknown scan position...');
            end
            
        end
    
        for sr = 1:height(T_SUPINE)
            if T_SUPINE{sr, air_measure} <= config.AIR_THRESHOLD
                T_SUPINE_BELOW_SUPINE_AIR_THRESH = cat(1, T_SUPINE_BELOW_SUPINE_AIR_THRESH, T_SUPINE(sr,:));
                T_PRONE_BELOW_SUPINE_AIR_THRESH  = cat(1, T_PRONE_BELOW_SUPINE_AIR_THRESH, T_PRONE(sr,:));
            end        
        end
    
        %{
        hf6a1 = figure();
        plot_scatter_and_correlate(hf6a1, T_SUPINE_ABOVE_SUPINE_AIR_THRESH{:,air_measure}, ...
                                          T_PRONE_ABOVE_SUPINE_AIR_THRESH{:,air_measure}, ...
                                          ['Pts. with SUPINE air > ' num2str(config.AIR_THRESHOLD) 'cm^3'], ...
                                          ['Volume of air within ' air_distance ' of prostate SUPINE [cm^3]'], ...
                                          ['Volume of air within ' air_distance ' of prostate PRONE  [cm^3]'], ...
                                          'air  ', ...
                                          'air  ', ...
                                          false);
    
        if SAVE_FIGS
            saveas(hf6a1, fullfile(config.PLOT_OUT_DIR, 'hf06a1.png'));
            savefig(hf6a1, fullfile(config.PLOT_OUT_DIR, 'hf06a1.fig'));
        end
        %}
    
        hf6a2 = figure();
        box_plot_and_ttest(hf6a2, cat(2,  T_SUPINE_BELOW_SUPINE_AIR_THRESH{:,air_measure}, ...
                                          T_PRONE_BELOW_SUPINE_AIR_THRESH{:,air_measure}), ...
                                          ['Pts. with SUPINE air <= ' num2str(config.AIR_THRESHOLD) 'cm^3'], ...
                                          ['Volume of air within ' air_distance ' of prostate [cm^3]']);
    
        if SAVE_FIGS
            saveas(hf6a2, fullfile(config.PLOT_OUT_DIR, 'hf06a2.png'));
            savefig(hf6a2, fullfile(config.PLOT_OUT_DIR, 'hf06a2.fig'));
        end
    
        %{
        hf6b1 = figure();
        plot_scatter_and_correlate(hf6b1, T_SUPINE_ABOVE_SUPINE_AIR_THRESH{:,dist_measure}, ...
                                          T_PRONE_ABOVE_SUPINE_AIR_THRESH{:,dist_measure}, ...
                                          [config.DIST_MEASURE_TO_USE ' OF GLAND: Pts. with SUPINE air > ' num2str(config.AIR_THRESHOLD) 'cm^3'], ...
                                          'Distortion measure [mm]', ...
                                          'Distortion measure [mm]', ...                                      
                                          'dist. SUPINE ', ...
                                          'dist. PRONE  ', ...
                                          false);
    
        if SAVE_FIGS
            saveas(hf6b1, fullfile(config.PLOT_OUT_DIR, 'hf06b1.png'));
            savefig(hf6b1, fullfile(config.PLOT_OUT_DIR, 'hf06b1.fig'));
        end
        %}
    
        hf6b2 = figure();
        box_plot_and_ttest(hf6b2, cat(2,  T_SUPINE_BELOW_SUPINE_AIR_THRESH{:,dist_measure}, ...
                                          T_PRONE_BELOW_SUPINE_AIR_THRESH{:,dist_measure}), ...
                                          [config.DIST_MEASURE_TO_USE ' OF GLAND: Pts. with SUPINE air <= ' num2str(config.AIR_THRESHOLD) 'cm^3'], ...
                                          'Distortion measure [mm]');
    
        if SAVE_FIGS
            saveas(hf6b2, fullfile(config.PLOT_OUT_DIR, 'hf06b2.png'));
            savefig(hf6b2, fullfile(config.PLOT_OUT_DIR, 'hf06b2.fig'));
        end
    
        hf7a1 = figure();
        plot_as_clusters(hf7a1, T_PRONE_BELOW_SUPINE_AIR_THRESH{:,air_measure}, ...
                             T_PRONE_BELOW_SUPINE_AIR_THRESH{:,dist_measure}, ...                           
                             T_SUPINE_BELOW_SUPINE_AIR_THRESH{:,air_measure}, ...
                             T_SUPINE_BELOW_SUPINE_AIR_THRESH{:,dist_measure}, ...  
                             [config.DIST_MEASURE_TO_USE ' OF GLAND: Pts. with SUPINE air <= ' num2str(config.AIR_THRESHOLD) 'cm^3'], ...
                             ['Volume of air within ' air_distance ' of prostate [cm^3]'], ...
                             'Distortion measure [mm]', ...
                             'Prone', 'Supine', ...
                             config.AIR_THRESHOLD ...
                             );
    
        if SAVE_FIGS
            saveas(hf7a1, fullfile(config.PLOT_OUT_DIR, 'hf07a1.png')); %#ok<*UNRCH>
            savefig(hf7a1, fullfile(config.PLOT_OUT_DIR, 'hf07a1.fig'));
        end
    
        return

    end


    function p_ttest = box_plot_and_ttest(hf, pts, f_title, y_label)
        
        if nargin < 4
            y_label = 'Distortion measure [mm]';
        end
        
        figure(hf);
        set(hf, 'Color','w');

        boxplot(pts, 'Labels',{'SUPINE','PRONE'});
        
        %ylim([0,20]);
        
        [x_pts_not_normal, p_x] = swtest(pts(:,1)); %#ok<ASGLU>
        [y_pts_not_normal, p_y] = swtest(pts(:,2)); %#ok<ASGLU>
        
        if ~x_pts_not_normal && ~y_pts_not_normal

            [~, p_ttest] = ttest(pts(:,1), pts(:,2));

            text(0.5, 0.8, ['t-test p = ' num2str(p_ttest, '%.3f')], 'Units','normalized');
        
        else
            
            %p_mwutest = ranksum(pts(:,1), pts(:,2));
            p_mwutest = signrank(pts(:,1), pts(:,2));  % Paired test here...
            
            text(0.5, 0.8, ['Wilcoxon p = ' num2str(p_mwutest, '%.3f')], 'Units','normalized');

        end
        
        xlabel('Patient position');
        ylabel(y_label); 
        
        title(f_title);
        
    end

    function [rho, p_val] = plot_scatter_and_correlate(hf, x_pts, y_pts, f_title, x_label, y_label, msg_mean_x, msg_mean_y, do_ttest, add_x_pts, add_y_pts) %#ok<DEFNU>
        
        if nargin < 4
            f_title = 'POSTERIOR OF GLAND';
        end
        
        if nargin < 5
            y_label = 'Distortion measure [mm]';
            x_label = 'Volume of air within --- of prostate [cm^3]';
            msg_mean_x = 'air  ';
            msg_mean_y = 'dist.';
            do_ttest = false;
        end
        
        if nargin < 10
            add_x_pts = [];
            add_y_pts = [];
        end

        figure(hf);
        set(hf, 'Color','w');
        
        scatterhist(x_pts, y_pts);
        
        if ~isempty(add_x_pts)
        
            hold on
            plot(gca, add_x_pts, add_y_pts, 'r+')
        
        end

        xlim([0,20]);
        ylim([0,20]);
        
        coeff = polyfit(x_pts, y_pts, 1);
        
        m = coeff(1);
        c = coeff(2);
        
        hold on
        bf = plot([0,20], [c,m*20+c], 'r.:');
        set(bf, 'LineWidth',2);
 
        [x_pts_not_normal, p_x] = swtest(x_pts); %#ok<ASGLU>
        [y_pts_not_normal, p_y] = swtest(y_pts); %#ok<ASGLU>
 
        if ~x_pts_not_normal && ~y_pts_not_normal
            [rho, p_val] = corr(x_pts, y_pts, 'Type','Pearson');
            text(1,15.0, ['Pearson r   = ' num2str(rho, '%.3f')]);
        else
            [rho, p_val] = corr(x_pts, y_pts, 'Type','Spearman');
            text(1,15.0, ['Spearman rho = ' num2str(rho, '%.3f')]);
        end
        
        text(1,13, ['p   = ' num2str(p_val, '%.3f')]);        
 
        text(1,11, ['(mean, median) ' msg_mean_y ' = (' num2str(mean(y_pts), '%.3f') ', ' num2str(median(y_pts), '%.3f') ')']);
        text(1, 9, ['(mean, median) ' msg_mean_x ' = (' num2str(mean(x_pts), '%.3f') ', ' num2str(median(x_pts), '%.3f') ')']);
   
        if do_ttest
            
            if ~x_pts_not_normal && ~y_pts_not_normal
                [~, p_ttest] = ttest(x_pts, y_pts);
                text(1,6.5, ['t-test p = ' num2str(p_ttest, '%.3f')]);
            else
                %p_mwutest = ranksum(x_pts, y_pts);
                p_mwutest = signrank(x_pts, y_pts);
                text(1,6.5, ['Wilcoxon p = ' num2str(p_mwutest, '%.3f')]);
            end
            
        end
           
        ylabel(y_label);
        xlabel(x_label);
        
        title(f_title);
        
    end

    function plot_as_clusters(hf, x_pts, y_pts, add_x_pts, add_y_pts, ...
                             f_title, x_label, y_label, desc_1, desc_2, threshold)
        
        figure(hf);
        set(hf, 'Color','w');
        
        plot(x_pts, y_pts, 'bo');
     
        xlim([0,10]);
        ylim([0,15]);
         
        if ~isempty(add_x_pts)
        
            hold on
            plot(gca, add_x_pts, add_y_pts, 'ro')
        
            plot(gca, [threshold,threshold],[0,20],'k:');
 
            legend(desc_1, desc_2, 'threshold');
        
        end

        ylabel(y_label);
        xlabel(x_label);
        
        title(f_title);
        
    end
    
end 
