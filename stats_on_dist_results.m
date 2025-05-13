

function stats_on_dist_results()

    addpath('C:\Users\...\my_data\UNI_SOURCE\Statistics\swtest');

    config.RESULTS_DIR = 'C:\Users\...\ANON_DICOM_DATA\CURRENT\ANON_PROS\anon-processing\results';
    
    config.RESULTS_XLS = 'Distortion-Results-ALL-With-Graphs-9.xlsx';
    
    config.PLOT_OUT_DIR = 'C:\Users\...\ANON_DICOM_DATA\CURRENT\ANON_PROS\anon-processing\results\stats';
    
    config.AIR_THRESHOLD = 4.0;
    
    % ------------------------------------------------------------------------------------
    
    if ~exist(config.PLOT_OUT_DIR, 'dir'), mkdir(config.PLOT_OUT_DIR); end
    
    T_RESULTS = readtable(fullfile(config.RESULTS_DIR, config.RESULTS_XLS), 'Sheet','ALL-Corrected-d-post');
   
    T_PRONE  = table();
    T_SUPINE = table();
    
    T_PRONE_ABOVE_THRESH  = table();
    T_SUPINE_ABOVE_THRESH = table();

    T_PRONE_BELOW_THRESH  = table();
    T_SUPINE_BELOW_THRESH = table();
    
    for r = 1:height(T_RESULTS)
       
        if     strcmpi(T_RESULTS{r, 'scan_position'}, 'PRONE')
            T_PRONE  = cat(1, T_PRONE,  T_RESULTS(r,:));
            if T_RESULTS{r, 'ref_prox_air_volume_3'} > config.AIR_THRESHOLD
                T_PRONE_ABOVE_THRESH = cat(1, T_PRONE_ABOVE_THRESH, T_RESULTS(r,:));
            else
                T_PRONE_BELOW_THRESH  = cat(1, T_PRONE_BELOW_THRESH, T_RESULTS(r,:));
            end
        elseif strcmpi(T_RESULTS{r, 'scan_position'}, 'SUPINE')
            T_SUPINE = cat(1, T_SUPINE, T_RESULTS(r,:));
            if T_RESULTS{r, 'ref_prox_air_volume_3'} > config.AIR_THRESHOLD
                T_SUPINE_ABOVE_THRESH = cat(1, T_SUPINE_ABOVE_THRESH, T_RESULTS(r,:));
            else
                T_SUPINE_BELOW_THRESH  = cat(1, T_SUPINE_BELOW_THRESH,  T_RESULTS(r,:));
            end
        else
            error('Unknown scan position...');
        end
        
    end

    hf1 = figure();
    plot_scatter_and_correlate(hf1, T_PRONE{:,'ref_prox_air_volume_3'}, T_PRONE{:,'d_res'}, 'WHOLE GLAND : prone');
    my_saveas(hf1, fullfile(config.PLOT_OUT_DIR, 'hf01.png'));
        
    hf2 = figure();
    plot_scatter_and_correlate(hf2, T_SUPINE{:,'ref_prox_air_volume_3'}, T_SUPINE{:,'d_res'}, 'WHOLE GLAND : supine');
    my_saveas(hf2, fullfile(config.PLOT_OUT_DIR, 'hf02.png'));

    hf3 = figure();
    plot_scatter_and_correlate(hf3, T_PRONE{:,'ref_prox_air_volume_3'}, T_PRONE{:,'d_post_res'}, 'POSTERIOR OF GLAND : prone');
    my_saveas(hf3, fullfile(config.PLOT_OUT_DIR, 'hf03.png'));
    
    hf4 = figure();
    plot_scatter_and_correlate(hf4, T_SUPINE{:,'ref_prox_air_volume_3'}, T_SUPINE{:,'d_post_res'}, 'POSTERIOR OF GLAND : supine');
    my_saveas(hf4, fullfile(config.PLOT_OUT_DIR, 'hf04.png'));
    
    hf5 = figure();
    plot_scatter_and_correlate(hf5, T_PRONE_ABOVE_THRESH{:,'ref_prox_air_volume_3'}, T_PRONE_ABOVE_THRESH{:,'d_post_res'}, 'POSTERIOR OF GLAND : prone : air > threshold');
    my_saveas(hf5, fullfile(config.PLOT_OUT_DIR, 'hf05.png'));

    hf6 = figure();
    plot_scatter_and_correlate(hf6, T_SUPINE_ABOVE_THRESH{:,'ref_prox_air_volume_3'}, T_SUPINE_ABOVE_THRESH{:,'d_post_res'}, 'POSTERIOR OF GLAND : supine : air > threshold');
    my_saveas(hf6, fullfile(config.PLOT_OUT_DIR, 'hf06.png'));

    hf7 = figure();
    plot_scatter_and_correlate(hf7, T_PRONE_BELOW_THRESH{:,'ref_prox_air_volume_3'}, T_PRONE_BELOW_THRESH{:,'d_post_res'}, 'POSTERIOR OF GLAND : prone : air < threshold');
    my_saveas(hf7, fullfile(config.PLOT_OUT_DIR, 'hf07.png'));

    hf8 = figure();
    plot_scatter_and_correlate(hf8, T_SUPINE_BELOW_THRESH{:,'ref_prox_air_volume_3'}, T_SUPINE_BELOW_THRESH{:,'d_post_res'}, 'POSTERIOR OF GLAND : supine : air < threshold');
    my_saveas(hf8, fullfile(config.PLOT_OUT_DIR, 'hf08.png'));

    hf9 = figure();
    plot_scatter_and_correlate(hf9, T_SUPINE{:,'d_res'}, T_PRONE{:,'d_res'}, 'WHOLE GLAND : prone vs supine distortion', ...
        'SUPINE distortion [mm]', 'PRONE distortion [mm]', 'supine dist.', 'prone dist.', ...
        true);
    my_saveas(hf9, fullfile(config.PLOT_OUT_DIR, 'hf09.png'));

    hf10 = figure();
    plot_scatter_and_correlate(hf10, T_SUPINE{:,'d_post_res'}, T_PRONE{:,'d_post_res'}, 'POSTERIOR OF GLAND : prone vs supine distortion', ...
        'SUPINE distortion [mm]', 'PRONE distortion [mm]', 'supine dist.', 'prone dist.', ...
        true);
    my_saveas(hf10, fullfile(config.PLOT_OUT_DIR, 'hf10.png'));
    
    hf11 = figure();
    box_plot_and_ttest(hf11, cat(2, T_SUPINE{:,'d_res'}, T_PRONE{:,'d_res'}), 'WHOLE GLAND : supine | prone distortion');
    my_saveas(hf11, fullfile(config.PLOT_OUT_DIR, 'hf11.png'));

    hf12 = figure();
    box_plot_and_ttest(hf12, cat(2, T_SUPINE{:,'d_post_res'}, T_PRONE{:,'d_post_res'}), 'POSTERIOR OF GLAND : supine | prone distortion');
    my_saveas(hf12, fullfile(config.PLOT_OUT_DIR, 'hf12.png'));
    
    hf13 = figure();
    plot_scatter_and_correlate(hf13, T_SUPINE{:,'ref_prox_air_volume_3'}, T_PRONE{:,'ref_prox_air_volume_3'}, 'AIR VOLUME : prone vs supine', ...
        'SUPINE air within 20 mm [cm^3]', 'PRONE air within 20 mm [cm^3]', 'supine air', 'prone air', ...
        true);
    my_saveas(hf13, fullfile(config.PLOT_OUT_DIR, 'hf13.png'));
       
    hf14 = figure();
    box_plot_and_ttest(hf14, ...
        cat(2, T_SUPINE{:,'ref_prox_air_volume_3'}, T_PRONE{:,'ref_prox_air_volume_3'}), 'AIR VOLUME : supine | prone', ...
        'Air volume within 20 mm [cm^3]');
    my_saveas(hf14, fullfile(config.PLOT_OUT_DIR, 'hf14.png'));
    
    return
    

    function my_saveas(hf, f_path_png)
        saveas(hf, f_path_png);
        savefig(hf, strrep(f_path_png, '.png', '.fig'));
    end

    function p_ttest = box_plot_and_ttest(hf, pts, f_title, y_label)
        
        if nargin < 4
            y_label = 'Distortion measure [mm]';
        end
        
        figure(hf);
        set(hf, 'Color','w');

        boxplot(pts, 'Labels',{'SUPINE','PRONE'});
        
        ylim([0,10]);
        
        [x_pts_not_normal, p_x] = swtest(pts(:,1)); %#ok<ASGLU>
        [y_pts_not_normal, p_y] = swtest(pts(:,2)); %#ok<ASGLU>
        
        if ~x_pts_not_normal && ~y_pts_not_normal

            [~, p_ttest] = ttest(pts(:,1), pts(:,2));

            text(1.5,8.5, ['t-test p = ' num2str(p_ttest, '%.3f')]);
        
        else
            
            %p_mwutest = ranksum(pts(:,1), pts(:,2));
            p_mwutest = signrank(pts(:,1), pts(:,2));  % Paired test here...
            
            text(1.5,8.5, ['Wilcoxon p = ' num2str(p_mwutest, '%.3f')]);

        end
        
        xlabel('Patient position');
        ylabel(y_label); 
        
        title(f_title);
        
    end

    function [rho, p_val] = plot_scatter_and_correlate(hf, x_pts, y_pts, f_title, x_label, y_label, msg_mean_x, msg_mean_y, do_ttest)
        
        if nargin < 4
            f_title = 'POSTERIOR OF GLAND';
        end
        
        if nargin < 5
            y_label = 'Distortion measure [mm]';
            x_label = 'Volume of air within 20 mm of prostate [cm^3]';
            msg_mean_x = 'air  ';
            msg_mean_y = 'dist.';
            do_ttest = false;
        end
        
        figure(hf);
        set(hf, 'Color','w');
        
        scatterhist(x_pts, y_pts);
        
        xlim([0,10]);
        ylim([0,10]);
        
        coeff = polyfit(x_pts, y_pts, 1);
        
        m = coeff(1);
        c = coeff(2);
        
        hold on
        bf = plot([0,10], [c,m*10+c], 'r.:');
        set(bf, 'LineWidth',2);
 
        [x_pts_not_normal, p_x] = swtest(x_pts); %#ok<ASGLU>
        [y_pts_not_normal, p_y] = swtest(y_pts); %#ok<ASGLU>
 
        if ~x_pts_not_normal && ~y_pts_not_normal
            [rho, p_val] = corr(x_pts, y_pts, 'Type','Pearson');
            text(1,9.0, ['Pearson r   = ' num2str(rho, '%.3f')]);
        else
            [rho, p_val] = corr(x_pts, y_pts, 'Type','Spearman');
            text(1,9.0, ['Spearman rho = ' num2str(rho, '%.3f')]);
        end
        
        text(1,8.5, ['p   = ' num2str(p_val, '%.3f')]);        
 
        text(1,7.5, ['(mean, median) ' msg_mean_y ' = (' num2str(mean(y_pts), '%.3f') ', ' num2str(median(y_pts), '%.3f') ')']);
        text(1,7.0, ['(mean, median) ' msg_mean_x ' = (' num2str(mean(x_pts), '%.3f') ', ' num2str(median(x_pts), '%.3f') ')']);
   
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
    
end