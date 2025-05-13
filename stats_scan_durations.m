

function stats_scan_durations()
% Compile a spreadsheet of all the scan durations for the prostate supine
% vs. prone study...

    ALL_RAW_DICOM_DIR = 'C:\Users\abg28\ANON_DICOM_DATA\CURRENT\ANON_TRISTAN\all-raw-data';

    OUT_DIR           = 'C:\Users\abg28\my_data\UNI_RESEARCH_PROJECTS\ANON_TRISTAN';

    OUT_FILE_NAME     = 'Prostate-Scan-Times.xlsx';

    % -----------------------------------------------------------------------------------------

    T_RESULTS = [];

    sd = dir(fullfile(ALL_RAW_DICOM_DIR, '2_Study_Imaging*'));

    nsd = length(sd);

    for s = 1:nsd

        SUB_DIR = fullfile(ALL_RAW_DICOM_DIR, sd(s).name);

        ed = dir(fullfile(SUB_DIR, '*'));

        ned = length(ed);

        for e = 3:ned

            exam_num = ed(e).name;

            EXAM_DIR = fullfile(SUB_DIR, ed(e).name);

            td = dir(fullfile(EXAM_DIR, '*'));
            
            ntd = length(td) - 2;

            assert(ntd == 1); 
            
            THIS_DIR = fullfile(EXAM_DIR, td(ntd+2).name);

            rd = dir(fullfile(THIS_DIR, '*'));

            nrd = length(rd);

            count = 0;

            for r = 3:nrd

                i = r-2;

                series_dir = fullfile(THIS_DIR, rd(r).name);

                series_num = str2double(rd(r).name(1:3));

                if series_num > 99
                    break
                end

                count = count + 1;

                d = dir(fullfile(series_dir, '*.dcm'));
  
                nd = length(d);

                H1 = dicominfo(fullfile(series_dir, d(1).name));
                Hn = dicominfo(fullfile(series_dir, d(nd).name));

                start_time{i} = H1.ContentTime; %#ok<*AGROW>
                end_time{i} = Hn.ContentTime;

                if i == 1
                    t0 = start_time{i};
                    t_elap_s(i) = 0;
                    t_int_s(i)  = 0;
                    t_ser_s(i)  = 0;
                else
                    t_elap_s(i)  = time_diff(t0, start_time{i});
                    t_int_s(i-1) = time_diff(start_time{i-1}, start_time{i});
                    t_ser_s(i)   = time_diff(start_time{i}, end_time{i});
                end

                series_name{i} = rd(r).name;
                exam_name{i}   = exam_num;

            end

            try
                series_dir = fullfile(THIS_DIR, series_name{count});
            catch
                disp('At error...');
            end

            d = dir(fullfile(series_dir, '*.dcm'));

            nd = length(d);

            Hn = dicominfo(fullfile(series_dir, d(nd).name));

            end_time = Hn.ContentTime;

            t_int_s(count) = time_diff(start_time{count}, end_time);
     

            for i = 1:count

                disp([exam_name{i} ' : ' series_name{i} ' : ' num2str(count) ' : ' start_time{i} ' : ' num2str(t_elap_s(i)) ' seconds' ' : ' num2str(t_int_s(i)) ' : ' num2str(t_ser_s(i)) ' seconds']);
               
            end

            exam_name   = exam_name';
            series_name = series_name';
            start_time  = start_time';
            t_elap_s    = t_elap_s';
            t_int_s     = t_int_s';
            t_ser_s     = t_ser_s';

            T_OUT = table(exam_name, series_name, start_time, t_elap_s, t_int_s, t_ser_s);

            clear exam_name series_name start_time end_time t_elap_s t_int_s t_ser_s

            if isempty(T_RESULTS)
                T_RESULTS = T_OUT;
            else
                T_RESULTS = cat(1, T_RESULTS, T_OUT);
            end

        end
     
    end

    writetable(T_RESULTS, fullfile(OUT_DIR, OUT_FILE_NAME));

    return


    function [t_secs, t_mins] = time_diff(t0, t1)


        t0 = t0(1:6);
        t1 = t1(1:6);

        dt0 = datetime(t0, 'InputFormat', 'HHmmss');
        dt1 = datetime(t1, 'InputFormat', 'HHmmss');

        tdiff = dt1 - dt0;

        t_mins = minutes(tdiff);
        t_secs = seconds(tdiff);

        return

    end
end
