%% this section makes time series and scatter plots for each site-year
% read in data
clear all; 
close all;
clc;
threshold_FCH4 = 0;
threshold_T = 0;
delta_T = 0.1;
degree = sprintf('%c', char(176));
id_name = ["(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)",...
    "(m)", "(n)", "(o)", "(p)", "(q)", "(r)"];
% variable names for analysis
for variable = 1:3
    if (variable==1)
        var_wanted = ["Year","DOY","FCH4","TS_1","TS_2","TS_3","TA",...
            "GPP_DT","GPP_NT","WTD","P_F","SWC"];
        var_dir = 'FCH4/';
    elseif (variable==2)
        var_wanted = ["Year","DOY","FCH4_F","TS_1","TS_2","TS_3","TA_F",...
            "GPP_DT","GPP_NT","WTD","P_F","SWC"];
        var_dir = 'FCH4_F/';
    elseif (variable==3)
        var_wanted = ["Year","DOY","FCH4_F_ANN","TS_1","TS_2","TS_3","TA_F",...
            "GPP_DT","GPP_NT","WTD","P_F","SWC"];
        var_dir = 'FCH4_F_ANN/';
    end
    % ecosystem type
    ecosystem_list = ["bog", "drained", "fen", "marsh", "rice", "wet_tundra",...
        "lake", "swamp", "peat_plateau", "salt_marsh", "upland"];
    for landtype = 1:length(ecosystem_list)
        ecosystem = char(erase(ecosystem_list(landtype),'"'));
        di=folderFiles(['../' ecosystem '/'],'*.csv');
        len=size(di);
        daily_data = [];
        for site=1:len(1)
            QC_data = [];
            daily_tmp = [];
            filename=['../' ecosystem '/' di(site,:)];
            % read in texts and values as a table
            T = readtable(filename);
            % read in the data header
            fileID = fopen(filename);
            tline = fgets(fileID);
            fclose(fileID);
            tline = strsplit(tline,',');
            header = string(erase(tline,'"'));
            for variables = 1:length(var_wanted) % loop through variables
                if (sum(contains(header,var_wanted(variables)))>=1) 
                    % if the requested variable exists
                    idx = find(strcmpi(header,var_wanted(variables)));
                    tmp = T{:,idx};
                    % variable will be stored as a cell array if it mixes texts &
                    % numerics
                    if (iscell(tmp))
                        % find nan values in the array
                        idx_NA=find(contains(tmp,'NA'));
                        idx_full = logical(1:length(tmp))';
                        idx_full(idx_NA) = 0;
                        xxx = find(idx_full==1);
                        yyy = find(idx_full==0);
                        % assign the numeric and sting extracted from the cell
                        % array to a double array with nan
                        dim = size(tmp);
                        data = zeros(dim);
                        data(xxx) = cellfun(@str2num,tmp(xxx));
                        data(yyy) = nan;
                    else
                        data = tmp;
                    end
                    daily_tmp(1:length(T{:,idx}), variables) = data;
                else % the requested data is not available
                    daily_tmp(1:length(T{:,idx}), variables) = nan;
                end 
            end
            % plot time series
            for T_type = 1:2
                clf;
                if (T_type==1) % Tair
                    col = 7;
                    yname = {'T_{air}'; ['(' degree 'C)']};
                    yname_sub = ['T_{air} (' degree 'C)'];
                    figname = '_Tair_ts';
                    dir = 'Tair/';
                elseif (T_type==2) % Tsoil
                    col = 4;
                    yname = {'T_{soil}'; ['(' degree 'C)']};
                    yname_sub = ['T_{soil} (' degree 'C)'];
                    figname = '_Tsoil_ts';
                    dir = 'Tsoil/';
                end
                yr_range = min(daily_tmp(:,1)):max(daily_tmp(:,1));
                if (di(site,1:5)=='JPBBY')
                    yr_range(end) = [];

                end
                if (sum(isfinite((daily_tmp(:,col))))>1)
                    % T
                    subplot(6,length(yr_range),[1:length(yr_range)])
                    doy_label = daily_tmp(:,1)+daily_tmp(:,2)/365;
                    plot(doy_label,daily_tmp(:,col),...
                        'linewidth',2)
                    axis([yr_range(1), yr_range(end)+1, min(daily_tmp(:,col)), max(10,max(daily_tmp(:,col)))]);
                    xticks([yr_range(1):yr_range(end)+1])
                    ylabel(yname)
                    t = text(0.02,0.98,'(a)','Units', 'Normalized',...
                        'VerticalAlignment', 'Top','FontSize',12);
                    % FCH4
                    subplot(6,length(yr_range),[length(yr_range)+1:2*length(yr_range)])
                    plot(doy_label,daily_tmp(:,3),...
                        'linewidth',2)
                    axis([yr_range(1), yr_range(end)+1, 0, max(10,max(daily_tmp(:,3)))]);
                    xticks([yr_range(1):yr_range(end)+1])
                    ylabel({'{\it F_{CH_4}}'; '(mg C m^{-2} d^{-1})'})
                    t = text(0.02,0.98,'(b)','Units', 'Normalized',...
                        'VerticalAlignment', 'Top','FontSize',12);
                    % WTD & Precip                        
                    subplot(6,length(yr_range),[2*length(yr_range)+1:3*length(yr_range)])
                    % precip
                    yyaxis left
                    [tmp, ia, ic] = unique(doy_label);
                    bar(doy_label(ia), daily_tmp(ia,11))
                    axis([yr_range(1), yr_range(end)+1, 0, max(daily_tmp(:,11))]);
                    xticks([yr_range(1):yr_range(end)+1])
                    ylabel({'Precip'; '(mm)'})
                    t = text(0.02,0.98,'(c)','Units', 'Normalized',...
                        'VerticalAlignment', 'Top','FontSize',12);
                    % WTD
                    if (di(site,1:5)=='JPBBY')
                        WTD = daily_tmp(:,10)*1000;
                    else
                        WTD = daily_tmp(:,10)*100;
                    end
                    yyaxis right
                    plot(doy_label, WTD,'x')
                    axis([yr_range(1), yr_range(end)+1, min(-10,min(WTD)), max(0,max(WTD))]);
                    xticks([yr_range(1):yr_range(end)+1])
                    ylabel({'WTD' '(cm)'})

                    % hysteresis pattern?
                    % loop through years
                    for yr = 1:length(yr_range)
                        idx = find(daily_tmp(:,1)==yr_range(yr));
                        data = daily_tmp(idx,:);
                        % Tair/Tsoil
                        xxx = data(:,col);
                        idx = find(xxx<=threshold_T);
                        xxx(idx) = nan;
                        % FCH4
                        yyy = data(:,3);
                        idx = find(yyy<=threshold_FCH4);
                        yyy(idx) = nan;
                        if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                            idx_valid = isfinite(xxx)&isfinite(yyy);
                            xxx = xxx(idx_valid);
                            yyy = yyy(idx_valid);
                            if (length(xxx)>1)
                                idx_Tmax = find(xxx==max(xxx));
                                idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
                                idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
                                x0 = mean(xxx(idx_Tmax));
                                y0 = mean(yyy(idx_Tmax));
                                if (length(idx_warming)>1 & length(idx_cooling)>1)
                                    colormap('redblue')
                                    subplot(6,length(yr_range),[3*length(yr_range)+yr:length(yr_range):5*length(yr_range)+yr])
                                    plot(xxx(1),yyy(1),'w.');
                                    hold on
                                    c = linspace(10,1,length(yyy));
                                    scatter(xxx,yyy,15,c,'filled')
                                    % warming branch
                                    % quadratic, window average
                                    xtmp = xxx(idx_warming);
                                    ytmp = yyy(idx_warming);
                                    T_range = 0:delta_T:ceil(max(xxx));                                        
                                    rg_coef = polyfix(xtmp, ytmp, 2, ...
                                        [0,x0],[0,y0]);
                                    % predicted values
                                    fch4_warming = polyval(rg_coef,T_range);
                                    fit_line_w = plot(T_range,fch4_warming,'r-','linewidth',2);
                                    % apparent activation energy
                                    [yfit_total, Ea_w] = YD2014_fit(xtmp(:), ytmp(:),T_range);
                                    % cooling branch
                                    xtmp = xxx(idx_cooling);
                                    ytmp = yyy(idx_cooling);                                        
                                    rg_coef = polyfix(xtmp, ytmp, 2, ...
                                        [0,x0],[0,y0]);
                                    % predicted values
                                    fch4_cooling = polyval(rg_coef,T_range);
                                    fit_line_c = plot(T_range,fch4_cooling,'b-','linewidth',2);
                                    % apparent activation energy
                                    [yfit_total, Ea_c] = YD2014_fit(xtmp(:), ytmp(:),T_range);
                                    % full-season
                                    xtmp = xxx(:);
                                    ytmp = yyy(:);                                        
                                    rg_coef = polyfix(xtmp, ytmp, 2, ...
                                        [0,x0],[0,y0]);
                                    fff = polyval(rg_coef,T_range);
                                    fit_line_w = plot(T_range,fff,'k-','linewidth',1);
                                    axis([min(xxx) max(xxx) 0 max(yyy)])
                                    id = char(erase(id_name(yr),'"'));
                                    t = text(0.02,0.98,[id],...
                                        'Units', 'Normalized', 'VerticalAlignment', 'Top');
                                    if (di(site,1:5)=='JPBBY')
                                        axis([0 x0 0 180])
                                        delTACH4 = [mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                                        norm_hys_area = sum(fch4_cooling-fch4_warming)*delta_T/...
                                            (T_range(end)-T_range(1));
                                        desp = '{\it H_{\mu}} = ';
                                        t = text(0.02,0.91,[desp num2str(sprintf('%.1f',delTACH4))],...
                                            'Units', 'Normalized', 'VerticalAlignment', 'Top');
                                        desp = '{\it H_A} = ';
                                        t = text(0.02,0.83,[desp num2str(sprintf('%.1f',norm_hys_area))],...
                                            'Units', 'Normalized', 'VerticalAlignment', 'Top');
                                    end
                                    xlabel(yname_sub,'FontSize',12)
                                end
                            end
                        end
                    end
                    % colorbar
                    c = colorbar('Position', [0.935, 0.1, 0.02, .36]);
                    set(c,'YTick',[])
                    set(c,'YDir','reverse')
                    ax1 = axes('Position',[0 0 1 1],'Visible','off');
                    axes(ax1) % sets ax1 to current axes
                    text(0.91,0.50,['End date'],'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                        'FontSize',12)
                    text(0.89,0.075,['Start date'],'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                        'FontSize',12)
                    % universal Title
                        text(0.48,0.98,[di(site,1:2) '-' di(site,3:5)],...
                            'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                            'FontSize',12,'FontWeight','bold')
                    t = text(0.02,0.20,{'CH_4 emission, {\it F_{CH_4}}'},'Units', 'Normalized',...
                        'VerticalAlignment', 'Top','FontSize',12);
                    t.Rotation = 90;
                    t = text(0.05,0.23,{'(mg C m^{-2} d^{-1})'},'Units', 'Normalized',...
                        'VerticalAlignment', 'Top','FontSize',12);
                    t.Rotation = 90;
                    print('-dpng','-r300',['../plots/time_series/' var_dir dir ecosystem '/' di(site,1:5) figname '.jpeg']);
                end
            end
        end
    end
end

%% this section makes time series and scatter plots for SH sites at each site-year
% use warming and cooling branch to represent NZKop 
% read in data
clear all; 
close all;
clc;
threshold_FCH4 = 0;
threshold_T = 0;
delta_T = 0.1;
degree = sprintf('%c', char(176));
id_name = ["(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)"];
% variable names for analysis
for variable = 1:3
    if (variable==1)
        var_wanted = ["Year","DOY","FCH4","TS_1","TS_2","TS_3","TA",...
            "GPP_DT","GPP_NT","WTD","P_F","SWC"];
        var_dir = 'FCH4/';
    elseif (variable==2)
        var_wanted = ["Year","DOY","FCH4_F","TS_1","TS_2","TS_3","TA_F",...
            "GPP_DT","GPP_NT","WTD","P_F","SWC"];
        var_dir = 'FCH4_F/';
    elseif (variable==3)
        var_wanted = ["Year","DOY","FCH4_F_ANN","TS_1","TS_2","TS_3","TA_F",...
            "GPP_DT","GPP_NT","WTD","P_F","SWC"];
        var_dir = 'FCH4_F_ANN/';
    end
    % ecosystem type
    ecosystem_list = ["bog", "drained", "fen", "marsh", "rice", "wet_tundra",...
        "lake", "swamp", "peat_plateau", "salt_marsh", "upland"];
    for landtype = 1:length(ecosystem_list);
        ecosystem = char(erase(ecosystem_list(landtype),'"'));
        di=folderFiles(['../' ecosystem '/'],'*.csv');
        len=size(di);
        daily_data = [];
        for site=1:len(1)
            if (di(site,1:5)=='NZKop')
                QC_data = [];
                daily_tmp = [];
                filename=['../' ecosystem '/' di(site,:)];
                % read in texts and values as a table
                T = readtable(filename);
                % read in the data header
                fileID = fopen(filename);
                tline = fgets(fileID);
                fclose(fileID);
                tline = strsplit(tline,',');
                header = string(erase(tline,'"'));
                for variables = 1:length(var_wanted) % loop through variables
                    if (sum(contains(header,var_wanted(variables)))>=1) 
                        % if the requested variable exists
                        idx = find(strcmpi(header,var_wanted(variables)));
                        tmp = T{:,idx};
                        % variable will be stored as a cell array if it mixes texts &
                        % numerics
                        if (iscell(tmp))
                            % find nan values in the array
                            idx_NA=find(contains(tmp,'NA'));
                            idx_full = logical(1:length(tmp))';
                            idx_full(idx_NA) = 0;
                            xxx = find(idx_full==1);
                            yyy = find(idx_full==0);
                            % assign the numeric and sting extracted from the cell
                            % array to a double array with nan
                            dim = size(tmp);
                            data = zeros(dim);
                            data(xxx) = cellfun(@str2num,tmp(xxx));
                            data(yyy) = nan;
                        else
                            data = tmp;
                        end
                        daily_tmp(1:length(T{:,idx}), variables) = data;
                    else % the requested data is not available
                        daily_tmp(1:length(T{:,idx}), variables) = nan;
                    end 
                end
                % sorting Southern Hemishpere measurements by warming-cooling branches
                    years = [min(daily_tmp(:,1)):max(daily_tmp(:,1))];
                    data_SH = []; 
                    for yr = 1:length(years)
                        idx = find(daily_tmp(:,1)==years(yr));
                        tmp = daily_tmp(idx,:);
                        % Tair axis
                        idx_Tmin = find(tmp(:,7)==min(tmp(:,7)));
                        idx_Tmax = find(tmp(1:180,7)==max(tmp(1:180,7)));
                        % create dummy year labels
                        if (yr==1)
                            data_SH = [data_SH; tmp(idx_Tmin:end,:)];
                        elseif (yr==length(years))
                            tmp(1:idx_Tmin,1) = years(yr)-1;
                            data_SH = [data_SH; tmp(1:idx_Tmin,:)];
                        else
                            tmp(1:idx_Tmin,1) = years(yr)-1;
                            data_SH = [data_SH; tmp];
                        end
                    end
                % plot time series
                for T_type = 1:2
                    clf;
                    if (T_type==1) % Tair
                        col = 7;
                        yname = {'T_{air}'; ['(' degree 'C)']};
                        yname_sub = ['T_{air} (' degree 'C)'];
                        figname = '_Tair_ts';
                        dir = 'Tair/';
                    elseif (T_type==2) % Tsoil
                        col = 4;
                        yname = {'T_{soil}'; ['(' degree 'C)']};
                        yname_sub = ['T_{soil} (' degree 'C)'];
                        figname = '_Tsoil_ts';
                        dir = 'Tsoil/';
                    end
                    yr_range = min(daily_tmp(:,1)):max(daily_tmp(:,1));                    
                    if (sum(isfinite((daily_tmp(:,col))))>1)
                        % T
                        subplot(6,length(yr_range),[1:length(yr_range)])
                        doy_label = daily_tmp(:,1)+daily_tmp(:,2)/365;
                        plot(doy_label,daily_tmp(:,col),...
                            'linewidth',2)
                        axis([yr_range(1), yr_range(end)+1, min(daily_tmp(:,col)), max(10,max(daily_tmp(:,col)))]);
                        xticks([yr_range(1):yr_range(end)+1])
                        ylabel(yname)
                        t = text(0.02,0.98,'(a)','Units', 'Normalized',...
                            'VerticalAlignment', 'Top','FontSize',12);
                        % FCH4
                        subplot(6,length(yr_range),[length(yr_range)+1:2*length(yr_range)])
                        plot(doy_label,daily_tmp(:,3),...
                            'linewidth',2)
                        axis([yr_range(1), yr_range(end)+1, 0, max(10,max(daily_tmp(:,3)))]);
                        xticks([yr_range(1):yr_range(end)+1])
                        ylabel({'{\it F_{CH_4}}'; '(mg C m^{-2} d^{-1})'})
                        t = text(0.02,0.98,'(b)','Units', 'Normalized',...
                            'VerticalAlignment', 'Top','FontSize',12);
                        % WTD & Precip                        
                        subplot(6,length(yr_range),[2*length(yr_range)+1:3*length(yr_range)])
                        % precip
                        yyaxis left
                        [tmp, ia, ic] = unique(doy_label);
                        bar(doy_label(ia), daily_tmp(ia,11))
                        axis([yr_range(1), yr_range(end)+1, 0, max(daily_tmp(:,11))]);
                        xticks([yr_range(1):yr_range(end)+1])
                        ylabel({'Precip'; '(mm)'})
                        t = text(0.02,0.98,'(c)','Units', 'Normalized',...
                            'VerticalAlignment', 'Top','FontSize',12);
                        % WTD
                        if (di(site,1:5)=='JPBBY')
                            WTD = daily_tmp(:,10)*1000;
                        else
                            WTD = daily_tmp(:,10)*100;
                        end
                        yyaxis right
                        plot(doy_label, WTD,'x')
                        axis([yr_range(1), yr_range(end)+1, min(-10,min(WTD)), max(0,max(WTD))]);
                        xticks([yr_range(1):yr_range(end)+1])
                        ylabel({'WTD' '(cm)'})
                        % hysteresis pattern?
                        % loop through years
                        if (di(site,1:5)=='JPBBY' | di(site,1:5)=='NZKop')
                            yr_range(end) = [];
                        end
                        for yr = 1:length(yr_range)
                            idx = find(data_SH(:,1)==yr_range(yr));
                            data = data_SH(idx,:);
                            % Tair/Tsoil
                            xxx = data(:,col);
                            idx = find(xxx<=threshold_T);
                            xxx(idx) = nan;
                            % FCH4
                            yyy = data(:,3);
                            idx = find(yyy<=threshold_FCH4);
                            yyy(idx) = nan;
                            if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                                idx_valid = isfinite(xxx)&isfinite(yyy);
                                xxx = xxx(idx_valid);
                                yyy = yyy(idx_valid);
                                if (length(xxx)>1)
                                    idx_Tmax = find(xxx==max(xxx));
                                    idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
                                    idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
                                    x0 = mean(xxx(idx_Tmax));
                                    y0 = mean(yyy(idx_Tmax));
                                    if (length(idx_warming)>1 & length(idx_cooling)>1)
                                        colormap('redblue')
                                        subplot(6,length(yr_range),[3*length(yr_range)+yr:length(yr_range):5*length(yr_range)+yr])
                                        plot(xxx(1),yyy(1),'w.');
                                        hold on
                                        c = linspace(10,1,length(yyy));
                                        scatter(xxx,yyy,15,c,'filled')
                                        % warming branch
                                        % quadratic, window average
                                        xtmp = xxx(idx_warming);
                                        ytmp = yyy(idx_warming);
                                        T_range = 0:delta_T:ceil(max(xxx));
                                        rg_coef = polyfix(xtmp, ytmp, 2, ...
                                            [0,x0],[0,y0]);
                                        % predicted values
                                        fch4_warming = polyval(rg_coef,T_range);
                                        fit_line_w = plot(T_range,fch4_warming,'r-','linewidth',2);
                                        % cooling branch
                                        xtmp = xxx(idx_cooling);
                                        ytmp = yyy(idx_cooling);
                                        T_range = 0:delta_T:ceil(max(xxx));
                                        rg_coef = polyfix(xtmp, ytmp, 2, ...
                                            [0,x0],[0,y0]);
                                        % predicted values
                                        fch4_cooling = polyval(rg_coef,T_range);
                                        fit_line_c = plot(T_range,fch4_cooling,'b-','linewidth',2);
                                        % full-season
                                        xtmp = xxx(:);
                                        ytmp = yyy(:);
                                        T_range = 0:delta_T:ceil(max(xxx));
                                        rg_coef = polyfix(xtmp, ytmp, 2, ...
                                            [0,x0],[0,y0]);
                                        fff = polyval(rg_coef,T_range);
                                        fit_line_w = plot(T_range,fff,'k-','linewidth',1);
                                        axis([min(T_range) max(T_range) 0 max(yyy)])
                                        id = char(erase(id_name(yr),'"'));
                                        t = text(0.02,0.98,[id],...
                                            'Units', 'Normalized', 'VerticalAlignment', 'Top');
                                        if (di(site,1:5)=='NZKop')
                                            axis([0 x0 0 160])
                                            delTACH4 = [mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                                            norm_hys_area = sum(fch4_cooling-fch4_warming)*delta_T/...
                                                (T_range(end)-T_range(1))/max(abs([fch4_cooling, fch4_warming]));
                                            desp = '{\it H_{\mu}} = ';
                                            t = text(0.02,0.83,[desp num2str(sprintf('%.1f',delTACH4))],...
                                                'Units', 'Normalized', 'VerticalAlignment', 'Top');
                                            desp = '{\it H_A} = ';
                                            t = text(0.02,0.91,[desp num2str(sprintf('%.1f',norm_hys_area))],...
                                                'Units', 'Normalized', 'VerticalAlignment', 'Top');
                                        end
                                        xlabel(yname_sub,'FontSize',12)
                                    end
                                end
                            end
                        end
                        % colorbar
                        c = colorbar('Position', [0.935, 0.1, 0.02, .36]);
                        set(c,'YTick',[])
                        set(c,'YDir','reverse')
                        ax1 = axes('Position',[0 0 1 1],'Visible','off');
                        axes(ax1) % sets ax1 to current axes
                        text(0.91,0.50,['End date'],'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                            'FontSize',12)
                        text(0.89,0.075,['Start date'],'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                            'FontSize',12)
                        % universal Title
                            text(0.48,0.98,[di(site,1:2) '-' di(site,3:5)],...
                                'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                                'FontSize',12,'FontWeight','bold')
                        % universal Y-axis
                        t = text(0.02,0.20,{'CH_4 emission, {\it F_{CH_4}}'},'Units', 'Normalized',...
                            'VerticalAlignment', 'Top','FontSize',12);
                        t.Rotation = 90;
                        t = text(0.05,0.23,{'(mg C m^{-2} d^{-1})'},'Units', 'Normalized',...
                            'VerticalAlignment', 'Top','FontSize',12);
                        t.Rotation = 90;
                        print('-dpng','-r300',['../plots/time_series/' var_dir dir ecosystem '/' di(site,1:5) figname '.jpeg']);
                    end
                end
            end
        end
    end
end

%% this section makes time series and scatter plots for JYBBY at each site-year
% read in data
clear all; 
close all;
clc;
threshold_FCH4 = 0;
threshold_T = 0;
delta_T = 0.1;
degree = sprintf('%c', char(176));
id_name = ["(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"];
% variable names for analysis
for variable = 1:3
    if (variable==1)
        var_wanted = ["Year","DOY","FCH4","TS_1","TS_2","TS_3","TA",...
            "GPP_DT","GPP_NT","WTD","P_F","SWC"];
        var_dir = 'FCH4/';
    elseif (variable==2)
        var_wanted = ["Year","DOY","FCH4_F","TS_1","TS_2","TS_3","TA_F",...
            "GPP_DT","GPP_NT","WTD","P_F","SWC"];
        var_dir = 'FCH4_F/';
    elseif (variable==3)
        var_wanted = ["Year","DOY","FCH4_F_ANN","TS_1","TS_2","TS_3","TA_F",...
            "GPP_DT","GPP_NT","WTD","P_F","SWC"];
        var_dir = 'FCH4_F_ANN/';
    end
    % ecosystem type
    ecosystem_list = ["bog", "drained", "fen", "marsh", "rice", "wet_tundra",...
        "lake", "swamp", "peat_plateau", "salt_marsh", "upland"];
    for landtype = 1
        ecosystem = char(erase(ecosystem_list(landtype),'"'));
        di=folderFiles(['../' ecosystem '/'],'*.csv');
        len=size(di);
        daily_data = [];
        for site=1:len(1)
            if (di(site,1:5)=='JPBBY')
                QC_data = [];
                daily_tmp = [];
                filename=['../' ecosystem '/' di(site,:)];
                % read in texts and values as a table
                T = readtable(filename);
                % read in the data header
                fileID = fopen(filename);
                tline = fgets(fileID);
                fclose(fileID);
                tline = strsplit(tline,',');
                header = string(erase(tline,'"'));
                for variables = 1:length(var_wanted) % loop through variables
                    if (sum(contains(header,var_wanted(variables)))>=1) 
                        % if the requested variable exists
                        idx = find(strcmpi(header,var_wanted(variables)));
                        tmp = T{:,idx};
                        % variable will be stored as a cell array if it mixes texts &
                        % numerics
                        if (iscell(tmp))
                            % find nan values in the array
                            idx_NA=find(contains(tmp,'NA'));
                            idx_full = logical(1:length(tmp))';
                            idx_full(idx_NA) = 0;
                            xxx = find(idx_full==1);
                            yyy = find(idx_full==0);
                            % assign the numeric and sting extracted from the cell
                            % array to a double array with nan
                            dim = size(tmp);
                            data = zeros(dim);
                            data(xxx) = cellfun(@str2num,tmp(xxx));
                            data(yyy) = nan;
                        else
                            data = tmp;
                        end
                        daily_tmp(1:length(T{:,idx}), variables) = data;
                    else % the requested data is not available
                        daily_tmp(1:length(T{:,idx}), variables) = nan;
                    end 
                end
                % plot time series
                for T_type = 1:2
                    clf;
                    if (T_type==1) % Tair
                        col = 7;
                        yname = {'T_{air}'; ['(' degree 'C)']};
                        yname_sub = ['T_{air} (' degree 'C)'];
                        figname = '_Tair_ts';
                        dir = 'Tair/';
                    elseif (T_type==2) % Tsoil
                        col = 4;
                        yname = {'T_{soil}'; ['(' degree 'C)']};
                        yname_sub = ['T_{soil} (' degree 'C)'];
                        figname = '_Tsoil_ts';
                        dir = 'Tsoil/';
                    end
                    yr_range = min(daily_tmp(:,1)):max(daily_tmp(:,1));
                    if (di(site,1:5)=='JPBBY')
                        yr_range(end) = [];

                    end
                    if (sum(isfinite((daily_tmp(:,col))))>1)
                        % T
                        subplot(6,length(yr_range),[1:length(yr_range)])
                        doy_label = daily_tmp(:,1)+daily_tmp(:,2)/365;
                        plot(doy_label,daily_tmp(:,col),...
                            'linewidth',2)
                        axis([yr_range(1), yr_range(end)+1, min(daily_tmp(:,col)), max(10,max(daily_tmp(:,col)))]);
                        xticks([yr_range(1):yr_range(end)+1])
                        ylabel(yname)
                        t = text(0.02,0.98,'(a)','Units', 'Normalized',...
                            'VerticalAlignment', 'Top','FontSize',12);
                        % FCH4
                        subplot(6,length(yr_range),[length(yr_range)+1:2*length(yr_range)])
                        plot(doy_label,daily_tmp(:,3),...
                            'linewidth',2)
                        axis([yr_range(1), yr_range(end)+1, 0, max(10,max(daily_tmp(:,3)))]);
                        xticks([yr_range(1):yr_range(end)+1])
                        ylabel({'{\it F_{CH_4}}'; '(mg C m^{-2} d^{-1})'})
                        t = text(0.02,0.98,'(b)','Units', 'Normalized',...
                            'VerticalAlignment', 'Top','FontSize',12);
                        % WTD & Precip                        
                        subplot(6,length(yr_range),[2*length(yr_range)+1:3*length(yr_range)])
                        % precip
                        yyaxis left
                        [tmp, ia, ic] = unique(doy_label);
                        bar(doy_label(ia), daily_tmp(ia,11))
                        axis([yr_range(1), yr_range(end)+1, 0, max(daily_tmp(:,11))]);
                        xticks([yr_range(1):yr_range(end)+1])
                        ylabel({'Precip'; '(mm)'})
                        t = text(0.02,0.98,'(c)','Units', 'Normalized',...
                            'VerticalAlignment', 'Top','FontSize',12);
                        % WTD
                        if (di(site,1:5)=='JPBBY')
                            WTD = daily_tmp(:,10)*1000;
                        else
                            WTD = daily_tmp(:,10)*100;
                        end
                        yyaxis right
                        plot(doy_label, WTD,'x')
                        axis([yr_range(1), yr_range(end)+1, min(-10,min(WTD)), max(0,max(WTD))]);
                        xticks([yr_range(1):yr_range(end)+1])
                        ylabel({'WTD' '(cm)'})                        
                        % hysteresis pattern?
                        % loop through years
                        for yr = 1:length(yr_range)
                            idx = find(daily_tmp(:,1)==yr_range(yr));
                            data = daily_tmp(idx,:);
                            % Tair/Tsoil
                            xxx = data(:,col);
                            idx = find(xxx<=threshold_T);
                            xxx(idx) = nan;
                            % FCH4
                            yyy = data(:,3);
                            idx = find(yyy<=threshold_FCH4);
                            yyy(idx) = nan;
                            if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                                idx_valid = isfinite(xxx)&isfinite(yyy);
                                xxx = xxx(idx_valid);
                                yyy = yyy(idx_valid);
                                if (length(xxx)>1)
                                    idx_Tmax = find(xxx==max(xxx));
                                    idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
                                    idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
                                    x0 = mean(xxx(idx_Tmax));
                                    y0 = mean(yyy(idx_Tmax));
                                    if (length(idx_warming)>1 & length(idx_cooling)>1)
                                        colormap('redblue')
                                        subplot(6,length(yr_range),[3*length(yr_range)+yr:length(yr_range):5*length(yr_range)+yr])
                                        plot(xxx(1),yyy(1),'w.');
                                        hold on
                                        c = linspace(10,1,length(yyy));
                                        scatter(xxx,yyy,15,c,'filled')
                                        % warming branch
                                        % quadratic, window average
                                        xtmp = xxx(idx_warming);
                                        ytmp = yyy(idx_warming);
                                        T_range = 0:delta_T:ceil(max(xxx));
                                        rg_coef = polyfix(xtmp, ytmp, 2, ...
                                            [0,x0],[0,y0]);
                                        % predicted values
                                        fch4_warming = polyval(rg_coef,T_range);
                                        fit_line_w = plot(T_range,fch4_warming,'r-','linewidth',2);
                                        % cooling branch
                                        xtmp = xxx(idx_cooling);
                                        ytmp = yyy(idx_cooling);
                                        T_range = 0:delta_T:ceil(max(xxx));
                                        rg_coef = polyfix(xtmp, ytmp, 2, ...
                                            [0,x0],[0,y0]);
                                        % predicted values
                                        fch4_cooling = polyval(rg_coef,T_range);
                                        fit_line_c = plot(T_range,fch4_cooling,'b-','linewidth',2);
                                        % full-season
                                        xtmp = xxx(:);
                                        ytmp = yyy(:);
                                        T_range = 0:delta_T:ceil(max(xxx));
                                        rg_coef = polyfix(xtmp, ytmp, 2, ...
                                            [0,x0],[0,y0]);
                                        fff = polyval(rg_coef,T_range);
                                        fit_line_w = plot(T_range,fff,'k-','linewidth',1);
                                        axis([min(xxx) max(xxx) 0 max(yyy)])
                                        id = char(erase(id_name(yr),'"'));
                                        t = text(0.02,0.98,[id],...
                                            'Units', 'Normalized', 'VerticalAlignment', 'Top');
                                        if (di(site,1:5)=='JPBBY')
                                            axis([0 x0 0 180])
                                            delTACH4 = [mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                                            norm_hys_area = sum(fch4_cooling-fch4_warming)*delta_T/...
                                                (T_range(end)-T_range(1))/max(abs([fch4_cooling, fch4_warming]));
                                            desp = '{\it H_{\mu}} = ';
                                            t = text(0.02,0.83,[desp num2str(sprintf('%.1f',delTACH4))],...
                                                'Units', 'Normalized', 'VerticalAlignment', 'Top');
                                            desp = '{\it H_A} = ';
                                            t = text(0.02,0.91,[desp num2str(sprintf('%.1f',norm_hys_area))],...
                                                'Units', 'Normalized', 'VerticalAlignment', 'Top');
                                        end
                                        xlabel(yname_sub,'FontSize',12)
                                    end
                                end
                            end
                        end
                        % colorbar
                        c = colorbar('Position', [0.935, 0.1, 0.02, .36]);
                        set(c,'YTick',[])
                        set(c,'YDir','reverse')
                        ax1 = axes('Position',[0 0 1 1],'Visible','off');
                        axes(ax1) % sets ax1 to current axes
                        text(0.91,0.50,['End date'],'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                            'FontSize',12)
                        text(0.89,0.075,['Start date'],'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                            'FontSize',12)
                        % universal Title
                            text(0.48,0.98,[di(site,1:2) '-' di(site,3:5)],...
                                'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                                'FontSize',12,'FontWeight','bold')
                        t = text(0.02,0.20,{'CH_4 emission, {\it F_{CH_4}}'},'Units', 'Normalized',...
                            'VerticalAlignment', 'Top','FontSize',12);
                        t.Rotation = 90;
                        t = text(0.05,0.23,{'(mg C m^{-2} d^{-1})'},'Units', 'Normalized',...
                            'VerticalAlignment', 'Top','FontSize',12);
                        t.Rotation = 90;
                        print('-dpng','-r300',['../plots/time_series/' var_dir dir ecosystem '/' di(site,1:5) figname '.jpeg']);
                    end
                end
            end
        end
    end
end
