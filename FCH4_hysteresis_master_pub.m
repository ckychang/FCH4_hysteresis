clear all;clc;close all;
mon = [31,28,31,30,31,30,31,31,30,31,30,31];
mon_sum = cumsum(mon);
delimiterIn = ',';
degree = sprintf('%c', char(176));
id_name = ["(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)",...
    "(m)", "(n)", "(o)", "(p)", "(q)", "(r)"];
% 6/1-10/31
ana_time = mon_sum(5)+1:mon_sum(10);
% count the number of site years
count_site_year = 1;
Tair_ts = [];
Tair_d_ts = [];
Tsoil1_ts = [];
Tsoil1_d_ts = [];
Tsoil2_ts = [];
Tsoil2_d_ts = [];
FCH4_ts = [];
FCH4_d_ts = [];
nt_ts = [];
nt_d_ts = [];
hys_ts = [];
hys_no_IAV_ts = [];
aggregate = [];

% variable names for analysis
var_wanted = ["Year","DOY","FCH4","TS_1","TS_3","TS_5","TA",...
    "GPP_DT","GPP_NT","LAT","NEE","NEE_F_ANN","WTD","WTD_F","P"];
% ecosystem type
ecosystem_list = ["bog", "fen", "marsh", "peat_plateau", "rice",...
        "salt_marsh", "swamp", "wet_tundra"];
for landtype = 1:length(ecosystem_list);
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
        % include site number
        daily_tmp(:,length(var_wanted)+1) = site;
        % include ecosystem number
        daily_tmp(:,length(var_wanted)+2) = landtype;
        % sorting Southern Hemishpere measurements by warming-cooling branches
        if (daily_tmp(1,10)<0)
            years = [min(daily_tmp(:,1)):max(daily_tmp(:,1))];
            data_SH = []; 
            for yr = 1:length(years)
                idx = find(daily_tmp(:,1)==years(yr));
                tmp = daily_tmp(idx,:);
                % Tair axis
                idx_Tmin = find(tmp(:,7)==min(tmp(:,7)));
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
            daily_data = [daily_data; data_SH];
        else
            daily_data = [daily_data; daily_tmp];
        end
    end
    eval([eval('ecosystem') '=daily_data;']);
    aggregate = [aggregate; daily_data];
end

% convert FCH4 from nmol/m2/s to mg/m2/d
aggregate(:,3) = aggregate(:,3)*(10^-9)*12*(10^3)*86400;

% assign columns
lat_col = 10; % LAT
NEE_col = 11; % NEE
NEE_F_col = 12; % NEE_F
WTD_col = 13; % WTD
WTD_F_col = 14; % WTD_F
P_col = 15; % precip
site_col = 16; % site
land_col = 17; % ecosystem type

% store site IDs
ecosystem_list = ["bog", "fen", "marsh", "peat_plateau", "rice",...
        "salt_marsh", "swamp", "wet_tundra"];
for landtype = 1:length(ecosystem_list);
    ecosystem = char(erase(ecosystem_list(landtype),'"'));
    di=folderFiles(['../' ecosystem '/'],'*.csv');
    len=size(di);
    daily_data = [];
    for site=1:len(1)
        site_ID(site, landtype) = {di(site,1:5)};
    end
end
            
pause

%% this section plots GPP-Tair distribution
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
window_size = [0.1, 0.05, 0.025]; % window sizes for the normalized T
clc;
close all
clear fit_line_x fit_line_y   
kkk = 8.62*10^-5;
xrange = linspace(-42.35,-39,4);
% define the colorcode for each landtype
mcolor = [.02 .50 .00; .90 .17 .31; .0 .0039 1.00; .00 .50 1.00; 1.00 .57 .69;...
    .11 .30 .24; .16 .32 .75; .00 .75 1.00; .76 .60 .42; .45 .63 .76; .53 .33 .04; 0 0 0];
% assign legend names
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
% normal axis
xrange = linspace(0.5,27,4);
xrange2 = linspace(0.5,27,3);
threshold_gs = 0;
close all
for T_type = 1:1
    hys_IAV = [];
    T_ts = [];
    FCH4_ts = [];
    h_plot = [];
    t_plot = [];
    clear Ea_daily_all
    clf;
    if (T_type==1) % air temperature
        col = 7;
        data = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = 'Gross Primary Productivity, GPP (g C m^{-2} d^{-1})';
        figname = 'TA_GPP_scatter';
    elseif (T_type==2) % soil temperature
        col = 4;
        data = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = 'CH_4 emission (mg C m^{-2} d^{-1})';
        figname = 'TS_Ea_scatter_normal_axes';
    end
    for landtype = 1:length(ecosystem_list)
        tmp_data = [];
        desp = char(erase(ecosystem_name(landtype),'"'));
        id = char(erase(id_name(landtype),'"'));
        % loop through all ecosystem types
        if (landtype<length(ecosystem_list)+1)
            idx = find(data(:,land_col)==landtype);
            tmp = data(idx,:);
        else
            tmp = data;
        end
        % include data during the period that T>0 C
        % loop through all ecosystem types
        number_site = max(tmp(:,site_col));
        for site=1:number_site
            idx_site = find(tmp(:,site_col)==site);
            data_site = tmp(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                % Tair/Tsoil
                xxx = data_site_year(:,col);
                % GPP
                yyy = data_site_year(:,8);
                idx = find(yyy<=0.05*max(yyy));
                yyy(idx) = nan;
                zzz = [xxx, yyy];
                tmp_data = [tmp_data; zzz];
            end
        end
        % nan filter
        idx = isfinite(tmp_data(:,1))&isfinite(tmp_data(:,2));
        xxx = tmp_data(idx,1);
        yyy = tmp_data(idx,2);
        % regression lines
        subplot(4, 2, landtype)
        % data points
        plot(0, 0,'w.');
        hold on
        sss = dscatter(xxx, yyy);
        colormap('Hot')
        axis tight
        t = text(0.02,1.2,[id ' ' desp],'Units', 'Normalized', 'VerticalAlignment', 'Top');
    end
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    axes(ax1) % sets ax1 to current axes
    t = text(0.03,0.16,'Gross Primary Productivity, GPP (\mumol C m^{-2} s^{-1})',...
        'Units', 'Normalized',...
        'VerticalAlignment', 'Top','FontSize',14);
    t.Rotation = 90;
    % universal X-axis
    t = text(0.41,0.05,xname,'Units', 'Normalized',...
        'VerticalAlignment', 'Top','FontSize',14);
    print('-dpng','-r300',['../plots/test/' figname '.jpeg']);
end

%% 
% this section compares the temperature dependence models at each landtype:
% (full-season quadratic) vs Chang2019 (intra-seasonal quadratic)
% only make the plot if there is a warming-cooling pair
clc;
kkk = 8.62*10^-5;
threshold_gpp = 0; % GPP > 0
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
delta_T = [0.1, 0.5, 1];

close all
for T_type = 1:1
    tmp_legend = [];
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figheader = 'TA_FCH4_comp';
    elseif (T_type==2) % soil temperature
        col = 4;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figheader = 'TS_FCH4_comp';
    end
    clf;
        for landtype = 1:length(ecosystem_list)
            desp = char(erase(ecosystem_name(landtype),'"'));
            id = char(erase(id_name(landtype),'"'));
            count = 0;
            ymax = [];
            full = [];
            warm = [];
            cool = [];
            Ea_ts = [];
            % loop through all ecosystem types
            idx = find(source(:,land_col)==landtype);
            data = source(idx,:);
            number_site = max(data(:,site_col));
            for site=1:number_site
                idx_site = find(data(:,site_col)==site);
                data_site = data(idx_site,:);
                range_year = min(data_site(:,1)):max(data_site(:,1));
                for yr = 1:length(range_year)
                    count = count+1;
                    idx_year = find(data_site(:,1)==range_year(yr));
                    data_site_year = data_site(idx_year,:);
                    % Tair/Tsoil
                    xxx = data_site_year(:,col);
                    idx = find(xxx<=threshold_T);
                    xxx(idx) = nan;
                    % FCH4
                    yyy = data_site_year(:,3);
                    idx = find(yyy<=threshold_FCH4);
                    yyy(idx) = nan;
                    if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                        idx_valid = isfinite(xxx)&isfinite(yyy);
                        xxx = xxx(idx_valid);
                        yyy = yyy(idx_valid);
                        idx_Tmax = find(xxx==max(xxx));
                        idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
                        idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
                        x0 = mean(xxx(idx_Tmax));
                        y0 = mean(yyy(idx_Tmax));
                        for branch = 1:3
                            if (branch==1) % full
                                idx = 1:length(xxx);
                                tmp = [xxx(idx), yyy(idx)];
                                full = [full; tmp];
                            elseif (branch==2) % warm
                                idx = idx_warming;
                                mcolor = '-r';
                                tmp = [xxx(idx), yyy(idx)];
                                warm = [warm; tmp];
                            elseif (branch==3) % cooling
                                idx = idx_cooling;
                                mcolor = '-b';
                                tmp = [xxx(idx), yyy(idx)];
                                cool = [cool; tmp];
                            end
                            xtmp = xxx(idx);
                            ytmp = yyy(idx);
                            rg_coef = polyfix(xtmp, ytmp, 2, ...
                                    [0,x0],[0,y0]);
                            if (branch==2) % warm
                                rg_coef_warm = rg_coef;
                            elseif (branch==3) % cooling
                                rg_coef_cool = rg_coef;
                            end
                        end
                        % plot the quadratic fits in warming and cooling
                        % branches for each site-year
                        if (isfinite(rg_coef_warm)==1 & isfinite(rg_coef_cool)==1)
                            plot_x = 0:0.1:x0;
                            plot_y_warming = polyval(rg_coef_warm,plot_x);
                            plot_y_cooling = polyval(rg_coef_cool,plot_x);
                            subplot(4,2,landtype)
                            ppp_warming = plot(plot_x,plot_y_warming,'-r','linewidth',1);
                            ppp_warming.Color(4) = 0.2;
                            hold on
                            ppp_cooling = plot(plot_x,plot_y_cooling,'-b','linewidth',1);
                            ppp_cooling.Color(4) = 0.2;
                            tmp = [plot_y_warming(:); plot_y_cooling(:)];
                            ymax = [ymax; max(tmp)];
                        end                            
                    end
                end
            end
            % superimpose Arrhenius fit
            xxx = full(:,1);
            yyy = full(:,2);
            xtmp = 1./(xxx+273.15);
            ytmp = log(yyy);
            % regression lines
            if (sum(isfinite(xtmp))>=2 & sum(isfinite(ytmp))>=2)
                rg_coef = polyfit(xtmp, ytmp, 1);
                fff = polyval(rg_coef,xtmp);
                Ea_daily = -rg_coef(1)*kkk;
                Ar_x = 1./xtmp-273.15;
                Ar_y = exp(fff);
                [Ar_x, idx] = sort(Ar_x);
                subplot(4,2,landtype)
                Ar_fit = plot(Ar_x,Ar_y(idx),'k-','linewidth',3);
                uistack(Ar_fit,'top');
                if (landtype==2)
                    ylim([0 300])
                elseif (landtype==3)
                    ylim([0 500])
                else
                    ylim([0 max(ymax)])
                end
                xlim([0 max(Ar_x)])
                t = text(0.02,0.98,[id ' ' desp],...
                    'Units', 'Normalized', ...
                    'VerticalAlignment', 'Top');
                t2 = text(0.02,0.75,['E_a = ' sprintf('%1.2f', Ea_daily)]...
                    ,'Units', 'Normalized', 'VerticalAlignment', 'Top');
            end
        end
        ax1 = axes('Position',[0 0 1 1],'Visible','off');
        axes(ax1) % sets ax1 to current axes
        tmp = [ppp_warming, ppp_cooling, Ar_fit];
            legend(tmp,{'Earlier','Later','Static'},...
            'location','northeastoutside','orientation','horizontal',...
            'box','off','FontSize',12)
        % universal Y-axis
        t = text(0.03,0.3,'CH_4 emission, {\it F_{CH_4}} (mg C m^{-2} d^{-1})',...
            'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',14);
        t.Rotation = 90;
        % universal X-axis
        t = text(0.41,0.08,xname,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',14);
        figname = [figheader];
        pause
        print('-dpng','-r300',['../plots/test/' figname '.jpeg']);
end

%% 
% this section estimates the FCH4 errors resulting from using T dependence inferred from FCH4-T 
% measurements collected from different periods
% aggregate all landtypes together
% evaluate bias for each branch in each site-year
clc;
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
delta_T = [0.1, 0.5, 1];
err_warm_agg = [];
err_cool_agg = [];
close all
for T_type = 1:1
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TA_FCH4_error_comp_TD_season_agg';
    elseif (T_type==2) % soil temperature
        col = 4;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TS_FCH4_error_comp_TD_season_agg';
    end
    err_warm = [];
    err_cool = [];
    err_mean = [];
    err_std = [];
    bias_agg = [];
    rmse_agg = [];
    r2_agg = [];
        clf;
        for landtype = 1:length(ecosystem_list)
            desp = char(erase(ecosystem_name(landtype),'"'));
            id = char(erase(id_name(landtype),'"'));
            ymax = [];
            full = [];
            warm = [];
            cool = [];
            Ea_ts = [];            
            % loop through all ecosystem types
            idx = find(source(:,land_col)==landtype);
            data = source(idx,:);
            number_site = max(data(:,site_col));
            for site=1:number_site
                idx_site = find(data(:,site_col)==site);
                data_site = data(idx_site,:);
                range_year = min(data_site(:,1)):max(data_site(:,1));
                for yr = 1:length(range_year)
                    idx_year = find(data_site(:,1)==range_year(yr));
                    data_site_year = data_site(idx_year,:);
                    % Tair/Tsoil
                    xxx = data_site_year(:,col);
                    idx = find(xxx<=threshold_T);
                    xxx(idx) = nan;
                    % FCH4
                    yyy = data_site_year(:,3);
                    idx = find(yyy<=threshold_FCH4);
                    yyy(idx) = nan;
                    if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                        idx_valid = isfinite(xxx)&isfinite(yyy);
                        xxx = xxx(idx_valid);
                        yyy = yyy(idx_valid);
                        idx_Tmax = find(xxx==max(xxx));
                        idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
                        idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
                        x0 = mean(xxx(idx_Tmax));
                        y0 = mean(yyy(idx_Tmax));
                        rmse_wam = [];
                        bias_wam = [];
                        r2_wam = [];
                        rmse_cool = [];
                        bias_cool = [];
                        r2_cool = [];
                        rmse_full = [];
                        bias_full = [];
                        r2_full = [];
                        rmse = [];
                        bias = [];
                        r2 = [];
                        if (length(idx_warming)>2 & length(idx_cooling)>2)
                            for branch = 1:3
                                if (branch==3) % full
                                    idx = 1:length(xxx);
                                    tmp = [xxx(idx), yyy(idx)];
                                elseif (branch==1) % warm
                                    idx = idx_warming;
                                    tmp = [xxx(idx), yyy(idx)];
                                elseif (branch==2) % cooling
                                    idx = idx_cooling;
                                    tmp = [xxx(idx), yyy(idx)];
                                end
                                % calculate TD for individual sets of
                                % measurements
                                xtmp = xxx(idx);
                                ytmp = yyy(idx);
                                rg_coef_full = polyfix(xtmp, ytmp, 2, ...
                                    [0,x0],[0,y0]);
                                % evaluate model performance for T
                                % dependence estimated from measurements
                                % collected from different periods
                                for i = 1:3
                                    if (i==3) % full
                                        idx = 1:length(xxx);
                                    elseif (i==1) % warm
                                        idx = idx_warming;
                                    elseif (i==2) % cooling
                                        idx = idx_cooling;
                                    end
                                    fff = polyval(rg_coef_full,xxx(idx));
                                    % abosolute bias in each season
                                    bias = [bias, (sum((fff-yyy(idx)))/sum(yyy(idx)))*100];
                                    stats = CalcPerf(fff, yyy(idx));
                                    rmse = [rmse, stats.RMSE];
                                end
                            end
                            bias_agg = [bias_agg; bias];
                            rmse_agg = [rmse_agg; rmse];
                            r2_agg = [r2_agg; r2];
                        end
                    end
                end
            end
        end
        
        for variable = 1:2
            clf;
            if (variable==1)
                data_tmp = bias_agg;
                yname = {'CH_4 emission bias'; '(%)'};
                fig_type = '_bias';
            elseif (variable==2)
                data_tmp = rmse_agg;
                yname = {'CH_4 emission RMSE'; '(mg C m^{-2} d^{-1})'};
                fig_type = '_rmse';
            end
            for season = 1:3
                if (season==1)
                    idx = 1:3:7;
                    yrange = [-60 230; 0 170];
                    desp = 'Earlier season measurements';
                elseif (season==2)
                    idx = 2:3:8;
                    yrange = [-80 90; 0 170];
                    desp = 'Later season measurements';
                elseif (season==3)
                    idx = 3:3:9;
                    yrange = [-80 90; 0 170];
                    desp = 'Full season measurements';
                end
                id = char(erase(id_name(season),'"'));
                subplot(3, 1, season)
                boxplot(data_tmp(:,idx), 'width',.7)
                % indicate mean
                hold on
                plot([1:3], mean(data_tmp(:,idx),'omitnan'), 'ro', 'linewidth', 2)
                set(findobj(gca,'type','line'),'linew',2)
                set(gca,'xticklabel',{'Earlier', 'Later', 'Full-season'},'FontSize',12)
                ylabel([yname],'FontSize',12);
                ylim(yrange(variable,:))
                xlim([.5 3.5])
                mean_value = mean(data_tmp(:,idx),'omitnan');
                median_value = median(data_tmp(:,idx));
                std_value = std(data_tmp(:,idx));
                x_pos = linspace(0.11, 0.78, 3);
                for j = 1:3
                    t = text(x_pos(j),0.99,[num2str(sprintf('%2.1f',mean_value(j))) char(177) ...
                            num2str(sprintf('%2.1f',std_value(j)))],...
                            'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
                end
                t = text(0.02,1.2,[id ' ' desp],'Units', 'Normalized',...
                    'VerticalAlignment', 'Top','FontSize',12);
            end
            print('-dpng','-r300',['../plots/test/' figname fig_type '.jpeg']);
        end        
end

%% this section analyzes the sensitivity of FCH4 hysteresis to env factors
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
    "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", ...
    "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)"];
TA_hys_percent = [];
TS_hys_percent = [];
delta_T = 0.1;
close all
clc;
for T_type = 1:1
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        fig_head = 'FCH4_hys_climate';
    elseif (T_type==2) % soil temperature
        col = 4;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
    end
    clf;
    delTACH4 = [];
    hys_area_ts = [];
    norm_hys_area_ts = [];
    f_hys_ts = [];
    climate_ts = [];
    intra_T_ts = [];
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));            
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                % Tair/Tsoil
                xxx = data_site_year(:,col);
                idx = find(xxx<=threshold_T);
                xxx(idx) = nan;
                zzz = data_site_year(:,4);
                idx = find(zzz<=threshold_T);
                data_site_year(idx,4) = nan;
                yyy = data_site_year(:,3);
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
                        % store env factos: MAT, Precip T range, changes in intra-seasonal mean T,
                        % mean GPP, Total GPP, GRS length, LAT, WTD,
                        % changes in intra-seasonal mean WTD
                        climate = [mean(data_site_year(idx_valid,7), 'omitnan'),...
                            mean(data_site_year(idx_valid,4), 'omitnan'),...
                            max(data_site_year(idx_valid,7))-min(data_site_year(idx_valid,7)),...
                            mean(data_site_year(idx_valid,8), 'omitnan')*(10^-6)*12*(10^3)*86400/1000,...
                            mean(xxx(idx_cooling), 'omitnan')-mean(xxx(idx_warming), 'omitnan'),...
                            mean(data_site_year(idx_cooling,4), 'omitnan')-mean(data_site_year(idx_warming,4), 'omitnan'),...
                            length(xxx), data_site_year(1,10),...
                            mean(data_site_year(idx_valid,13), 'omitnan')*100,...
                            mean(data_site_year(idx_cooling,13), 'omitnan')-mean(data_site_year(idx_warming,13), 'omitnan')*100];
                        % temperature range recorded in the site-year
                        xfit = 0:delta_T:ceil(max(xxx));
                        if (length(idx_warming)>2 & length(idx_cooling)>2)
                            clear FCH4_store fch4_warming fch4_cooling
                            for branch = 1:2
                                if (branch==1) % warming
                                    idx = idx_warming;
                                elseif (branch==2) % cooling
                                    idx = idx_cooling;
                                end
                                % quadratic, window average
                                xtmp = xxx(idx);
                                ytmp = yyy(idx);
                                if (length(idx)>2)
                                    rg_coef = polyfix(xtmp, ytmp, 2, ...
                                        [0,x0],[0,y0]);
                                    if (branch==1) % warming
                                        fch4_warming = polyval(rg_coef, xfit);
                                    elseif (branch==2) % cooling
                                        fch4_cooling = polyval(rg_coef, xfit);
                                    end
                                end
                            end
                            hys_area = sum(fch4_cooling-fch4_warming)*delta_T;
                            hys_area_ts = [hys_area_ts; hys_area];
                            norm_hys_area_ts = [norm_hys_area_ts; ...
                                hys_area/(xfit(end)-xfit(1))/max(abs([fch4_cooling, fch4_warming]))];
                            delTACH4 = [delTACH4; mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                            f_hys = max(fch4_cooling-fch4_warming)/max(yyy)*100;
                            f_hys_ts = [f_hys_ts; f_hys];
                            climate_ts = [climate_ts; climate];
                        end
                    end
                end
            end
        end
    end
    % report the number of site-years in each metric
    num_siteyear(T_type) = length(hys_area_ts);
    for hys = 1:2
        if (hys==1)
            yyy = delTACH4;
            yname = '{\it H_\mu}';
            yrange = [-200 200];
            class = '_env_Hmu_v2';
        elseif (hys==2)
            yyy = norm_hys_area_ts;
            yname = '{\it H_A}';
            yrange = [-1 1];
            class = '_env_HA_v2';
        end
        clf;
        for variable = 1:8%length(climate_ts(1,:))
            if (variable==1)
                xname = ['Seasonal mean T_{air} (' degree 'C)'];
            elseif (variable==2)
                xname = ['Seasonal mean T_{soil} (' degree 'C)'];
%                 xname = ['Seasonal total precipitation (mm)'];
            elseif (variable==3)
                xname = ['Seasonal T_{air} range (' degree 'C)'];
            elseif (variable==4)
                xname = ['Seasonal mean GPP (g C m^{-2} d^{-1})'];
            elseif (variable==5)
                xname = ['Intra-seasonal changes in mean T_{air} (' degree 'C)'];
            elseif (variable==6)
                xname = ['Intra-seasonal changes in mean T_{soil} (' degree 'C)'];
            elseif (variable==7)
                xname = ['Frost-free season length (day)'];
            elseif (variable==8)
                xname = ['Latitude (' degree ')'];
            elseif (variable==9)
                xname = ['Seasonal mean WTD (cm)'];
            elseif (variable==10)
                xname = ['Intra-seasonal changes in mean WTD (cm)'];
            end
            subplot(4, 2, variable)
            idx = isfinite(climate_ts(:,variable))&isfinite(yyy);
            R2 = corr(climate_ts(idx,variable), yyy(idx))^2;
            % line fit
            rg_fit_tmp = polyfit(climate_ts(idx,variable), yyy(idx), 1);
            xxx = linspace(min(climate_ts(:,variable)),max(climate_ts(:,variable)),10);
            fff = polyval(rg_fit_tmp,xxx);
            plot(xxx,fff,'b-','Linewidth',2)
            hold on
            ppp = dscatter(climate_ts(idx,variable), yyy(idx));
            colormap('Hot')
            axis tight
            ylim(yrange)
            xlabel(xname, 'FontSize', 12);
            ylabel(yname, 'FontSize', 12);
            id = char(erase(id_name(variable),'"'));
            t = text(0.02,1.4,[id ' R^2 = ' num2str(sprintf('%1.2f',R2))],...
                'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                    'FontSize',12);
        end
        print('-dpng','-r300',['../plots/test/' fig_head class '.jpeg']);
    end
end

%% this section integrate GPP hysteresis using different metrics
% this section calculates mean GPP differences at each site-year and plot
% its distribution
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
TA_hys_percent = [];
TS_hys_percent = [];
delta_T = 0.1;
close all
for T_type = 1:2
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'GPP'; '(\mumol C m^{-2} s^{-1})'};
        figname = 'TA_GPP_hys_dist_no_win_agg';
    elseif (T_type==2) % soil temperature
        col = 4;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'GPP'; '(\mumol C m^{-2} s^{-1})'};
        figname = 'TS_GPP_hys_dist_no_win_agg';
    end
    clf;
    delTACH4 = [];
    hys_area_ts = [];
    norm_hys_area_ts = [];
    f_hys_ts = [];
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));            
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                % Tair/Tsoil
                xxx = data_site_year(:,col);
                idx = find(xxx<=threshold_T);
                xxx(idx) = nan;
                % GPP_DT
                yyy = data_site_year(:,8);
                % conver umol C m^{-2} s^{-1} to mgC/m2/d
                yyy = yyy*(10^-6)*12*(10^3)*86400;
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
                        % temperature range recorded in the site-year
                        xfit = 0:delta_T:ceil(max(xxx));
                        if (length(idx_warming)>2 & length(idx_cooling)>2)
                            clear FCH4_store fch4_warming fch4_cooling
                            for branch = 1:2
                                if (branch==1) % warming
                                    idx = idx_warming;
                                elseif (branch==2) % cooling
                                    idx = idx_cooling;
                                end
                                % quadratic, window average
                                xtmp = xxx(idx);
                                ytmp = yyy(idx);
                                if (length(idx)>2)
                                    rg_coef = polyfix(xtmp, ytmp, 2, ...
                                        [0,x0],[0,y0]);
                                    if (branch==1) % warming
                                        fch4_warming = polyval(rg_coef, xfit);
                                    elseif (branch==2) % cooling
                                        fch4_cooling = polyval(rg_coef, xfit);
                                    end
                                end
                            end
                            hys_area = sum(fch4_cooling-fch4_warming)*delta_T;
                            hys_area_ts = [hys_area_ts; hys_area];
                            norm_hys_area_ts = [norm_hys_area_ts; ...
                                hys_area/(xfit(end)-xfit(1))/max(abs([fch4_cooling, fch4_warming]))];
                            delTACH4 = [delTACH4; mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                            f_hys = max(fch4_cooling-fch4_warming)/max(yyy)*100;
                            f_hys_ts = [f_hys_ts; f_hys];
                        end
                    end
                end
            end
        end
    end
    % report the number of site-years in each metric
    num_siteyear(T_type) = length(hys_area_ts);
    clf;
    % make plot- delta FCH4/hys_are distribution
    for variable = 1:2
        if (variable==2)
            tmp = delTACH4;
            window = linspace(-6100, 6100, 100);
            xname = 'Mean seasonal GPP Hysteresis, {\it H_{\mu}} (mg C m^{-2} d^{-1})';
            id1 = '(c)';
            id2 = '(d)';
            plot_id1 = [6:8];
            plot_id2 = [9:10];
            
            fig_header = 'H_mu';
            symbol = '{\it H_{\mu}}';
        elseif (variable==1)
            tmp = norm_hys_area_ts;
            window = linspace(-1, 1, 100);
            id1 = '(a)';
            id2 = '(b)';
            plot_id1 = [1:3];
            plot_id2 = [4:5];
            xname = 'Normalized area of seasonal GPP Hysteresis, {\it H_A}';
            fig_header = 'H_A';
            symbol = '{\it H_A}';
        elseif (variable==3)
            tmp = f_hys_ts;
            window = linspace(-50, 100, 50);
            id = '(c)';
            xname = 'Hysteresis fraction (%)';
        end
        clear frequency
        for i = 2:length(window) % loop through Ahys windows
            idx = find(tmp>window(i-1) & tmp<=window(i));
            frequency(i) = length(idx);
        end
        pos_percent = length(find((tmp)>0))/length(tmp);
        if (T_type==1)
            TA_hys_percent = [TA_hys_percent; pos_percent];
        elseif (T_type==2)
            TS_hys_percent = [TS_hys_percent; pos_percent];
        end
        idx = isfinite(tmp);
        tmp = tmp(idx);
        % numerical mean and std
        [mu] = [mean(tmp)];
        [mu_std] = [std(tmp)];
        subplot(10, 1, plot_id1)
        bbb = bar(window,frequency,'k');
        bbb.FaceAlpha = 0.3;
        bbb.LineStyle = 'none';
        hold on
        t = text(0.01,0.98,id1,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
        % statistics
        text(0.01,0.85,['Positive' symbol ' in ' num2str(sprintf('%2.0f',pos_percent*100)) '% site-years' ],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        text(0.01,0.70,['Mean ' symbol ' = ' num2str(sprintf('%2.2f',mu))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % skewness
        skew = skewness(tmp);
        text(0.01,0.55,['Skewness ' ' = ' num2str(sprintf('%2.2f',skew))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % separate the positive and negative regimes
        if (variable==1)
            text(0.6,1.25,'Positive hysteresis','Units', 'Normalized', ...
                'VerticalAlignment', 'Top','Fontsize',12)
            text(0.12,1.25,'Negative hysteresis','Units', 'Normalized', ...
                'VerticalAlignment', 'Top','Fontsize',12)
            % arrow 1
            y1 = [.94 .94];
            x1 = [.13 .52];
            annotation('doublearrow',x1,y1,'Color','k');
            % arrow 2
            y1 = [.94 .94];
            x1 = [.52 .91];
            annotation('doublearrow',x1,y1,'Color','k');
        end
        % separate negative and positive hysteresis regimes
        plot([0,0],[0 max(frequency)],'r:','linewidth',2)
        ax = gca;
        ax.XLim = [min(window) max(window)];
        ax.YLim = [0, max(frequency)];
        ylabel('Frequency','Fontsize',12)
        set(gca,'xticklabel',[])
        plot_pos = subplot(10, 1, plot_id2);
        h= subplot('position', plot_pos.Position+[0.01, 0.08, -0.01, -0.08] ); 
        bbb = boxplot(tmp, 'width',.8, 'orientation', 'horizontal');
        set(findobj(gca,'type','line'),'linew',2)
        xlim([min(window) max(window)])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        ylabel([symbol '    '],'FontSize',12)
        xlh = xlabel([xname],'Fontsize',12);
        xlh.Position(2) = xlh.Position(2) + 0.1;
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle')
        t = text(0.02,0.95,id2,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);        
    end
    print('-dpng','-r300',['../plots/test/' figname 'agg' '.jpeg']);
end

%% 
% this section calculates mean FCH4 differences at each site-year and plot
% its relation to intra-seasonal changes in WTD
% stack bar plots for pdf distribution, with positive and negative WTD
clc;
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
TA_hys_percent = [];
TS_hys_percent = [];
delta_T = 0.1;
close all
for T_type = 1:1
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TA_hys_dist_WTD2';
    end
    clf;
    delTACH4 = [];
    hys_area_ts = [];
    norm_hys_area_ts = [];
    f_hys_ts = [];
    delta_WTD_ts = [];
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                % Tair/Tsoil
                xxx = data_site_year(:,col);
                idx = find(xxx<=threshold_T);
                xxx(idx) = nan;
                % FCH4
                yyy = data_site_year(:,3);
                idx = find(yyy<=threshold_FCH4);
                yyy(idx) = nan;
                % WTD in cm
                wtd = data_site_year(:,13)*100;
                if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                    idx_valid = isfinite(xxx)&isfinite(yyy);
                    xxx = xxx(idx_valid);
                    yyy = yyy(idx_valid);
                    wtd = wtd(idx_valid);
                    if (length(xxx)>1)
                        idx_Tmax = find(xxx==max(xxx));
                        idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
                        idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
                        x0 = mean(xxx(idx_Tmax));
                        y0 = mean(yyy(idx_Tmax));
                        % temperature range recorded in the site-year
                        xfit = 0:delta_T:ceil(max(xxx));
                        if (length(idx_warming)>2 & length(idx_cooling)>2)
                            clear FCH4_store fch4_warming fch4_cooling
                            for branch = 1:2
                                if (branch==1) % warming
                                    idx = idx_warming;
                                elseif (branch==2) % cooling
                                    idx = idx_cooling;
                                end
                                % quadratic, window average
                                xtmp = xxx(idx);
                                ytmp = yyy(idx);
%                                     ztmp = wtd(idx);
                                if (length(idx)>2)
                                    rg_coef = polyfix(xtmp, ytmp, 2, ...
                                        [0,x0],[0,y0]);
                                    if (branch==1) % warming
                                        fch4_warming = polyval(rg_coef, xfit);
                                        wtd_warming = mean(wtd(idx),'omitnan');
                                    elseif (branch==2) % cooling
                                        fch4_cooling = polyval(rg_coef, xfit);
                                        wtd_cooling = mean(wtd(idx),'omitnan');
                                    end
                                end
                            end
                            hys_area = sum(fch4_cooling-fch4_warming)*delta_T;
                            hys_area_ts = [hys_area_ts; hys_area];
                            norm_hys_area_ts = [norm_hys_area_ts; ...
                                hys_area/(xfit(end)-xfit(1))/max(abs([fch4_cooling, fch4_warming]))];
                            delTACH4 = [delTACH4; mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                            delta_WTD_ts = [delta_WTD_ts; (wtd_cooling-wtd_warming)];
                            f_hys = max(fch4_cooling-fch4_warming)/max(yyy)*100;
                            f_hys_ts = [f_hys_ts; f_hys];
                        end
                    end
                end
            end
        end
    end
    % report the number of site-years in each metric
    num_siteyear(T_type) = length(hys_area_ts);
    % separate delta WTD
    xname = 'Normalized area of seasonal CH_4 emission hysteresis, {\it H_A}';
    fig_header = 'H_A';
    symbol = '{\it H_A}';
    window = linspace(-1, 1, 100);
    clf;
    for variable = 1:2
        if (variable==1) % delta WTD>0
            id = '(a)';
            id2 = '(b)';
            desp = '\Delta WTD > 0';
            data = delta_WTD_ts;
            data(data<0) = nan;
            plot_a = 1:3;
            plot_b = 4:5;
        elseif (variable==2)
            id = '(c)';
            id2 = '(d)';
            desp = '\Delta WTD < 0';
            data = delta_WTD_ts;
            data(data>0) = nan;
            plot_a = 6:8;
            plot_b = 9:10;
        end
        % when WTD info is provided
        idx = isfinite(data)&isfinite(norm_hys_area_ts);
        xtmp = data(idx);
        ytmp = norm_hys_area_ts(idx);
        % when WTD info is provided
        tmp = squeeze(frequency_wtd(:,variable));
        pos_percent = length(find((ytmp)>0))/length(ytmp);
        idx = isfinite(ytmp);
        ytmp = ytmp(idx);
        mu = mean(ytmp);
        mu_std = std(ytmp);
        subplot(10,1,plot_a)
        bbb = bar(window,tmp);
        bbb.FaceColor = 'k';
        bbb.FaceAlpha = 0.3;
        bbb.LineStyle = 'none';
        ax = gca;
        ax.YLim = [0 max(tmp)];
        hold on
        % statistics
        % percent of postive WTD in negative hysteresis
        text(0.01,0.85,['Positive' symbol ' in ' num2str(sprintf('%2.0f',pos_percent*100)) '% site-years' ],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        text(0.01,0.70,['Mean ' symbol ' = ' num2str(sprintf('%2.2f',mu))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % skewness
        skew = skewness(ytmp);
        text(0.01,0.55,['Skewness ' ' = ' num2str(sprintf('%2.2f',skew))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)        
        t = text(0.01,0.98,[id desp],'Units', 'Normalized', 'VerticalAlignment', 'Top');
        % separate negative and positive hysteresis regimes
        plot([0,0],[0 max(tmp)],'r:','linewidth',2)
        ylabel('Frequency','Fontsize',12)
        set(gca,'xticklabel',[])
        
        plot_pos = subplot(10, 1, plot_b);
        h= subplot('position', plot_pos.Position+[0.01, 0.08, -0.01, -0.08] ); 
        bbb = boxplot(ytmp, 'width',.8, 'orientation', 'horizontal');
        set(findobj(gca,'type','line'),'linew',2)
        xlim([min(window) max(window)])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        ylabel([symbol '    '],'FontSize',12)
        xlh = xlabel([xname],'Fontsize',12);
        xlh.Position(2) = xlh.Position(2) + 0.1;
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle')
        t = text(0.02,0.95,[id2 desp],'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
    end
    % description axes
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    axes(ax1) % sets ax1 to current axes
    text(0.63,0.980,'Positive hysteresis','Units', 'Normalized', ...
            'VerticalAlignment', 'Top','Fontsize',12)
    text(0.24,0.980,'Negative hysteresis','Units', 'Normalized', ...
        'VerticalAlignment', 'Top','Fontsize',12)
    % arrow 1
    y1 = [.94 .94];
    x1 = [.13 .52];
    annotation('doublearrow',x1,y1,'Color','k');
    % arrow 2
    y1 = [.94 .94];
    x1 = [.52 .91];
    annotation('doublearrow',x1,y1,'Color','k');
    print('-dpng','-r300',['../plots/test/' figname '_sub.jpeg']);
end


%%
% this section calculates mean FCH4 differences at each site-year and plot
% its distribution
% use mean FCH4 during -0.5~0.5C as the end point, instead of 0 FCH4 at 0C
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
TA_hys_percent = [];
TS_hys_percent = [];
delta_T = 0.1;
close all
for T_type = 1:3
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TA_hys_dist_meanFCH4';
    elseif (T_type==2) % soil temperature
        col = 4;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TS_hys_dist_meanFCH4';
    elseif (T_type==3) % soil temperature
        col = 6;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TS3_hys_dist_meanFCH4';
    end
    clf;
    delTACH4 = [];
    hys_area_ts = [];
    norm_hys_area_ts = [];
    f_hys_ts = [];
    FCH4_min_ts = [];
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));
            
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                % FCH4 during -0.5~0.5 C
                xxx = data_site_year(:,col);
                yyy = data_site_year(:,3);
                idx = find(xxx>=-0.5&xxx<=0.5);
                FCH4_min = max(0, mean(yyy(idx),'omitnan'));
                FCH4_min_ts = [FCH4_min_ts; FCH4_min];
                % Tair/Tsoil
                xxx = data_site_year(:,col);
                idx = find(xxx<=threshold_T);
                xxx(idx) = nan;
                % FCH4
                yyy = data_site_year(:,3);
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
                        % temperature range recorded in the site-year
                        xfit = 0:delta_T:ceil(max(xxx));
                        if (length(idx_warming)>2 & length(idx_cooling)>2)
                            clear FCH4_store fch4_warming fch4_cooling
                            for branch = 1:2
                                if (branch==1) % warming
                                    idx = idx_warming;
                                elseif (branch==2) % cooling
                                    idx = idx_cooling;
                                end
                                % quadratic, window average
                                xtmp = xxx(idx);
                                ytmp = yyy(idx);
                                if (length(idx)>2)
                                    rg_coef = polyfix(xtmp, ytmp, 2, ...
                                        [0,x0],[FCH4_min,y0]);
                                    if (branch==1) % warming
                                        fch4_warming = polyval(rg_coef, xfit);
                                    elseif (branch==2) % cooling
                                        fch4_cooling = polyval(rg_coef, xfit);
                                    end
                                end
                            end
                            hys_area = sum(fch4_cooling-fch4_warming)*delta_T;
                            hys_area_ts = [hys_area_ts; hys_area];
                            norm_hys_area_ts = [norm_hys_area_ts; ...
                                hys_area/(xfit(end)-xfit(1))/max(abs([fch4_cooling, fch4_warming]))];
                            delTACH4 = [delTACH4; mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                            f_hys = max(fch4_cooling-fch4_warming)/max(yyy)*100;
                            f_hys_ts = [f_hys_ts; f_hys];
                        end
                    end
                end
            end
        end
    end
    % report the number of site-years in each metric
    num_siteyear(T_type) = length(hys_area_ts);
    % make plot- delta FCH4/hys_are distribution
    for variable = 1:2
        if (variable==1)
            tmp = delTACH4;
            window = linspace(-210, 211, 100);
            id = '(a)';
            fig_header = 'H_mu';
            symbol = '{\it H_{\mu}}';
            xname = 'Mean seasonal CH_4 emission hysteresis, {\it H_{\mu}} (mg C m^{-2} d^{-1})';
        elseif (variable==2)
            tmp = norm_hys_area_ts;
            window = linspace(-1, 1, 100);
            id = '(b)';
            xname = 'Normalized area of seasonal CH_4 emission hysteresis, {\it H_A}';
            fig_header = 'H_A';
            symbol = '{\it H_A}';
        elseif (variable==3)
            tmp = f_hys_ts;
            window = linspace(-50, 100, 50);
            id = '(c)';
            xname = 'Hysteresis fraction (%)';
        end
        clear frequency
        for i = 2:length(window) % loop through Ahys windows
            idx = find(tmp>window(i-1) & tmp<=window(i));
            frequency(i) = length(idx);
        end
        pos_percent = length(find((tmp)>0))/length(tmp);
        if (T_type==1)
            TA_hys_percent = [TA_hys_percent; pos_percent];
        elseif (T_type==2)
            TS_hys_percent = [TS_hys_percent; pos_percent];
        end
        idx = isfinite(tmp);
        tmp = tmp(idx);
        % numerical mean and std
        [mu] = [mean(tmp)];
        [mu_std] = [std(tmp)];
        clf;
        subplot(5, 1, 1:4)
        bbb = bar(window,frequency,'k');
        bbb.FaceAlpha = 0.3;
        bbb.LineStyle = 'none';
        hold on
        t = text(0.01,0.95,'(a)','Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
        % statistics
        text(0.01,0.90,['Positive' symbol ' in ' num2str(sprintf('%2.0f',pos_percent*100)) '% site-years' ],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        text(0.01,0.845,['Mean ' symbol ' = ' num2str(sprintf('%2.2f',mu))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % skewness
        skew = skewness(tmp);
        text(0.01,0.78,['Skewness ' ' = ' num2str(sprintf('%2.2f',skew))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % separate the positive and negative regimes
        text(0.63,1.1,'Positive hysteresis','Units', 'Normalized', ...
            'VerticalAlignment', 'Top','Fontsize',12)
        text(0.12,1.1,'Negative hysteresis','Units', 'Normalized', ...
            'VerticalAlignment', 'Top','Fontsize',12)
        % arrow 1
        y1 = [.94 .94];
        x1 = [.13 .52];
        annotation('doublearrow',x1,y1,'Color','k');
        % arrow 2
        y1 = [.94 .94];
        x1 = [.52 .91];
        annotation('doublearrow',x1,y1,'Color','k');
        % separate negative and positive hysteresis regimes
        plot([0,0],[0 max(frequency)],'r:','linewidth',3)
        ax = gca;
        ax.XLim = [min(window) max(window)];
        ax.YLim = [0, max(frequency)];
        ylabel('Frequency','Fontsize',12)
        subplot(5, 1, 5)
        boxplot(tmp, 'width',0.8, 'orientation', 'horizontal')
        set(findobj(gca,'type','line'),'linew',2)
        xlim([min(window) max(window)])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        ylabel([symbol '    '],'FontSize',12)
        xlh = xlabel([xname],'Fontsize',12);
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle')
        t = text(0.02,0.95,'(b)','Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
        set(findobj(gca,'type','line'),'linew',2)
        print('-dpng','-r300',['../plots/test/' figname fig_header '.jpeg']);
    end
end

%%
% this section calculates apparent activation energy for FCH4 differences
% include Ea vs GPP
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
TA_hys_percent = [];
TS_hys_percent = [];
delta_T = 0.1;
close all
for T_type = 1:1
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TA_FCH4_Ea_dist_';
    end
    clf;
    FCH4_Ea_ts = [];
    gpp_ts = [];
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));            
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                % Tair/Tsoil
                xxx = data_site_year(:,col);
                idx = find(xxx<=threshold_T);
                xxx(idx) = nan;
                % FCH4
                yyy = data_site_year(:,3);
                idx = find(yyy<=threshold_FCH4);
                yyy(idx) = nan;
                % GPP_DT
                gpp = data_site_year(:,8);
                % conver umol C m^{-2} s^{-1} to mgC/m2/d
                gpp = gpp*(10^-6)*12*86400;
                if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                    idx_valid = isfinite(xxx)&isfinite(yyy);
                    xxx = xxx(idx_valid);
                    yyy = yyy(idx_valid);
                    gpp = gpp(idx_valid);
                    if (length(xxx)>1)
                        idx_Tmax = find(xxx==max(xxx));
                        idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
                        idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
                        xfit = floor(min(xxx)):delta_T:ceil(max(xxx));
                        if (length(idx_warming)>2 & length(idx_cooling)>2)
                            Ea_tmp = [];
                            gpp_tmp = [];
                            for branch = 1:3
                                if (branch==1) % warming
                                    idx = idx_warming;
                                elseif (branch==2) % cooling
                                    idx = idx_cooling;
                                elseif (branch==3) % full-season
                                    idx = 1:length(xxx);
                                end
                                % Boltzman function for Ea
                                xtmp = xxx(idx);
                                ytmp = yyy(idx);                                    
                                if (length(idx)>2)
                                    [yfit_total, Ea] = YD2014_fit(xtmp(:), ytmp(:),xfit);
                                    Ea_tmp = [Ea_tmp, Ea];
                                    gpp_tmp = [gpp_tmp, mean(gpp(idx), 'omitnan')];
                                end
                            end
                            FCH4_Ea_ts = [FCH4_Ea_ts; Ea_tmp];
                            gpp_ts = [gpp_ts; gpp_tmp];
                        end
                    end
                end
            end
        end
    end
    % integrate fig
    clf;
    % Ea at individual periods
    median_value = median(FCH4_Ea_ts);
    std_value = std(FCH4_Ea_ts);
    subplot(2, 3, 1:3)
    boxplot(FCH4_Ea_ts, 'width',0.7, 'orientation', 'vertical')
    set(findobj(gca,'type','line'),'linew',2)
    ylim([-1 2.2])
    set(gca,'xticklabel',{'Earlier', 'Later', 'Full-season'}, ...
        'FontSize',12)
    ylh = ylabel({'CH_4 emission activation energy, \it{E_a}'; '(eV)'},'FontSize',12);
    x_pos = linspace(0.17, 0.76, 3);
    t = text(0.02,0.98,['(a)'],...
                'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
    for j = 1:3
        t = text(x_pos(j),0.99,[num2str(sprintf('%2.1f',median_value(j))) char(177) ...
                num2str(sprintf('%2.1f',std_value(j)))],...
                'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
    end
    % Ea vs GPP
    for i = 1:3
        if (i==1)
            id = '(b)';
            desp = 'Earlier';
        elseif (i==2)
            id = '(c)';
            desp = 'Later';
        elseif (i==3)
            id = '(d)';
            desp = 'Full-season';
        end
        idx = isfinite(gpp_ts(:,i))&isfinite(FCH4_Ea_ts(:,i));
        subplot(2, 3, 3+i)
        ppp = dscatter(gpp_ts(idx,i), FCH4_Ea_ts(idx,i));
        colormap('Hot')
        box on
        xlabel('GPP (g C m^{-2} d^{-1})', 'FontSize', 12);
        ylabel({'\it{E_a} (eV)'}, 'FontSize', 12);
        axis([0 15 -.2 3])
        R2 = corr(gpp_ts(idx,i), FCH4_Ea_ts(idx,i))^2;
        t = text(0.02,0.98,[id ' ' desp],...
                'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
        t = text(0.05,0.86,['R^2 = ', num2str(sprintf('%2.1f',R2))],...
                'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
    end
    print('-dpng','-r300',['../plots/test/TA_FCH4_Ea_dist_agg_v2.jpeg']);
end

%%
% this section calculates mean FCH4 differences at each site-year and plot
% its distribution
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
TA_hys_percent = [];
TS_hys_percent = [];
delta_T = 0.1;
close all
for T_type = 1:3
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TA_hys_dist_no_win';
    elseif (T_type==2) % soil temperature
        col = 4;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TS_hys_dist_no_win';
    elseif (T_type==3) % soil temperature
        col = 6;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TS3_hys_dist_no_win';
    end
    clf;
    delTACH4 = [];
    hys_area_ts = [];
    norm_hys_area_ts = [];
    f_hys_ts = [];
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));            
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                % Tair/Tsoil
                xxx = data_site_year(:,col);
                idx = find(xxx<=threshold_T);
                xxx(idx) = nan;
                % FCH4
                yyy = data_site_year(:,3);
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
                        % temperature range recorded in the site-year
                        xfit = 0:delta_T:ceil(max(xxx))
                        if (length(idx_warming)>2 & length(idx_cooling)>2)
                            % mean FCH4@cooling - mean FCH4@warming
                            delTACH4 = [delTACH4; mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                            clear FCH4_store fch4_warming fch4_cooling
                            for branch = 1:2
                                if (branch==1) % warming
                                    idx = idx_warming;
                                elseif (branch==2) % cooling
                                    idx = idx_cooling;
                                end
                                % quadratic fits
                                xtmp = xxx(idx); % T
                                ytmp = yyy(idx); % FCH4
                                if (length(idx)>2)
                                    rg_coef = polyfix(xtmp, ytmp, 2, ...
                                        [0,x0],[0,y0]);
                                    if (branch==1) % warming
                                        fch4_warming = polyval(rg_coef, xfit);
                                    elseif (branch==2) % cooling
                                        fch4_cooling = polyval(rg_coef, xfit);
                                    end
                                end
                            end
                            hys_area = sum(fch4_cooling-fch4_warming)*delta_T;
                            hys_area_ts = [hys_area_ts; hys_area];
                            norm_hys_area_ts = [norm_hys_area_ts; ...
                                hys_area/(xfit(end)-xfit(1))/max(abs([fch4_cooling, fch4_warming]))];
                            f_hys = max(fch4_cooling-fch4_warming)/max(yyy)*100;
                            f_hys_ts = [f_hys_ts; f_hys];
                        end
                    end
                end
            end
        end
    end
    % report the number of site-years in each metric
    num_siteyear(T_type) = length(hys_area_ts);
    % make plot- delta FCH4/hys_are distribution
    for variable = 1:2
        if (variable==1)
            tmp = delTACH4;
            window = linspace(-210, 210, 100);
            id = '(a)';
            fig_header = 'H_mu';
            symbol = '{\it H_{\mu}}';
            xname = 'Mean seasonal CH_4 emission hysteresis, {\it H_{\mu}} (mg C m^{-2} d^{-1})';
        elseif (variable==2)
            tmp = norm_hys_area_ts;
            window = linspace(-1, 1, 100);
            id = '(b)';
            xname = 'Normalized area of seasonal CH_4 emission hysteresis, {\it H_A}';
            fig_header = 'H_A';
            symbol = '{\it H_A}';
        end
        clear frequency
        for i = 2:length(window) % loop through Ahys windows
            idx = find(tmp>window(i-1) & tmp<=window(i));
            frequency(i) = length(idx);
        end
        pos_percent = length(find((tmp)>0))/length(tmp);
        if (T_type==1)
            TA_hys_percent = [TA_hys_percent; pos_percent];
        elseif (T_type==2)
            TS_hys_percent = [TS_hys_percent; pos_percent];
        end
        idx = isfinite(tmp);
        tmp = tmp(idx);
        % numerical mean and std
        [mu] = [mean(tmp)];
        [mu_std] = [std(tmp)];
        clf;
        subplot(5, 1, 1:4)
        bbb = bar(window,frequency,'k');
        bbb.FaceAlpha = 0.3;
        bbb.LineStyle = 'none';
        hold on
        t = text(0.01,0.95,'(a)','Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
        % statistics
        text(0.01,0.90,['Positive' symbol ' in ' num2str(sprintf('%2.0f',pos_percent*100)) '% site-years' ],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        text(0.01,0.845,['Mean ' symbol ' = ' num2str(sprintf('%2.2f',mu))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % skewness
        skew = skewness(tmp);
        text(0.01,0.78,['Skewness ' ' = ' num2str(sprintf('%2.2f',skew))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % separate the positive and negative regimes
        text(0.63,1.1,'Positive hysteresis','Units', 'Normalized', ...
            'VerticalAlignment', 'Top','Fontsize',12)
        text(0.12,1.1,'Negative hysteresis','Units', 'Normalized', ...
            'VerticalAlignment', 'Top','Fontsize',12)
        % arrow 1
        y1 = [.94 .94];
        x1 = [.13 .52];
        annotation('doublearrow',x1,y1,'Color','k');
        % arrow 2
        y1 = [.94 .94];
        x1 = [.52 .91];
        annotation('doublearrow',x1,y1,'Color','k');
        % separate negative and positive hysteresis regimes
        plot([0,0],[0 max(frequency)],'r:','linewidth',3)
        ax = gca;
        ax.XLim = [min(window) max(window)];
        ax.YLim = [0, max(frequency)];
        ylabel('Frequency','Fontsize',12)
        subplot(5, 1, 5)
        boxplot(tmp, 'width',0.8, 'orientation', 'horizontal')
        set(findobj(gca,'type','line'),'linew',2)
        xlim([min(window) max(window)])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        ylabel([symbol '    '],'FontSize',12)
        xlh = xlabel([xname],'Fontsize',12);
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle')
        t = text(0.02,0.95,'(b)','Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
        set(findobj(gca,'type','line'),'linew',2)
        print('-dpng','-r300',['../plots/test/' figname fig_header '.jpeg']);
    end  
end

%%
% this section calculates mean FCH4 differences at each site-year and plot
% its distribution
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
TA_hys_percent = [];
TS_hys_percent = [];
delta_T = 0.1;
close all
for T_type = 1:3
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TA_hys_dist_no_win';
    elseif (T_type==2) % soil temperature
        col = 4;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TS_hys_dist_no_win';
    elseif (T_type==3) % soil temperature
        col = 6;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TS3_hys_dist_no_win';
    end
    clf;
    delTACH4 = [];
    hys_area_ts = [];
    norm_hys_area_ts = [];
    f_hys_ts = [];
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));            
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                % Tair/Tsoil
                xxx = data_site_year(:,col);
                idx = find(xxx<=threshold_T);
                xxx(idx) = nan;
                % FCH4
                yyy = data_site_year(:,3);
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
                        % temperature range recorded in the site-year
                        xfit = 0:delta_T:ceil(max(xxx));
                        if (length(idx_warming)>2 & length(idx_cooling)>2)
                            % mean FCH4@cooling - mean FCH4@warming
                            delTACH4 = [delTACH4; mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                            clear FCH4_store fch4_warming fch4_cooling
                            for branch = 1:2
                                if (branch==1) % warming
                                    idx = idx_warming;
                                elseif (branch==2) % cooling
                                    idx = idx_cooling;
                                end
                                % quadratic fits
                                xtmp = xxx(idx); % T
                                ytmp = yyy(idx); % FCH4
                                if (length(idx)>2)
                                    rg_coef = polyfix(xtmp, ytmp, 2, ...
                                        [0,x0],[0,y0]);
                                    if (branch==1) % warming
                                        fch4_warming = polyval(rg_coef, xfit);
                                    elseif (branch==2) % cooling
                                        fch4_cooling = polyval(rg_coef, xfit);
                                    end
                                end
                            end
                            hys_area = sum(fch4_cooling-fch4_warming)*delta_T;
                            hys_area_ts = [hys_area_ts; hys_area];
                            norm_hys_area_ts = [norm_hys_area_ts; ...
                                hys_area/(xfit(end)-xfit(1))/max(abs([fch4_cooling, fch4_warming]))];
                            f_hys = max(fch4_cooling-fch4_warming)/max(yyy)*100;
                            f_hys_ts = [f_hys_ts; f_hys];
                        end
                    end
                end
            end
        end
    end
    % report the number of site-years in each metric
    num_siteyear(T_type) = length(hys_area_ts);
    % make plot- delta FCH4/hys_are distribution
    for variable = 1:2
        if (variable==2)
            tmp = delTACH4;
            window = linspace(-210, 210, 100);
            id1 = '(c)';
            id2 = '(d)';
            plot_id1 = [6:8];
            plot_id2 = [9:10];            
            fig_header = 'H_mu';
            symbol = '{\it H_{\mu}}';
            xname = 'Mean seasonal CH_4 emission hysteresis, {\it H_{\mu}} (mg C m^{-2} d^{-1})';
        elseif (variable==1)
            tmp = norm_hys_area_ts;
            window = linspace(-1, 1, 100);
            id1 = '(a)';
            id2 = '(b)';
            plot_id1 = [1:3];
            plot_id2 = [4:5];
            xname = 'Normalized area of seasonal CH_4 emission hysteresis, {\it H_A}';
            fig_header = 'H_A';
            symbol = '{\it H_A}';
        end
        clear frequency
        for i = 2:length(window) % loop through Ahys windows
            idx = find(tmp>window(i-1) & tmp<=window(i));
            frequency(i) = length(idx);
        end
        pos_percent = length(find((tmp)>0))/length(tmp);
        if (T_type==1)
            TA_hys_percent = [TA_hys_percent; pos_percent];
        elseif (T_type==2)
            TS_hys_percent = [TS_hys_percent; pos_percent];
        end
        idx = isfinite(tmp);
        tmp = tmp(idx);
        % numerical mean and std
        [mu] = [mean(tmp)];
        [mu_std] = [std(tmp)];
        plot_pos = subplot(10, 1, plot_id1);
        h= subplot('position', plot_pos.Position+[0, 0, 0, 0] ); 
        bbb = bar(window,frequency,'k');
        bbb.FaceAlpha = 0.3;
        bbb.LineStyle = 'none';
        hold on
        t = text(0.01,0.98,id1,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
        % statistics
        text(0.01,0.85,['Positive' symbol ' in ' num2str(sprintf('%2.0f',pos_percent*100)) '% site-years' ],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        text(0.01,0.70,['Mean ' symbol ' = ' num2str(sprintf('%2.2f',mu))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % skewness
        skew = skewness(tmp);
        text(0.01,0.55,['Skewness ' ' = ' num2str(sprintf('%2.2f',skew))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % separate the positive and negative regimes
        if (variable==1)
            text(0.6,1.25,'Positive hysteresis','Units', 'Normalized', ...
                'VerticalAlignment', 'Top','Fontsize',12)
            text(0.12,1.25,'Negative hysteresis','Units', 'Normalized', ...
                'VerticalAlignment', 'Top','Fontsize',12)
            % arrow 1
            y1 = [.94 .94];
            x1 = [.13 .52];
            annotation('doublearrow',x1,y1,'Color','k');
            % arrow 2
            y1 = [.94 .94];
            x1 = [.52 .91];
            annotation('doublearrow',x1,y1,'Color','k');
        end
        % separate negative and positive hysteresis regimes
        plot([0,0],[0 max(frequency)],'r:','linewidth',2)
        ax = gca;
        ax.XLim = [min(window) max(window)];
        ax.YLim = [0, max(frequency)];
        ylabel('Frequency','Fontsize',12)
        set(gca,'xticklabel',[])
        plot_pos = subplot(10, 1, plot_id2);
        h= subplot('position', plot_pos.Position+[0.01, 0.08, -0.01, -0.08] ); 
        bbb = boxplot(tmp, 'width',.8, 'orientation', 'horizontal');
        set(findobj(gca,'type','line'),'linew',2)
        xlim([min(window) max(window)])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        ylabel([symbol '    '],'FontSize',12)
        xlh = xlabel([xname],'Fontsize',12);
        xlh.Position(2) = xlh.Position(2) + 0.1;
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle')
        t = text(0.02,0.95,id2,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
    end
    print('-dpng','-r300',['../plots/test/' figname '_agg.jpeg']);  
end

%%
% this section calculates mean FCH4 differences at each site-year and plot
% its distribution **based on monthly mean estimates**
clc;
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
TA_hys_percent = [];
TS_hys_percent = [];
delta_T = 0.1;
close all
for T_type = 1:1
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TA_hys_dist_no_win';
    end
    clf;
    delTACH4 = [];
    hys_area_ts = [];
    norm_hys_area_ts = [];
    f_hys_ts = [];
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));            
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                % Tair/Tsoil
                xxx = data_site_year(:,col);
                idx = find(xxx<=threshold_T);
                xxx(idx) = nan;
                % FCH4
                yyy = data_site_year(:,3);
                idx = find(yyy<=threshold_FCH4);
                yyy(idx) = nan;
                % DOY
                zzz = data_site_year(:,2);
                if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                    idx_valid = isfinite(xxx)&isfinite(yyy);
                    xxx = xxx(idx_valid);
                    yyy = yyy(idx_valid);
                    zzz = zzz(idx_valid);
                    if (length(xxx)>1)
                        idx_Tmax = find(xxx==max(xxx));
                        idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
                        idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
                        x0 = mean(xxx(idx_Tmax));
                        y0 = mean(yyy(idx_Tmax));
                        xfit = 0:delta_T:ceil(max(xxx));
                        if (length(idx_warming)>2 & length(idx_cooling)>2)
                            clear FCH4_store fch4_warming fch4_cooling
                            for branch = 1:2
                                if (branch==1) % warming
                                    idx = idx_warming;
                                elseif (branch==2) % cooling
                                    idx = idx_cooling;
                                end
                                % collect data
                                xtmp_mon = [];
                                ytmp_mon = [];
                                xtmp = xxx(idx); % T
                                ytmp = yyy(idx); % FCH4
                                ztmp = zzz(idx); % DOY
                                % monthly mean
                                for m = 1:12
                                    if (m==1)
                                        idx_doy = find(ztmp<=mon_sum(1));
                                    else
                                        idx_doy = find(ztmp<=mon_sum(m)&ztmp>mon_sum(m-1));
                                    end
                                    if (length(idx_doy)>0)
                                        xtmp_mon = [xtmp_mon; mean(xtmp(idx_doy),'omitnan')];
                                        ytmp_mon = [ytmp_mon; mean(ytmp(idx_doy),'omitnan')];
                                    end
                                end
                                % quadratic fits
                                if (length(xtmp_mon)>=2)
                                    rg_coef = polyfix(xtmp_mon, ytmp_mon, 2, ...
                                        [0,x0],[0,y0]);
                                    if (branch==1) % warming
                                        fch4_warming = polyval(rg_coef, xfit);
                                        fch4_warming_mean = mean(ytmp_mon,'omitnan');
                                    elseif (branch==2) % cooling
                                        fch4_cooling = polyval(rg_coef, xfit);
                                        fch4_cooling_mean = mean(ytmp_mon,'omitnan');
                                    end
                                else
                                    if (branch==1) % warming
                                        fch4_warming = NaN;
                                        fch4_warming_mean = NaN;
                                    elseif (branch==2) % cooling
                                        fch4_cooling = NaN;
                                        fch4_cooling_mean = NaN;
                                    end
                                end
                            end
                            delTACH4 = [delTACH4; fch4_cooling_mean-fch4_warming_mean];
                            hys_area = sum(fch4_cooling-fch4_warming)*delta_T;
                            hys_area_ts = [hys_area_ts; hys_area];
                            norm_hys_area_ts = [norm_hys_area_ts; ...
                                hys_area/(xfit(end)-xfit(1))/max(abs([fch4_cooling, fch4_warming]))];
                            f_hys = max(fch4_cooling-fch4_warming)/max(yyy)*100;
                            f_hys_ts = [f_hys_ts; f_hys];
                        end
                    end
                end
            end
        end
    end
    % report the number of site-years in each metric
    num_siteyear(T_type) = length(hys_area_ts);
    for variable = 1:2
        if (variable==2)
            tmp = delTACH4;
            window = linspace(-210, 210, 100);
            id1 = '(c)';
            id2 = '(d)';
            plot_id1 = [6:8];
            plot_id2 = [9:10];            
            fig_header = 'H_mu';
            symbol = '{\it H_{\mu}}';
            xname = 'Mean seasonal CH_4 emission hysteresis, {\it H_{\mu}} (mg C m^{-2} d^{-1})';
        elseif (variable==1)
            tmp = norm_hys_area_ts;
            window = linspace(-1, 1, 100);
            id1 = '(a)';
            id2 = '(b)';
            plot_id1 = [1:3];
            plot_id2 = [4:5];
            xname = 'Normalized area of seasonal CH_4 emission hysteresis, {\it H_A}';
            fig_header = 'H_A';
            symbol = '{\it H_A}';
        end
        clear frequency
        for i = 2:length(window) % loop through Ahys windows
            idx = find(tmp>window(i-1) & tmp<=window(i));
            frequency(i) = length(idx);
        end
        pos_percent = length(find((tmp)>0))/length(tmp);
        if (T_type==1)
            TA_hys_percent = [TA_hys_percent; pos_percent];
        elseif (T_type==2)
            TS_hys_percent = [TS_hys_percent; pos_percent];
        end
        idx = isfinite(tmp);
        tmp = tmp(idx);
        % numerical mean and std
        [mu] = [mean(tmp)];
        [mu_std] = [std(tmp)];
        plot_pos = subplot(10, 1, plot_id1);
        h= subplot('position', plot_pos.Position+[0, 0, 0, 0] ); 
        bbb = bar(window,frequency,'k');
        bbb.FaceAlpha = 0.3;
        bbb.LineStyle = 'none';
        hold on
        t = text(0.01,0.98,id1,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
        % statistics
        text(0.01,0.85,['Positive' symbol ' in ' num2str(sprintf('%2.0f',pos_percent*100)) '% site-years' ],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        text(0.01,0.70,['Mean ' symbol ' = ' num2str(sprintf('%2.2f',mu))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % skewness
        skew = skewness(tmp);
        text(0.01,0.55,['Skewness ' ' = ' num2str(sprintf('%2.2f',skew))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % separate the positive and negative regimes
        if (variable==1)
            text(0.6,1.25,'Positive hysteresis','Units', 'Normalized', ...
                'VerticalAlignment', 'Top','Fontsize',12)
            text(0.12,1.25,'Negative hysteresis','Units', 'Normalized', ...
                'VerticalAlignment', 'Top','Fontsize',12)
            % arrow 1
            y1 = [.94 .94];
            x1 = [.13 .52];
            annotation('doublearrow',x1,y1,'Color','k');
            % arrow 2
            y1 = [.94 .94];
            x1 = [.52 .91];
            annotation('doublearrow',x1,y1,'Color','k');
        end
        % separate negative and positive hysteresis regimes
        plot([0,0],[0 max(frequency)],'r:','linewidth',2)
        ax = gca;
        ax.XLim = [min(window) max(window)];
        ax.YLim = [0, max(frequency)];
        ylabel('Frequency','Fontsize',12)
        set(gca,'xticklabel',[])
        plot_pos = subplot(10, 1, plot_id2);
        h= subplot('position', plot_pos.Position+[0.01, 0.08, -0.01, -0.08] ); 
        bbb = boxplot(tmp, 'width',.8, 'orientation', 'horizontal');
        set(findobj(gca,'type','line'),'linew',2)
        xlim([min(window) max(window)])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        ylabel([symbol '    '],'FontSize',12)
        xlh = xlabel([xname],'Fontsize',12);
        xlh.Position(2) = xlh.Position(2) + 0.1;
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle')
        t = text(0.02,0.95,id2,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
    end
    print('-dpng','-r300',['../plots/test/' figname '_mon_agg.jpeg']);
end

%% This section integrates hysteresis distribution for Tsoil metrics
% this section calculates mean FCH4 differences at each site-year and plot
% its distribution
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
TA_hys_percent = [];
TS_hys_percent = [];
delta_T = 0.1;
close all
for T_type = 2:3
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TA_hys_dist_no_win';
    elseif (T_type==2) % soil temperature
        col = 4;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TS_hys_dist_no_win_agg';
    elseif (T_type==3) % soil temperature
        col = 5;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TS3_hys_dist_no_win_agg';
    end
    clf;
    delTACH4 = [];
    hys_area_ts = [];
    norm_hys_area_ts = [];
    f_hys_ts = [];
    count = 0;
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                % Tair/Tsoil
                xxx = data_site_year(:,col);
                idx = find(xxx<=threshold_T);
                xxx(idx) = nan;
                % FCH4
                yyy = data_site_year(:,3);
                idx = find(yyy<=threshold_FCH4);
                yyy(idx) = nan;
                if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                    count = count+1;
                    idx_valid = isfinite(xxx)&isfinite(yyy);
                    xxx = xxx(idx_valid);
                    yyy = yyy(idx_valid);
                    if (length(xxx)>1)
                        idx_Tmax = find(xxx==max(xxx));
                            idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
                            idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
                        x0 = mean(xxx(idx_Tmax));
                        y0 = mean(yyy(idx_Tmax));
                        % temperature range recorded in the site-year
                        xfit = 0:delta_T:ceil(max(xxx));
                        if (length(idx_warming)>2 & length(idx_cooling)>2)
                            % mean FCH4@cooling - mean FCH4@warming
                            delTACH4 = [delTACH4; mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                            clear FCH4_store fch4_warming fch4_cooling
                            for branch = 1:2
                                if (branch==1) % warming
                                    idx = idx_warming;
                                elseif (branch==2) % cooling
                                    idx = idx_cooling;
                                end
                                % quadratic fits
                                xtmp = xxx(idx); % T
                                ytmp = yyy(idx); % FCH4
                                if (length(idx)>2)
                                    rg_coef = polyfix(xtmp, ytmp, 2, ...
                                        [0,x0],[0,y0]);
                                    if (branch==1) % warming
                                        fch4_warming = polyval(rg_coef, xfit);
                                    elseif (branch==2) % cooling
                                        fch4_cooling = polyval(rg_coef, xfit);
                                    end
                                end
                            end
                            hys_area = sum(fch4_cooling-fch4_warming)*delta_T;
                            hys_area_ts = [hys_area_ts; hys_area];
                            norm_hys_area_ts = [norm_hys_area_ts; ...
                                hys_area/(xfit(end)-xfit(1))/max(abs([fch4_cooling, fch4_warming]))];
                            f_hys = max(fch4_cooling-fch4_warming)/max(yyy)*100;
                            f_hys_ts = [f_hys_ts; f_hys];
                        end
                    end
                end
            end
        end
    end
    site_year_num(T_type) = count;
    % report the number of site-years in each metric
    num_siteyear(T_type) = length(hys_area_ts);
    clf;
    % make plot- delta FCH4/hys_are distribution
    for variable = 1:2
        if (variable==2)
            tmp = delTACH4;
            window = linspace(-210, 210, 100);
            id1 = '(c)';
            id2 = '(d)';
            plot_id1 = [6:8];
            plot_id2 = [9:10];            
            fig_header = 'H_mu';
            symbol = '{\it H_{\mu}}';
            xname = 'Mean seasonal CH_4 emission hysteresis, {\it H_{\mu}} (mg C m^{-2} d^{-1})';
        elseif (variable==1)
            tmp = norm_hys_area_ts;
            window = linspace(-1, 1, 100);
            id1 = '(a)';
            id2 = '(b)';
            plot_id1 = [1:3];
            plot_id2 = [4:5];
            xname = 'Normalized area of seasonal CH_4 emission hysteresis, {\it H_A}';
            fig_header = 'H_A';
            symbol = '{\it H_A}';
        end
        clear frequency
        for i = 2:length(window) % loop through Ahys windows
            idx = find(tmp>window(i-1) & tmp<=window(i));
            frequency(i) = length(idx);
        end
        pos_percent = length(find((tmp)>0))/length(tmp);
        if (T_type==1)
            TA_hys_percent = [TA_hys_percent; pos_percent];
        elseif (T_type==2)
            TS_hys_percent = [TS_hys_percent; pos_percent];
        end
        idx = isfinite(tmp);
        tmp = tmp(idx);
        % numerical mean and std
        [mu] = [mean(tmp)];
        [mu_std] = [std(tmp)];
        plot_pos = subplot(10, 1, plot_id1);
        h= subplot('position', plot_pos.Position+[0, 0, 0, 0] ); 
        bbb = bar(window,frequency,'k');
        bbb.FaceAlpha = 0.3;
        bbb.LineStyle = 'none';
        hold on
        t = text(0.01,0.98,id1,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
        % statistics
        text(0.01,0.85,['Positive' symbol ' in ' num2str(sprintf('%2.0f',pos_percent*100)) '% site-years' ],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        text(0.01,0.70,['Mean ' symbol ' = ' num2str(sprintf('%2.2f',mu))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % skewness
        skew = skewness(tmp);
        text(0.01,0.55,['Skewness ' ' = ' num2str(sprintf('%2.2f',skew))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % separate the positive and negative regimes
        if (variable==1)
            text(0.6,1.25,'Positive hysteresis','Units', 'Normalized', ...
                'VerticalAlignment', 'Top','Fontsize',12)
            text(0.12,1.25,'Negative hysteresis','Units', 'Normalized', ...
                'VerticalAlignment', 'Top','Fontsize',12)
            % arrow 1
            y1 = [.94 .94];
            x1 = [.13 .52];
            annotation('doublearrow',x1,y1,'Color','k');
            % arrow 2
            y1 = [.94 .94];
            x1 = [.52 .91];
            annotation('doublearrow',x1,y1,'Color','k');
        end
        % separate negative and positive hysteresis regimes
        plot([0,0],[0 max(frequency)],'r:','linewidth',2)
        ax = gca;
        ax.XLim = [min(window) max(window)];
        ax.YLim = [0, max(frequency)];
        ylabel('Frequency','Fontsize',12)
        set(gca,'xticklabel',[])
        plot_pos = subplot(10, 1, plot_id2);
        h= subplot('position', plot_pos.Position+[0.01, 0.08, -0.01, -0.08] ); 
        bbb = boxplot(tmp, 'width',.8, 'orientation', 'horizontal');
        set(findobj(gca,'type','line'),'linew',2)
        xlim([min(window) max(window)])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        ylabel([symbol '    '],'FontSize',12)
        xlh = xlabel([xname],'Fontsize',12);
        xlh.Position(2) = xlh.Position(2) + 0.1;
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle')
        t = text(0.02,0.95,id2,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
    end
    print('-dpng','-r300',['../plots/test/' figname 'agg' '.jpeg']);    
end

%%
% this section calculates mean FCH4 differences at each site-year and plot
% its distribution
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
TA_hys_percent = [];
TS_hys_percent = [];
delta_T = 0.1;
close all
for T_type = 1:1
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TA_hys_dist_no_win_gpp0_';
    elseif (T_type==2) % soil temperature
        col = 4;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TS_hys_dist_no_win_gpp0_';
    elseif (T_type==3) % soil temperature
        col = 6;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TS3_hys_dist_no_win_gpp0_';
    end
    clf;
    delTACH4 = [];
    hys_area_ts = [];
    norm_hys_area_ts = [];
    f_hys_ts = [];
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));            
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                threshold_gpp = 0;
                idx_gs = find(data_site_year(:,8)>threshold_gpp);
                if (length(idx_gs)>threshold_gs)
                    % Tair/Tsoil
                    xxx = data_site_year(idx_gs,col);
                    % FCH4
                    yyy = data_site_year(idx_gs,3);
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
                            % temperature range recorded in the site-year
                            xfit = 0:delta_T:ceil(max(xxx));
                            if (length(idx_warming)>2 & length(idx_cooling)>2)
                                clear FCH4_store fch4_warming fch4_cooling
                                for branch = 1:2
                                    if (branch==1) % warming
                                        idx = idx_warming;
                                    elseif (branch==2) % cooling
                                        idx = idx_cooling;
                                    end
                                    % quadratic, window average
                                    xtmp = xxx(idx);
                                    ytmp = yyy(idx);
                                    if (length(idx)>2)
                                        rg_coef = polyfix(xtmp, ytmp, 2, ...
                                            [0,x0],[0,y0]);
                                        if (branch==1) % warming
                                            fch4_warming = polyval(rg_coef, xfit);
                                        elseif (branch==2) % cooling
                                            fch4_cooling = polyval(rg_coef, xfit);
                                        end
                                    end
                                end
                                hys_area = sum(fch4_cooling-fch4_warming)*delta_T;
                                hys_area_ts = [hys_area_ts; hys_area];
                                norm_hys_area_ts = [norm_hys_area_ts; ...
                                    hys_area/(xfit(end)-xfit(1))/max(abs([fch4_cooling, fch4_warming]))];
                                delTACH4 = [delTACH4; mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                                f_hys = max(fch4_cooling-fch4_warming)/max(yyy)*100;
                                f_hys_ts = [f_hys_ts; f_hys];
                            end
                        end
                    end
                end
            end
        end
    end
    % report the number of site-years in each metric
    num_siteyear(T_type) = length(hys_area_ts);
    clf;
    % make plot- delta FCH4/hys_are distribution
    for variable = 1:2
        if (variable==2)
            tmp = delTACH4;
            window = linspace(-210, 210, 100);
            id1 = '(c)';
            id2 = '(d)';
            plot_id1 = [6:8];
            plot_id2 = [9:10];            
            fig_header = 'H_mu';
            symbol = '{\it H_{\mu}}';
            xname = 'Mean seasonal CH_4 emission hysteresis, {\it H_{\mu}} (mg C m^{-2} d^{-1})';
        elseif (variable==1)
            tmp = norm_hys_area_ts;
            window = linspace(-1, 1, 100);
            id1 = '(a)';
            id2 = '(b)';
            plot_id1 = [1:3];
            plot_id2 = [4:5];
            xname = 'Normalized area of seasonal CH_4 emission hysteresis, {\it H_A}';
            fig_header = 'H_A';
            symbol = '{\it H_A}';
        end
        clear frequency
        for i = 2:length(window) % loop through Ahys windows
            idx = find(tmp>window(i-1) & tmp<=window(i));
            frequency(i) = length(idx);
        end
        pos_percent = length(find((tmp)>0))/length(tmp);
        if (T_type==1)
            TA_hys_percent = [TA_hys_percent; pos_percent];
        elseif (T_type==2)
            TS_hys_percent = [TS_hys_percent; pos_percent];
        end
        idx = isfinite(tmp);
        tmp = tmp(idx);
        % numerical mean and std
        [mu] = [mean(tmp)];
        [mu_std] = [std(tmp)];
        plot_pos = subplot(10, 1, plot_id1);
        h= subplot('position', plot_pos.Position+[0, 0, 0, 0] ); 
        bbb = bar(window,frequency,'k');
        bbb.FaceAlpha = 0.3;
        bbb.LineStyle = 'none';
        hold on
        t = text(0.01,0.98,id1,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
        % statistics
        text(0.01,0.85,['Positive' symbol ' in ' num2str(sprintf('%2.0f',pos_percent*100)) '% site-years' ],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        text(0.01,0.70,['Mean ' symbol ' = ' num2str(sprintf('%2.2f',mu))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % skewness
        skew = skewness(tmp);
        text(0.01,0.55,['Skewness ' ' = ' num2str(sprintf('%2.2f',skew))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % separate the positive and negative regimes
        if (variable==1)
            text(0.6,1.25,'Positive hysteresis','Units', 'Normalized', ...
                'VerticalAlignment', 'Top','Fontsize',12)
            text(0.12,1.25,'Negative hysteresis','Units', 'Normalized', ...
                'VerticalAlignment', 'Top','Fontsize',12)
            % arrow 1
            y1 = [.94 .94];
            x1 = [.13 .52];
            annotation('doublearrow',x1,y1,'Color','k');
            % arrow 2
            y1 = [.94 .94];
            x1 = [.52 .91];
            annotation('doublearrow',x1,y1,'Color','k');
        end
        % separate negative and positive hysteresis regimes
        plot([0,0],[0 max(frequency)],'r:','linewidth',2)
        ax = gca;
        ax.XLim = [min(window) max(window)];
        ax.YLim = [0, max(frequency)];
        ylabel('Frequency','Fontsize',12)
        set(gca,'xticklabel',[])
        plot_pos = subplot(10, 1, plot_id2);
        h= subplot('position', plot_pos.Position+[0.01, 0.08, -0.01, -0.08] ); 
        bbb = boxplot(tmp, 'width',.8, 'orientation', 'horizontal');
        set(findobj(gca,'type','line'),'linew',2)
        xlim([min(window) max(window)])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        ylabel([symbol '    '],'FontSize',12)
        xlh = xlabel([xname],'Fontsize',12);
        xlh.Position(2) = xlh.Position(2) + 0.1;
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle')
        t = text(0.02,0.95,id2,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
    end
    print('-dpng','-r300',['../plots/test/' figname 'agg' '.jpeg']);
end

%%
% this section calculates mean FCH4 differences at each site-year and plot
% its distribution
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
TA_hys_percent = [];
TS_hys_percent = [];
delta_T = 0.1;
close all
for T_type = 1:1
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TA_hys_dist_no_win_gpp5_';
    end
    clf;
    delTACH4 = [];
    hys_area_ts = [];
    norm_hys_area_ts = [];
    f_hys_ts = [];
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));            
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                threshold_gpp = 0.05*max(data_site_year(:,8));
                idx_gs = find(data_site_year(:,8)>threshold_gpp);
                if (length(idx_gs)>threshold_gs & threshold_gpp>0)
                    % Tair/Tsoil
                    xxx = data_site_year(idx_gs,col);
                    % FCH4
                    yyy = data_site_year(idx_gs,3);
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
                            % temperature range recorded in the site-year
                            xfit = 0:delta_T:ceil(max(xxx));
                            if (length(idx_warming)>2 & length(idx_cooling)>2)
                                clear FCH4_store fch4_warming fch4_cooling
                                for branch = 1:2
                                    if (branch==1) % warming
                                        idx = idx_warming;
                                    elseif (branch==2) % cooling
                                        idx = idx_cooling;
                                    end
                                    % quadratic, window average
                                    xtmp = xxx(idx);
                                    ytmp = yyy(idx);
                                    if (length(idx)>2)
                                        rg_coef = polyfix(xtmp, ytmp, 2, ...
                                            [0,x0],[0,y0]);
                                        if (branch==1) % warming
                                            fch4_warming = polyval(rg_coef, xfit);
                                        elseif (branch==2) % cooling
                                            fch4_cooling = polyval(rg_coef, xfit);
                                        end
                                    end
                                end
                                hys_area = sum(fch4_cooling-fch4_warming)*delta_T;
                                hys_area_ts = [hys_area_ts; hys_area];
                                norm_hys_area_ts = [norm_hys_area_ts; ...
                                    hys_area/(xfit(end)-xfit(1))/max(abs([fch4_cooling, fch4_warming]))];
                                delTACH4 = [delTACH4; mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                                f_hys = max(fch4_cooling-fch4_warming)/max(yyy)*100;
                                f_hys_ts = [f_hys_ts; f_hys];
                            end
                        end
                    end
                end
            end
        end
    end
    % report the number of site-years in each metric
    num_siteyear(T_type) = length(hys_area_ts);
    clf;
    % make plot- delta FCH4/hys_are distribution
    for variable = 1:2
        if (variable==2)
            tmp = delTACH4;
            window = linspace(-210, 210, 100);
            id1 = '(c)';
            id2 = '(d)';
            plot_id1 = [6:8];
            plot_id2 = [9:10];            
            fig_header = 'H_mu';
            symbol = '{\it H_{\mu}}';
            xname = 'Mean seasonal CH_4 emission hysteresis, {\it H_{\mu}} (mg C m^{-2} d^{-1})';
        elseif (variable==1)
            tmp = norm_hys_area_ts;
            window = linspace(-1, 1, 100);
            id1 = '(a)';
            id2 = '(b)';
            plot_id1 = [1:3];
            plot_id2 = [4:5];
            xname = 'Normalized area of seasonal CH_4 emission hysteresis, {\it H_A}';
            fig_header = 'H_A';
            symbol = '{\it H_A}';
        end
        clear frequency
        for i = 2:length(window) % loop through Ahys windows
            idx = find(tmp>window(i-1) & tmp<=window(i));
            frequency(i) = length(idx);
        end
        pos_percent = length(find((tmp)>0))/length(tmp);
        if (T_type==1)
            TA_hys_percent = [TA_hys_percent; pos_percent];
        elseif (T_type==2)
            TS_hys_percent = [TS_hys_percent; pos_percent];
        end
        idx = isfinite(tmp);
        tmp = tmp(idx);
        % numerical mean and std
        [mu] = [mean(tmp)];
        [mu_std] = [std(tmp)];
        plot_pos = subplot(10, 1, plot_id1);
        h= subplot('position', plot_pos.Position+[0, 0, 0, 0] ); 
        bbb = bar(window,frequency,'k');
        bbb.FaceAlpha = 0.3;
        bbb.LineStyle = 'none';
        hold on
        t = text(0.01,0.98,id1,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
        % statistics
        text(0.01,0.85,['Positive' symbol ' in ' num2str(sprintf('%2.0f',pos_percent*100)) '% site-years' ],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        text(0.01,0.70,['Mean ' symbol ' = ' num2str(sprintf('%2.2f',mu))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % skewness
        skew = skewness(tmp);
        text(0.01,0.55,['Skewness ' ' = ' num2str(sprintf('%2.2f',skew))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontsize',12)
        % separate the positive and negative regimes
        if (variable==1)
            text(0.6,1.25,'Positive hysteresis','Units', 'Normalized', ...
                'VerticalAlignment', 'Top','Fontsize',12)
            text(0.12,1.25,'Negative hysteresis','Units', 'Normalized', ...
                'VerticalAlignment', 'Top','Fontsize',12)
            % arrow 1
            y1 = [.94 .94];
            x1 = [.13 .52];
            annotation('doublearrow',x1,y1,'Color','k');
            % arrow 2
            y1 = [.94 .94];
            x1 = [.52 .91];
            annotation('doublearrow',x1,y1,'Color','k');
        end
        % separate negative and positive hysteresis regimes
        plot([0,0],[0 max(frequency)],'r:','linewidth',2)
        ax = gca;
        ax.XLim = [min(window) max(window)];
        ax.YLim = [0, max(frequency)];
        ylabel('Frequency','Fontsize',12)
        set(gca,'xticklabel',[])
        plot_pos = subplot(10, 1, plot_id2);
        h= subplot('position', plot_pos.Position+[0.01, 0.08, -0.01, -0.08] ); 
        bbb = boxplot(tmp, 'width',.8, 'orientation', 'horizontal');
        set(findobj(gca,'type','line'),'linew',2)
        xlim([min(window) max(window)])        
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        ylabel([symbol '    '],'FontSize',12)
        xlh = xlabel([xname],'Fontsize',12);
        xlh.Position(2) = xlh.Position(2) + 0.1;
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle')
        t = text(0.02,0.95,id2,'Units', 'Normalized',...
            'VerticalAlignment', 'Top','FontSize',12);
    end
    print('-dpng','-r300',['../plots/test/' figname 'agg' '.jpeg']);
end

%% 
% this section calculates mean FCH4 differences at each site-year and plot
% its distribution with mean seasonal T/LAT/wetland type
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
delta_T = 0.1;
close all
for T_type = 1:2
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TA_hys_dist_no_win';
    elseif (T_type==2) % soil temperature
        col = 4;
        source = aggregate;
        xname = ['Soil temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TS_hys_dist_no_win';
    end
    clf;
    delTACH4 = [];
    hys_area_ts = [];
    norm_hys_area_ts = [];
    f_hys_ts = [];
    mst_ts = [];
    lat_ts = [];
    type_ts = [];
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                if (length(data_site_year)>1)
                    lat_site_year = data_site_year(1,lat_col);
                end
                % Tair/Tsoil
                xxx = data_site_year(:,col);
                idx = find(xxx<=threshold_T);
                xxx(idx) = nan;
                % FCH4
                yyy = data_site_year(:,3);
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
                        % temperature range recorded in the site-year
                        xfit = 0:delta_T:ceil(max(xxx));
                        if (length(idx_warming)>2 & length(idx_cooling)>2)
                            clear FCH4_store fch4_warming fch4_cooling
                            for branch = 1:2
                                if (branch==1) % warming
                                    idx = idx_warming;
                                elseif (branch==2) % cooling
                                    idx = idx_cooling;
                                end
                                % quadratic, window average
                                xtmp = xxx(idx);
                                ytmp = yyy(idx);
                                if (length(idx)>2)
                                    rg_coef = polyfix(xtmp, ytmp, 2, ...
                                        [0,x0],[0,y0]);
                                    if (branch==1) % warming
                                        fch4_warming = polyval(rg_coef, xfit);
                                    elseif (branch==2) % cooling
                                        fch4_cooling = polyval(rg_coef, xfit);
                                    end
                                end
                            end
                            hys_area = sum(fch4_cooling-fch4_warming)*delta_T;
                            hys_area_ts = [hys_area_ts; hys_area];
                            norm_hys_area_ts = [norm_hys_area_ts; ...
                                hys_area/(xfit(end)-xfit(1))/max(abs([fch4_cooling, fch4_warming]))];
                            delTACH4 = [delTACH4; mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                            f_hys = max(fch4_cooling-fch4_warming)/max(yyy)*100;
                            f_hys_ts = [f_hys_ts; f_hys];
                            mst_ts = [mst_ts; mean(xxx)];
                            lat_ts = [lat_ts; lat_site_year];
                            type_ts = [type_ts; landtype];
                        end
                    end
                end
            end
        end
    end
    % make plot- delta FCH4/hys_are distribution
    x_window = linspace(-1, 1, 100);
    xname = ['Normalized hysteresis area, {\it H_A}'];
    % assign categories
    mst_window = [1:5:40, 40];
    lat_window = [1:10:80, 80];
    for variable = 1:3
        close all
        if (variable==1)
            group = mst_ts;
            group_window = mst_window;
            group_name = ["1-5", "6-10", "11-15", "16-20", "21-25", "26-30", ...
                "31-35", "36-40"];
            group_unit = [degree 'C'];
            fig_header = '_MST';
        elseif (variable==2)
            group = lat_ts;
            group_window = lat_window;
            group_name = ["1~10", "11~20", "21~30", "31~40", "41~50", "51~60", ...
                "61~70", "71~80"];
            group_unit = [degree];
            fig_header = '_LAT';
        elseif (variable==3)
            group = type_ts;
            group_window = 1:9;
            group_name = ecosystem_name;
            group_unit = [' '];
            fig_header = '_TYPE';
        end
        for subgroup = 1:length(group_window)-1
            idx = find(group>=group_window(subgroup) & group<group_window(subgroup+1));
            tmp = norm_hys_area_ts(idx);
            clear frequency
            for i = 2:length(x_window) % loop through Ahys windows
                idx = find(tmp>x_window(i-1) & tmp<=x_window(i));
                frequency(i) = length(idx);
            end
            if (max(frequency)>=1) % distribution exists
                subplot(4,2,subgroup)
                bbb = bar(x_window,frequency,'k');
                bbb.FaceAlpha = 0.3;
                bbb.LineStyle = 'none';
                hold on
                % separate negative and positive hysteresis regimes
                plot([0,0],[0 max(frequency)],'r:','linewidth',2)
                idx = isfinite(tmp);
                tmp = tmp(idx);
                mu = mean(tmp);
                mu_std = std(tmp);
                ax = gca;
                ax.XLim = [min(x_window) max(x_window)];
                ax.YLim = [0, max(frequency)];
                desp = char(erase(group_name(subgroup),'"'));
                id = char(erase(id_name(subgroup),'"'));
                t = text(0.02,0.98,[id],...
                    'Units', 'Normalized', 'VerticalAlignment', 'Top');
                t = text(0.02,0.80,[desp group_unit],...
                    'Units', 'Normalized', 'VerticalAlignment', 'Top');
                ylabel('Frequency','Fontsize',12)
            end
        end
        ax1 = axes('Position',[0 0 1 1],'Visible','off');
        axes(ax1) % sets ax1 to current axes
        % universal X-axis
        if (T_type==1 & variable~=1)
            t = text(0.15,0.07,[xname],'Units', 'Normalized',...
                'VerticalAlignment', 'Top','FontSize',12);
            t = text(0.59,0.07,[xname],'Units', 'Normalized',...
                'VerticalAlignment', 'Top','FontSize',12);
        else
            t = text(0.15,0.27,[xname],'Units', 'Normalized',...
                'VerticalAlignment', 'Top','FontSize',12);
            t = text(0.59,0.27,[xname],'Units', 'Normalized',...
                'VerticalAlignment', 'Top','FontSize',12);
        end
        % description axes
        axes(ax1) % sets ax1 to current axes
        text(0.31,0.98,'Positive hysteresis','Units', 'Normalized', 'VerticalAlignment', 'Top')
        text(0.12,0.98,'Negative hysteresis','Units', 'Normalized', 'VerticalAlignment', 'Top')
        % arrow 1
        y1 = [.94 .94];
        x1 = [.13 .30];
        annotation('doublearrow',x1,y1,'Color','k');
        % arrow 2
        y2 = [.94 .94];
        x2 = [.30 .47];
        annotation('doublearrow',x2,y2,'Color','k');
        text(0.75,0.98,'Positive hysteresis','Units', 'Normalized', 'VerticalAlignment', 'Top')
        text(0.57,0.98,'Negative hysteresis','Units', 'Normalized', 'VerticalAlignment', 'Top')
        % arrow 1
        y3 = [.94 .94];
        x3 = [.57 .74];
        annotation('doublearrow',x3,y3,'Color','k');
        % arrow 2
        y4 = [.94 .94];
        x4 = [.74 .91];
        annotation('doublearrow',x4,y4,'Color','k');
        print('-dpng','-r300',['../plots/test/' figname fig_header '.jpeg']);
    end  
end

%%
% this section calculate the NRMSE with different temperature dependence
% scaling and compare NRMSE sensitivity to ecosystem type
% NRMSE comparason- 
% site-year specific TD vs ecosystem-specific TD vs aggregated TD
% use quadratic fit instead of Arrhenius
clc;
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];

id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
close all
for T_type = 1:1
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Ecosystem type'];
        yname = {'NRMSE'};
        figname_head = 'TA_NRMSE_comp_quadratic';
    elseif (T_type==2) % soil temperature
        col = 4;
        source = aggregate;
        xname = ['Ecosystem type'];
        yname = {'NRMSE'};
        figname_head = 'TS_NRMSE_comp_quadratic';
    end
    clf;
    clear rmse nrmse nrmse2 nrmse3 r2_value
    fch4_obs = [];
    fch4_est = [];
    data_QC = [];
%     % GPP threshold, only include data when gpp>5%of annual max GPP
%     for landtype = 1:length(ecosystem_name)
%         tmp_data = [];
%         
%         desp = char(erase(ecosystem_name(landtype),'"'));
%         id = char(erase(id_name(landtype),'"'));
%         % loop through all ecosystem types
%         if (landtype<length(ecosystem_name)+1)
%             idx = find(source(:,land_col)==landtype);
%             data = source(idx,:);
%         else
%             data = source;
%         end
%         % include data during the period that GPP gt 5% of annual max GPP
%         % loop through all ecosystem types
%         number_site = max(data(:,site_col));
%         for site=1:number_site
%             idx_site = find(data(:,site_col)==site);
%             data_site = data(idx_site,:);
%             range_year = min(data_site(:,1)):max(data_site(:,1));
%             for yr = 1:length(range_year)
%                 idx_year = find(data_site(:,1)==range_year(yr));
%                 data_site_year = data_site(idx_year,:);
%                 threshold_gpp = 0.05*max(data_site_year(:,8));
%                 idx_gs = find(data_site_year(:,8)>threshold_gpp);
%                 data_QC = [data_QC; data_site_year(idx_gs,:)];
%                 if (length(idx_gs)>threshold_gs & threshold_gpp>0)
%                     % Tair/Tsoil
%                     xxx = data_site_year(idx_gs,col);
% %                     idx = find(xxx<=threshold_T);
% %                     xxx(idx) = nan;
%                     % FCH4
%                     yyy = data_site_year(idx_gs,3);
%                     idx = find(yyy<=threshold_FCH4);
%                     yyy(idx) = nan;
%                     zzz = [xxx, yyy];
%                     tmp_data = [tmp_data; zzz];
%                 end
%             end
%         end
%     end
%         % nan filter
%         idx = isfinite(tmp_data(:,1))&isfinite(tmp_data(:,2));
%         xxx = tmp_data(idx,1);
%         yyy = tmp_data(idx,2);
% %     tmp = source(:,8);
% %     idx_gs = find(tmp<=0);
    idx = isfinite(source(:,col))&isfinite(source(:,3));
    % aggregate_Ea_filterd TD
    xxx = source(idx,col);
    yyy = source(idx,3);
%     zzz = source(idx,8);
% %     % growing season
% %     xxx(idx_gs) = nan;
% %     yyy(idx_gs) = nan;
%     valid T & FCH4
    idx = find(xxx<=threshold_T | yyy<=threshold_FCH4);
    xxx(idx) = [];
    yyy(idx) = [];
    xtmp = 1./(xxx+273.15);
    ytmp = log(yyy);
    rg_fit = polyfit(xtmp, ytmp, 1);
    for landtype = 1:length(ecosystem_list)
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(source(:,site_col));
        for scale = 1:6
            % decreased spatial and temporal description with increased scale#
            store = [];
            fff = [];
            if (scale==1)
                % scale=1, site-specific temperature dependence, 
                %           each site, each year, each season
                for site=1:number_site
                    idx_site = find(data(:,site_col)==site);
                    data_site = data(idx_site,:);
                    range_year = min(data_site(:,1)):max(data_site(:,1));
                    for yr = 1:length(range_year)
                        idx_year = find(data_site(:,1)==range_year(yr));
                        data_site_year = data_site(idx_year,:);
                        % Tair/Tsoil
                        xxx = data_site_year(:,col);
                        idx = find(xxx<=threshold_T);
                        xxx(idx) = nan;
                        % FCH4
                        yyy = data_site_year(:,3);
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
                                    % warming branch
                                    rg_fit_tmp = polyfix(xxx(1:idx_Tmax), yyy(1:idx_Tmax), 2, ...
                                        [0,x0],[0,y0]);
                                    % predicted values
                                    fff = polyval(rg_fit_tmp,xxx(1:idx_Tmax));
                                    % cooling branch
                                    rg_fit_tmp = polyfix(xxx(idx_Tmax:end), yyy(idx_Tmax:end), 2, ...
                                        [0,x0],[0,y0]);
                                    % predicted values
                                    fff = [fff; polyval(rg_fit_tmp,xxx(idx_Tmax+1:end))];
                                    store_tmp = [xxx, yyy, fff];
                                    store = [store; store_tmp];
                                    fch4_obs = [fch4_obs; yyy];
                                    fch4_est = [fch4_est; fff];
                                end
                            end
                        end
                    end
                end
            elseif (scale==2)
                % scale=2, site-specific temperature dependence, 
                %           each site, each year, all season
                for site=1:number_site
                    idx_site = find(data(:,site_col)==site);
                    data_site = data(idx_site,:);
                    range_year = min(data_site(:,1)):max(data_site(:,1));
                    for yr = 1:length(range_year)
                        idx_year = find(data_site(:,1)==range_year(yr));
                        data_site_year = data_site(idx_year,:);
                        % Tair/Tsoil
                        xxx = data_site_year(:,col);
                        idx = find(xxx<=threshold_T);
                        xxx(idx) = nan;
                        % FCH4
                        yyy = data_site_year(:,3);
                        idx = find(yyy<=threshold_FCH4);
                        yyy(idx) = nan;
                        if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                            idx_valid = isfinite(xxx)&isfinite(yyy);
                            xxx = xxx(idx_valid);
                            yyy = yyy(idx_valid);
                            if (length(xxx)>1)
                                rg_fit_tmp = polyfix(xxx(:), yyy(:), 2, ...
                                    [0],[0]);
                                % predicted values
                                fff = polyval(rg_fit_tmp,xxx(:));
                                store_tmp = [xxx, yyy, fff];
                                store = [store; store_tmp];
                            end
                        end
                    end
                end
            elseif (scale==3)
                % scale=3, site-specific temperature dependence , 
                %           each site, all year, all season
                for site=1:number_site
                    idx_site = find(data(:,site_col)==site);
                    data_site = data(idx_site,:);
                    range_year = min(data_site(:,1)):max(data_site(:,1));
                    % Tair/Tsoil
                    xxx = data_site(:,col);
                    idx = find(xxx<=threshold_T);
                    xxx(idx) = nan;
                    % FCH4
                    yyy = data_site(:,3);
                    idx = find(yyy<=threshold_FCH4);
                    yyy(idx) = nan;
                    if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                        idx_valid = isfinite(xxx)&isfinite(yyy);
                        xxx = xxx(idx_valid);
                        yyy = yyy(idx_valid);
                        if (length(xxx)>1)
                            rg_fit_tmp = polyfix(xxx(:), yyy(:), 2, 0, 0);
                            % predicted values
                            fff = polyval(rg_fit_tmp,xxx(:));
                            store_tmp = [xxx, yyy, fff];
                            store = [store; store_tmp];
                        end
                    end
                end
            elseif (scale==4)
                % scale=4, ecosystem-specific temperature dependence , 
                %           all sites within a ecosystem type, all year, all season
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
                        rg_fit_tmp = polyfix(xxx(:), yyy(:), 2, 0, 0);
                        % predicted values
                        fff = polyval(rg_fit_tmp,xxx(:));
                        store_tmp = [xxx, yyy, fff];
                        store = [store; store_tmp];
                    end
                end
            elseif (scale==5)
                % scale=5, aggregated temperature dependence , 
                %           all ecosystem types, all year, all season
                % quadratic fit
                idx = isfinite(data(:,col))&isfinite(data(:,3));
                xxx = data(idx,col);
                yyy = data(idx,3);
                % valid T & FCH4
                idx = find(xxx<=threshold_T | yyy<=threshold_FCH4);
                xxx(idx) = [];
                yyy(idx) = [];
                if (length(xxx)>1)
                    rg_fit_tmp = polyfix(xxx(:), yyy(:), 2, 0, 0);
                    % predicted values
                    fff = polyval(rg_fit_tmp,xxx(:));
                    store_tmp = [xxx, yyy, fff];
                    store = [store; store_tmp];
                end
            elseif (scale==6)
                % scale=6, aggregated temperature dependence , 
                %           all ecosystem types, all year, all season
                idx = isfinite(data(:,col))&isfinite(data(:,3));
                xxx = data(idx,col);
                yyy = data(idx,3);
                % valid T & FCH4
                idx = find(xxx<=threshold_T | yyy<=threshold_FCH4);
                xxx(idx) = [];
                yyy(idx) = [];
                if (length(xxx)>1)
                    % convert T and FCH4 to right units for Ea calculation
                    xxx = 1./(xxx+273.15);
                    yyy = log(yyy);
                    % predicted values
                    fff = polyval(rg_fit,xxx(:));
                    store_tmp = [xxx, exp(yyy), exp(fff)];
                    store = [store; store_tmp];
                end
            end
            if (length(store)>1) % the requested data exists
                % NRMSE for CH4 flux predicted by ecosystem specific temperature dependence
                stats = CalcPerf((store(:,2)), (store(:,3)));
                rmse(landtype,scale,1) = stats.RMSE;
                nrmse(landtype,scale,1) = stats.NRMSE*100;
                nrmse2(landtype,scale,1) = stats.NRMSE2*100;
                nrmse3(landtype,scale,1) = stats.NRMSE3*100;
                % R2 values
                yresid = store(:,2) - store(:,3);
                SSresid = sum(yresid.^2);
                SStotal = (length(store(:,2))-1) * var(store(:,2));
                r2 = 1 - SSresid/SStotal;
                r2_value(landtype,scale) = r2;
                % calculate RMSE in individual temperature windows
                % >10; >20;
                for i = 2:3
                    if (i==2)
                        idx = find(1./store(:,1)-273.15>=10);
                    elseif (i==3)
                        idx = find(1./store(:,1)-273.15>=20);
                    end
                    if (length(idx)>1)
                        stats = CalcPerf((store(idx,2)), (store(idx,3)));
                        rmse(landtype,scale,i) = stats.RMSE;
                        nrmse(landtype,scale,i) = stats.NRMSE*100;
                        nrmse2(landtype,scale,i) = stats.NRMSE2*100;
                        nrmse3(landtype,scale,i) = stats.NRMSE3*100;
                    else
                        rmse(landtype,scale,i) = nan;
                        nrmse(landtype,scale,i) = nan;
                        nrmse2(landtype,scale,i) = nan;
                        nrmse3(landtype,scale,i) = nan;
                    end
                end
            else 
                rmse(landtype,scale,1:3) = nan;
                nrmse(landtype,scale,1:3) = nan;
                nrmse2(landtype,scale,1:3) = nan;
                nrmse3(landtype,scale,1:3) = nan;
                r2_value(landtype,scale) = nan;
            end
        end
    end
end
% error statistics for the best model that includes all variability
NRMSE = sqrt(sum((fch4_obs-fch4_est).^2)/sum((fch4_obs-mean(fch4_obs)).^2));
% NRMSD, normalized root-mean-square deviation/error
NRMSD = sqrt(mean((fch4_obs-fch4_est).^2))/mean(fch4_obs);
% mean error percentage
MBE_best = (sum((fch4_est-fch4_obs))/sum(fch4_obs))*100;
% mean absolute error percentage (MAE%)
NMBE_best = (sum(abs(fch4_est-fch4_obs))/sum(fch4_obs))*100;
MAPE_best = mean(abs((fch4_obs-fch4_est)./fch4_obs));
% R2 values
yresid = fch4_obs - fch4_est;
SSresid = sum(yresid.^2);
SStotal = (length(fch4_obs)-1) * var(fch4_obs);
r2_best = 1 - SSresid/SStotal;

%% 
% this section estimates the FCH4 errors resulting from using mean
% full-season TD derived from all site-years, i.e., f(T)
% aggregate all landtypes together
% evaluate bias for each branch in each site-year
clc;
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
threshold_gpp = 0;
threshold_gs = 0;
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
delta_T = [0.1, 0.5, 1];
err_warm_agg = [];
err_cool_agg = [];
close all
for T_type = 1:1
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Air temperature (' degree 'C)'];
        yname = {'CH_4 emission'; '(mg C m^{-2} d^{-1})'};
        figname = 'TA_FCH4_error_comp_no_win_agg_v3';
    end
    err_warm = [];
    err_cool = [];
    err_warm_agg = [];
    err_cool_agg = [];
    err_full = [];
    err_full_agg = [];
    full_obs = [];
    full_est = [];
    err_mean = [];
    err_std = [];
    % aggregated TD
    xxx = source(:,col);
    idx = find(xxx<=threshold_T);
    xxx(idx) = nan;
    % FCH4
    yyy = source(:,3);
    idx = find(yyy<=threshold_FCH4);
    yyy(idx) = nan;
    if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
        idx_valid = isfinite(xxx)&isfinite(yyy);
        xxx = xxx(idx_valid);
        yyy = yyy(idx_valid);
        if (length(xxx)>1)
            idx = isfinite(xxx) & isfinite(yyy);
            rg_fit = polyfix(xxx(idx), yyy(idx), 2,...
                [0],[0]);
        end
    end
    for landtype = 1:length(ecosystem_list)
        desp = char(erase(ecosystem_name(landtype),'"'));
        id = char(erase(id_name(landtype),'"'));
        ymax = [];
        full = [];
        warm = [];
        cool = [];
        Ea_ts = [];            
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        for site=1:number_site
            idx_site = find(data(:,site_col)==site);
            data_site = data(idx_site,:);
            range_year = min(data_site(:,1)):max(data_site(:,1));
            for yr = 1:length(range_year)
                idx_year = find(data_site(:,1)==range_year(yr));
                data_site_year = data_site(idx_year,:);
                % Tair/Tsoil
                xxx = data_site_year(:,col);
                idx = find(xxx<=threshold_T);
                xxx(idx) = nan;
%                         FCH4
                yyy = data_site_year(:,3);
                idx = find(yyy<=threshold_FCH4);
                yyy(idx) = nan;
                if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                    idx_valid = isfinite(xxx)&isfinite(yyy);
                    xxx = xxx(idx_valid);
                    yyy = yyy(idx_valid);
                    idx_Tmax = find(xxx==max(xxx));
                    idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
                    idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
                    x0 = mean(xxx(idx_Tmax));
                    y0 = mean(yyy(idx_Tmax));
                    for branch = 1:3
                        if (branch==1) % full
                            idx = 1:length(xxx);
                            tmp = [xxx(idx), yyy(idx)];
                            full = [full; tmp];
                        elseif (branch==2) % warm
                            idx = idx_warming;
                            mcolor = '-r';
                            tmp = [xxx(idx), yyy(idx)];
                            warm = [warm; tmp];
                        elseif (branch==3) % cooling
                            idx = idx_cooling;
                            mcolor = '-b';
                            tmp = [xxx(idx), yyy(idx)];
                            cool = [cool; tmp];
                        end
                        % calculate full-season TD in each site-year
                        if (branch==1 & length(idx)>2)
                            xtmp = xxx(idx);
                            ytmp = yyy(idx);
                            rg_coef_full = polyfix(xtmp, ytmp, 2, ...
                                [0,x0],[0,y0]);
                        end
                        % plot the quadratic fits in warming and cooling
                        % branches for each site-year
                        if (length(idx)>2)
                            % average in each T window
                            xtmp = xxx(idx);
                            ytmp = yyy(idx);
                            % use FCH4 T dependence derived across all
                            % site-years
                            fff = polyval(rg_fit,xtmp);
                            err = (fff-ytmp)./mean(ytmp)*100;
                        end
                        if (branch==1)
                            err_full = [err_full; mean(err)];
                            err_full_agg = [err_full_agg; err];
                            full_obs = [full_obs; ytmp];
                            full_est = [full_est; fff];
                        elseif (branch==2)
                            err_warm = [err_warm; mean(err)];
                            err_warm_agg = [err_warm_agg; err];
                        elseif (branch==3)
                            err_cool = [err_cool; mean(err)];
                            err_cool_agg = [err_cool_agg; err];
                        end
                    end
                end
            end
        end
    end
end

% error statistics
NRMSE = sqrt(sum((full_obs-full_est).^2)/sum((full_obs-mean(full_obs)).^2));
% NRMSD, normalized root-mean-square deviation/error
NRMSD = sqrt(mean((full_obs-full_est).^2))/mean(full_obs);
% mean error percentage
MBE_T_alone = (sum((full_est-full_obs))/sum(full_obs))*100;
% mean absolute error percentage (MAE%)
NMBE_T_alone = (sum(abs(full_est-full_obs))/sum(full_obs))*100;
MAPE_T_alone = mean(abs((full_obs-full_est)./full_obs));
% R2 values
yresid = full_obs - full_est;
SSresid = sum(yresid.^2);
SStotal = (length(full_obs)-1) * var(full_obs);
r2_T_alone = 1 - SSresid/SStotal;

%% make scatter plots for the best and worst FCH4 models, embedded RMSE distribution
clc;
close all
for i = 1:2
    if (i==2)
        xtmp = fch4_obs;
        ytmp = fch4_est;
        id = '(b)';
        R2 = 0.76;
        run_name = '{\it F_{CH_4} = f (T,site,IAV,ISV)}';
    elseif (i==1)
        xtmp = full_obs;
        ytmp = full_est;
        id = '(a)';
        R2 = 0.21;
        run_name = '{\it F_{CH_4} = f (T)}';
    end
    idx = find(xtmp>1600);
    xtmp(idx) = [];
    ytmp(idx) = [];
    % mean absolute error percentage (MAE%)
    NMBE = (sum(abs(ytmp-xtmp))/sum(xtmp))*100;
    subplot(2,2,i)
    % one to one
    plot([min(xtmp),max(xtmp)],[min(xtmp),max(xtmp)],'k--')
    hold on
    % all FCH4
    rg_fit_tmp = polyfit(xtmp, ytmp, 1);
    xxx = linspace(min(xtmp),max(xtmp),10);
    fff = polyval(rg_fit_tmp,xxx);
    plot(xxx,fff,'b-','Linewidth',2)
    ppp = dscatter(xtmp, ytmp);
    colormap('Hot')
    axis([min(xtmp) max(xtmp) min(xtmp) max(xtmp)])
    xticks([0 500 1000 1500])
    yticks([0 500 1000 1500])
    xlabel('Observed CH_4 emission (mg C m^{-2} d^{-1})', 'FontSize', 12);
    ylabel({'Modeled CH_4 emission'; '(mg C m^{-2} d^{-1})'}, 'FontSize', 12);
    t = text(0.02,0.98,[id ' ' run_name],...
        'Units', 'Normalized', 'VerticalAlignment', 'Top',...
            'FontSize',12);
    t = text(0.02,0.85,['R^2 = ' num2str(sprintf('%1.2f',R2))],...
        'Units', 'Normalized', 'VerticalAlignment', 'Top',...
            'FontSize',12);
    t = text(0.02,0.71,['|bias| = ' num2str(sprintf('%2.0f',NMBE)) '%'],...
        'Units', 'Normalized', 'VerticalAlignment', 'Top',...
            'FontSize',12);
end
% error distribution
fch4_bin = 0:50:1600;
clear NMBE_bin num_data rmse_bin
for i = 1:2
    if (i==2)
        xtmp = fch4_obs;
        ytmp = fch4_est;
    elseif (i==1)
        xtmp = full_obs;
        ytmp = full_est;
    end
    idx = find(xtmp>1600);
    xtmp(idx) = [];
    ytmp(idx) = [];
    % calculate the bias and number of data in each bin
    for bin = 1:length(fch4_bin)-1
        idx = find(xtmp>=fch4_bin(bin)&xtmp<fch4_bin(bin+1));
        % mean absolute error percentage (MAE%)
        NMBE_bin(bin,i) = (sum(abs(ytmp(idx)-xtmp(idx)))/sum(xtmp(idx)))*100;
        num_data(bin,i) = length(idx);
        stats = CalcPerf(xtmp(idx), ytmp(idx));
        rmse_bin(bin,i) = stats.RMSE;
    end
end
% plot in stacked form 
subplot(2, 2, 3:4)
yyaxis left
bbb = bar(fch4_bin(1:end-1), rmse_bin, 'BarWidth', 1.2, 'EdgeColor', 'none');
ylabel({'Root Mean Square Error'; '(mg C m^{-2} d^{-1})'}, 'FontSize', 12);
xlabel('Observed CH_4 emission bin (mg C m^{-2} d^{-1})', 'FontSize', 12);
ax1 = gca;
ax1.YColor = 'k';
bbb(2).FaceColor = [0.8500 0.3250 0.0980];
t = text(0.02,0.98,['(c)'],...
    'Units', 'Normalized', 'VerticalAlignment', 'Top',...
    'FontSize',12);
axis tight

yyaxis right
plot(fch4_bin(1:end-1), num_data(:,2), '-','Color',[5 128 0]/255,...
    'linewidth', 2)
ylabel({'Number of data points'}, 'FontSize', 12);
ax1 = gca;
ax1.YColor = [5 128 0]/255;
ax1.YTick = [500 5000 10000 15000];
xlim([fch4_bin(1) fch4_bin(end-1)])
lgd = legend(bbb,{'{\it F_{CH_4} = f (T)}','{\it F_{CH_4} = f (T,site,IAV,ISV)}'},...
    'box','off','orientation','horizontal','FontSize',12);
lgd.Position = [0.24 0.38 0.3643 0.0452];
pause
figname=['../plots/test/' 'TA_FCH4_mip_agg_rmse3'];
print('-dpng','-r300',[figname '.jpeg']);
