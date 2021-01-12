clear all;clc;close all;
mon = [31,28,31,30,31,30,31,31,30,31,30,31];
mon_sum = cumsum(mon);
delimiterIn = ',';
degree = sprintf('%c', char(176));
make_plot = 1; %'YES';
% QC threshold
QC = [0];
% 6/1-10/31
ana_time = mon_sum(5)+1:mon_sum(10);
% count the number of site years
count = 1;
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
var_wanted = ["Year","DOY","FCH4","TS_1","TS_2","TS_3","TA",...
    "GPP_DT","GPP_NT","LAT","NEE","NEE_F_ANN","WTD","WTD_F",...
    "P","WS","PA","LE"];
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
            if (variables==length(var_wanted))
                % convert LE to ET
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
                    daily_tmp(1:length(T{:,idx}), variables) = data...
                        /(2.5*10^6)*86400;
                else % the requested data is not available
                    daily_tmp(1:length(T{:,idx}), variables) = nan;
                end
            else
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
                idx_Tmax = find(tmp(1:180,7)==max(tmp(1:180,7)));
                % create dummy year labels
                if (yr==1)
                    data_SH = [data_SH; tmp(idx_Tmin:end,:)];
                elseif (yr==length(years))
                    tmp(1:idx_Tmin) = years(yr)-1;
                    data_SH = [data_SH; tmp(1:idx_Tmin,:)];
                else
                    tmp(1:idx_Tmin) = years(yr)-1;
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
% convert nmol/m2/s to mg/m2/d
aggregate(:,3) = aggregate(:,3)*(10^-9)*12*(10^3)*86400;
% assign columns
lat_col = 10; % LAT
NEE_col = 11; % NEE
NEE_F_col = 12; % NEE_F
WTD_col = 13; % WTD
WTD_F_col = 14; % WTD_F
AP_col = 17; % air pressure
site_col = 19; % site
land_col = 20; % ecosystem type

% filter out questionable data
tmp = aggregate(:,AP_col);
tmp(tmp==0) = nan;
aggregate(:,AP_col) = tmp;
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

% Inter Annual Variability
kkk = 8.62*10^-5;
threshold_T = 0; % T > 0 C
threshold_FCH4 = 0; % FCH4 > 0 mgC m^-2 d^-1
ecosystem_name = ["Bog", "Fen", "Marsh", "Peat plateau", "Rice", ...
    "Salt marsh", "Swamp", "Wet tundra"];
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"];
delta_T = 0.1;
close all
for T_type = 1:1
    site_year_count = 1;
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Ecosystem type'];
        yname = {'NRMSE'};
        figname_head = 'TA_coef_comp_';
    elseif (T_type==2) % soil temperature
        col = 4;
        source = aggregate;
        xname = ['Ecosystem type'];
        yname = {'NRMSE'};
        figname_head = 'TS_coef_comp_';
    end
    clf;
    clear rmse nrmse nrmse2 nrmse3 smaple_number r2_value_a r2_value_b
    idx = isfinite(source(:,col))&isfinite(source(:,3));
    % aggregated TD
    xxx = source(idx,col);
    yyy = source(idx,3);
    % valid T & FCH4
    idx = find(xxx<=threshold_T | yyy<=threshold_FCH4);
    xxx(idx) = [];
    yyy(idx) = [];
    xtmp = 1./(xxx+273.15);
    ytmp = log(yyy);
    rg_fit = polyfit(xtmp, ytmp, 1);
    % set up the matrix for IAV analysis
    IAV_ana = [];
    r2_ecosystem = [];
    r2_site = [];
    r2_season = [];
    r2_season_acc = [];
    site_num = 1;
    for landtype = 1:length(ecosystem_name)
        clf;
        count = 1;
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        land_tmp = [];
        season_acc_tmp = [];
        for scale = 1:1
            % decreased spatial and temporal description with increased scale#
            store = [];
            fff = [];
            rg_coef_warming = [];
            rg_coef_cooling = [];
            rg_par = [];
            if (scale==1)
                % scale=1, site-specific temperature dependence, 
                %           each site, each year, each season
                for site=1:number_site
                    idx_site = find(data(:,site_col)==site);
                    data_site = data(idx_site,:);
                    range_year = min(data_site(:,1)):max(data_site(:,1));
                    site_tmp = [];
                    for yr = 1:length(range_year)
                        season_tmp = [];
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
                        % environmental conditions
                        zzz_mean = data_site_year(:,[4, 7, 14, 16, 17]);
                        zzz_sum = data_site_year(:,[9, 15, 18]);
                        gpp = data_site_year(:,[9]);
                        date = data_site_year(:,[2]);
                        if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                            idx_valid = isfinite(xxx)&isfinite(yyy);
                            xxx = xxx(idx_valid);
                            yyy = yyy(idx_valid);
                            zzz_mean = zzz_mean(idx_valid,:);
                            zzz_sum = zzz_sum(idx_valid,:);
                            gpp = gpp(idx_valid);
                            date = date(idx_valid);
                            if (length(xxx)>1)
                                idx_Tmax = find(xxx==max(xxx));
                                idx_GPP_max = find(gpp==max(gpp));
                                idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
                                idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
                                x0 = mean(xxx(idx_Tmax));
                                y0 = mean(yyy(idx_Tmax));
                                if (length(idx_warming)>1 & length(idx_cooling)>1)
                                    env_mean = mean(zzz_mean,1,'omitnan');
                                    env_sum = sum(zzz_sum,1,'omitnan');
                                    if (length(idx_GPP_max)>0)
                                        lag_T_GPP = date(idx_GPP_max)-date(idx_Tmax);
                                    else
                                        lag_T_GPP = nan;
                                    end
                                    tmp = [env_mean, env_sum];
                                    store = [store; tmp];
                                    tmp = [x0, y0];
                                    rg_par = [rg_par; tmp];
                                    for branch = 1:2
                                        if (branch==1) % warming
                                            idx = idx_warming;
                                            mcolor = '^r';
                                        elseif (branch==2) % cooling
                                            idx = idx_cooling;
                                            mcolor = 'ob';
                                        end
                                        % quadratic, window average
                                        xtmp = xxx(idx);
                                        ytmp = yyy(idx);
                                        rg_coef = polyfix(xtmp, ytmp, 2, ...
                                            [0,x0],[0,y0]);
                                        val = polyval(rg_coef, xtmp);
                                        % R2 values
                                        yresid = ytmp - val;
                                        SSresid = sum(yresid.^2);
                                        SStotal = (length(ytmp)-1) * var(ytmp);
                                        r2 = 1 - SSresid/SStotal;
                                        % store the obs and modeled FCH4 pairs
                                        tmp = [ytmp, val];
                                        land_tmp = [land_tmp; tmp];
                                        site_tmp = [site_tmp; tmp];
                                        season_tmp = [season_tmp; tmp];
                                        season_acc_tmp = [season_acc_tmp; tmp];
                                        % IAV analysis, j=1 year, j=2 r2,
                                        % j=3 a, j=4 b, j=5 x0, j=6 y0, j=7
                                        % branch, j=8 site, j=9 landtype,
                                        % j=10 number of data points, j=11
                                        % mean(Tsoil), j=12 mean(Tair), j=13 WTD
                                        % j=14 wind speed, j=15 air pressure
                                        % j= 16 GPP, j=17 precip, j=18 ET
                                        % j=19 site number
                                        if (branch==1) % warming
                                            idx = idx_warming;
                                        elseif (branch==2) % cooling
                                            idx = idx_cooling;
                                        end
                                        tmp = [data_site_year(1,1), r2, rg_coef(1:2), ...
                                            x0, y0, branch, data_site_year(1,site_col), ...
                                            data_site_year(1,land_col), length(idx),...
                                            env_mean,env_sum, site_num, ...
                                            data_site_year(1,lat_col), site_year_count, lag_T_GPP];
                                        IAV_ana = [IAV_ana; tmp];
                                        if (branch==1) % warming
                                            rg_coef_warming = [rg_coef_warming; rg_coef];
                                        elseif (branch==2) % cooling
                                            rg_coef_cooling = [rg_coef_cooling; rg_coef];
                                        end
                                    end
                                    % R2 values, seasonal
                                    yresid = season_tmp(:,1) - season_tmp(:,2);
                                    SSresid = sum(yresid.^2);
                                    SStotal = (length(season_tmp(:,1))-1) * var(season_tmp(:,1));
                                    r2 = 1 - SSresid/SStotal;
                                    % store the obs and modeled FCH4 pairs
                                    tmp = [r2, site, landtype];
                                    r2_season = [r2_season; tmp];
                                    % R2 values, seasonal accumulated
                                    yresid = season_acc_tmp(:,1) - season_acc_tmp(:,2);
                                    SSresid = sum(yresid.^2);
                                    SStotal = (length(season_acc_tmp(:,1))-1) * var(season_acc_tmp(:,1));
                                    r2 = 1 - SSresid/SStotal;
                                    % store the obs and modeled FCH4 pairs
                                    tmp = [r2, site, landtype];
                                    r2_season_acc = [r2_season_acc; tmp];
                                end
                                site_year_count = site_year_count+1;
                            end                            
                        end                        
                    end
                    site_num = site_num+1;
                    count = count+1;
                    if (length(site_tmp)>0)
                        % R2 values, site
                        yresid = site_tmp(:,1) - site_tmp(:,2);
                        SSresid = sum(yresid.^2);
                        SStotal = (length(site_tmp(:,1))-1) * var(site_tmp(:,1));
                        r2 = 1 - SSresid/SStotal;
                        % store the obs and modeled FCH4 pairs
                        tmp = [r2, site, landtype];
                        r2_site = [r2_site; tmp];
                    end
                end
            end
        end
        for i = 1:2
            if (i==1)
                tmp =rg_coef_warming;
            elseif (i==2)
                tmp =rg_coef_cooling;
            end
            for j = 1:4
                idx = isfinite(store(:,j)) & isfinite(tmp(:,1));
                coef = corrcoef(store(idx,j), tmp(idx,1));
                r2_value_a(i, j, landtype) = coef(1,2).^2;
                smaple_number(i, j, landtype) = sum(idx);
                idx = isfinite(store(:,j)) & isfinite(tmp(:,2));
                coef = corrcoef(store(idx,j), tmp(idx,2));
                r2_value_b(i, j, landtype) = coef(1,2).^2;
            end
            
        end
        % R2 values
        yresid = land_tmp(:,1) - land_tmp(:,2);
        SSresid = sum(yresid.^2);
        SStotal = (length(land_tmp(:,1))-1) * var(land_tmp(:,1));
        r2 = 1 - SSresid/SStotal;
        % store the obs and modeled FCH4 pairs
        r2_ecosystem = [r2_ecosystem, r2];
    end
end

pause
%% use the same set of RF to predict ahys during warming and cooling
% include time lag between GPP&T
% RF model without year, precip, GPP, WTD, AP, WS
close all
clc;
Year = (IAV_ana(:,1));
Branch = categorical(IAV_ana(:,7));
Landtype = categorical(IAV_ana(:,9));
Site = categorical(IAV_ana(:,19));
LAT = IAV_ana(:,20);
T_max = IAV_ana(:,5);
Site_year = categorical(IAV_ana(:,21));
% observed
FCH4_T_max = IAV_ana(:,6);
% RF modeled
% FCH4_T_max = oobPredict(model_fch4_Tmax);
T_mean = IAV_ana(:,12);
GPP = IAV_ana(:,16);
GPP_T_lag = IAV_ana(:,22);
WTD = IAV_ana(:,13);
Precip = IAV_ana(:,17);
WS = IAV_ana(:,14);
AP = IAV_ana(:,15);
% ET = IAV_ana(:,18);
hys = IAV_ana(:,3);
X = table(FCH4_T_max, Site, T_max, Site_year, T_mean, LAT, GPP, Precip, ...
    Landtype, Branch,   ...
    hys);
xname = {'FCH4_T_max'; 'Site'; 'T_max'; 'Site-year'; 'T_mean'; 'LAT'; ...
    'GPP'; 'Precip'; 'Type'; 'Branch'};
rng('default'); % For reproducibility
NumTrees = 150;
t = templateTree('NumPredictorsToSample','all', ...
    'PredictorSelection', 'curvature', 'Surrogate', 'on');
rng('default'); % For reproducibility
Mdl = fitrensemble(X, 'hys', 'Method', 'bag', 'NumLearningCycles', 200, 'Learners', t);
model_ahys_obs = Mdl;
% % predictor importance 
imp = oobPermutedPredictorImportance(Mdl);
[imp, idx] = sort(imp,'descend');
clf;
subplot(2,1,1)
bar(imp);
ylabel({'Unbiased Predictor'; 'Importance Estimates'}, 'FontSize', 12);
xlabel('Predictors', 'FontSize', 12);
h = gca;
h.XTickLabel = xname(idx);
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
t = text(0.02,0.98,['(a)'],...
    'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12);
subplot(2,1,2)
% R2
yHat = oobPredict(Mdl);
R2 = corr(Mdl.Y,yHat)^2;
% one to one
plot([min(Mdl.Y),max(Mdl.Y)],[min(Mdl.Y),max(Mdl.Y)],'k--')
hold on
% line fit
rg_fit_tmp = polyfit(Mdl.Y, yHat, 1);
xxx = linspace(min(Mdl.Y),max(Mdl.Y),10);
fff = polyval(rg_fit_tmp,xxx);
plot(xxx,fff,'b-','Linewidth',2)
ppp = dscatter(Mdl.Y, yHat);
colormap('Hot')
c = colorbar('Position', [0.92, 0.17, 0.02, 0.3]);
axis([min(Mdl.Y) max(Mdl.Y) min(Mdl.Y) max(Mdl.Y)])
xlabel('Hysteresis parameter', 'FontSize', 12);
ylabel({'Estimated'; 'hysteresis parameter'}, 'FontSize', 12);
t = text(0.02,0.98,['(b) R^2 = ' num2str(sprintf('%1.2f',R2))],...
    'Units', 'Normalized', 'VerticalAlignment', 'Top',...
        'FontSize',12);
figname=['../plots/test/' 'RF_ahys'];
print('-dpng','-r300',[figname '.jpeg']);

%% compare the estimated FCH4 across all models, without window avg
% hysteretic scaling is a function of end_T, mgs_T, landtype, GPP, WTD, Precip
% build up the array that define the hysteretic scaling for all sites
% j=1 site, j=2 landtype, j=3 hysteresis parameter, j=4 FCH4@Tmax, 
% j=5 branch, j=6 year
clc;
clear hys_scaling
hys_scaling(:,1) = IAV_ana(:,8);
hys_scaling(:,2) = IAV_ana(:,9);
hys_scaling(:,3) = oobPredict(model_ahys_obs);
hys_scaling(:,4) = IAV_ana(:,6);
hys_scaling(:,5) = IAV_ana(:,7);
hys_scaling(:,6) = IAV_ana(:,1);
% predict FCH4 with all methods and compare their performance
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
for T_type = 1:1
    if (T_type==1) % air temperature
        col = 7;
        source = aggregate;
        xname = ['Ecosystem type'];
        yname = {'NRMSE'};
        figname_head = 'TA_NRMSE_';
    end
    clf;
    clear rmse nrmse nrmse2 nrmse3 r2_value
    idx = find(source(:,8)<=threshold_gpp);
    source(idx,:) = nan;
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
        % loop through all ecosystem types
        idx = find(source(:,land_col)==landtype);
        data = source(idx,:);
        number_site = max(data(:,site_col));
        idx = find(hys_scaling(:,2)==landtype);
        hys_tmp = hys_scaling(idx,:);
        for scale = 1:7
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
                                    for branch = 1:2
                                        if (branch==1) % warming
                                            idx = idx_warming;
                                        elseif (branch==2) % cooling
                                            idx = idx_cooling;
                                        end
                                        % quadratic, window average
                                        xtmp = xxx(idx);
                                        ytmp = yyy(idx);
                                        T_range = floor(min(xxx)):delta_T:ceil(max(xxx));
                                        idx = isfinite(xtmp) & isfinite(ytmp);
                                        rg_coef = polyfix(xtmp(idx), ytmp(idx), 2, ...
                                            [0,x0],[0,y0]);
                                        if (branch==1) % warming
                                            rg_coef_warming = rg_coef;
                                        elseif (branch==2) % cooling
                                            rg_coef_cooling = rg_coef;
                                        end
                                    end
                                    % predicted values, warming branch
                                    fff = polyval(rg_coef_warming,xxx(1:idx_Tmax));
                                    % predicted values, cooling branch
                                    fff = [fff; polyval(rg_coef_cooling,xxx(idx_Tmax+1:end))];
                                    store_tmp = [xxx, yyy, fff];
                                    store = [store; store_tmp];
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
                                xtmp = xxx(:);
                                ytmp = yyy(:);
                                T_range = floor(min(xxx)):delta_T:ceil(max(xxx));
                                idx = isfinite(xtmp) & isfinite(ytmp);
                                rg_fit_tmp = polyfix(xtmp(idx), ytmp(idx), 2,...
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
                            xtmp = xxx(:);
                            ytmp = yyy(:);
                            T_range = floor(min(xxx)):delta_T:ceil(max(xxx));
                            idx = isfinite(xtmp) & isfinite(ytmp);
                            rg_fit_tmp = polyfix(xtmp(idx), ytmp(idx), 2,...
                                [0],[0]);
                            % predicted values
                            fff = polyval(rg_fit_tmp,xxx(:));
                            store_tmp = [xxx, yyy, fff];
                            store = [store; store_tmp];
                        end
                    end
                end
            elseif (scale==4)
                % scale=4, using the hysteretic relations from RF
                for site=1:number_site
                    idx_site = find(hys_tmp(:,1)==site);
                    hys_site = hys_tmp(idx_site,:);
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
                        % hysteretic scaling
                        if (sum(isfinite(xxx))>=2 & sum(isfinite(yyy))>=2) % the requested data is available
                            idx_valid = isfinite(xxx)&isfinite(yyy);
                            xxx = xxx(idx_valid);
                            yyy = yyy(idx_valid);
                            if (length(xxx)>1)
                                idx_Tmax = find(xxx==max(xxx));
                                idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
                                idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
                                % hysteretic scaling
                                idx = find(hys_site(:,6)==range_year(yr));
                                c = mean(hys_site(idx,4));
                                idx = find(hys_site(:,5)==1 & hys_site(:,6)==range_year(yr));
                                a_warming = hys_site(idx,3);
                                b_warming = c/mean(xxx(idx_Tmax))-a_warming*mean(xxx(idx_Tmax));
                                idx = find(hys_site(:,5)==2 & hys_site(:,6)==range_year(yr));
                                a_cooling = hys_site(idx,3);
                                b_cooling = c/mean(xxx(idx_Tmax))-a_cooling*mean(xxx(idx_Tmax));
                                if (length(idx_warming)>1 & length(idx_cooling)>1)
                                    for branch = 1:2
                                        if (branch==1) % warming
                                            idx = idx_warming;
                                            rg_coef_warming = [a_warming,b_warming,0];
                                        elseif (branch==2) % cooling
                                            idx = idx_cooling;
                                            rg_coef_cooling = [a_cooling,b_cooling,0];
                                        end
                                    end
                                    % predicted values, warming branch
                                    fff = polyval(rg_coef_warming,xxx(1:idx_Tmax));
                                    % predicted values, cooling branch
                                    fff = [fff; polyval(rg_coef_cooling,xxx(idx_Tmax+1:end))];
                                    store_tmp = [xxx, yyy, fff];
                                    store = [store; store_tmp];
                                end
                            end
                        end
                    end
                end
            elseif (scale==6)
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
                        xtmp = xxx(:);
                        ytmp = yyy(:);
                        T_range = floor(min(xxx)):delta_T:ceil(max(xxx));
                        idx = isfinite(xtmp) & isfinite(ytmp);
                        rg_fit_tmp = polyfix(xtmp(idx), ytmp(idx), 2,...
                            [0],[0]);
                        % predicted values
                        fff = polyval(rg_fit_tmp,xxx(:));
                        store_tmp = [xxx, yyy, fff];
                        store = [store; store_tmp];
                    end
                end
            elseif (scale==7)
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
                    % predicted values
                    fff = polyval(rg_fit,xxx(:));
                    store_tmp = [xxx, yyy, fff];
                    store = [store; store_tmp];
                end
            elseif (scale==5)
                % scale=7, mean temperature dependence across site & years, 
                %           each season, each landtype
                warming = [];
                cooling = [];
                for site=1:number_site
                    idx_site = find(data(:,site_col)==site);
                    data_site = data(idx_site,:);
                    range_year = min(data_site(:,1)):max(data_site(:,1));
                    for yr = 1:length(range_year)
                        idx_year = find(data_site(:,1)==range_year(yr));
                        data_site_year = data_site(idx_year,:);
                        idx_gs = find(data_site_year(:,8)>threshold_gpp);
                        if (length(idx_gs)>threshold_gs)
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
                                    if (length(idx_warming)>1 & length(idx_cooling)>1)
                                        % warming branch
                                        tmp = [xxx(idx_warming), yyy(idx_warming)];
                                        warming = [warming; tmp];
                                        % cooling branch
                                        tmp = [xxx(idx_cooling(2:end)), yyy(idx_cooling(2:end))];
                                        cooling = [cooling; tmp];
                                    end
                                end
                            end
                        end
                    end
                end
                % calculate mean FCH4-T dependence at the warming and
                % cooling branches
                if (length(warming)>0)
                    tmp = [warming(:,1); cooling(:,1)];
                    T_window = floor(min(tmp)):.5:ceil(max(tmp));
                    idx = isfinite(warming(:,1)) & isfinite(warming(:,2));
                    T = warming(idx,1);
                    FCH4 = warming(idx,2);
                    idx = find(T==max(T));
                    x0 = T(idx);
                    y0 = FCH4(idx);
                    rg_fit_tmp = polyfix(T, FCH4, 2, [0,x0],[0,y0]);
                    % predicted values
                    fff = polyval(rg_fit_tmp,warming(:, 1));
                    store_tmp = [warming(:, 1), warming(:, 2), fff];
                end
                if (length(cooling)>0)
                    clear T FCH4
                    idx = isfinite(cooling(:,1)) & isfinite(cooling(:,2));
                    T = cooling(idx,1);
                    FCH4 = cooling(idx,2);
                    idx = find(T==max(T));
                    x0 = T(idx);
                    y0 = FCH4(idx);
                    rg_fit_tmp = polyfix(T, FCH4, 2, [0,x0],[0,y0]);
                    % predicted values
                    fff = polyval(rg_fit_tmp,cooling(:, 1));
                    tmp = [cooling(:, 1), cooling(:, 2), fff];
                end
                if (length(warming)>0 | length(cooling)>0)
                    store = [store; store_tmp; tmp];
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
                % absolute bias
                bias = (sum((store(:,3)-store(:,2)))/sum(store(:,2)))*100;
                bias_value(landtype,scale) = bias;
                bias = (sum(abs(store(:,3)-store(:,2)))/sum(store(:,2)))*100;
                abs_bias_value(landtype,scale) = bias;
                % Mean Absolute Percentage Error
                MAPE(landtype,scale,1) = stats.Mape;
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
                        MAPE(landtype,scale,i) = stats.Mape;
                    else
                        rmse(landtype,scale,i) = nan;
                        nrmse(landtype,scale,i) = nan;
                        nrmse2(landtype,scale,i) = nan;
                        nrmse3(landtype,scale,i) = nan;
                        MAPE(landtype,scale,i) = nan;
                    end
                end
            else 
                rmse(landtype,scale,1:3) = nan;
                nrmse(landtype,scale,1:3) = nan;
                nrmse2(landtype,scale,1:3) = nan;
                nrmse3(landtype,scale,1:3) = nan;
                r2_value(landtype,scale) = nan;
                bias_value(landtype,scale) = nan;
                abs_bias_value(landtype,scale) = nan;
                MAPE(landtype,scale,1:3) = nan;
            end
        end
        % data point distribution
        idx = find(source(:,land_col)==landtype);
        if (length(idx)>1)
            fraction(landtype,1) = max(source(idx,site_col));
            xxx = source(idx,col);
            yyy = source(idx,3);
            % valid T & FCH4
            idx_QC = find(xxx>threshold_T & yyy>threshold_FCH4);
            % number of QCed data points
            fraction(landtype,2) = length(idx_QC);
            % number of total data points
            fraction(landtype,3) = length(idx);
            mean_fch4(landtype) = mean(source(idx,3),'omitnan');
        end
    end
    % make plot- 
    % data distribution
    close all;
    subplot(5,1,5)
    yyaxis right
    xtmp = 1.25:8.25;
    bbb = bar(xtmp-.07,fraction(:,2));
    bbb.BarWidth = 0.3;
    ylabel('Data points', 'FontSize', 12);
    t = text(0.02,0.95,'(b)',...
            'Units', 'Normalized', 'VerticalAlignment', 'Top');
    yyaxis left
    xtmp = 0.75:7.75;
    sss = bar(xtmp+.07,fraction(:,1));
    sss.BarWidth = 0.3;
    ylabel('Sites', 'FontSize', 12);
    set(gca, 'SortMethod','depth')
    ax = gca;
    ax.XLim = [.5 8.5];
    subplot(5,1,1:4)
    err_array = abs_bias_value;
    tmp = flip(err_array, 2);
    hhh = image(tmp','CDataMapping','scaled');
    % set nan as blank
    set(hhh,'AlphaData',~isnan(tmp'))
    t = text(0.02,1.05,'(a)',...
            'Units', 'Normalized', 'VerticalAlignment', 'Top');
    colormap(jet(256))
    caxis([20 100])
    c = colorbar('NorthOutside', 'Position', [0.18, 0.93, 0.7, 0.03]);
    % overlay grid lines
    [rows, columns, numberOfColorChannels] = size(tmp');
    hold on;
    lineSpacing = 1;
    for row = .5 : lineSpacing : rows+.5
        line([.5, columns+.5], [row, row], 'Color', 'k', 'LineWidth', 1);
    end
    for col = .5 : lineSpacing : columns+.5
        line([col, col], [.5, rows+.5], 'Color', 'k', 'LineWidth', 1);
    end
    ax = ancestor(hhh, 'axes');
    yrule = ax.YAxis;
    % Change properties of the axes
    ax.YTick = [1:7];
    ax.YTick = [1:7];
    ax.YTickLabel = {'{\it f (T)}','{\it f (T,type)}','{\it f (T,type,ISV)}',...
        '{\it f (hybrid)}',...
        '{\it f (T,site)}','{\it f (T,site,IAV)}','{\it f (T,site,IAV,ISV)}',...
        };
    ax.XTick = [ ];
    % Change properties of the ruler
    yrule.FontSize = 11;
    yrule.TickLabelRotation = 30;
    ax_bar = ancestor(bbb, 'axes');
    xrule = ax_bar.XAxis;
    % Change properties of the axes
    ax_bar.XTick = [1:8];
    ax_bar.XTickLabel = {'Bog', 'Fen', 'Marsh', 'Peat plateau', 'Rice', ...
            'Salt marsh', 'Swamp', 'Wet tundra'};
    % Change properties of the ruler
    xrule.FontSize = 12;
    % set colorbar description
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    axes(ax1) % sets ax1 to current axes
%         text(0.915,0.975,['(%)'],'Units', 'Normalized', 'VerticalAlignment', 'Top')
    % arrow 1
    x1 = [.92 .92];
    y1 = [.83 .65];
    annotation('doublearrow',x1,y1,'Color','k');
    t1 = text(0.98,.85,'Ecosystem-type','Units', 'Normalized', ...
        'VerticalAlignment', 'Top', 'FontSize', 12);
    t1.Rotation = 270;
    t1 = text(0.95,.80,'variability','Units', 'Normalized', ...
        'VerticalAlignment', 'Top', 'FontSize', 12);
    t1.Rotation = 270;
%     % arrow 2
    x2 = [.92 .92];
    y2 = [.65 .28];
    annotation('doublearrow',x2,y2,'Color','k');
    t2 = text(0.96,.63,'Ecosystem-site variability','Units', 'Normalized', ...
        'VerticalAlignment', 'Top', 'FontSize', 12);
    t2.Rotation = 270;
    % set colorbar style
    text(0.89,0.97,['|bias| (%)'],'Units', 'Normalized', 'VerticalAlignment', 'Top')
    figname=['../plots/test/' figname_head 'bias_RF_obsFCH4'];
    print('-dpng','-r300',[figname '.jpeg']);
end
