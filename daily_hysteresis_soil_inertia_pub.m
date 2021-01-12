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
var_wanted = ["Year","DOY","FCH4","TS_1","TS_2","TS_3","TS_4","TS_5","TA",...
    "GPP_DT","GPP_NT","LAT","NEE","NEE_F_ANN","WTD","WTD_F"];
% read in data
filename=['../marsh/USMyb.csv'];
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

% convert FCH4 from nmol/m2/s to mg/m2/d
usmyb = daily_tmp;
usmyb(:,3) = usmyb(:,3)*(10^-9)*12*(10^3)*86400;
                    
pause

%% plot the FCH4-T relationships from air to 5 layers of soil T
% individual hysteresis plots
% j=4 TS1, j=5 TS2, j=6 TS3, j=7 TS4, j=8 TS5, j=9 TA 
% USMyb only have j=6-9 in 2017
clc; close all
id_name = ["(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)",...
    "(m)", "(n)", "(o)", "(p)", "(q)", "(r)"];
mcolor = [.02 .50 .00; .90 .17 .31; .0 .0039 1.00; .00 .50 1.00; 1.00 .57 .69;...
    .11 .30 .24; .16 .32 .75; .00 .75 1.00; .76 .60 .42; .45 .63 .76; .53 .33 .04; 0 0 0];
threshold_T = 0;
threshold_FCH4 = 0;
delta_T = 0.1;
obs_yr = 2017:2017;
xname = ['Temperature (' degree 'C)'];
yname = {'CH_4 emission, {\it F_{CH_4}}'; '(mg C m^{-2} d^{-1})'};
figname = 'FCH4_T_hysteresis_layered_v2';
% aggregate plot
plot_id = [];
for yr = 1:length(obs_yr)
    idx_year = find(usmyb(:,1)==obs_yr(yr));
    data_site_year = usmyb(idx_year,:);
    for T_type = 6:9
        col = T_type;
        % set marker
        if (T_type==4)
            style = '.';
        elseif (T_type==5)
            style = '^';
        elseif (T_type==6)
            style = 's';
            ii = 1;
        elseif (T_type==7)
            style = 'o';
            ii = 2;
        elseif (T_type==8)
            style = 'x';
            ii = 3;
        elseif (T_type==9)
            style = '.';
            ii = 12;
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
            subplot(2, 4, 1:4)
            ppp = plot(xxx, yyy, style, 'MarkerFaceColor', mcolor(ii,:),...
                'MarkerEdgeColor', mcolor(ii,:), 'MarkerSize', 7);
            hold on
            plot_id = [plot_id, ppp];
        end
        if (T_type==9)
            lgd = legend(plot_id,{'T_{soil}, 8 cm','T_{soil}, 16 cm',...
                'T_{soil}, 32 cm', 'T_{air}'},'orientation','vertical',...
                    'box','off');
            lgd.FontSize = 12;
            lgd.Position = [0.05 0.75 0.3643 0.0452];
            xlabel(xname, 'FontSize', 12);
            ylabel(yname, 'FontSize', 12);
            title('Apparent temperature sensitivity inferred from US-Myb in 2017',...
                'FontSize', 12)
            t = text(0.02,0.95,['(a)'],...
                    'Units', 'Normalized', 'VerticalAlignment', 'Top');
        end
    end
end
plot_id = [];
for yr = 1:length(obs_yr)
    idx_year = find(usmyb(:,1)==obs_yr(yr));
    data_site_year = usmyb(idx_year,:);
    for T_type = 6:9
        % set marker
        if (T_type==4)
            style = '.';
        elseif (T_type==5)
            style = '^';
        elseif (T_type==6)
            id = '(b) T_{soil}, 8 cm';
            plot_id = 5;
            xname = 'T_{soil}';
        elseif (T_type==7)
            id = '(c) T_{soil}, 16 cm';
            plot_id = 6;
            xname = 'T_{soil}';
        elseif (T_type==8)
            id = '(d) T_{soil}, 32 cm';
            plot_id = 7;
            xname = 'T_{soil}';
        elseif (T_type==9)
            id = '(e) T_{air}';
            plot_id = 8;
            xname = 'T_{air}';
        end
        col = T_type;
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
            xfit = floor(min(xxx)):delta_T:ceil(max(xxx));
            idx_Tmax = find(xxx==max(xxx));
            idx_warming = 1:idx_Tmax; % warming branch is from the beginning of the year to Tmax
            idx_cooling = idx_Tmax:length(xxx); % cooling branch is from Tmax to the end of the year
            x0 = mean(xxx(idx_Tmax));
            y0 = mean(yyy(idx_Tmax));
            if (length(idx_warming)>1 & length(idx_cooling)>1)
                colormap('redblue')
                subplot(2, 4, plot_id)
                plot(xxx(1),yyy(1),'w.');
                hold on
                c = linspace(10,1,length(yyy));
                scatter(xxx,yyy,15,c,'filled')
                xlabel([xname ' (' degree 'C)'], 'FontSize', 12)
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
                [yfit_total, Ea_w] = YD2014_fit(xtmp(:), ytmp(:),xfit);
%                 % cooling branch
                xtmp = xxx(idx_cooling);
                ytmp = yyy(idx_cooling);
                rg_coef = polyfix(xtmp, ytmp, 2, ...
                    [0,x0],[0,y0]);
%                 % predicted values
                fch4_cooling = polyval(rg_coef,T_range);
                fit_line_c = plot(T_range,fch4_cooling,'b-','linewidth',2);
                % apparent activation energy
                [yfit_total, Ea_c] = YD2014_fit(xtmp(:), ytmp(:),xfit);
                axis([min(xxx) max(xxx) 0 max(yyy)])
                t = text(0.02,0.98,[id],...
                    'Units', 'Normalized', 'VerticalAlignment', 'Top');
                delTACH4 = [mean(yyy(idx_cooling))-mean(yyy(idx_warming))];
                norm_hys_area = sum(fch4_cooling-fch4_warming)*delta_T/...
                     (T_range(end)-T_range(1))/max(abs([fch4_cooling, fch4_warming]));
                desp = '{\it H_{\mu}} = ';
                t = text(0.02,0.78,[desp num2str(sprintf('%.1f',delTACH4))],...
                    'Units', 'Normalized', 'VerticalAlignment', 'Top');
                desp = '{\it H_A} = ';
                t = text(0.02,0.88,[desp num2str(sprintf('%.1f',norm_hys_area))],...
                    'Units', 'Normalized', 'VerticalAlignment', 'Top');
                if (T_type==6)
                    ylabel(yname, 'FontSize', 12)
                end
            end
        end
    end
end
c = colorbar('Position', [0.935, 0.1, 0.02, .36]);
set(c,'YTick',[])
set(c,'YDir','reverse')
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1) % sets ax1 to current axes
text(0.91,0.50,['End date'],'Units', 'Normalized', 'VerticalAlignment', 'Top',...
    'FontSize',12)
text(0.89,0.075,['Start date'],'Units', 'Normalized', 'VerticalAlignment', 'Top',...
    'FontSize',12)
axes(ax1) % sets ax1 to current axes
lgd = legend([fit_line_w, fit_line_c],{'Earlier', 'Later'},'orientation','horizontal',...
        'box','off');
lgd.FontSize = 12;
lgd.Position = [0.60 0.45 0.3643 0.0452];
print('-dpng','-r300',['../plots/test/' figname '.jpeg']);
