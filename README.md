# FCH4_hysteresis

This repository contains all code needed to reproduce the results and analysis presented in Chang et al. 2021: Chang, K.-Y., Riley, W.J., Knox, S.H., Jackson, R.B., McNicol, G., Poulter, B., Aurela, M., Baldocchi, D., Bansal, S., Bohrer, G., Campbell, D.I., Cescatti, A. Chu, H., Delwiche, K.B., Desai, A.R., Euskirchen, E., Friborg, T., Goeckede, M., Helbig, M., Hemes, K.S., Hirano, T., Iwata, H., Kang, M., Keenan, T., Krauss, K.W., Lohila, A., Mammarella, I., Mitra, B., Miyata, A., Nilsson, M.B., Noomets, A., Oechel, W.C., Papale, D., Peichl, M., Reba, M.L., Rinne, J., Runkle, B.R.K., Ryu, Y., Sachs, T., Schäfer, K.V.R., Schimd, H.P., Shurpali, N., Sonnentag, O., Tang, A.C.I., Torn, M.S., Trotta, C., Tuittila, E.-S., Ueyama, M., Vargas, R., Vesala, T., Windham-Myers, L., Zhang, S., Zona, D., 2020. Substantial hysteresis in emergent temperature sensitivity of global wetland CH4 emissions, in revision in Nature Communications.

All time series plots presented in Chang et al. 2021 were produced by Hysteresis_time_series_plots_pub.m. The random-forest model and the associated plots presented in Chang et al. 2021 were generated by EC_daily_hysteresis_RF_model_pub.m. All other figures and analysis presented in Chang et al. 2021 were produced by FCH4_hysteresis_master_pub.m. Scripts were built and examined in in Matlab (MathWorks Inc., 2019, version 9.7.0). 

The CalcPerf and polyfix MATLAB functions were used to calculate error estimates and line fits, respectively. These functions can be downloaded from https://www.mathworks.com/matlabcentral/fileexchange/60833-psnr-mse-r-rmse-nrmse-mape-calculating and https://www.mathworks.com/matlabcentral/fileexchange/54207-polyfix-x-y-n-xfix-yfix-xder-dydx, respectively. 

Note that due to data policy restrictions, the original FLUXNET-CH4 database used in Chang et al. 2021 is not redistributed here, but is freely available for download at https://fluxnet.org/. 
