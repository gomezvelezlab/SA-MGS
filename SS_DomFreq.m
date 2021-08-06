
% Selection of Dominant Frequencies to be used on the Groundwater
% Spectral Solutions presented by Perez & Gomez-Velez 2020 (WRR submitted).

% This code is adapted from the 2DSpecTools: A 2D spectral analysis
% toolkit for MATLAB by Taylor Perrorn (http://web.mit.edu/perron/www/downloads.html)
% Reference: Perron, J. T., Kirchner, J. W., & Dietrich, W. E. (2008). 
% Spectral signatures of characteristic spatial scales and nonfractal structure in landscapes. 
% Journal of Geophysical Research: Earth Surface, 113(4), 1–14. https://doi.org/10.1029/2007JF000866

%   Copyright:  Gabriel Perez
%   Gomez-Velez Lab Group. Vanderbilt University

%   email:   gabriel.perez.mesa@vanderbilt.edu
%	06  Aug	2021,   MATLAB	2019b  version
%	IF  YOU	PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%	Perez. G., et al., (2021)  Identification of Characteristic Spatial Scales to Improve the Performance 
%   of Analytical Spectral Solutions to the Groundwater Flow Equation (Submitted to WRR)

function [fx_dom,fy_dom,P_Value_Results]=SS_DomFreq(Z,Step,N_fxy,p_value_t,Method_mean_spectrum,f_cut,f_max,Prop,Iter,Path_export_Fig)

    % Input arguments

    % Z:                        A NxM Matrix containing the DEM in meters
    % Step:                     DEM resolution in meters
    % N_fxy:                    Number of pairs (fx and fy) of dominat frequencies to be sampled
    % p_value_t:                P-value threshold (e.g. 0.95) to select the dominant frequencies. 
    % Method_mean_spectrum:     1: Uses Random DEMs to get the average spectrum, 0:uses the average spectrum from the observed DEM
    % f_cut:                    Optional parameter to increment the
    % freqeuncy resolution for  f<f_cut. If f_cut=0 there is not cut and freqeuncy resolution is the same provided by the Discrete Fourier Transform
    % f_max:                    This is the largest frequency to be sampled. In practice this should be a couple of times larger than the DEM resolution
    % Prop:                     This is the proportion of frequencies to be sampled in the range [0 f_cut]
    % Iter:                     Number of iterations in the estimation of the average power spectrum
    % Path_export_Fig:          Path to save the figures

    % Output arguments

    % fx_dom:                   Dominant frequencies in x-direction
    % fy_dom:                   Dominant frequencies in y-direction
    % P_values_Results:         P-values matrix for the 2D power spectrum
 
%% Load Perron's code
addpath([pwd '\Perron_Codes'])
rng(10);
%% Visualize Top Boundary Pressure Head
figure
[x,y]=meshgrid(1:length(Z(1,:)),1:length(Z(:,1)));
mesh(x.*Step,y.*Step,Z); shading interp; camlight left; lighting gouraud;   
saveas(gcf,[Path_export_Fig '\TopBoundary_PressureHead.bmp']); close

%% Power Spectrum
Z = Detrend(Z);                                                             % The DEM must be detrended in order to get the right signals. 
pad = 1;                                                                    % Padding the DEM with zeros (more effcient for FFT)
window = 1;                                                                 % Windowing
dx = Step;                                                                  % grid spacing in the x-direction
dy = Step;                                                                  % grid spacing in the y-direction

[Pm, fm, Pv, fv] = fft2D(Z,dx,dy,pad,window);                               %2D FFT
figure; subplot(2,1,1); [axf, ~] = SpecPlot1D(fv,Pv);                       % Plot 1D power Spectrum
hold(axf,'on')
nbin = 20;                                                                  % Number of bins to get the binned 1D power spectrum
B = bin(log10(fv),log10(Pv),nbin,0);                                        % Bin the log-transformed data
plot(axf,10.^B(:,1),10.^B(:,2),'ok','markerfacecolor','w');                 % Plot the binned values
fit = robustfit(B(:,1),B(:,2));                                             % Fit a trend with the form P ~ 1/f^n, and plot it
plot(axf,10.^B(:,1),10^fit(1)*(10.^B(:,1)).^fit(2),'k');
xlabel(axf,[]);

% Get the roll-off frequency 
% The frequency at which the spectral roll-off occurs was identified for 
% each spectrum by fitting  least squares regression lines to the 
% log-transformed, binned data. Beggining with the lowest frequency point, 
% and subsequently moving through all the points,  the lines was fit to 
% all points at freqeuncies greater than or equal to the  present point.
% The roll-off frequency was identified as the point at which the magnitude
% of the regression line slope (beta) reached a maximum (Perron et al., 2008)

fv_mean=10.^B(:,1); Pv_mean=10.^B(:,2);
Log_fv_mean=log10(fv_mean); Log_Pv_mean=log10(Pv_mean);

for i=1:size(fv_mean,1)-1
    X=Log_fv_mean(i:end);  Y=Log_Pv_mean(i:end);
    tbl_reg=table(Y,X);
    mdl = fitlm(tbl_reg,'Y ~ X');
    R2(i)=mdl.Rsquared.Ordinary;
    Beta(i)=abs(mdl.Coefficients.Estimate(2));                              % This is the magnitude of the regression line, then it is in absolute value
end

subplot(2,1,2)
yyaxis left
plot(fv_mean(1:end-1),R2,'--o'); hold on
ylim([min(R2)-0.1 1]); ylabel('R^2');
yyaxis right
plot(fv_mean(1:end-1),Beta,'--*');
set(gca,'Xscale','log');  ylabel('\beta');
ylim([0 max(Beta)+1]);
[~, pos_Temp]=max(Beta);                                                    % Locate the roll-off frequency

if pos_Temp<=2                                                              % Take at least 3 points for the regression
    pos_Temp=3;
end

if pos_Temp==1                                                              % There is not f_roll_off, then select the lowest R2 as f_roll_off
   pos_Temp=size(fv_mean,1);
   f_roll_off=fv_mean(end); 
else
   f_roll_off=fv_mean(pos_Temp(1)); 
end

yyaxis left
plot([f_roll_off f_roll_off],[0.8 1],'k');
xlabel('Radial frequency (m^{-1})');

f_seg_roll_off=fv_mean(1:pos_Temp(1));                                      % Get the segment of the binned values below the roll-off frequency

% Use this fit to normalize the 2D spectrum. Note that this fit can be used
% as the average power spectrum (which is the binned power spectrum), 
% however these values are biased by the peaks.
% A better alternative is to use Random DEMs with similar statistics than the
% observed DEM in order to get the average mean spectrum

Pmn = Pm./(10^fit(1)*fm.^fit(2));                                           % Normalized 2D spectrum
figure; SpecPlot2D(fm,log10(Pmn));                                          % Plot the normalized 2D spectrum.  The labels is an approximation

%% Get the average power spectrum from "Iter" iteration of random topographic surfaces using the diamond-square algorithm 
tsize=max(size(Z));                                                         % Get the higher dim in x and y
if floor(log(tsize-1)/log(2))==log(tsize-1)/log(2)
else 
    tsize=2^(ceil(log(tsize)/log(2)))+1;                                    % Go to the closer power of 2 +1
end
range=abs(max(max(Z))-min(min(Z)));                                         % Get the range in elevation

H_Vector=0:0.1:1; 
Log_p2=log10(10^fit(1)*(10.^f_seg_roll_off).^fit(2));                       % Using line fitted to the binned data
    
poolobj = parpool;                                                          % Use parallel computing for better performance
parfor j=1:size(H_Vector,2)   
    H=H_Vector(j);                                                          % Define roughness. This is the parameter to be found
    Pv_r_mean=0; fv_r_mean=0;
    for i=1:Iter
        Z_r=createFractalTerrain(tsize, range, H);                          % Get the random DEM with diamond-square algorithm
        Z_r=Z_r(1:size(Z,1),1:size(Z,2));                                   % Get the same dimensions of the observed DEM        
        [~, ~, Pv_r, fv_r] = fft2D(Z_r,dx,dy,pad,window);                   % Get the 1D power spectrum
        Pv_r_mean=Pv_r_mean+Pv_r; fv_r_mean=fv_r_mean+fv_r;
    end
    Pv_r_mean=Pv_r_mean./Iter;
    fv_r_mean=fv_r_mean./Iter;
    nbin = 20; % Number of bins
    B_r_mean = bin(log10(fv_r_mean),log10(Pv_r_mean),nbin,0);               % Bin the log-transformed data
    fit_r_mean = robustfit(B_r_mean(:,1),B_r_mean(:,2));                    % Get the response just below of the roll-off frequency
    Log_p1=log10(10^fit_r_mean(1)*(10.^f_seg_roll_off).^fit_r_mean(2));     % Using line fitted to the binned data. Use the same f than Log_p2
    MSE_H(j)=immse(Log_p1,Log_p2);
end
delete(poolobj)

[~, Pos_MSE_H_min]=min(MSE_H);                                              % Find the lowest mse (best H)
H_best=H_Vector(Pos_MSE_H_min);

figure
plot(H_Vector,MSE_H,'--ko'); xlabel('H (Roughness)'); ylabel('MSE (m^4)');  % Plot H vs R2
saveas(gcf,[Path_export_Fig '\H_MSE.bmp']); close


%% Get 1D power spectrum for the best H

Pv_r_mean=0; fv_r_mean=0;
for i=1:Iter
    Z_r=createFractalTerrain(tsize, range, H_best);                         % Get the random DEM with diamond-square algorithm  
    Z_r=Z_r(1:size(Z,1),1:size(Z,2));                                       % Get the same dimensions of the observed DEM
    
    [~, ~, Pv_r, fv_r] = fft2D(Z_r,dx,dy,pad,window);                       % Get the 1D power spectrum
    Pv_r_mean=Pv_r_mean+Pv_r; fv_r_mean=fv_r_mean+fv_r;
end
Pv_r_mean=Pv_r_mean./Iter;
fv_r_mean=fv_r_mean./Iter;
figure
[axf_r_mean, ~] = SpecPlot1D(fv_r_mean,Pv_r_mean);
nbin = 20;                                                                  % Number of bins
B_r_mean = bin(log10(fv_r_mean),log10(Pv_r_mean),nbin,0);                   % Bin the log-transformed data

hold(axf_r_mean,'on');                                                      % Plot the binned values
plot(axf_r_mean,10.^B_r_mean(:,1),10.^B_r_mean(:,2),'ok','markerfacecolor','w'); 
fit_r_mean = robustfit(B_r_mean(:,1),B_r_mean(:,2));
plot(axf_r_mean,10.^B_r_mean(:,1),10^fit_r_mean(1)*(10.^B_r_mean(:,1)).^fit_r_mean(2),'k');
saveas(gcf,[Path_export_Fig '\Average_Power_Spectrum.bmp']); close;
plot(axf,10.^B_r_mean(:,1),10.^B_r_mean(:,2),'--bo');
saveas(gcf,[Path_export_Fig '\Power_Spectrum_2D.bmp']); close
saveas(gcf,[Path_export_Fig '\Power_Spectrum_1D.bmp']); close

%% Plot the normalized 2D spectrum and get the p-values

% Get the p-values, with the normalized power spectrum using the mean spectrum from the random DEMs with the H_best
if Method_mean_spectrum==1
    Pmn_r = Pm./(10^fit_r_mean(1)*fm.^fit_r_mean(2));                       % p-values from random DEMs
else
    Pmn_r=Pmn;                                                              % Get the values from the mean spectrum
end


figure;
scatter(fm(:),Pmn_r(:),5,[.5 .5 .5]); hold on;                              % Plot Normalized Power
xlabel('Radial frequency (m^{-1})'); ylabel('Normalized power');

V=2;                                                                        % Degrees of freedom
chi_threshold = chi2inv(p_value_t,V)./2;                                    % Get the chi square threshold for the significance levels
plot([10^-4 1],[chi_threshold chi_threshold],'k');
set(gca,'Xscale','log');
saveas(gcf,[Path_export_Fig '\Normalized_Power_1D.bmp']); close
% Plot
p_matrix = chi2cdf(2.*Pmn_r,V);                                             % Get the significance levels
p_matrix=p_matrix(2:end,2:end);                                             % Remove the first column and row in order to get the right position
[nfy, nfx] = size(fm);
nyq = fm(nfy/2+1,1);                                                        % The Nyquist frequency in the x direction. We pad the data, so it is a square matrix, then nyq=nxq
Vector_freqX=-nyq:2*nyq/(nfx-2):nyq;
Vector_freqY=nyq:-2*nyq/(nfx-2):-nyq;
[x_contour,y_contour] = meshgrid(Vector_freqX,Vector_freqY);
figure; contour(x_contour,y_contour,p_matrix,5); colorbar;                  % Plot the p_values based on the origincal DFT resolution
xlabel('x frequency'); ylabel('y frequency'); grid on
saveas(gcf,[Path_export_Fig '\Raw_P_values_2D.bmp']); close
%% Select the dominant frequencies based on the p-value

if f_cut~=0   
    [~, pos_cut]=min(abs(abs(x_contour(1,1:floor(size(x_contour,2)/2)))-f_cut)); % Find the position of the closest freqeuncy in the fx or fy matrix (x_contour or y_contour)
    p_Matrix_cut=p_matrix(pos_cut:end-pos_cut+1,pos_cut:end-pos_cut+1);
    x_contour_cut=x_contour(pos_cut:end-pos_cut+1,pos_cut:end-pos_cut+1);
    y_contour_cut=y_contour(pos_cut:end-pos_cut+1,pos_cut:end-pos_cut+1);
    P_value_cut.p_Matrix_cut=p_Matrix_cut; 
    P_value_cut.x_contour_cut=x_contour_cut; 
    P_value_cut.y_contour_cut=y_contour_cut;
    figure
    contourf(x_contour_cut,y_contour_cut,p_Matrix_cut,5); colorbar
    saveas(gcf,[Path_export_Fig '\Raw_P_values_2D_Range1.bmp']); close
    
    Vector_freqX_interp2=-f_cut:(abs(2*f_cut))/1000:f_cut;                  % Refine matrix with 1000 intervals between -f_cut and f_cut using spline interpolation
    Vector_freqY_interp2=-f_cut:(abs(2*f_cut))/1000:f_cut;
    [x_interp2,y_interp2] = meshgrid(Vector_freqX_interp2,Vector_freqY_interp2);
    p_matrix_f=interp2(x_contour_cut,y_contour_cut,p_Matrix_cut,x_interp2,y_interp2,'spline');
    p_matrix_f(p_matrix_f<p_value_t)=NaN;                                   % Remove the undesire frequencies
    p_matrix_f(p_matrix_f>1)=0.9999;                                        % This is an artifact caused by the interpolation method
    
    figure
    contourf(x_interp2,y_interp2,p_matrix_f,5); colorbar;                   % Plot results in specific range
    xlabel('x frequency'); ylabel('y frequency'); grid on; hold on
    % Selection of Frequencies
    Loc_N=(size(p_matrix_f,2)+1)/2;
    p_matrix_df=p_matrix_f(:,Loc_N+1:end);                                  % Look for dominant frequencies in the I and IV quadrant, we can ignore II and III quadrants because of Symmetry
    x_contour_df=x_interp2(:,Loc_N+1:end);
    y_contour_df=y_interp2(:,Loc_N+1:end);
    fx_dom_cut=[]; fy_dom_cut=[]; fxy_dom_cut=[]; 
    P_values_cut=[];
    N_cut=floor(Prop*N_fxy);                                                % This will be the Propoportion of frequencies of the total frequencies sampled between [0 f_cut]
        
    [row, col] = find(~isnan(p_matrix_df));
    Random_Selection=randperm(size(row,1),N_cut);                           %Selection based on a random selection for values learger than the p-value threshold
    for i=1:N_cut
        fx_dom_cut(i,1)=x_contour_df(row(Random_Selection(i)), col(Random_Selection(i)));
        fy_dom_cut(i,1)=y_contour_df(row(Random_Selection(i)), col(Random_Selection(i)));
        P_values_cut(i,1)=p_matrix_df(row(Random_Selection(i)), col(Random_Selection(i))); 
    end

    hold on
    scatter3(fx_dom_cut, fy_dom_cut, ones(size(fx_dom_cut,1),1),'filled');
    saveas(gcf,[Path_export_Fig '\Selected_Frequencies_Range1.bmp']); close
     
    p_matrix_f=p_matrix;
    p_matrix_f_median=p_matrix_f;                                           % Do the selection of frequencies for the rest of frequencies within the interval [f_cut f_max] 
    
    Vector_freqX_interp2=x_contour(1,1):(-x_contour(1,1)-x_contour(1,1))/2000:-x_contour(1,1); % Refine matrix with 2000 intervals between -f_cut and f_cut using spline interpolation
    Vector_freqY_interp2=x_contour(1,1):(-x_contour(1,1)-x_contour(1,1))/2000:-x_contour(1,1);
    [x_interp2,y_interp2] = meshgrid(Vector_freqX_interp2,Vector_freqY_interp2);
    p_matrix_f=interp2(x_contour,y_contour,p_matrix_f,x_interp2,y_interp2,'spline');
    p_matrix_f(p_matrix_f<p_value_t)=NaN;
    % Plot results with just the median filter 
    figure
    p_matrix_f_median(p_matrix_f_median<p_value_t)=NaN;
    contourf(x_contour,y_contour,p_matrix_f_median,5); colorbar
    xlabel('x frequency'); ylabel('y frequency'); grid on; hold on
    Loc_N=(size(p_matrix_f,2)+1)/2;
    p_matrix_df=p_matrix_f(:,Loc_N+1:end);                                  % Look for dominant frequencies in the I and IV quadrant, we can ignore II and III quadrants because of Symmetry
    x_contour_df=x_interp2(:,Loc_N+1:end);
    y_contour_df=y_interp2(:,Loc_N+1:end);
    fx_dom=[]; fy_dom=[]; fxy_dom=[]; 
    P_values=[];
    N_cut=floor((1-Prop)*N_fxy);                                            % This will be the (1-Prop)% of frequencies of the total frequencies

    
    [row, col] = find(~isnan(p_matrix_df) & x_contour_df>-f_max & x_contour_df<f_max & y_contour_df>-f_max & y_contour_df<f_max);
    Random_Selection=randperm(size(row,1),N_cut);
    for i=1:N_cut
    fx_dom(i,1)=x_contour_df(row(Random_Selection(i)), col(Random_Selection(i)));
    fy_dom(i,1)=y_contour_df(row(Random_Selection(i)), col(Random_Selection(i)));
    P_values(i,1)=p_matrix_df(row(Random_Selection(i)), col(Random_Selection(i))); % Check the selected p_values
    end

    hold on
    scatter3(fx_dom, fy_dom, ones(size(fx_dom,1),1),'filled');
    saveas(gcf,[Path_export_Fig '\Selected_Frequencies_Range2.bmp']); close
     
    fx_dom=[fx_dom_cut; fx_dom]; % Merge the final selected frequencies and p_values
    fy_dom=[fy_dom_cut; fy_dom];
    P_values=[P_values_cut; P_values];  
    
elseif f_cut==0
    p_matrix_f=p_matrix;
    p_matrix_f_median=p_matrix_f;
    
    Vector_freqX_interp2=x_contour(1,1):(-x_contour(1,1)-x_contour(1,1))/2000:-x_contour(1,1);  % Refine mesh
    Vector_freqY_interp2=x_contour(1,1):(-x_contour(1,1)-x_contour(1,1))/2000:-x_contour(1,1);
    [x_interp2,y_interp2] = meshgrid(Vector_freqX_interp2,Vector_freqY_interp2);
    p_matrix_f=interp2(x_contour,y_contour,p_matrix_f,x_interp2,y_interp2,'spline');
    p_matrix_f(p_matrix_f<p_value_t)=NaN;
    % Plot results with just the median filter 
    figure
    p_matrix_f_median(p_matrix_f_median<p_value_t)=NaN;
    contourf(x_contour,y_contour,p_matrix_f_median,5); colorbar
    xlabel('x frequency'); ylabel('y frequency'); grid on; hold on
    Loc_N=(size(p_matrix_f,2)+1)/2;
    p_matrix_df=p_matrix_f(:,Loc_N+1:end);                                  % Look for dominant frequencies in the I and IV quadrant, we can ignore II and III quadrants because of Symmetry
    x_contour_df=x_interp2(:,Loc_N+1:end);
    y_contour_df=y_interp2(:,Loc_N+1:end);
    fx_dom=[]; fy_dom=[];
    P_values=[];
    N_cut=N_fxy;
    
    [row, col] = find(~isnan(p_matrix_df) & x_contour_df>-f_max & x_contour_df<f_max & y_contour_df>-f_max & y_contour_df<f_max);
    Random_Selection=randperm(size(row,1),N_cut);
    for i=1:N_cut
        fx_dom(i,1)=x_contour_df(row(Random_Selection(i)), col(Random_Selection(i)));
        fy_dom(i,1)=y_contour_df(row(Random_Selection(i)), col(Random_Selection(i)));
        P_values(i,1)=p_matrix_df(row(Random_Selection(i)), col(Random_Selection(i))); % Check the selected p_values
    end
        
    hold on
    scatter3(fx_dom, fy_dom, ones(size(fx_dom,1),1),'filled');
    saveas(gcf,[Path_export_Fig '\Selected_Frequencies.bmp']); close   
end

% Save results
P_Value_Results.P_values=P_values;
P_Value_Results.Pm=Pm;
P_Value_Results.fm=fm;
P_Value_Results.x_contour=x_contour;
P_Value_Results.y_contour=y_contour;
P_Value_Results.p_matrix=p_matrix;
    if f_cut~=0   
        P_Value_Results.P_value_cut=P_value_cut;
    end


end



















