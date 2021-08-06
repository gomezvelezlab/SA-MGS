%% Examples to execute the SS_DomFreq, SS_Top ans SS_Vel for the Groundwater 
% Spectral Solutions presented by Perez & Gomez-Velez 2020 (WRR submitted).

%   Copyright:  Gabriel Perez
%   Gomez-Velez Lab Group. Vanderbilt University

%   email:   gabriel.perez.mesa@vanderbilt.edu
%	06  Aug	2021,   MATLAB	2019b  version
%	IF  YOU	PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%	Perez. G., et al., (2021)  Identification of Characteristic Spatial Scales to Improve the Performance 
%   of Analytical Spectral Solutions to the Groundwater Flow Equation (Submitted to WRR)

%   Input
%   Test=1 --> Tothian Basin
%   Test=2 --> Sandy Bedforms
%   Test=3 --> Rio Hondo Basin
%   Test=4 --> Random DEM

%% Load surface Examples with specific parameters
% The parameters are conditioned by the spatial scale of the DEM. 
% Note: The parameters p_value_t, f_cut, f_max, Prop, and N_fxy are defined by the
% User

Test=3;

if Test==1
    load([pwd '\Examples\Example_Tothian\Surf_Tothian']); 
    p_value_t=0.95; f_cut=0.0004; f_max=0.1; Prop=0.999; N_fxy=392; 
    Path_export_Fig=[pwd '\Examples\Example_Tothian'];                     % Path to save the figures
    
elseif Test==2
    load([pwd '\Examples\Example_Bedforms\Surf_Bedforms']);
    p_value_t=0.95; f_cut=30; f_max=50; Prop=0.999; N_fxy=392;
    Path_export_Fig=[pwd '\Examples\Example_Bedforms'];                    % Path to save the figures
elseif Test==3
    load([pwd '\Examples\Example_RioHondo\Surf_RioHondo']);
    p_value_t=0.5; f_cut=0.0004; f_max=0.004; Prop=0.7; N_fxy=2450;
    Path_export_Fig=[pwd '\Examples\Example_RioHondo'];                    % Path to save the figures
elseif Test==4
    load([pwd '\Examples\Example_RandomDEM\Surf_RandomDEM']);
    p_value_t=0.5; f_cut=0.01; f_max=0.1; Prop=0.7; N_fxy=392;
    Path_export_Fig=[pwd '\Examples\Example_RandomDEM'];                   % Path to save the figures
    zm=zm+50;
end

% General Parameters
                                                                           % Number of pairs (fx and fy) of dominat frequencies to be sampled
Method_mean_spectrum=1;                                                    % 1: Uses Random DEMs to get the average spectrum, 0:uses the average spectrum from the observed DEM
Iter=10;                                                                   % Number of iterations in the estimation of the average power spectrum

%% Selection of Dominant Frequencies
Z=zm;                                                                      % A NxM Matrix containing the DEM in meters
Step=xm(1,2)-xm(1,1);                                                      % DEM resolution in meters
[fx_dom,fy_dom,P_Value_Results]=SS_DomFreq(Z,Step,N_fxy,p_value_t,Method_mean_spectrum,f_cut,f_max,Prop,Iter,Path_export_Fig);

% Include dominant frequencies in order to capture the regional slopes in x
% and y direction. This can be aproximated by wavelenghts couple of 5 and
% 10 times the domain lenght. This step is only needed if there is an
% important regional slope on the domain.

X_lenght=size(Z,1)*Step;
Y_lenght=size(Z,2)*Step;
fx_dom(end-1,1)=1/(5*X_lenght); fy_dom(end-1,1)=1/(10*X_lenght);
fx_dom(end,1)=1/(10*Y_lenght); fy_dom(end,1)=1/(5*Y_lenght);

%% Estimation of the water table using Dominant Frequencies
X=xm(1,1:end); Y=ym(1:end,1)'; Z=zm;
Vector_wavelenghts=[1./fx_dom 1./fy_dom]; 
StepWise_Flag=1;

[paramhat,kx,ky,h_top,mean_error,av]=SS_Top(X,Y,Z,Vector_wavelenghts,StepWise_Flag);

%% Estimation of the Pressure head and velocity field
fact=1;                                                                     % Multiples of grid step used in X and Y for computational efficiency
N_Z=3;                                                                      % Number of grid steps with depth
Cond=1e-4;                                                                  % Hydraulic Conductivity m/s, for Dunes use 9.8067e-4 m/s, for Tot use 3.2*10^-6 m/s, For Rio Hondo use 1e-4;
poro=0.3;                                                                   % Porosity
z_min=0;
z_max=-100;                                                                 % Max Depth in meters in the output
dp=3*max(X_lenght,Y_lenght);                                                % Depth flow domain in meters. Recommended to be 2 or 3 times the max domain lenght. For Dunes use 2m, for Rio Hondo use 1600m. For others use 3*max(X_lenght,Y_lenght)

[U,V,W,x,y,z,h]=SS_Vel(X,Y,fact,N_Z,av,kx,ky,paramhat,Cond,poro,z_min,z_max,dp);

%% Plot some results

% Power spectrum
figure; set(gcf,'color','white') 
f_rad=sqrt((kx./(2*pi)).^2+(ky./(2*pi)).^2); f_rad=[f_rad f_rad];
pdf_2=histogram(abs(paramhat)./(1./f_rad),100);
scatter((pdf_2.BinEdges(1:end-1)+pdf_2.BinEdges(2:end))./2,pdf_2.BinCounts./sum(pdf_2.BinCounts),30,'MarkerEdgeColor',[0 0 0],...
    'MarkerFaceColor',[.7 .7 .7],'LineWidth',.5);
set(gca,'yscale','log'); set(gca,'xscale','log');
xlabel('|h_m|/\lambda (m)'); ylabel('Probability Density'); hold on; box on; grid on; set(gca,'FontSize',12);

% Reconstructed Water Table
figure; set(gcf,'color','white') 
subplot(1,2,1)
surf(xm, ym, zm); shading interp;  camlight left;  lighting gouraud;    hold on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z Data (m) ');
subplot(1,2,2)
surf(xm, ym, h_top); shading interp;  camlight left;  lighting gouraud;    hold on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z Fourier (m)');

% Darcy Velocities on Top Layer
figure; set(gcf,'color','white') 
ColorScale_R=[0,109,44;49,163,84;116,196,118;186,228,179]./256; % 
ColorScale_D=[239,243,255;189,215,231;107,174,214;8,81,156]./256; %  
ColorScale=[ColorScale_R;ColorScale_D];
subplot(1,3,1)
contourf(xm, ym,poro.*U(:,:,1),'edgecolor','none'); xlabel('X (m)'); ylabel('Y (m)'); ax1=gca; colormap(ax1,ColorScale);
c=colorbar('northoutside'); c.Label.String = 'U [m/s]';  caxis([-max(max(abs(poro.*U(:,:,1)))) max(max(abs(poro.*U(:,:,1))))]);
subplot(1,3,2)
contourf(xm, ym,poro.*V(:,:,1),'edgecolor','none'); xlabel('X (m)'); ylabel('Y (m)'); ax1=gca; colormap(ax1,ColorScale);
c=colorbar('northoutside'); c.Label.String = 'V [m/s]';   caxis([-max(max(abs(poro.*V(:,:,1)))) max(max(abs(poro.*V(:,:,1))))]);
subplot(1,3,3)
contourf(xm, ym,poro.*W(:,:,1),'edgecolor','none'); xlabel('X (m)'); ylabel('Y (m)'); ax1=gca; colormap(ax1,ColorScale);
c=colorbar('northoutside'); c.Label.String = 'W [m/s]';   caxis([-max(max(abs(poro.*W(:,:,1)))) max(max(abs(poro.*W(:,:,1))))]);





