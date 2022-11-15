% Comparison Between Numerical and Analytical solutions

%% Define parameters
Add_Dispersion=0;

if Add_Dispersion==1
    Dispersion_scaling=1; % scaling factor to explore the effect of dispersion [-] Use 1
    Stepwise=1; % Do the PDE solution in two segments
else
    Dispersion_scaling=1e-3; % scaling factor to explore the effect of dispersion [-] Use 1
    Stepwise=0;
end

% Porous media parameters
tau_max=10; % Maximum residence time within the streamline [d] 
t_max=20; % Maximum simulation time [d]

% Parameters to calculate disperssion coefficient
t_max_dispersion=10; % Maximum simulation time [d]
porosity = 0.3; % Sediment porosity [-]
q = 0.0097; % Darcy flux [m/d]. For sandy bedforms the average q in COMSOL is q=0.0097 m/d
L=(q/porosity)*t_max_dispersion; % Characteristic length [m]
alpha_L=L/20; % Longitudinal dispersivity [m]
DM = 1e-9*24*60*60; % molecular difussion [m2/d]
DL=Dispersion_scaling*(alpha_L*(q/porosity)+DM);
u=q/porosity; % velocity [m/d]

% Retardation coefficients
R_0 = 1; % Retardation factor for dissolved oxygen [-]
R_1 = 1; % Retardation factor for NH4 (Ammonium) [-]
R_2 = 1; % Retardation factor for NO3 (Nitrate) [-]
R_3 = 1; % Retardation factor for N2gas (Nitrogen) [-]

% Stream concentrations
C_0_s = 10; % Stream concentration of dissolved oxygen in the stream [mg/L]. We assumed a constant value for all the t
C_1_s = 5.46; %5.46; % Stream concentration of NH4 (Ammonium) in the stream [mg/L]
C_2_s = 1.325; % Stream concentration of NO3 (Nitrate) in the stream [mg/L]
C_3_s = 0; % Stream concentration of N2gas (Nitrogen) in the stream [mg/L]

C_0_lim= 3; % DO trheshold where microbial respiration is inhibited [mg/L]

% Initial concentrations within sediment
C_0_0 = C_0_lim; % Initial concentration of dissolved oxygen [mg/L]. We assumed a constant value for all the t
C_1_0 = 0; % Initial concentration of NH4 (Ammonium)  [mg/L]
C_2_0 = 0; % Initial concentration of NO3 (Nitrate)  [mg/L]
C_3_0 = 0; % Initial concentration of N2gas (Nitrogen)  [mg/L]


%% Correct reaction rates for temperature  
% Rate coeffcients
Arrhenius=1;
if(Arrhenius == 1)
    T=6; % Water temperature [C]
    KR_20=0.1; % Rate coefficient for R at 20C [1/d];
    KN_20=3.46; % Rate coefficient for N at 20C [1/d];
    KD_20=1.65; % Rate coefficient for D at 20C [1/d];
    KC_20=1.0; % Rate coefficient for C at 20C [1/d];
    Phi_R=1.047; % Dimesnionless temperature coefficient for R;
    Phi_N=1.040; % Dimesnionless temperature coefficient for N;
    Phi_D=1.045; % Dimesnionless temperature coefficient for D;
    Phi_C=1.047; % Dimesnionless temperature coefficient for C;
    KR= KR_20*Phi_R^(T-20); % Rate coeffcient of biomass respiration [1/d]
    KN= KN_20*Phi_N^(T-20); % Rate coefficient of nitrification  [1/d]
    KD= KD_20*Phi_D^(T-20); % Rate coefficient of denitrification  [1/d]
    KC= KC_20*Phi_C^(T-20); % Rate coefficient of biomass uptake due to microbial assimilatory reduction of nitrate  [1/d]
    KRN=KR+KN;
else
    % Values taken from Marzadri 2011, Table 3 Site A1
    KR= 0.053; % Rate coeffcient of biomass respiration [1/d]
    KN= 9.903; % Rate coefficient of nitrification  [1/d]
    KD= 2.922; % Rate coefficient of denitrification  [1/d]
    KC= 0.523; % Rate coefficient of biomass uptake due to microbial assimilatory reduction of nitrate  [1/d]
    KRN=KR+KN;
end


%% Analytical solution
if Add_Dispersion==0
    DL=0; 
end
Pulse=0;
tau = [0, 10.^(-4:0.005:log10(tau_max))]; % [d]
t = [0, 10.^(-4:0.005:log10(t_max))];
[C_0_path, tau_lim]=compute_DO_Concentration(C_0_s,C_0_lim,KR,KN,t,tau,DL,u);
[C_1_path]=compute_NH4_Concentration(C_0_s,C_1_s,C_0_lim,tau_lim,KR,KN,t,tau,Pulse,DL,u);
[C_2_path]=compute_NO3_Concentration(C_0_s,C_1_s,C_2_s,C_0_lim,tau_lim,KR,KN,KC,KD,t,tau,Pulse,DL,u);
[C_3_path]=compute_N2gas_Concentration(C_0_s,C_1_s,C_2_s,C_3_s,C_0_lim,tau_lim,KR,KN,KC,KD,t,tau,Pulse,DL,u);

%% Get the tau_lim from the analytical solution. For dispersion case, tau_lim depends on time, however, here we assumed that tau_lim is constant, and approaximated as the
% tau_lim obtained from the analytical solution at the last time t.
pos_tau_lim= find(C_3_path(:,end)==0);
tau_lim=tau(pos_tau_lim(end));

%% Numerical solution
[sol]=compute_PDE_Concentration(Dispersion_scaling,alpha_L,q,porosity,DM,KRN, KR, KN, KC, KD,C_0_lim, C_0_s, C_1_s, C_2_s, C_3_s, C_0_0, C_1_0, C_2_0, C_3_0, R_0, R_1, R_2, R_3,tau_max,t_max,tau_lim,Stepwise);
tau_an= 10.^(log10(tau_lim):0.005:log10(tau_max)); % [d]

%% Some Figures

Time_plot=[0.1 1 10 50]; % Days
for i=1:length(Time_plot)
    [val,idx]=min(abs(t-Time_plot(i)));
    idx_plot(i)=idx;
end

%
figure; set(gcf,'color','white');
for i=1:4
    Time=idx_plot(i);
    subplot(2,2,i);
    scatter(tau(1:20:end), (sol(Time,(1:20:end),1)/C_0_s)', 'LineWidth',1, 'MarkerFaceColor','k', 'MarkerEdgeColor','k'); hold on;
    scatter(tau(1:20:end), (sol(Time,(1:20:end),2)/C_1_s)', 'LineWidth',1, 'MarkerFaceColor',[0.00,0.45,0.74], 'MarkerEdgeColor','k'); 
    scatter(tau(1:20:end), (sol(Time,(1:20:end),3)/C_2_s)', 'LineWidth',1, 'MarkerFaceColor',[0.85,0.33,0.10], 'MarkerEdgeColor','k');
    scatter(tau(1:20:end), (sol(Time,(1:20:end),4)), 'LineWidth',1, 'MarkerFaceColor',[0.93,0.69,0.13], 'MarkerEdgeColor','k');
    %scatter(tau_an(1:20:end), (sol_an(Time,(1:20:end))), 'LineWidth',1, 'MarkerFaceColor',[0.93,0.69,0.13], 'MarkerEdgeColor','k');
    plot(tau', (C_0_path(:,Time)/C_0_s)', 'LineWidth',2,'Color','k'); 
    plot(tau', (C_1_path(:,Time)/C_1_s)', 'LineWidth',2,'Color',[0.00,0.45,0.74]); 
    plot(tau', (C_2_path(:,Time)/C_2_s)', 'LineWidth',2,'Color',[0.85,0.33,0.10]);
    plot(tau', (C_3_path(:,Time)), 'LineWidth',2,'Color',[0.93,0.69,0.13]);
    plot(tau',ones(length(tau),1).*C_0_lim/C_0_s,'--k');
    xlabel('\tau [days]'); ylabel('C/{C_0} [-] and N_2 [mg/l]'); set(gca,'Xscale','log')
    xlim([1e-3 10]);  set(gca,'Xscale','log'); grid on;
    title(['Solution at t=' num2str(t(Time),2) 'days']); set(gca,'FontSize',12)
    if i==1 
        legend({'Numerical DO'; 'Numerical NH$^+_4$';'Numerical NO$^-_3$';'Numerical N$_2$'  ; 'Analytical DO';  'Analytical NH$^+_4$';'Analytical NO$^-_3$';'Analytical N$_2$'},'interpreter','latex')
    end
end







