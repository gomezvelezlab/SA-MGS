%% Function to calculate the NH3 concentration. Taken from Marzadri 2011
% Coded by Gabriel Perez
% gabriel.perez.mesa@vanderbilt.edu
% Date: 11/01/2022

%% Note: The model presented in Marzadri 2011 and Marzadri 2012 is different.
% There are discrepancies between the equations!

% Inputs:
% C_1_s: Initial concentration of NH4 in the stream [mg/L]
% C_2_s: Initial concentration of NO3 in the stream [mg/L]
% KN: % Rate coefficient of nitrification  [1/d]
% KC: Rate coefficient of biomass uptake due to microbial assimilatory reduction of nitrate  [1/d]
% KD: Rate coefficient of denitrification  [1/d]
% t: time [d]. It has lenght [N]
% tau: Residence time [d]. It has lenght [M]
% tau_lim: time location where there is a change of aerobic and anaerobic conditions [d]
% Pulse: 1=Uses an instantaneous pulse of C_0_s at t=0, and zero for t>0
% Pulse: 0=Uses a constant and steady C_O_s suring the entire simulation
% DL: Disperssion coefficient [m2/d]
% u: darcy velocity [m/d]

% Output
% C_2_matrix: NO3 concentration at a given t and tau [mg/L] NxM matrix
% N: number of t steps
% M: number of tau steps

function [C_2_matrix]=compute_NO3_Concentration(C_0_s,C_1_s,C_2_s,C_0_lim,tau_lim,KR,KN,KC,KD,t,tau,Pulse,DL,u)

[t_matrix,tau_matrix]=meshgrid(t,tau);

if Pulse==1
    C_1_s=C_1_s.*((t_matrix-tau_matrix)==0);
    C_2_s=C_2_s.*((t_matrix-tau_matrix)==0);
else
    C_1_s=C_1_s.*ones(size(t_matrix));
    C_2_s=C_2_s.*ones(size(t_matrix));
end

if DL==0
    % For t>tau and tau<tau_lim
    C_2_matrix=C_2_s.*exp(-KC.*tau_matrix)+C_1_s.*(KN/(KC-KN)).*(exp(-KN.*tau_matrix)-exp(-KC.*tau_matrix));
    % For t>tau and tau>tau_lim
    %C_2_matrix_an=C_2_s.*exp(-KC.*tau_lim-KD.*(tau_matrix-tau_lim))+C_1_s.*(KN/(KC-KN)).*(exp(-KN.*tau_lim)-exp(-KC.*tau_lim-KD.*(tau_matrix-tau_lim))); %This is the original equation by Marzadri 2011, This is incorrect
    C_2_matrix_an=C_2_s.*exp(-KC.*tau_lim-KD.*(tau_matrix-tau_lim))+C_1_s.*(KN/(KC-KN)).*(exp(-KN.*tau_lim-KD.*(tau_matrix-tau_lim))-exp(-KC.*tau_lim-KD.*(tau_matrix-tau_lim))); % Correction!
    C_2_matrix(t_matrix>=tau_matrix & tau_matrix>tau_lim)=C_2_matrix_an(t_matrix>=tau_matrix & tau_matrix>tau_lim);
    % For t<tau
    C_2_matrix(t_matrix<tau_matrix)=0;

else
    % For tau<tau_lim

    % Approach using tau_lim as a constant

    %     C_1_ox=exp((u^2*tau_matrix./(2*DL)).*(1+sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_matrix+u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))+...
    %         exp((u^2*tau_matrix./(2*DL)).*(1-sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_matrix-u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)));
    %
    %     C_1_matrix=0.5.*C_1_s.*C_1_ox;
    %
    %     C_2_ox=exp((u^2*tau_matrix./(2*DL)).*(1+sqrt(1+4*KC*DL./(u^2)))).*erfc((u.*tau_matrix+u.*t_matrix.*sqrt(1+4*KC*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))+...
    %         exp((u^2*tau_matrix./(2*DL)).*(1-sqrt(1+4*KC*DL./(u^2)))).*erfc((u.*tau_matrix-u.*t_matrix.*sqrt(1+4*KC*DL/(u^2)))./(2*sqrt(DL.*t_matrix)));
    %
    %     C_2_matrix=(KN/(KC-KN)).*C_1_matrix+0.5.*(C_2_s-C_1_s.*(KN/(KC-KN))).*C_2_ox; % Corrected. The 0.5 was removed from Equation B2 Marzadri 2012
    %
    %     % For tau>tau_lim
    %
    %    C_2_aox= exp((u^2*tau_matrix./(2*DL)).*(1+sqrt(1+4*KD*DL./(u^2)))).*erfc((u.*tau_matrix+u.*t_matrix.*sqrt(1+4*KD*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))+...
    %             exp((u^2*tau_matrix./(2*DL)).*(1-sqrt(1+4*KD*DL./(u^2)))).*erfc((u.*tau_matrix-u.*t_matrix.*sqrt(1+4*KD*DL/(u^2)))./(2*sqrt(DL.*t_matrix)));
    %
    %     C_1_ox_lim=exp((u^2*tau_lim./(2*DL)).*(1+sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_lim+u.*tau_lim.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*tau_lim)))+...
    %         exp((u^2*tau_lim./(2*DL)).*(1-sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_lim-u.*tau_lim.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*tau_lim)));
    %
    %     C_2_ox_lim=exp((u^2*tau_lim./(2*DL)).*(1+sqrt(1+4*KC*DL./(u^2)))).*erfc((u.*tau_lim+u.*tau_lim.*sqrt(1+4*KC*DL/(u^2)))./(2*sqrt(DL.*tau_lim)))+...
    %         exp((u^2*tau_lim./(2*DL)).*(1-sqrt(1+4*KC*DL./(u^2)))).*erfc((u.*tau_lim-u.*tau_lim.*sqrt(1+4*KC*DL/(u^2)))./(2*sqrt(DL.*tau_lim)));
    %
    %     C_2_lim=0.5.*(KN/(KC-KN)).*C_1_ox_lim+0.5.*(C_2_s-C_1_s.*(KN/(KC-KN))).*C_2_ox_lim;
    %
    %     C_2_matrix_an=0.5.*C_2_aox.*C_2_lim;
    %     C_2_matrix(tau_matrix>tau_lim)=C_2_matrix_an(tau_matrix>tau_lim);


    % Approach using tau_lim as function of time

    KRN=KR+KN;
    C_0_matrix=C_0_lim.*exp(-KRN.*t_matrix).*(1-0.5.*erfc((u.*tau_matrix-u.*t_matrix)./(2.*sqrt(DL.*t_matrix)))-0.5.*exp(u.^2.*tau_matrix./DL).*erfc((u.*tau_matrix+u.*t_matrix)./(2.*sqrt(DL.*t_matrix))))+...
        0.5.*C_0_s.*(exp((u.^2.*tau_matrix./(2*DL)).*(1+sqrt(1+(4*KN*DL/(u^2))))).*erfc((u.*tau_matrix+u.*t_matrix.*sqrt(1+(4.*KRN*DL/(u^2))))./(2*sqrt(DL.*t_matrix)))+...
        exp((u^2*tau_matrix./(2*DL)).*(1-sqrt(1+4*KRN*DL./(u^2)))).*erfc((u.*tau_matrix-u.*t_matrix.*sqrt(1+4*KRN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))     );

    [~, Pos]=min((C_0_matrix-C_0_lim).^2,[],1,'linear');
    tau_lim_vector=tau_matrix(Pos);
    [tau_lim_matrix,tau_matrix]=meshgrid(tau_lim_vector,tau);

    % For tau<tau_lim

    C_1_ox=exp((u^2*tau_matrix./(2*DL)).*(1+sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_matrix+u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))+...
        exp((u^2*tau_matrix./(2*DL)).*(1-sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_matrix-u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)));

    C_1_matrix=0.5.*C_1_s.*C_1_ox;

    C_2_ox=exp((u^2*tau_matrix./(2*DL)).*(1+sqrt(1+4*KC*DL./(u^2)))).*erfc((u.*tau_matrix+u.*t_matrix.*sqrt(1+4*KC*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))+...
        exp((u^2*tau_matrix./(2*DL)).*(1-sqrt(1+4*KC*DL./(u^2)))).*erfc((u.*tau_matrix-u.*t_matrix.*sqrt(1+4*KC*DL/(u^2)))./(2*sqrt(DL.*t_matrix)));

    C_2_matrix=(KN/(KC-KN)).*C_1_matrix+0.5.*(C_2_s-C_1_s.*(KN/(KC-KN))).*C_2_ox; % Corrected. The 0.5 was removed from Equation B2 Marzadri 2012

    % For tau>tau_lim

    C_2_aox= exp((u^2*(tau_matrix-tau_lim_matrix)./(2*DL)).*(1+sqrt(1+4*KD*DL./(u^2)))).*erfc((u.*(tau_matrix-tau_lim_matrix)+u.*t_matrix.*sqrt(1+4*KD*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))+...
        exp((u^2*(tau_matrix-tau_lim_matrix)./(2*DL)).*(1-sqrt(1+4*KD*DL./(u^2)))).*erfc((u.*(tau_matrix-tau_lim_matrix)-u.*t_matrix.*sqrt(1+4*KD*DL/(u^2)))./(2*sqrt(DL.*t_matrix)));

    C_1_ox_lim=exp((u^2*tau_lim_matrix./(2*DL)).*(1+sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_lim_matrix+u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))+...
        exp((u^2*tau_lim_matrix./(2*DL)).*(1-sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_lim_matrix-u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)));

    C_1_matrix_lim=0.5.*C_1_s.*C_1_ox_lim;

    C_2_ox_lim=exp((u^2*tau_lim_matrix./(2*DL)).*(1+sqrt(1+4*KC*DL./(u^2)))).*erfc((u.*tau_lim_matrix+u.*t_matrix.*sqrt(1+4*KC*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))+...
        exp((u^2*tau_lim_matrix./(2*DL)).*(1-sqrt(1+4*KC*DL./(u^2)))).*erfc((u.*tau_lim_matrix-u.*t_matrix.*sqrt(1+4*KC*DL/(u^2)))./(2*sqrt(DL.*t_matrix)));

    C_2_lim=(KN/(KC-KN)).*C_1_matrix_lim+0.5.*(C_2_s-C_1_s.*(KN/(KC-KN))).*C_2_ox_lim;

    C_2_matrix_an=0.5.*C_2_lim.*C_2_aox;

    C_2_matrix(tau_matrix>tau_lim_matrix)=C_2_matrix_an(tau_matrix>tau_lim_matrix);




end


end
