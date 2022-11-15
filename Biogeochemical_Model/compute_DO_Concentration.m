%% Function to calculate the DO concentration. Taken from Marzadri 2011
% Coded by Gabriel Perez
% gabriel.perez.mesa@vanderbilt.edu
% Date: 11/01/2022

% Inputs:
% C_0_s: Initial concentration of dissolved oxygen in the stream [mg/L]
% C_0_lim: DO trheshold where microbial respiration is inhibited [mg/L]
% KR: Rate coeffcient of biomass respiration [1/d]
% KN: Rate coefficient of nitrification  [1/d]
% t: time [d]. It has lenght [N]
% tau: Residence time [d]. It has lenght [M]
% DL: Disperssion coefficient [m2/d]
% u: darcy velocity [m/d]

% Output
% C_0_matrix: DO concentration at a given t and tau. NxM matrix
% tau_lim: time location where there is a change of aerobic and anaerobic conditions [d]

function [C_0_matrix, tau_lim]=compute_DO_Concentration(C_0_s,C_0_lim,KR,KN,t,tau,DL,u)
sympref('HeavisideAtOrigin',0); % Set Heaviside function equal to zero for heaviside(0)=0;
[t_matrix,tau_matrix]=meshgrid(t,tau);

KRN=KR+KN;
tau_lim=(1/KRN)*log(C_0_s/C_0_lim);

if DL==0

    % For t>tau
    C_0_matrix=C_0_s.*exp(-KRN.*tau_matrix);
    % For t<tau
    C_0_matrix_2=C_0_lim.*exp(-KRN.*t_matrix);
    C_0_matrix(t_matrix<tau_matrix)=C_0_matrix_2(t_matrix<tau_matrix);
%     C_0_matrix(tau_matrix>tau_lim)=C_0_lim; % DO is not longer depleted below the threshold value
%     C_0_matrix(C_0_matrix<C_0_lim)=C_0_lim; % DO is not longer depleted below the threshold value

else
    C_0_matrix=C_0_lim.*exp(-KRN.*t_matrix).*(1-0.5.*erfc((u.*tau_matrix-u.*t_matrix)./(2.*sqrt(DL.*t_matrix)))-0.5.*exp(u.^2.*tau_matrix./DL).*erfc((u.*tau_matrix+u.*t_matrix)./(2.*sqrt(DL.*t_matrix))))+...
        0.5.*C_0_s.*(exp((u.^2.*tau_matrix./(2*DL)).*(1+sqrt(1+(4*KN*DL/(u^2))))).*erfc((u.*tau_matrix+u.*t_matrix.*sqrt(1+(4.*KRN*DL/(u^2))))./(2*sqrt(DL.*t_matrix)))+...
        exp((u^2*tau_matrix./(2*DL)).*(1-sqrt(1+4*KRN*DL./(u^2)))).*erfc((u.*tau_matrix-u.*t_matrix.*sqrt(1+4*KRN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))     );
end


end







