%% Function to calculate the NH4 concentration. Taken from Marzadri 2011
% Coded by Gabriel Perez
% gabriel.perez.mesa@vanderbilt.edu
% Date: 11/01/2022

% Inputs:
% C_1_s: Initial concentration of NH4 in the stream [mg/L]
% R1= Retardation coefficient for NH4
% KN: % Rate coefficient of nitrification  [1/d]
% t: time [d]. It has lenght [N]
% tau: Residence time [d]. It has lenght [M]
% tau_lim: time location where there is a change of aerobic and anaerobic conditions [d]
% Pulse: 1=Uses an instantaneous pulse of C_0_s at t=0, and zero for t>0
% Pulse: 0=Uses a constant and steady C_O_s suring the entire simulation
% DL: Disperssion coefficient [m2/d]
% u: darcy velocity [m/d]

% Output
% C_1_matrix: NH4 concentration at a given t and tau [mg/L] NxM matrix
% N: number of t steps
% M: number of tau steps

function [C_1_matrix]=compute_NH4_Concentration(C_0_s,C_1_s,C_0_lim,tau_lim,KR,KN,t,tau,Pulse,DL,u)

[t_matrix,tau_matrix]=meshgrid(t,tau);
if Pulse==1
    C_1_s=C_1_s.*((t_matrix-tau_matrix)==0);
else
    C_1_s=C_1_s.*ones(size(t_matrix));
end

if DL==0

    % For t>=tau and tau<tau_lim
    C_1_matrix=C_1_s.*exp(-KN.*tau_matrix);
    % For t>=tau and tau>tau_lim
    C_1_matrix_an=C_1_s.*exp(-KN.*tau_lim);
    Flag=(t_matrix>=tau_matrix) & (tau_matrix>tau_lim);
    C_1_matrix(Flag)=C_1_matrix_an(Flag);
    % For t<tau
    C_1_matrix(t_matrix<tau_matrix)=0;

else

    % Approach using tau_lim as a constant

    %     % For tau<tau_lim
    %     C_1_ox=exp((u^2*tau_matrix./(2*DL)).*(1+sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_matrix+u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))+...
    %         exp((u^2*tau_matrix./(2*DL)).*(1-sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_matrix-u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)));
    %
    %     C_1_matrix=0.5.*C_1_s.*C_1_ox;
    %     % For tau>tau_lim
    %
    %     C_1_matrix_an=0.5.*C_1_s.*(exp((u^2*tau_lim./(2*DL)).*(1+sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_lim+u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))+...
    %         exp((u^2*tau_lim./(2*DL)).*(1-sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_lim-u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))); % This is corrected!
    %     C_1_matrix(tau_matrix>tau_lim)=C_1_matrix_an(tau_matrix>tau_lim);


    % Approach using tau_lim as function of time
    KRN=KR+KN;
    C_0_matrix=C_0_lim.*exp(-KRN.*t_matrix).*(1-0.5.*erfc((u.*tau_matrix-u.*t_matrix)./(2.*sqrt(DL.*t_matrix)))-0.5.*exp(u.^2.*tau_matrix./DL).*erfc((u.*tau_matrix+u.*t_matrix)./(2.*sqrt(DL.*t_matrix))))+...
        0.5.*C_0_s.*(exp((u.^2.*tau_matrix./(2*DL)).*(1+sqrt(1+(4*KN*DL/(u^2))))).*erfc((u.*tau_matrix+u.*t_matrix.*sqrt(1+(4.*KRN*DL/(u^2))))./(2*sqrt(DL.*t_matrix)))+...
        exp((u^2*tau_matrix./(2*DL)).*(1-sqrt(1+4*KRN*DL./(u^2)))).*erfc((u.*tau_matrix-u.*t_matrix.*sqrt(1+4*KRN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))     );

    [~, Pos]=min((C_0_matrix-C_0_lim).^2,[],1,'linear');
    tau_lim_vector=tau_matrix(Pos);
    [tau_lim_matrix,tau_matrix]=meshgrid(tau_lim_vector,tau);


    % For tau<tau_lim
    C_1_ox=exp((u.^2.*tau_matrix./(2*DL)).*(1+sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_matrix+u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))+...
        exp((u^2*tau_matrix./(2*DL)).*(1-sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_matrix-u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)));

    C_1_matrix=0.5.*C_1_s.*C_1_ox;
    % For tau>tau_lim

    C_1_matrix_an=0.5.*C_1_s.*(exp((u^2*tau_lim_matrix./(2*DL)).*(1+sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_lim_matrix+u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))+...
        exp((u^2*tau_lim_matrix./(2*DL)).*(1-sqrt(1+4*KN*DL./(u^2)))).*erfc((u.*tau_lim_matrix-u.*t_matrix.*sqrt(1+4*KN*DL/(u^2)))./(2*sqrt(DL.*t_matrix)))); % This is corrected!

    C_1_matrix(tau_matrix>tau_lim_matrix)=C_1_matrix_an(tau_matrix>tau_lim_matrix);








end

end
