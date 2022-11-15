%% Function to calculate the DO, NH4, NO3, and N2gas using numerical solution of the PDE
% Coded by Gabriel Perez and Jesus Gomez-Velez
% gabriel.perez.mesa@vanderbilt.edu
% Date: 11/01/2022

% This version calculates the nitrogen N2 by performing the individual PDE solutions in the aerobic and anaerobic region.
% This is required since that the numerical solution does not capture the
% "abrupt" change of N2 in the interface between aerobic and anaerobic
% condition


%Stepwise: 1=Solve the PDE in two intervals, interval (1) 0<tau<tau_lim,
%interval (2) tau_lim<tau<tau_tau_max
% 0: Solve the PDE as a continous interval 0<tau<tau_tau_max

function [sol]=compute_PDE_Concentration(Dispersion_scaling,alpha_L,q,porosity,DM,KRN, KR, KN, KC, KD,C_0_lim, C_0_s, C_1_s, C_2_s, C_3_s, C_0_0, C_1_0, C_2_0, C_3_0, R_0, R_1, R_2, R_3,tau_max,t_max,tau_lim,Stepwise)


%% Prepare solver parameters
p = struct('DL',Dispersion_scaling*(alpha_L*(q/porosity)+DM), ...
    'V',q/porosity, ...
    'tau_lim',(1/KRN)*log(C_0_s/C_0_lim), ...
    'K', [KRN, KR, KN, KC, KD], ...
    'cs', [C_0_s; C_1_s; C_2_s; C_3_s], ...
    'c0', [C_0_0; C_1_0; C_2_0; C_3_0], ...
    'C_0_lim', C_0_lim, ...
    'R', [R_0; R_1; R_2; R_3]);

%% Numerical solution for the aerobic region.
% A good approximation for the location in the transition between aerobic
% and anaerobic is based on the purely advective case, where
%tau_lim=(1/KRN)*log(C_0_s/C_0_lim); % [days]

% Numerical parameterization
tau = [0, 10.^(-4:0.005:log10(tau_max))]; % [d]
t = [0, 10.^(-4:0.005:log10(t_max))];

options = odeset('RelTol',1e-3,'AbsTol',1e-3); % Use this when the solution does not converge
m = 0;
sympref('HeavisideAtOrigin',0);
sol = pdepe(m,@(x,t,u,dudx)pdefun(x,t,u,dudx,p), @(x)icfun(x,p.c0), @(xL,uL,xR,uR,t)bcfun(xL,uL,xR,uR,t,p.cs),tau,t,options);

if Stepwise==0
    disp('Solution is completed')
    return  
end
disp('Solution for aerobic region is completed')

%% Numerical solution for the anaerobic region
[tau_matrix,t_matrix]=meshgrid(tau,t); % New order to agree with the NO3_sol dimensions
NO3_sol=sol(:,:,3); % These are the results for Nitrate, and these used as input to solve the anaerobic part for Nitrogen.
if tau_lim<tau_max

    tau = 10.^(log10(tau_lim):0.005:log10(tau_max)); % [d]
    t = [0, 10.^(-4:0.005:log10(t_max))];

    options = odeset('RelTol',1e-3,'AbsTol',1e-3); % Use this when the solution does not converge
    m = 0;
    sympref('HeavisideAtOrigin',0);
    sol_an = pdepe(m,@(x,t,u,dudx)pdefun_an(x,t,u,dudx,p,NO3_sol,t_matrix,tau_matrix), @(x)icfun_an(x,p.c0), @(xL,uL,xR,uR,t)bcfun_an(xL,uL,xR,uR,t,p.cs),tau,t,options);
    sol(:,:,4)=0; 
    sol((size(sol,1)-size(sol_an,1)+1):size(sol,1),(size(sol,2)-size(sol_an,2)+1):size(sol,2),4)=sol_an;
    disp('Solution for anaerobic region is completed')

end

%% Define functions
% Functions for the aerobic region
    function val = get_f(x,C,K,tau_lim,C_0_lim)
        val = [-K(1)*C(1); ...
            -K(3)*C(2)*(1-heaviside(C_0_lim-C(1))); ...
            -(K(4)*(1-heaviside(C_0_lim-C(1))) + K(5)*heaviside(C_0_lim-C(1)))*C(3) + K(3)*C(2)*(1-heaviside(C_0_lim-C(1)));
            K(5)*C(3)*heaviside(C_0_lim-C(1))];
    end

    function [c,f,s] = pdefun(x,t,u,dudx,p)
        c = p.R;
        f = p.DL/(p.V^2)*(dudx);
        s = -dudx + get_f(x,u,p.K,p.tau_lim,p.C_0_lim);
    end

    function u0 = icfun(x,u0_vals)
        u0 = u0_vals;
    end

    function [pL,qL,pR,qR] = bcfun(xL,uL,xR,uR,t,ub_vals)
        pL = uL-ub_vals;
        qL = [0; 0; 0; 0];
        pR = [0; 0; 0; 0];
        qR =[1; 1; 1; 1];
    end

% Function for the anaerobic region
    function val = get_f_an(x,t,C,K,tau_lim,C_0_lim,NO3_sol,t_matrix,tau_matrix)
        C_NO3=interp2(tau_matrix,t_matrix,NO3_sol,x,t);
        val = K(5)*C_NO3; 
    end

    function [c,f,s] = pdefun_an(x,t,u,dudx,p,NO3_sol,t_matrix,tau_matrix)
        c = p.R(4);
        f = p.DL/(p.V^2)*(dudx);
        s = -dudx + get_f_an(x,t,u,p.K,p.tau_lim,p.C_0_lim,NO3_sol,t_matrix,tau_matrix);
    end

    function u0 = icfun_an(x,u0_vals)
        u0 = u0_vals(4);
    end

    function [pL,qL,pR,qR] = bcfun_an(xL,uL,xR,uR,t,ub_vals)
        pL = uL-ub_vals(4);
        qL = 0;
        pR = 0;
        qR =1;
    end

end