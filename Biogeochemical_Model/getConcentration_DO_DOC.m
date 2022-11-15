function [F_do,F_doc] = getConcentration_DO_DOC(t, t_max,k_do, c_do_i, c_doc_i, poc_0, k_d, rho_b, alpha, poro)

    s = k_do*t;
    s_max = k_do*t_max;
    
    M_O = 32*1e-3; % [kg/mol]
    M_C = 30*1e-3; % [kg/mol]
    N_doc_i = c_doc_i/M_C; % [mol/m3]
    N_do_i = c_do_i/M_O;  % [mol/m3]

    Gamma = (rho_b*k_d*alpha)/(k_do);
    Pi = poc_0/(k_d*M_C*N_do_i);
    F_doc_i = N_doc_i/N_do_i;

    F_do =  exp(-s);
    
    
    if(Gamma == 1)
        
        % For Gamma == 1
        fun = @(x)exp(-x)*(Pi*(exp(x)-1) - x + F_doc_i); 
        s1 = fzero(fun,0);
        if(isnan(s1)) s1 = inf; end
        s2 = (s1 - 1) + exp(-s1)/(Gamma*Pi);
        s3 = s2 + 1;
        F_do_2 = -Gamma*Pi*(s2 - s1)+exp(-s1);

        inII = s> s1 & s<= s2;
        inIII = s> s2;
        
        F_do(inII) = -Pi*(s(inII) - s1) + exp(-s1);
        F_do(inIII) = F_do_2*exp(-(s(inIII)-s2));

        F_doc =  exp(-s)*(F_doc_i - Pi*(1-exp(-s)) - s);
        F_doc(inII) = 0;
        F_doc(inIII) = Pi*(1 - exp(-(s(inIII)-s2))) - F_do_2*(s(inIII) - s2)*exp(-(s(inIII) - s2));
        
        
    else
        
        % For Gamma != 1    
        fun = @(x)(((exp(-x) - exp(-Gamma*x))/(1-Gamma)) + Pi*(1-exp(-Gamma*x)) + F_doc_i*exp(-Gamma*x));    
        s1 = fzero(fun,0);
        if(isnan(s1)) s1 = inf; end
        s2 = (s1 - 1) + exp(-s1)/(Gamma*Pi);
        s3 = s2 + 1;
        F_do_2 = -Gamma*Pi*(s2 - s1)+exp(-s1);
        
        inII = s> s1 & s<= s2;
        inIII = s> s2;
        
        F_do(inII) = -Gamma*Pi*(s(inII)-s1)+exp(-s1);
        F_do(inIII) = F_do_2*exp(-(s(inIII)-s2));

        F_doc =  ((exp(-s) - exp(Gamma*s))/(1-Gamma)) + Pi* (1-exp(-Gamma*s)) + F_doc_i*exp(-Gamma*s);
        F_doc(inII) = 0;
        F_doc(inIII) = Gamma*Pi*(1 - exp(-(s(inIII)-s2))) - F_do_2*(s(inIII) - s2).*exp(-(s(inIII) - s2));
        
    end


end