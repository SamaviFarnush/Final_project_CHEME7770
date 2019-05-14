function dxdt=positive_feedback(t,x,L)
    D = 33.5*60*60;                         %u-m2/hr
    W = 20;                                 %um
    R = 50;                                 %um
    h = 2;                                  %um
    effective_time = 3.1416*R^2*L/(D*W);    %hr
    
    compartment_volume  = 3.1416*R^2*h*10^-15;                      %L
    RNAP                = 28.57/2;                                  %uM 
    Gp                  = 10/2;                                     %uM
    Ribosome            = 0.035/2;                                  %uM
    protein_halflife    = 110/60;                                   %hr
    
    e_x     = 60*60*60;                                             %nt/hr
    e_L     = 16.5*60*60;                                           %aa/hr
    K_IX    = 3600/42;                                              %McClure hr^-1
    K_X     = 0.05;                                                 %uM                                               
    K_IL    = 3600/15;                                              %BNID:109525 1/min
    K_L     = 454.64;
    
    %read length
    LX_1    = 300;                                                 %nt
    LT_1    = 300/3;                                               %aa
    kd_m    = 60/10;                                               %1/hr, Karzbrun
    kd_p    = log(2)/(protein_halflife);                           %1/hr

    KE_1    = e_x/LX_1;                                            %1/hr
    tau_1   = KE_1/K_IX; 
    KE_4    = e_L/LT_1;                                            %1/hr
    tau_4   = KE_4/K_IL;
    
    S       = eye(4,4);

    A       = zeros(4,4);                                          %1/hr
    A(1,1)  = -(kd_m+1/effective_time);
    A(2,2)  = -(kd_p+1/effective_time);         
    A(3,3)  = -(kd_m+D*W/3.1416/R^2/L);
    A(4,4)  = -(kd_p+D*W/3.1416/R^2/L);
    
        %araC
        u1      = 1;
        T_X1    = u1*KE_1*RNAP*Gp/(K_X*tau_1+(tau_1+1)*Gp);         %uM/hr
        T_L1    = KE_4*Ribosome*x(1)/(K_L*tau_4+(tau_4+1)*x(1));    %uM/hr
        
        %GFP
        W       = 0;
        K       = 0.3;                                                    
        n       = 0.5;
        W1      = 50;
        f       = x(2)^n/(K^n+x(2)^n);
        u3      = (W+W1*f)/(1+W+W1*f);                                       
        
        T_X2    = u3*KE_1*RNAP*Gp/(K_X*tau_1+(tau_1+1)*Gp);         %uM/hr                               
        T_L2    = KE_4*Ribosome*x(3)/(K_L*tau_4+(tau_4+1)*x(3));    %uM/hr

    if (t>=0.005*L+0.75)
        r = [ T_X1;
              T_L1;
              T_X2;
              T_L2];
    else
        r = [T_X1;
             T_L1;
             0;
             0];
    end
    
    dxdt=(A*x+S*r);
end