function dxdt=communication(t,x,L)
    D = 33.5*60*60;                             %u-m2/hr
    W = 20;                                     %um
    R = 50;                                 	%um
    h = 2;                                      %um
    effective_time1 = 3.1416*R^2*(L+200)/(D*W);     %hr
    effective_time2 = 3.1416*R^2*(L+300)/(D*W); %hr
    effective_time3 = 3.1416*R^2*(L)/(D*W);     %hr
    
    compartment_volume  = 3.1416*R^2*h*10^-15;                      %L
    RNAP                = 28.57;                                    %uM 
    Gp                  = 10;                                       %uM
    Ribosome            = 0.035;                                    %uM
    protein_halflife    = 110/60;                                   %hr
    
    e_x     = 60*60*60;                                                %nt/hr
    e_L     = 60*16.5*60;                                              %aa/hr
    K_IX    = 60*60/42;                                                %McClure hr^-1
    K_X     = 0.05;                                                    %uM                                               
    K_IL    = 60*60/15;                                                %BNID:109525 1/hr
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
    
    S       = eye(6,6);

    A       = zeros(6,6);                                          %1/hr
    A(1,1)  = -(+1/effective_time1);
    A(2,2)  = -(+1/effective_time1);         
    A(3,3)  = -(+1/effective_time1);
    A(4,4)  = -(+1/effective_time1);
    A(5,5)  = -(+1/effective_time2);
    A(6,6)  = -(+1/effective_time2);
    
       %GFP
        W       = 0.9;
        K       = 0.3;                                                    
        n       = 1;
        W1      = 100;
        f       = (x(6)/effective_time3*t)^n/(K^n+(x(6)/effective_time3*t)^n);
        u1      = W/(1+W+W1*f);                                       
        
        T_X1    = 0.25*u1*KE_1*RNAP*Gp/(K_X*tau_1+(tau_1+1)*0.5*Gp);         %uM/hr                               
        T_L1    = KE_4*0.5*Ribosome*x(1)/(K_L*tau_4+(tau_4+1)*x(1));    %uM/hr
        
        %sigma38
        W       = 10;
        K       = 0.3;                                                    
        n       = 0.5;
        W1      = 50;
        f       = (x(6)/effective_time3*t)^n/(K^n+(x(6)/effective_time3*t)^n);
        u2      = W/(1+W+W1*f);    
        T_X2    = 0.25*u2*KE_1*RNAP*Gp/(K_X*tau_1+(tau_1+1)*0.5*Gp);    %uM/hr
        T_L2    = KE_4*0.5*Ribosome*x(3)/(K_L*tau_4+(tau_4+1)*x(3));    %uM/hr
        
        %CI 
        W       = 0.05;
        K       = 0.3;                                                    
        n       = 9;
        W1      = 500;
        f       = (x(4)/effective_time3*t)^n/(K^n+(x(4)/effective_time3*t)^n);
        u3      = (W+W1*f)/(1+W+W1*f);                                       
        
        T_X3    = u3*KE_1*RNAP*Gp/(K_X*tau_1+(tau_1+1)*Gp);    %uM/hr                               
        T_L3    = KE_4*Ribosome*x(5)/(K_L*tau_4+(tau_4+1)*x(5));    %uM/hr
        
        r = [ T_X1;
              T_L1;
              T_X2;
              T_L2;
              T_X3;
              T_L3];
       
       dxdt = (A*x+S*r);
end