function dxdt=GM1(t,x)
    D1 = 15*60*60;                             %u-m2/hr
    W  = 20;                                     %um
    R  = 50;                                     %um
    h  = 2;                                      %um
    L  = 200;   
    effective_time  = 3.1416*R^2*L/(D1*W);        %hr
    effective_time1 = 3.1416*R^2*100/(D1*W);
    
   k1   = 0.5;
   km_1 = 10;
   E1   = 0.1;
   k2   = 5;
   km_2 = 1;
   E2   = 2.5;
   
   dxdt = [-k1*E1*x(1)/(km_1+x(1))-x(1)*(1/effective_time+1/effective_time1);
            k1*E1*x(1)/(km_1+x(1))-x(2)*(1/effective_time+1/effective_time1);
            k2*E2*x(2)/effective_time*t/(km_2+x(2)/effective_time*t)-x(3)*2/effective_time1];
  end