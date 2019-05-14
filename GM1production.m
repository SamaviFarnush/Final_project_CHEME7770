%Miya Bidon & Samavi Farnush Bint E Naser
%CHEME 7770
%Final Project
%8 May 2019
%--------------------------------------------------------------------------%

clc
clear all

x0 = [1;                                    %GM3
      0.0;                                  %GM2
      0.0;                                  %GM1
      ];

end_sim = 6;
time    = [0:1/120:end_sim];                  
[t,X]   = ode45('GM1',time,x0);
X       = X*10^3; 

figure(1)
plot(t,X(:,(2:3)))
legend('GM2','GM1')
xlabel("time (hr)")
ylabel("concentration, nM")
title("GM1 production")