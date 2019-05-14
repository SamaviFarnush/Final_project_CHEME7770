%Miya Bidon & Samavi Farnush Bint E Naser
%CHEME 7770
%Final Project
%8 May 2019
%--------------------------------------------------------------------------%

clc
clear all

%unregulated
x0 = [0.0;                                 %mRNA1
      0.0;                                 %protein1                
     ];

end_sim = 4;

time    = [0:1/60:end_sim];
L       = [50:50:300];

for i = 1:length(L)
    Length  = L(i);
    [t,X]   = ode45(@(t,x) model(t,x,Length),time,x0);
    X       = X(:,2);
    F       = 12.3*X;
    C(:,i)  = F;
end 

figure(1)
plot(t,C);
legend('50 um','100 um','150 um','200 um','250 um','300 um')
xlabel("time (hr)")
ylabel("F (AU)")
xlim([0 4])
title("unregulated")
hold on


%positive feedback
x0      = [0.0;                                %mRNA1
           0.0;                                %araC                
           0.0;                                %mRNA1
           0.0;];                              %GFP

end_sim     = 4;

time        = [0:1/60:end_sim];
L           = [50:50:300];

for i=1:length(L)
    Length = L(i);
    [t,Y]   = ode45(@(t,x) positive_feedback(t,x,Length),time,x0);
    Y       = Y(:,4);
    F2      = 12.3*Y;
    C2(:,i) = F2;
end 


figure(2)
plot(t,C2);
legend('50 um','100 um','150 um','200 um','250 um','300 um')
xlabel("time (hr)")
ylabel("F (AU)")
ylim([0 8])
title("positive feedback")
hold on


%Negative feedback
x0      = [0.0;                                %mRNA1
           0.0;                                %cro                
           0.0;                                %mRNA1
           0.0;];                              %GFP

end_sim     = 4;

time        = [0:1/60:end_sim];
L           = [50:50:300];

for i=1:length(L)
    Length  = L(i);
    [t,Z]   = ode45(@(t,x) negative_feedback(t,x,Length),time,x0);
    Z       = Z(:,2);
    F3      = 12.3*Z;
    C3(:,i) = F3;
end 

figure(3)
plot(t,C3);
legend('50 um','100 um','150 um','200 um','250 um','300 um')
xlabel("time (hr)")
ylabel("F (AU)")
xlim([0 4])
ylim([0 2])
title("negative feedback")
hold on


%Activator_repressor 1
x0      = [0.0;                                %mRNA1
           0.0;                                %sigma28                
           0.0;                                %mRNA2
           0.0;                                %GFP
           0.0;                                %mRNA3
           0.0;];                              %CI

end_sim     = 15;

time        = [0:1/60:end_sim];
L           = [50:50:300];

for i=1:length(L)
    Length = L(i);
    [t,V]   = ode45(@(t,x) activator_repressor1(t,x,Length),time,x0);
    V       = V(:,4);
    F4      = 12.3*V;
    C4(:,i) = F4;
end 

figure(4)
plot(t,C4(:,(2:6)));
legend('100 um','150 um','200 um','250 um','300 um')
xlabel("time (hr)")
ylabel("F (AU)")
title("activator-repressor network-1")
hold on


%Activator_repressor 2
x0      = [0.0;                                %mRNA1
           0.0;                                %GFP                
           0.0;                                %mRNA2
           0.0;                                %sigma38
           0.0;                                %mRNA3
           0.0;];                              %CI

end_sim     = 6;

time        = [0:1/60:end_sim];
L           = [50:50:300];

for i=1:length(L)
    Length = L(i);
    [t,U]   = ode45(@(t,x) activator_repressor2(t,x,Length),time,x0);
    U       = U(:,2);
    F5      = 12.3*U;
    C5(:,i) = F5;
end 

figure(5)
plot(t,C5(:,(2:6)));
legend('100 um','150 um','200 um','250 um','300 um')
xlabel("time (hr)")
ylabel("F (AU)")
title("activator-repressor network-2")
hold on


%Communication
x0      = [0.0;                                %mRNA1
           0.0;                                %GFP                
           0.0;                                %mRNA2
           0.0;                                %sigma38
           0.0;                                %mRNA3
           0.0;];                              %CI

end_sim     = 8;

time        = [0:1/60:end_sim];
L           = [50:50:300];

for i=1:length(L)
    Length = L(i);
    [t,S]   = ode45(@(t,x) communication(t,x,Length),time,x0);
    S       = S(:,2);
    F6      = 12.3*S;
    C6(:,i) = F6;
end 

figure(6)
plot(t,C6);
legend('50 um','100 um','150 um','200 um','250 um','300 um')
xlabel("time (hr)")
ylabel("F (AU)")
ylim([0 10])
title("communication between DNA compartments")
hold on