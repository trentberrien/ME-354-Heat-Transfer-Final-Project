clc;
clear;

Ti = 15; %assumed starting temp of nodes
Tinf = 15; %deg C, temp of cooling water
qdot = 100000; %W/m^2
d = .0025; %m

W = .0025; %Max 17mm
H = .0125; %Max 19mm
a = .010; %9mm to 12mm
b = .01; %10mm to 14mm
L = .09; %90mm
k = 190; %190 W/m-K
Vf = 7e-5;
nu = 1.139e-6;
Pr = 8.09;

Dh = 4*W*H/(2*H+W);
u = Vf/(W*H);
Re = u*Dh/nu;
f = (.79*log(Re)-1.64)^-2;
Nu = .125*f*Re*(Pr^(1/3));
h = Nu*k/Dh;
        
Rtot = .02*(k+h*(.02-H))/((.02*k*h*W*L+(k+h*(.02-H)))*(2*k*L*(a-(W/2))));
As = 2*a*L;
        
U=1/(Rtot*As);

%tstep calcualtion
alpha = 8.3e-5; %google diffusivity of alumminum
deltmax = d^2/(4*alpha*(1+h*d/k));
tstep = .01*deltmax;
tau = (alpha*tstep)/(d^2);
t = 0;
i = 1;

%init node temps to Ti
T1 = Ti;
T2 = Ti;
T3 = Ti;
T4 = Ti;
T5 = Ti;
T6 = Ti;
T7 = Ti;
T8 = Ti;
T9 = Ti;
T10 = Ti;
T11 = Ti;
T12 = Ti;
T13 = Ti;
T14 = Ti;
T15 = Ti;
T16 = Ti;
T17 = Ti;
T18 = Ti;
T19 = Ti;
T20 = Ti;
T21 = Ti;
T22 = Ti;
T23 = Ti;
T24 = Ti;
T25 = Ti;
T26 = Ti;
T27 = Ti;
T28 = Ti;
T29 = Ti;
T30 = Ti;
T31 = Ti;
T32 = Ti;

plotData = [0 T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13 T14 T15 T16 T17 T18 T19 T20 T21 T22 T23 T24 T25 T26 T27 T28 T29 T30 T31 T32];

%Formulas
%type 4: TX*(1-4*tau)+tau*(2*qdot*d/k+TX+2*TX+TX);
%type 5: TX*(1-4*tau)+tau*(TX+TX+TX+TX);
%type 6: TX*(1-4*tau-4*tau*h*d/(3*k))+(tau/3)*(4*h*d*Tinf/k+2*TX+4*TX+2*TX+4*TX);
%type 7: TX*(1-4*tau-2*tau*h*d/k)+tau*(TX+TX+2*TX+2*h*d*Tinf/k);
%type 8: TX + 2*tau*(-3*TX+.5*TX+.5*TX+TX);

while (t <= .5)

    T1new = T1 + 2*tau*(-2*T1+.5*T2+.5*T2+T5);
    T2new = T2 + 2*tau*(-2*T2+.5*T1+.5*T3+T6);
    T3new = T3 + 2*tau*(-2*T3+.5*T2+.5*T4+T7);
    T4new = T4 + 2*tau*(-2*T4+.5*T3+.5*T3+T8);
    T5new = T5*(1-4*tau)+tau*(T1+T6+T6+T9);
    T6new = T6*(1-4*tau)+tau*(T2+T5+T7+T10);
    T7new = T7*(1-4*tau-4*tau*h*d/(3*k))+(tau/3)*(4*h*d*Tinf/k+2*T8+4*T3+2*T11+4*T6);
    T8new = T8*(1-4*tau-2*tau*h*d/k)+tau*(T7+T7+2*T4+2*h*d*Tinf/k);
    T9new = T9*(1-4*tau)+tau*(T5+T10+T10+T13);
    T10new = T10*(1-4*tau)+tau*(T6+T9+T11+T14);
    T11new = T11*(1-4*tau-2*tau*h*d/k)+tau*(T7+T15+2*T10+2*h*d*Tinf/k);
    %T12new = 15
    T13new = T13*(1-4*tau)+tau*(T9+T14+T14+T17);
    T14new = T14*(1-4*tau)+tau*(T10+T13+T15+T18);
    T15new = T15*(1-4*tau-2*tau*h*d/k)+tau*(T11+T19+2*T14+2*h*d*Tinf/k);
    %T16new = 15
    T17new = T17*(1-4*tau)+tau*(T13+T18+T18+T21);
    T18new = T18*(1-4*tau)+tau*(T14+T17+T19+T22);
    T19new = T19*(1-4*tau-2*tau*h*d/k)+tau*(T15+T23+2*T18+2*h*d*Tinf/k);
    %T20new = 15
    T21new = T21*(1-4*tau)+tau*(T17+T22+T22+T25);
    T22new = T22*(1-4*tau)+tau*(T18+T21+T23+T25);
    T23new = T23*(1-4*tau-2*tau*h*d/k)+tau*(T19+T23+2*T22+2*h*d*Tinf/k);
    %T24new = 15
    T25new = T25*(1-4*tau)+tau*(T21+T26+T26+T29);
    T26new = T26*(1-4*tau)+tau*(T22+T25+T27+T30);
    T27new = T27*(1-4*tau-4*tau*h*d/(3*k))+(tau/3)*(4*h*d*Tinf/k+2*T23+4*T26+2*T28+4*T31);
    T28new = T28*(1-4*tau-2*tau*h*d/k)+tau*(T27+T27+2*T32+2*h*d*Tinf/k);
    T29new = T29*(1-4*tau)+tau*(2*qdot*d/k+T30+2*T25+T30);
    T30new = T30*(1-4*tau)+tau*(2*qdot*d/k+T29+2*T26+T31);
    T31new = T31*(1-4*tau)+tau*(2*qdot*d/k+T30+2*T27+T32);
    T32new = T32*(1-4*tau)+tau*(2*qdot*d/k+T31+2*T28+T31);


    i = i +1;
    t = t+tstep;

    plotData(i,1) = t;
    plotData(i,2) = T1new; %
    plotData(i,3) = T2new; %
    plotData(i,4) = T3new; %
    plotData(i,5) = T4new; %
    plotData(i,6) = T5new;
    plotData(i,7) = T6new;
    plotData(i,8) = T7new;
    plotData(i,9) = T8new;
    plotData(i,10) = T9new;
    plotData(i,11) = T10new;
    plotData(i,12) = T11new;
    plotData(i,13) = T12; %
    plotData(i,14) = T13new;
    plotData(i,15) = T14new;
    plotData(i,16) = T15new;
    plotData(i,17) = T16; %
    plotData(i,18) = T17new;
    plotData(i,19) = T18new;
    plotData(i,20) = T19new;
    plotData(i,21) = T20; %
    plotData(i,22) = T21new;
    plotData(i,23) = T22new;
    plotData(i,24) = T23new;
    plotData(i,25) = T24; %
    plotData(i,26) = T25new;
    plotData(i,27) = T26new;
    plotData(i,28) = T27new;
    plotData(i,29) = T28new;
    plotData(i,30) = T29new;
    plotData(i,31) = T30new;
    plotData(i,32) = T31new;
    plotData(i,33) = T32new;

    T1 = T1new;
    T2 = T2new;
    T3 = T3new;
    T4 = T4new;
    T5 = T5new;
    T6 = T6new;
    T7 = T7new;
    T8 = T8new;
    T9 = T9new;
    T10 = T10new;
    T11 = T11new;
    %T12 = T12new;
    T13 = T13new;
    T14 = T14new;
    T15 = T15new;
    %T16 = T16new;
    T17 = T17new;
    T18 = T18new;
    T19 = T19new;
    %T20 = T20new;
    T21 = T21new;
    T22 = T22new;
    T23 = T23new;
    %T24 = T24new;
    T25 = T25new;
    T26 = T26new;
    T27 = T27new;
    T28 = T28new;
    T29 = T29new;
    T30 = T30new;
    T31 = T31new;
    T32 = T32new;
end

if max([T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13 T14 T15 T16 T17 T18 T19 T20 T21 T22 T23 T24 T25 T26 T27 T28 T29 T30 T31 T32]) > 40
    fprintf("TOO HOT");
else
    fprintf("Did not exceed 40C. Max temp was: " + max([T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13 T14 T15 T16 T17 T18 T19 T20 T21 T22 T23 T24 T25 T26 T27 T28 T29 T30 T31 T32]));
end
for i=2:1:32
    plot(plotData(:,1),plotData(:,i));
    hold on
end
xlabel('Time (s)');
ylabel('Temp (deg C)');
legend({'T1','T2','T3','T4','T5','T6','T7','T8','T10','T11','T12','T13','T14','T15','T16','T17','T18','T19','T20','T21','T22','T23','T24','T25','T26','T27','T28','T29','T30','T31','T32'}, 'Location', 'southeast');