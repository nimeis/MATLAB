clear
Flc = 230;
C1pC3 = 50;
pickBeta = false;
acBeta = 160;
C2 = 10*10^-6;
Va = 40;
Veb = 0.7;
Av = 52;
RE1 = 100;
Rs = 100;
Rin = 6000;
vt=0.026;
Vcc = 15;
desiredSwing = 13.25;
VecSat = 0.5;
VecSatThreshold = 0.25;
numberOfPointsPerGraph = 101;
numberOfGraphsPerGraph = 10;
numberGraohToUse4Analysis = 4;
interpolationPointGap = 0.1;

VecLimit = VecSat + VecSatThreshold + desiredSwing;
Vecstart = VecLimit-desiredSwing;
VecQ = ((VecLimit - Vecstart)/2) + Vecstart;

VecMax = Vcc;


file = "transitorCurveData.csv";
data = readmatrix(file);
Ib = data(1,:);
Ic = data(2,:);
Vec = data(3,:);
shape = [numberOfPointsPerGraph,numberOfGraphsPerGraph];
Ib = transpose(reshape(Ib,shape));
Ic = transpose(reshape(Ic,shape));
Vec = transpose(reshape(Vec,shape));
Ib = Ib./1000000;
Ic = Ic./1000;
plot(Vec(1,:),Ic(1,:),Vec(2,:),Ic(2,:),Vec(3,:),Ic(3,:),Vec(4,:),Ic(4,:),Vec(5,:),Ic(5,:),Vec(6,:),Ic(6,:),Vec(7,:),Ic(7,:),Vec(8,:),Ic(8,:),Vec(9,:),Ic(9,:),Vec(10,:),Ic(10,:))
Vec_cutoffmid = (Vec(numberGraohToUse4Analysis,:) < (VecQ+interpolationPointGap)) & (Vec(numberGraohToUse4Analysis,:) > (VecQ-interpolationPointGap));
Vec_cutofflow = (Vec(numberGraohToUse4Analysis-1,:) < (VecQ+interpolationPointGap)) & (Vec(numberGraohToUse4Analysis-1,:) > (VecQ-interpolationPointGap));
Vec_cutoffhigh = (Vec(numberGraohToUse4Analysis+1,:) < (VecQ+interpolationPointGap)) & (Vec(numberGraohToUse4Analysis+1,:) > (VecQ-interpolationPointGap));
Vecmid = Vec(numberGraohToUse4Analysis,Vec_cutoffmid);
Ibmid = Ib(numberGraohToUse4Analysis,Vec_cutoffmid);
Icmid = Ic(numberGraohToUse4Analysis,Vec_cutoffmid);

Veclow = Vec(numberGraohToUse4Analysis-1,Vec_cutofflow);
Iblow = Ib(numberGraohToUse4Analysis-1,Vec_cutofflow);
Iclow = Ic(numberGraohToUse4Analysis-1,Vec_cutofflow);

Vechigh = Vec(numberGraohToUse4Analysis+1,Vec_cutoffhigh);
Ibhigh = Ib(numberGraohToUse4Analysis+1,Vec_cutoffhigh);
Ichigh = Ic(numberGraohToUse4Analysis+1,Vec_cutoffhigh);

%plot(Vec,Ic)
p1mid = [Vecmid(1),Icmid(1)];
p2mid = [Vecmid(end),Icmid(end)];
p1low = [Veclow(1),Iclow(1)];
p2low = [Veclow(end),Iclow(end)];
p1high = [Vechigh(1),Ichigh(1)];
p2high = [Vechigh(end),Ichigh(end)];
dymid = (p2mid(2)-p1mid(2));
dxmid = (p2mid(1)-p1mid(1));
dylow = (p2low(2)-p1low(2));
dxlow = (p2low(1)-p1low(1));
dyhigh = (p2high(2)-p1high(2));
dxhigh = (p2high(1)-p1high(1));
Mmid = dymid/dxmid;
Mlow = dylow/dxlow;
Mhigh = dyhigh/dxhigh;
p1s = [p1low;p1mid;p1high];
Ms = [Mlow,Mmid,Mhigh];

Iblow = mean(Iblow);
Ibhigh = mean(Ibhigh);
Ibmid = mean(Ibmid);

Ibs = [Iblow,Ibmid,Ibhigh];

VecQ;
Icq = calcIc(2,VecQ,Ms,p1s);
IcLowq = calcIc(1,VecQ,Ms,p1s);
IcHighq = calcIc(3,VecQ,Ms,p1s);

Ibq = Ibs(2);

dcB = (Icq)/(Ibs(2))*-1
dIbu = Ibs(3)-Ibs(2);
dIcu = IcHighq-Icq;
dIbd = Ibs(2)-Ibs(1);
dIcd = Icq-IcLowq;

acBu = (dIcu/dIbu)*-1;
acBd = (dIcd/dIbd)*-1;


if (pickBeta)
    acB = acBeta;
else
    acB = mean([acBd,acBd]);
end

gm = Icq/(vt);

Rc = 52*((1/gm)+RE1);

Rpi = acB/gm;

Rstuff = Rpi+(RE1*acB);

Rth = ((Rin-Rs)*Rstuff)/(Rstuff-((Rin-Rs)));

alpha = (acB+1)/acB;

Vcq = Icq*Rc;
Veq = Vcq+VecQ;

Ie = alpha*Icq;

RE2 = ((Vcc-Veq)-(Ie*RE1))/Ie;

Vth = Vcc-(Ie*(RE1+RE2))-Veb-((-1*Ibq)*Rth);

R1 = (Vcc*(Rth))/Vth;

R2 = 1/((1/Rth)-(1/R1));






Rs
R1
R2
RE1
RE2
Rc
Icq
VecQ



reqC1 = Rin+50
reqC2 = 1/((1/((1/gm)+((1/((1/Rth)+(1/Rs)))/(acB+1))+RE1))+(1/RE2));
re1rpi = 1/((1/RE1)+(1/Rpi));
Ro = Va/Icq;
Ro = Ro*(1+(gm*(re1rpi)));
reqC3 = 1/((1/Rc)+(1/Ro));

f2 = 1/(2*pi*C2*reqC2);
f1 = (Flc-f2)*(C1pC3/100);
f3 = (Flc-f2)*((100-C1pC3)/100);
C1 = (1/(f1))/(2*pi*reqC1);
C3 = (1/(f3))/(2*pi*reqC3);





Rs
R1
R2
RE1
RE2
Rc
C1
C2
C3
VecQ
































function Ic4Vce = calcIc(lowmidhigh,Vec,Ms,Ps)
x1 = Ps(lowmidhigh,1);
y1 = Ps(lowmidhigh,2);
M = Ms(lowmidhigh);
if lowmidhigh == 1
        Ic4Vce = (M*(Vec-x1))+y1;
elseif lowmidhigh == 2
        Ic4Vce = (M*(Vec-x1))+y1;
elseif lowmidhigh == 3
        Ic4Vce = (M*(Vec-x1))+y1;
else
    Ic4Vce = -1;
end
end
