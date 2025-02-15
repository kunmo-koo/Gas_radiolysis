clc
clear

global numberOfSpecies
global D
global generation


numberOfSpecies = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setting up initial variables (Dose, D, generation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order: [H H2 O O2 OH HO2 H2O2 H2O];
% PV=nRT n/V = P/RT
% P = 1atm, T=298K, R=0.082 initial concentration = 4.09E-2 Mol/L

molarMass = [1.008 2.016 15.999 31.999 17.008 33.007 34.015 18.015];
molecularVolume = [2.31 6.12 6.11 16.3 8.42 18.61 19.21 13.1];

for i = 1:numberOfSpecies
    Mab(i) = 2/(1/molarMass(i) + 1/molarMass(8));
    D(i) = 0.00143 * 1e14 * 298.^1.75 / (sqrt(Mab(i)) * ((molecularVolume(i).^1/3) + molecularVolume(8).^1/3).^2 * 0.013); %1atm, 0.5atm, 0.013 atm
end


disp(D);
% Initial concentrations (micromolar)%

% Order: [electrons H+ OH- H2O2 HO2- H OH O- HO2 O2- O2 H2 O3- O3 HO3 H2O];
initialConcentration = [0 0 0 0 0 0 0 4.09e4*0.013];


name = ["H" "H2" "O" "O2" "OH" "HO2" "H2O2" "H2O"];
name2 = ["H" "H_{2}" "O" "O_{2}" "OH" "HO_{2}" "H_{2}O_{2}" "H_{2}O"];

doseRate = 1e5 * 2.814 * 6.2e-11 / (3.141592 * 2.5e-21) ; % in Gy/s 10 torr water 200keV * 5* math.pow(10,-6)

gValues = [7.4/100 0.5/100 1.1/100 0 6.3/100 0 0 -7.4/100];  % molecules / eV

density = 8.04e-4* 0.013; %g/cm3

generation = 1e6 * density * doseRate * gValues / (6.022e23 * 1.6e-19); % generation per second in mM

disp(generation);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setting up initial concentration and geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry description = 100 nm 1D mesh, index = 1001
% probe size = 0.1 nm, probe interval = 1 nm, 10 point (10nm) scan

mesh = zeros(numberOfSpecies,1001);

for i = 1:numberOfSpecies
    mesh(i,:) = initialConcentration(i);
end


x=linspace (0,100,1001); % 0 to 100 nm, 1001 index

time = linspace(1e-9,1e-4,400000);


arg=mesh;
m=0;

bfun=@(x)icfun(x,arg);

sol = pdepe(m,@pdefun,bfun,@bcfun,x,time);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot timewise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
file1 = fopen('time.txt','w');
vout = VideoWriter('output.avi');
vout.Quality=95;
vout.FrameRate = 30;
open(vout);

for i = 1:400000
    if i<=399 || rem(i,100)==0
        u = sol(i,:,:);
        plt=squeeze(u).';

        semilogy(x,[plt(1,:);plt(2,:);plt(3,:);plt(4,:);plt(5,:);plt(6,:);plt(7,:);plt(8,:)],"LineWidth",1);

        xlabel('X (\mum)');
        ylabel('C (\muMol)');

        xlim([30 70]);
        ylim([1e-10 1e-1]);
        title(sprintf('Time: %.2d sec',time(i)));
        legend(name2(1),name2(2),name2(3),name2(4),name2(5),name2(6),name2(7),name2(8),'Location','southeastoutside');
        set(gcf,'position',[500 500 900 600]);
        pause(0.001);
        frame=getframe(gcf);
        writeVideo(vout,frame);
        drawnow;
    end
  

    if i==1 || i==37 || i==397 || i==3997 || i==39997 || i==399997
        fprintf(file1,'\n\n\nt = %.2d sec\n',time(i));
        fprintf(file1,'x (um)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8));
        for j=1:1000
            fprintf(file1,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',x(j),plt(1,j),plt(2,j),plt(3,j),plt(4,j),plt(5,j),plt(6,j),plt(7,j),plt(8,j));
        end
    end

end
close(vout);
fclose(file1);


figure(2);
cplot = squeeze(sol(:,456,:)).';
loglog(time,[cplot(1,:);cplot(2,:);cplot(3,:);cplot(4,:);cplot(5,:);cplot(6,:);cplot(7,:);cplot(8,:)]);
title('Center');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file2=fopen('center.txt','w');
fprintf(file2,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8));
for j=1:400000
    if rem(j,100)==0 || j<=300
        fprintf(file2,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),cplot(1,j),cplot(2,j),cplot(3,j),cplot(4,j),cplot(5,j),cplot(6,j),cplot(7,j),cplot(8,j));
    end
end
fclose(file2);


figure(3);
oplot = squeeze(sol(:,465,:)).';
loglog(time,[oplot(1,:);oplot(2,:);oplot(3,:);oplot(4,:);oplot(5,:);oplot(6,:);oplot(7,:);oplot(8,:)]);
title('Between Beam');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file3=fopen('undersample.txt','w');
fprintf(file3,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8));
for j=1:400000
    if rem(j,100)==0 || j<=300
        fprintf(file3,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),oplot(1,j),oplot(2,j),oplot(3,j),oplot(4,j),oplot(5,j),oplot(6,j),oplot(7,j),oplot(8,j));
    end
end
fclose(file3);


figure(4);
eplot = squeeze(sol(:,551,:)).';
loglog(time,[eplot(1,:);eplot(2,:);eplot(3,:);eplot(4,:);eplot(5,:);eplot(6,:);eplot(7,:);eplot(8,:)]);
title('Edge');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file4=fopen('edge.txt','w');
fprintf(file4,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8));
for j=1:400000
    if rem(j,100)==0 || j<=300
        fprintf(file4,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),eplot(1,j),eplot(2,j),eplot(3,j),eplot(4,j),eplot(5,j),eplot(6,j),eplot(7,j),eplot(8,j));
    end
end
fclose(file4);


figure(5);
dplot = squeeze(sol(:,751,:)).';
loglog(time,[dplot(1,:);dplot(2,:);dplot(3,:);dplot(4,:);dplot(5,:);dplot(6,:);dplot(7,:);dplot(8,:)]);
title('Distant');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file5=fopen('distant.txt','w');
fprintf(file5,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8));
for j=1:400000
    if rem(j,100)==0 || j<=300
        fprintf(file5,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),dplot(1,j),dplot(2,j),dplot(3,j),dplot(4,j),dplot(5,j),dplot(6,j),dplot(7,j),dplot(8,j));
    end
end
fclose(file5);



function [c,f,s] = pdefun(x,t,C,dCdx)

global D;
global generation;

c = [1/D(1);1/D(2);1/D(3);1/D(4);1/D(5);1/D(6);1/D(7);1/D(8)];
f = [dCdx(1);dCdx(2);dCdx(3);dCdx(4);dCdx(5);dCdx(6);dCdx(7);dCdx(8)];

%%%%%% Probe Generation if the t<=dwell time (continuous probes interlacing) %%%%%%

n = floor(t/3e-6);
n = rem(n,10);

if x>=(45.45+n) && x<=(45.55+n)
    s = [generation(1)/D(1);generation(2)/D(2);generation(3)/D(3);generation(4)/D(4);generation(5)/D(5);generation(6)/D(6);generation(7)/D(7);generation(8)/D(8)];

else
    s = [0;0;0;0;0;0;0;0];

end

%%%%%% Reaction between the spiecies %%%%%%

molarToMicromolar = 10^6;

k = zeros(18,1);
k(1) = 3.36*10^10;
k(2) = 1.58*10^12;
k(3) = 1.44*10^12;
k(4) = 2.07*10^10;
k(5) = 1.64*10^10;
k(6) = 3.68*10^9;
k(7) = 3.86*10^10;
k(8) = 2.40*10^9;
k(9) = 1.80*10^9;
k(10) = 1.58*10^10; 
k(11) = 2.43*10^7;
k(12) = 9.00*10^8;
k(13) = 1.30*10^4;
k(14) = 8.39*10^8;
k(15) = 3.87*10^8;
k(16) = 5.11*10^3;
k(17) = 1.58*10^-20;
k(18) = 1.17*10^-4;

%%convert from molar to micromolar
k = k./(molarToMicromolar);

% Order: [H H2 O O2 OH HO2 H2O2 H2O];

%Reactant concentrations (micromolar)%
H_mono   = C(1); 
H2       = C(2);
O_mono   = C(3);
O2       = C(4);
OH       = C(5);
HO2      = C(6);
H2O2     = C(7);
H2O      = C(8);


%Reaction Set%
% Note: H2O is divided out as a reactant
r = zeros(18,1);
r(1)  = k(1)*H_mono*H_mono;
r(2)  = k(2)*H_mono*OH;
r(3)  = k(3)*OH*OH;
r(4)  = k(4)*O_mono*OH;
r(5)  = k(5)*H_mono*O2;
r(6)  = k(6)*H_mono*HO2;
r(7)  = k(7)*H_mono*HO2;
r(8)  = k(8)*H_mono*HO2;
r(9)  = k(9)*HO2*HO2;
r(10) = k(10)*H_mono*O_mono;
r(11) = k(11)*H_mono*H2O2;
r(12) = k(12)*OH*H2O2;
r(13) = k(13)*OH*H2;
r(14) = k(14)*OH*OH;
r(15) = k(15)*O_mono*O_mono;
r(16) = k(16)*O_mono*H2;
r(17) = k(17)*H2O2;
r(18) = k(18)*H_mono*H2O;


products = zeros(16,1);


%%H%%
products(1) = -2*r(1) - r(2) - r(5) - r(6) -r(7) -r(8) -r(10) -r(11) -r(18) + r(4) +r(13) +r(16);

%%H2%%
products(2) = -r(13) -r(16) + r(1) +r(6) +r(18);

%%O%%
products(3) = -r(4) -r(10) -2*r(15) -r(16) +r(8) +r(14);

%%O2%%
products(4) = -r(5) +r(4) +r(6) +r(9) +r(15);

%%OH%%
products(5) = -r(2) -2*r(3) -r(4) -r(12) -r(13) -2*r(14) +2*r(7) +r(10) +r(11) +r(16) +2*r(17) +r(18);

%%HO2%%
products(6) = -r(6) -r(7) -r(8) -2*r(9) +r(5) +r(12);

%%H2O2%%
products(7) = -r(11) - r(12) - r(17) +r(3) +r(9);

%%H2O%%
products(8) = -r(18) +r(2) +r(8) +r(11) +r(12) +r(13) +r(14);


s(1) = s(1) + products(1)/D(1);
s(2) = s(2) + products(2)/D(2);
s(3) = s(3) + products(3)/D(3);
s(4) = s(4) + products(4)/D(4);
s(5) = s(5) + products(5)/D(5);
s(6) = s(6) + products(6)/D(6);
s(7) = s(7) + products(7)/D(7);
s(8) = s(8) + products(8)/D(8);


end

function C0 = icfun(x,arg)

C0=[arg(1,round(x*10) + 1);arg(2,round(x*10) + 1);arg(3,round(x*10) + 1);arg(4,round(x*10) + 1);arg(5,round(x*10) + 1);arg(6,round(x*10) + 1);arg(7,round(x*10) + 1);arg(8,round(x*10) + 1)];

end

function [pL,qL,pR,qR] = bcfun(xL,CL,xR,CR,t)

pR = CR;
pR(8) = CR(8)-4.09e4*0.013;
qR = [0;0;0;0;0;0;0;0];

pL = CL;
pL(8) = CL(8)-4.09e4*0.013;
qL = [0;0;0;0;0;0;0;0];

end

