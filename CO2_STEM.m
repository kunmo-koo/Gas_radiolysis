clc
clear

global numberOfSpecies
global D
global generation


numberOfSpecies = 18;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setting up initial variables (Dose, D, generation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order: [e-, C, C+, O, O+, O2, O2+, CO, CO+, CO2, CO2+,CO3+, CO4+,CO4-, C2O2+, C2O3+, C2O4+, C3O4+]
% PV=nRT n/V = P/RT
% P = 1atm, T=298K, R=0.082 initial concentration = 4.09E-2 Mol/L

molarMass = [5.49e-4 12.011 12.011 15.999 15.999 31.998 31.998 28.010 28.010 44.009 44.009 60.008 76.007 76.007 56.020 72.019 88.018 100.029];
molecularVolume = [0 15.9 15.9 6.11 6.11 16.3 16.3 22.01 22.01 32.2 32.2 38.31 48.5 48.5 48.1 54.21 64.4 80.3];

for i = 1:numberOfSpecies
    Mab(i) = 2/(1/molarMass(i) + 1/molarMass(10));
    D(i) = 0.00143 * 1e14 * 298.^1.75 / (sqrt(Mab(i)) * ((molecularVolume(i).^1/3) + molecularVolume(10).^1/3).^2 * 0.013); %1atm, 0.0133322 atm
end
disp(D);

initialConcentration = [0 0 0 0 0 0 0 0 0 4.09e4*0.013 0 0 0 0 0 0 0 0];  % uMol

name = ["e-" "C" "C+" "O" "O+" "O2" "O2+" "CO" "CO+" "CO2" "CO2+" "CO3+" "CO4+" "CO4-" "C2O2+" "C2O3+" "C2O4+" "C3O4+"];
name2 = ["e^{-}" "C" "C^{+}" "O" "O^{+}" "O_{2}" "O_{2}^{+}" "CO" "CO^{+}" "CO_{2}" "CO_{2}^{+}" "CO_{3}^{+}" "CO_{4}^{+}" "CO_{4}^{-}" "C_{2}O_{2}^{+}" "C_{2}O_{3}^{+}" "C_{2}O_{4}^{+}" "C_{3}O_{4}^{+}"];

doseRate = 2.480e5 * 6.2e-11 / (3.141592 * 2.5e-21) ; % in CO2
disp(doseRate);

gValues = [2.96/100 0 0.07/100 0.51/100 0.21/100 0 0 0.21/100 0.51/100 -3.03/100 2.24/100 0 0 0 0 0 0 0];  % molecules / eV

density=1.96e-3 * 0.013; %g/cm3

generation = 1e6 * density * doseRate * gValues / (6.022e23 *1.6e-19); % generation per second in mM

%disp(generation);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setting up initial concentration and geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry description = 100 nm 1D mesh, index = 1001
% probe size = 0.1 nm, probe interval = 1 nm, 10 point (10nm) scan

mesh = zeros(numberOfSpecies,1001);

for i = 1:numberOfSpecies
    mesh(i,:) = initialConcentration(i);
end


time = linspace(1e-9,1e-4,200000);
x=linspace (0,100,1001);

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

for i = 1:200000
    if i<=399 || rem(i,100)==0
        u = sol(i,:,:);
        plt=squeeze(u).';
        semilogy(x,[plt(10,:);plt(11,:);plt(8,:)],'-',"LineWidth",1);
        hold on;
        semilogy(x,[plt(4,:)],'--',"LineWidth",1);
        hold on;
        semilogy(x,[plt(5,:);plt(3,:);plt(1,:);plt(2,:);plt(7,:)],'-',"LineWidth",1);
        hold off;

        line_color = ["#37AD6B";"#8E8E00";"#7D4E4E";"#00CECE";"#B177DE";"#1A6FDF";"#525252";"#CC9900";"#FF6300"];

        colororder(line_color);

        xlabel('X (nm)');
        ylabel('C (uMol)')
        xlim([30 70]);
        ylim([1e-10 1e5]);
        title(sprintf('Time: %.2d sec',time(i)));
        legend(name2(10),name2(11),name2(8),name2(4),name2(5),name2(3),name2(1),name2(2),name2(7),'Location','southeastoutside');
        set(gcf,'position',[500 500 900 600]);
        pause(0.001);
        frame=getframe(gcf);
        writeVideo(vout,frame);
        drawnow;
    end

    if i==1 || i==17 || i==197 || i==1997 || i==19997 || i==199997
        u = sol(i,:,:);
        plt=squeeze(u).';
        fprintf(file1,'\n\n\nt = %.2d sec\n',time(i));
        fprintf(file1,'x (nm)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16),name(17),name(18));
        for j=1:1000
            fprintf(file1,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',x(j),plt(1,j),plt(2,j),plt(3,j),plt(4,j),plt(5,j),plt(6,j),plt(7,j),plt(8,j),plt(9,j),plt(10,j),plt(11,j),plt(12,j),plt(13,j),plt(14,j),plt(15,j),plt(16,j),plt(17,j),plt(18,j));
        end
    end

end
close(vout);
fclose(file1);

figure(2);
cplot = squeeze(sol(:,466,:)).';
loglog(time,[cplot(1,:);cplot(2,:);cplot(3,:);cplot(4,:);cplot(5,:);cplot(6,:);cplot(7,:);cplot(8,:);cplot(9,:);cplot(10,:);cplot(11,:);cplot(12,:);cplot(13,:);cplot(14,:);cplot(15,:);cplot(16,:);cplot(17,:);cplot(18,:)]);
title('Center');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16),name(17),name(18));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file2=fopen('center.txt','w');
fprintf(file2,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16),name(17),name(18));
for j=1:200000
    if rem(j,100)==0 || j<=300
        fprintf(file2,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),cplot(1,j),cplot(2,j),cplot(3,j),cplot(4,j),cplot(5,j),cplot(6,j),cplot(7,j),cplot(8,j),cplot(9,j),cplot(10,j),cplot(11,j),cplot(12,j),cplot(13,j),cplot(14,j),cplot(15,j),cplot(16,j),cplot(17,j),cplot(18,j));
    end
end
fclose(file2);

figure(3);
oplot = squeeze(sol(:,510,:)).';
loglog(time,[oplot(1,:);oplot(2,:);oplot(3,:);oplot(4,:);oplot(5,:);oplot(6,:);oplot(7,:);oplot(8,:);oplot(9,:);oplot(10,:);oplot(11,:);oplot(12,:);oplot(13,:);oplot(14,:);oplot(15,:);oplot(16,:);oplot(17,:);oplot(18,:)]);
title('Between Beam');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16),name(17),name(18));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file3=fopen('oversample.txt','w');
fprintf(file3,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16),name(17),name(18));
for j=1:200000
    if rem(j,100)==0 || j<=300
        fprintf(file3,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),oplot(1,j),oplot(2,j),oplot(3,j),oplot(4,j),oplot(5,j),oplot(6,j),oplot(7,j),oplot(8,j),oplot(9,j),oplot(10,j),oplot(11,j),oplot(12,j),oplot(13,j),oplot(14,j),oplot(15,j),oplot(16,j),oplot(17,j),oplot(18,j));
    end
end
fclose(file3);

figure(4);
eplot = squeeze(sol(:,546,:)).';
loglog(time,[eplot(1,:);eplot(2,:);eplot(3,:);eplot(4,:);eplot(5,:);eplot(6,:);eplot(7,:);eplot(8,:);eplot(9,:);eplot(10,:);eplot(11,:);eplot(12,:);eplot(13,:);eplot(14,:);eplot(15,:);eplot(16,:);eplot(17,:);eplot(18,:)]);
title('Edge');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16),name(17),name(18));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file4=fopen('edge.txt','w');
fprintf(file4,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16),name(17),name(18));
for j=1:200000
    if rem(j,100)==0 || j<=300
        fprintf(file4,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),eplot(1,j),eplot(2,j),eplot(3,j),eplot(4,j),eplot(5,j),eplot(6,j),eplot(7,j),eplot(8,j),eplot(9,j),eplot(10,j),eplot(11,j),eplot(12,j),eplot(13,j),eplot(14,j),eplot(15,j),eplot(16,j),eplot(17,j),eplot(18,j));
    end
end
fclose(file4);


figure(5);
dplot = squeeze(sol(:,750,:)).';
loglog(time,[dplot(1,:);dplot(2,:);dplot(3,:);dplot(4,:);dplot(5,:);dplot(6,:);dplot(7,:);dplot(8,:);dplot(9,:);dplot(10,:);dplot(11,:);dplot(12,:);dplot(13,:);dplot(14,:);dplot(15,:);dplot(16,:);dplot(17,:);dplot(18,:)]);
title('Distant');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16),name(17),name(18));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file5=fopen('distant.txt','w');
fprintf(file5,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16),name(17),name(18));
for j=1:200000
    if rem(j,100)==0 || j<=300
        fprintf(file5,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),dplot(1,j),dplot(2,j),dplot(3,j),dplot(4,j),dplot(5,j),dplot(6,j),dplot(7,j),dplot(8,j),dplot(9,j),dplot(10,j),dplot(11,j),dplot(12,j),dplot(13,j),dplot(14,j),dplot(15,j),dplot(16,j),dplot(17,j),dplot(18,j));
    end
end
fclose(file5);
close all;


function [c,f,s] = pdefun(x,t,C,dCdx)

global D;
global generation;

c = [1/D(1);1/D(2);1/D(3);1/D(4);1/D(5);1/D(6);1/D(7);1/D(8);1/D(9);1/D(10);1/D(11);1/D(12);1/D(13);1/D(14);1/D(15);1/D(16);1/D(17);1/D(18)];
f = [dCdx(1);dCdx(2);dCdx(3);dCdx(4);dCdx(5);dCdx(6);dCdx(7);dCdx(8);dCdx(9);dCdx(10);dCdx(11);dCdx(12);dCdx(13);dCdx(14);dCdx(15);dCdx(16);dCdx(17);dCdx(18)];

n = ceil(t/3e-6);
n = rem(n,10);

if x>=(45.45+n) && x<=(45.55+n)
    s = [generation(1)/D(1);generation(2)/D(2);generation(3)/D(3);generation(4)/D(4);generation(5)/D(5);generation(6)/D(6);generation(7)/D(7);generation(8)/D(8);generation(9)/D(9);generation(10)/D(10);generation(11)/D(11);generation(12)/D(12);generation(13)/D(13);generation(14)/D(14);generation(15)/D(15);generation(16)/D(16);generation(17)/D(17);generation(18)/D(18)];

else
    s = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

end


for i = 1:18
    if C(i)<0
        C(i)=0;
    end
end

%%%%%% Reaction between the spiecies %%%%%%

k = 1e-6 .* [8.91e-8; 9.03e-8; 9.03e-8; 1.32e-7; 1.81e-7; 1.69e-7; 1.81e-7; 6.02e10; 8.43e11; 1.21e10; 1.20e11; 6.02e7; 9.03e10; 6.02e11; 2.71e-8; 1.20e11; 6.02e8; 1.26e9; 8.43e-11; 1.81e11; 9.64e11; 6.62e11; 1.32e11; 1.20e11; 6.02e9; 1.81e-7; 1.81e-7; 1.99e-9; 1.20e15; 3.87e14; 1.21e14; 1.27e14; 6.08e14; 6.02e14; 6.08e12; 6.08e14; 6.08e12; 3.03e14; 1.21e15; 1.21e15; 1.86e2; 1.24e1; 6.20e0; 1.86e2; 6.20e0; 1.86e2; 1.86e0; 1.86e0; 1.81e2; 1.81e-9; 1.81e14; 1.81e14; 0; 6.02e1; 1.77e-12; 1.81e10; 1.81e-8; 1.81e-9; 3.01e12; 1.75e-1; 1.75e2; 1.81e2; 1.81e2; 1.75e2; 1.81e14; 1.81e14; 1.81e12; 1.81e12; 1.81e14; 1.81e14; 6.02e13; 6.02e13; 6.02e13; 1.81e-9; 1.81e-9; 1.81e-8]; %6.02e29

% Reaction Set%
r = zeros (76,1);

r(1) = k(1) * C(9) * C(8);
r(2) = k(2) * C(9);
r(3) = k(3) * C(9) * C(6);
r(4) = k(4) * C(11);
r(5) = k(5) * C(11) * C(8);
r(6) = k(6) * C(7);
r(7) = k(7) * C(7) * C(8);
r(8) = k(8) * C(11) * C(6);
r(9) = k(9) * C(9);
r(10) = k(10) * C(5) * C(6);
r(11) = k(11) * C(9) * C(6);
r(12) = k(12) * C(13) * C(8);
r(13) = k(13) * C(17) * C(6);
r(14) = k(14) * C(12) * C(8);
r(15) = k(15) * C(17) * C(6);
r(16) = k(16) * C(12);
r(17) = k(17) * C(15) * C(6);
r(18) = k(18) * C(15);
r(19) = k(19) * C(15);
r(20) = k(20) * C(5);
r(21) = k(21) * C(3);
r(22) = k(22) * C(3) * C(6);
r(23) = k(23) * C(17) * C(8);
r(24) = k(24) * C(16) * C(8);
r(25) = k(25) * C(16) * C(6);
r(26) = k(26) * C(17) * C(8);
r(27) = k(27) * C(16) * C(8);
r(28) = k(28) * C(1) * C(6);
r(29) = k(29) * C(18) * C(1);
r(30) = k(30) * C(9) * C(1);
r(31) = k(31) * C(11) * C(1);
r(32) = k(32) * C(7) * C(1);
r(33) = k(33) * C(12) * C(1);
r(34) = k(34) * C(12) * C(1);
r(35) = k(35) * C(13) * C(1);
r(36) = k(36) * C(13) * C(1);
r(37) = k(37) * C(13) * C(1);
r(38) = k(38) * C(15) * C(1);
r(39) = k(39) * C(17) * C(1);
r(40) = k(40) * C(16) * C(1);
r(41) = k(41) * C(11) * C(1);
r(42) = k(42) * C(9) * C(1);
r(43) = k(43) * C(5) * C(1);
r(44) = k(44) * C(7) * C(1);
r(45) = k(45) * C(3) * C(1);
r(46) = k(46) * C(12) * C(1);
r(47) = k(47) * C(12) * C(1);
r(48) = k(48) * C(13) * C(1);
r(49) = k(49) * C(15) * C(1);
r(50) = k(50) * C(16) * C(1);
r(51) = k(51) * C(9) * C(14);
r(52) = k(52) * C(7) * C(14);
r(53) = k(53) * C(10);
r(54) = k(54) * C(2);
r(55) = k(55) * C(4) * C(4);
r(56) = k(56) * C(2) * C(6);
r(57) = k(57) * C(13) * C(8);
r(58) = k(58) * C(15) * C(1);
r(59) = k(59) * C(15) * C(1);
r(60) = k(60) * C(13) * C(1);
r(61) = k(61) * C(13) * C(1);
r(62) = k(62) * C(18) * C(1);
r(63) = k(63) * C(16) * C(1);
r(64) = k(64) * C(17) * C(1);
r(65) = k(65) * C(18) * C(14);
r(66) = k(66) * C(11) * C(14);
r(67) = k(67) * C(7) * C(14);
r(68) = k(68) * C(12) * C(14);
r(69) = k(69) * C(12) * C(14);
r(70) = k(70) * C(17) * C(14);
r(71) = k(71) * C(15) * C(14);
r(72) = k(72) * C(16) * C(14);
r(73) = k(73) * C(13) * C(14);
r(74) = k(74) * C(16) * C(6);
r(75) = k(75) * C(15) * C(6);
r(76) = k(76) * C(12) * C(8);


products = zeros(18,1);
products(1) = -r(28) - r(29) - r(30) - r(31) - r(32) - r(33) - r(34) - r(35) - r(36) - r(37) - r(38) - r(39) - r(40) - r(41) - r(42) - r(43) - r(44) - r(45) - r(46) - r(47) - r(48) - r(49) - r(50) - r(58) - r(59) - r(60) - r(61) - r(62) - r(63) - r(64);
products(2) = r(30) + r(42) + r(45) + r(58) + r(59) - r(54) - r(56);
products(3) = -r(21) - r(22) - r(45);
products(4) = r(10) + r(22) + r(30) + r(31) + 2 * r(32) + r(34) + 2 * r(36) + r(37) + r(39) + r(40) + r(41) + r(42) + r(43) + 2 * r(44) + r(47) + r(53) + r(56) + r(59) + 2 * r(61) + r(64) + r(69) - 2 * r(55);
products(5) = -r(10) - r(20) - r(43);
products(6) = r(33) + r(35) + r(37) + r(46) + r(48) + r(51) + 2 * r(52) + r(55) + r(60) + r(66) + r(67) + 2 * r(68) + r(69) + r(70) + r(71) + r(72) + 2 * r(73) - r(3) - r(8) - r(10) - r(11) - r(13) - r(15) - r(17) - r(22) - r(25) - r(28) - r(56) - r(74) - r(75);
products(7) = r(8) + r(10) + r(11) + r(13) + r(15) + r(16) + r(17) + r(20) + r(25) + r(74) + r(75) - r(6) - r(7) - r(32) - r(44) - r(52) - r(67);
products(8) = r(9) + r(11) + r(16) + 2 * r(17) + r(18) + r(20) + r(21) + r(25) + 2 * r(29) + r(31) + r(33) + r(37) + 2 * r(38) + r(39) + 2 * r(40) + r(41) + r(46) + 2 * r(49) + r(50) + r(51) + r(53) + 2 * r(54) + r(56) + r(59) + r(60) + 2 * r(62) + 2 * r(63) + r(64) + 2 * r(65) + r(67) + r(68) + 2 * r(71) + r(72) + r(74) + 2 * r(75) - r(1) - r(5) - r(7) - r(12) - r(14) - r(23) - r(24) - r(26) - r(27) - r(57) - r(76);
products(9) = r(18) + r(21) + r(22) - r(1) - r(2) - r(3) - r(9) - r(11) - r(30) - r(42) - r(51);
products(10) = r(8) + r(12) + 2 * r(13) + r(14) + r(15) + r(23) + r(24) + r(25) + r(26) + r(27) + r(29) + r(34) + r(35) + r(36) + r(39) + r(47) + r(48) + r(50) + r(51) + r(52) + r(57) + r(58) + r(65) + 2 * r(66) + r(69) + 2 * r(70) + r(72) + r(73) + r(76) - r(2) - r(4) - r(6) - r(9) - r(19) - r(20) - r(21) - r(28) - r(53)- r(54) - r(75);
products(11) = r(9) + r(14) + r(76) - r(4) - r(5) - r(8) - r(31) - r(41) - r(66);
products(12) = r(3) + r(7) + r(12) + r(57) - r(14) - r(16) - r(33) - r(34) - r(46) - r(47) - r(68) - r(69) - r(76);
products(13) = r(6) - r(12) - r(35) - r(36) - r(37) - r(48) - r(57) - r(60) - r(61) - r(73);
products(14) = r(28) - r(51) - r(52) -r(65) - r(66) - r(67) - r(68) - r(69) - r(70) - r(71) - r(72) - r(73);
products(15) = r(1) + r(24) + r(27) - r(17) - r(18) - r(19) - r(38) - r(49) - r(58) - r(59) - r(71) - r(75);
products(16) = r(2) + r(5) + r(23) + r(26) - r(24) - r(25) - r(27) - r(40) - r(50) - r(63) - r(72) - r(74);
products(17) = r(4) - r(13) - r(15) - r(23) - r(26) - r(39) - r(64) - r(70);
products(18) = r(19) - r(29) - r(62) - r(65);


s(1) = s(1) + products(1)/D(1);
s(2) = s(2) + products(2)/D(2);
s(3) = s(3) + products(3)/D(3);
s(4) = s(4) + products(4)/D(4);
s(5) = s(5) + products(5)/D(5);
s(6) = s(6) + products(6)/D(6);
s(7) = s(7) + products(7)/D(7);
s(8) = s(8) + products(8)/D(8);
s(9) = s(9) + products(9)/D(9);
s(10) = s(10) + products(10)/D(10);
s(11) = s(11) + products(11)/D(11);
s(12) = s(12) + products(12)/D(12);
s(13) = s(13) + products(13)/D(13);
s(14) = s(14) + products(14)/D(14);
s(15) = s(15) + products(15)/D(15);
s(16) = s(16) + products(16)/D(16);
s(17) = s(17) + products(17)/D(17);
s(18) = s(18) + products(18)/D(18);


end

function C0 = icfun(x,arg)

C0=[arg(1,round(x*10) + 1);arg(2,round(x*10) + 1);arg(3,round(x*10) + 1);arg(4,round(x*10) + 1);arg(5,round(x*10) + 1);arg(6,round(x*10) + 1);arg(7,round(x*10) + 1);arg(8,round(x*10) + 1);arg(9,round(x*10) + 1);arg(10,round(x*10) + 1);arg(11,round(x*10) + 1);arg(12,round(x*10) + 1);arg(13,round(x*10) + 1);arg(14,round(x*10) + 1);arg(15,round(x*10) + 1);arg(16,round(x*10) + 1);arg(17,round(x*10) + 1);arg(18,round(x*10) + 1)];

end

function [pL,qL,pR,qR] = bcfun(xL,CL,xR,CR,t)

pR = CR;
pR(10) = CR(10)-4.09e4 *0.013;
qR = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

pL = CL;
pL(10) = CL(10)-4.09e4 *0.013;
qL = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

end

