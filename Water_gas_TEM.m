clc
clear

global numberOfSpecies
global D
global generation


numberOfSpecies = 16;

% Order: [electrons H+ OH- H2O2 HO2- H OH O- HO2 O2- O2 H2 O3- O3 HO3 H2O];

molarMass = [5.49e-4 1.008 17.008 34.015 33.007 1.008 17.008 15.999 33.007 31.999 31.999 2.016 47.998 47.998 49.006 18.015];
molecularVolume = [0 0 8.42 19.21 18.61 2.31 8.42 6.11 18.61 16.3 16.3 6.12 22.41 22.41 24.72 13.1];

for i = 1:numberOfSpecies
    Mab(i) = 2/(1/molarMass(i) + 1/molarMass(16));
    D(i) = 0.00143 * 1e14 * 298.^1.75 / (sqrt(Mab(i)) * ((molecularVolume(i).^1/3) + molecularVolume(16).^1/3).^2 * 0.013); %1atm, 0.5atm, 0.013 atm
end

D = D .* 1e-6; % um2/s

% Order: [electrons H+ OH- H2O2 HO2- H OH O- HO2 O2- O2 H2 O3- O3 HO3 H2O];
initialConcentration = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.09e4*0.013];


name = ["e-" "H+" "OH-" "H2O2" "HO2-" "H" "OH" "O-" "HO2" "O2-" "O2" "H2" "O3-" "O3" "HO3" "H2O"];
name2 = ["e^{-}" "H^{+}" "OH^{-}" "H_{2}O_{2}" "HO_{2}^{-}" "H" "OH" "O^{-}" "HO_{2}" "O_{2}^{-}" "O_{2}" "H_{2}" "O_{3}^{-}" "O_{3}" "HO_{3}" "H_{2}O"];

r = 2e-6; 
doseRate = 1e5 * 2.814 * 4.2e-9 / (3.141592 * r^2) ; % screen current 4.2e-9A

gValues = [3.47/100 4.42/100 0.95/100 0.47/100 0 1/100 3.63/100 0 0.08/100 0 0 0.17/100 0 0 0 -5.68/100];  % molecules / eV

density = 8.04e-4* 0.013; %g/cm3

generation = 1e6 * density * doseRate * gValues / (6.022e23 *1.6e-19); % generation in uM


mesh = zeros(numberOfSpecies,1001);

for i = 1:numberOfSpecies
    mesh(i,:) = initialConcentration(i);
end


x=linspace (0,50,1001);


time = linspace(1e-9,1e-4,400001);


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

for i = 1:400001
    if i<=400 || rem(i,100)==0
        u = sol(i,:,:);
        plt=squeeze(u).';

        semilogy(x,[plt(7,:);plt(4,:);plt(6,:);plt(12,:);plt(1,:);plt(10,:)],"LineWidth",1);
        line_color = ["#CC9900";"#242424";"#F14040";"#1A6FDF";"#AFAFAF";"#FB6501"];
        colororder(line_color);

        xlabel('X (\mum)');
        ylabel('C (\muMol)');

        xlim([17 33]);
        ylim([1e-10 1e-1]);
        title(sprintf('Time: %.2d sec',time(i)));
        legend(name2(7),name2(4),name2(6),name2(12),name2(1),name2(10),'Location','southeastoutside');
        set(gcf,'position',[500 500 900 600]);
        pause(0.001);
        frame=getframe(gcf);
        writeVideo(vout,frame);
        drawnow;
    end

    if i==1 || i==37 || i==397 || i==3997 || i==39997 || i==399997
        u = sol(i,:,:);
        plt=squeeze(u).';
        fprintf(file1,'\n\n\nt = %.2d sec\n',time(i));
        fprintf(file1,'x (um)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16));
        for j=1:1000
            fprintf(file1,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',x(j),plt(1,j),plt(2,j),plt(3,j),plt(4,j),plt(5,j),plt(6,j),plt(7,j),plt(8,j),plt(9,j),plt(10,j),plt(11,j),plt(12,j),plt(13,j),plt(14,j),plt(15,j),plt(16,j));
        end
    end

end
close(vout);
fclose(file1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot placewise (center x=505, edge x=546, distance x= 750)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
cplot = squeeze(sol(:,500,:)).';
loglog(time,[cplot(1,:);cplot(2,:);cplot(3,:);cplot(4,:);cplot(5,:);cplot(6,:);cplot(7,:);cplot(8,:);cplot(9,:);cplot(10,:);cplot(11,:);cplot(12,:);cplot(13,:);cplot(14,:);cplot(15,:);cplot(16,:)]);
title('Center');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file2=fopen('center.txt','w');
fprintf(file2,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16));
for j=1:400001
    if rem(j,100)==0 || j<=300
        fprintf(file2,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),cplot(1,j),cplot(2,j),cplot(3,j),cplot(4,j),cplot(5,j),cplot(6,j),cplot(7,j),cplot(8,j),cplot(9,j),cplot(10,j),cplot(11,j),cplot(12,j),cplot(13,j),cplot(14,j),cplot(15,j),cplot(16,j));
    end
end
fclose(file2);

figure(4);
eplot = squeeze(sol(:,541,:)).';
loglog(time,[eplot(1,:);eplot(2,:);eplot(3,:);eplot(4,:);eplot(5,:);eplot(6,:);eplot(7,:);eplot(8,:);eplot(9,:);eplot(10,:);eplot(11,:);eplot(12,:);eplot(13,:);eplot(14,:);eplot(15,:);eplot(16,:)]);
title('Edge');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file4=fopen('edge.txt','w');
fprintf(file4,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16));
for j=1:400001
    if rem(j,100)==0 || j<=300
        fprintf(file4,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),eplot(1,j),eplot(2,j),eplot(3,j),eplot(4,j),eplot(5,j),eplot(6,j),eplot(7,j),eplot(8,j),eplot(9,j),eplot(10,j),eplot(11,j),eplot(12,j),eplot(13,j),eplot(14,j),eplot(15,j),eplot(16,j));
    end
end
fclose(file4);


figure(5);
dplot = squeeze(sol(:,750,:)).';
loglog(time,[dplot(1,:);dplot(2,:);dplot(3,:);dplot(4,:);dplot(5,:);dplot(6,:);dplot(7,:);dplot(8,:);dplot(9,:);dplot(10,:);dplot(11,:);dplot(12,:);dplot(13,:);dplot(14,:);dplot(15,:);dplot(16,:)]);
title('Distant');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file5=fopen('distant.txt','w');
fprintf(file5,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10),name(11),name(12),name(13),name(14),name(15),name(16));
for j=1:400001
    if rem(j,100)==0 || j<=300
        fprintf(file5,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),dplot(1,j),dplot(2,j),dplot(3,j),dplot(4,j),dplot(5,j),dplot(6,j),dplot(7,j),dplot(8,j),dplot(9,j),dplot(10,j),dplot(11,j),dplot(12,j),dplot(13,j),dplot(14,j),dplot(15,j),dplot(16,j));
    end
end
fclose(file5);



function [c,f,s] = pdefun(x,t,C,dCdx)

global D;
global generation;

c = [1/D(1);1/D(2);1/D(3);1/D(4);1/D(5);1/D(6);1/D(7);1/D(8);1/D(9);1/D(10);1/D(11);1/D(12);1/D(13);1/D(14);1/D(15);1/D(16)];
f = [dCdx(1);dCdx(2);dCdx(3);dCdx(4);dCdx(5);dCdx(6);dCdx(7);dCdx(8);dCdx(9);dCdx(10);dCdx(11);dCdx(12);dCdx(13);dCdx(14);dCdx(15);dCdx(16)];


if x>=23 && x<=27
    s = [generation(1)/D(1);generation(2)/D(2);generation(3)/D(3);generation(4)/D(4);generation(5)/D(5);generation(6)/D(6);generation(7)/D(7);generation(8)/D(8);generation(9)/D(9);generation(10)/D(10);generation(11)/D(11);generation(12)/D(12);generation(13)/D(13);generation(14)/D(14);generation(15)/D(15);generation(16)/D(16)];

else
    s = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

end

%%%%%% Override zero when C is negative %%%%%%

for i = 1:16
    if C(i)<0
        C(i)=0;
    end
end

%equilibria%
K = zeros(5,1);

%H2O <=> H+ + OH-&
K(1) = 10^-13.999; %Molar^2

%H2O2 <=> H+ + HO2-%
K(2) = 10^-11.65;%Molar

%OH <=> H+ + O-%
K(3) = 10^-11.9;%Molar

%HO2 <=> H+ + O2-%
K(4) = 10^-4.57;%Molar

%H <=> H+ + e-%
K(5) = 10^-9.77;%Molar

%%%%%% Reaction between the spiecies %%%%%%

molarToMicromolar = 10^6;

% Note, some are out of numerical order so they are defined when they are
% used in other rate constants. To convert numbering scheme from what is
% presented here to what it presented in the paper [ref], add 6 to the
% index. (i.e., Reaction 25's rate constant has an index of 19)
k = zeros(73,1);
k(1) = 1.4*10^11;
k(2) = k(1)*K(1)*molarToMicromolar^2;%M/s
k(4) = 5*10^10;
k(3) = k(4)*K(2)*molarToMicromolar;%/s
k(5) = 1.3*10^10;
k(6) = k(5)*(K(1)/K(2))*molarToMicromolar;%%/s
k(7) = 1.9*10;
k(8) = 2.2*10^7;
k(10) = 2.3*10^10;
k(9) = k(10)*K(5)*molarToMicromolar;%/s
k(11) = 1.3*10^10;
k(12) = k(11)*(K(1)/K(3))*molarToMicromolar;%/s
k(14) = 10^11;
k(13) = k(14)*K(3)*molarToMicromolar;%/s
k(16) = 5*10^10;
k(15) = k(16)*K(4)*molarToMicromolar;%/s
k(17) = 5*10^10;
k(18) = k(17)*(K(1)/K(4))*molarToMicromolar;%/s%
k(19) = 3*10^10;
k(20) = 1.1*10^10;
k(21) = 1.3*10^10 ;
k(22) = 2*10^10;
k(23) = 1.9*10^10;
k(24) = 5.5*10^9; 
k(25) = 2.5*10^10 ;
k(26) = 3.5*10^9;
k(27) = 2.2*10^10 ;
k(28) = 1.6*10^10; 
k(29) = 3.6*10^10;
k(30) = 1.1*10;
k(31) = 10^10;
k(32) = 9*10^7;
k(33) = 10^10;
k(34) = 7.8*10^9;
k(35) = 7.0*10^9;
k(36) = 9*10^7;
k(37) = 2.1*10^10;
k(38) = 1.8*10^10;
k(39) = 1.8*10^10;
k(40) = 3.8*10^10;
k(41) = 3.6*10^9;
k(42) = 6*10^9; 
k(43) = 8.2*10^9;
k(44) = 4.3*10^7;
k(45) = 2.7*10^7;
k(46) = 2.5*10^10;
k(47) = 7.5*10^9;
k(48) = 2.6*10^9;
k(49) = 6*10^9;
k(50) = 1.1*10^8;
k(51) = 8*10^7;
k(52) = 7*10^5;
k(53) = 6*10^9;
k(54) = 5*10^-1;
k(55) = 5*10^-1;
k(56) = 6*10^9;
k(57) = 5*10^8;
k(58) = 10^2;
k(59) = 6*10^8;
k(60) = 1.3*10^-1;
k(61) = 1.3*10^-1;
k(62) = 10^4 ;
k(63) = 1.5*10^9;
k(64) = 10^9 ;
k(65) = 3.6*10^9;
k(66) = 8*10^7;
k(67) = 5*10^8;
k(68) = 4*10^8;
k(69) = 7*10^8;
k(70) = 5*10^9;
k(71) = 3.3*10^3*molarToMicromolar;%/s
k(72) = 9*10^10;
k(73) = 1.1*10^5*molarToMicromolar;%/s


%%convert from molar to micromolar
k = k./(molarToMicromolar);

%Reactant concentrations (micromolar)%
electrons = C(1); 
H_ions    = C(2);
OH_ions   = C(3);
H2O2      = C(4);
HO2_ions  = C(5);
H_mono    = C(6);
OH        = C(7);
O_ions    = C(8);
HO2       = C(9);
O2_ions   = C(10);
O2        = C(11);
H2        = C(12);
O3_ions   = C(13);
O3        = C(14);
HO3       = C(15);
H2O       = C(16);


%Reaction Set%
% Note: H2O is divided out as a reactant
r = zeros(73,1);
r(1)  = k(1)*H_ions*OH_ions;
r(2)  = k(2);
r(3)  = k(3)*H2O2;
r(4)  = k(4)*H_ions*HO2_ions;
r(5)  = k(5)*H2O2*OH_ions;
r(6)  = k(6)*HO2_ions;
r(7)  = k(7)*electrons*H2O;
r(8)  = k(8)*H_mono*OH_ions;
r(9)  = k(9)*H_mono;
r(10) = k(10)*electrons*H_ions;
r(11) = k(11)*OH*OH_ions;
r(12) = k(12)*O_ions;
r(13) = k(13)*OH;
r(14) = k(14)*O_ions*H_ions;
r(15) = k(15)*HO2;
r(16) = k(16)*O2_ions*H_ions;
r(17) = k(17)*HO2*OH_ions;
r(18) = k(18)*O2_ions;
r(19) = k(19)*electrons*OH;
r(20) = k(20)*electrons*H2O2;
r(21) = k(21)*electrons*O2_ions;
r(22) = k(22)*electrons*HO2;
r(23) = k(23)*electrons*O2;
r(24) = k(24)*electrons^2;
r(25) = k(25)*electrons*H_mono;
r(26) = k(26)*electrons*HO2_ions;
r(27) = k(27)*electrons*O_ions;
r(28) = k(28)*electrons*O3_ions;
r(29) = k(29)*electrons*O3;
r(30) = k(30)*H_mono*H2O; 
r(31) = k(31)*H_mono*O_ions;
r(32) = k(32)*H_mono*HO2_ions;
r(33) = k(33)*H_mono*O3_ions;
r(34) = k(34)*H_mono^2;
r(35) = k(35)*H_mono*OH;
r(36) = k(36)*H_mono*H2O2;
r(37) = k(37)*H_mono*O2;
r(38) = k(38)*H_mono*HO2;
r(39) = k(39)*H_mono*O2_ions;
r(40) = k(40)*H_mono*O3;
r(41) = k(41)*OH^2;
r(42) = k(42)*OH*HO2; 
r(43) = k(43)*OH*O2_ions;
r(44) = k(44)*OH*H2;
r(45) = k(45)*OH*H2O2;
r(46) = k(46)*OH*O_ions;
r(47) = k(47)*OH*HO2_ions;
r(48) = k(48)*OH*O3_ions;
r(49) = k(49)*OH*O3_ions;
r(50) = k(50)*OH*O3;
r(51) = k(51)*HO2*O2_ions;
r(52) = k(52)*HO2^2;
r(53) = k(53)*HO2*O_ions;
r(54) = k(54)*HO2*H2O2;
r(55) = k(55)*HO2*HO2_ions;
r(56) = k(56)*HO2*O3_ions;
r(57) = k(57)*HO2*O3;
r(58) = k(58)*O2_ions^2;
r(59) = k(59)*O2_ions*O_ions;
r(60) = k(60)*O2_ions*H2O2;
r(61) = k(61)*O2_ions*HO2_ions;
r(62) = k(62)*O2_ions*O3_ions;
r(63) = k(63)*O2_ions*O3;
r(64) = k(64)*O_ions^2;
r(65) = k(65)*O_ions*O2;
r(66) = k(66)*O_ions*H2;
r(67) = k(67)*O_ions*H2O2;
r(68) = k(68)*O_ions*HO2_ions;
r(69) = k(69)*O_ions*O3_ions;
r(70) = k(70)*O_ions*O3;
r(71) = k(71)*O3_ions;
r(72) = k(72)*O3_ions*H_ions;
r(73) = k(73)*HO3;

products = zeros(16,1);


%%electrons%%
products(1) = -r(7) + r(8) + r(9) - r(10) - r(19) - r(20) - r(21) - r(22)...
    -r(23) - 2*r(24) - r(25) - r(26) - r(27) - r(28) - r(29);

%%H+%%
products(2) = -r(1) + r(2) + r(3) - r(4) + r(9) - r(10) + r(13) - r(14)...
    + r(15) - r(16) + r(49) - r(72);

%%OH-%%
products(3) = -r(1) + r(2) - r(5) + r(6) + r(7) - r(8) - r(11) + r(12) -...
    r(17)+ r(18) + r(19) + r(20) + r(21) + 2*r(24) + r(25) + r(26) + ...
    2*r(27) + 2*r(28) + r(31) + r(32) + r(33) + r(43) + r(47) + r(48) + ...
    r(53) + r(55) + r(56) + 2*r(58) + 2*r(59) + r(60) + r(61) + 2*r(62)...
    + r(64) + r(66) + r(68);

%%H2O2%%
products(4) = - r(3) + r(4) -r(5) + r(6) - r(20) - r(36) + r(38) + r(41)...
    - r(45) + r(52) - r(54) + r(58) - r(60) - r(67)   ;

%%HO2-%%
products(5) = r(3) - r(4) + r(5) - r(6) + r(21) + r(22) - r(26) - r(32)...
    + r(39) + r(46) - r(47) + r(51) - r(55) - r(61) + r(64) - r(68);

%%H%%
products(6) = r(7) - r(8) - r(9) + r(10) - r(25) - r(30) - r(31) - r(32)...
    - r(33) - 2*r(34) - r(35) - r(36) - r(37) - r(38) - r(39) - r(40) + r(44)...
    + r(66);

%%OH%%
products(7) = -r(11) + r(12) - r(13) + r(14) - r(19) + r(20) + r(30) ...
    + r(32) - r(35) + r(36) - 2*r(41) - r(42) - r(43) - r(44) - r(45) -...
    r(46) - r(47) - r(48) - r(49) - r(50) + r(54) + r(55) + r(60) + r(72)...
    + r(73);

%%O-%%
products(8) = r(11) - r(12) + r(13) - r(14) + r(26) - r(27) - r(31) - r(46)...
    - r(53) - r(59) + r(61) - 2*r(64) - r(65) - r(66) - r(67) - r(68) -...
    r(69) - r(70) + r(71) ;

%%HO2%%
products(9) = -r(15) + r(16) - r(17) + r(18) - r(22) + r(37) - r(38) - ...
    r(42) + r(45) + r(47) + r(50) - r(51) - 2*r(52) - r(53) - r(54) - r(55)...
    - r(56) - r(57);

%%O2-%%
products(10) = r(15) - r(16) + r(17) - r(18) - r(21) + r(23) - r(39) - ...
    r(43) + 2*r(49) - r(51) - 2*r(58) - r(59) - r(60) - r(61) - r(62) - r(63)...
    + r(67) + r(68) + 2*r(69) + r(70);

%%O2%%
products(11) = -r(23) + r(28) + r(33) - r(37) + r(42) + r(43) + r(50) ...
    + r(51) + r(52) + r(53) + r(54) + r(55) + 2*r(56) + r(57) + r(58) + r(59)...
    + r(60) + r(61) + 2*r(62) + r(63) - r(65) + r(70) + r(71) + r(72) + r(73);

%%H2%%
products(12) = r(24) + r(25) + r(30) + r(34) - r(44) - r(66);

%%O3-%%
products(13) = -r(28) + r(29) - r(33) - r(48) - r(49) - r(56) - r(62) ...
    + r(63) + r(65) - r(69) - r(71) - r(72);

%%O3%%
products(14) = -r(29) - r(40) + r(48) - r(50) - r(57) - r(63) - r(70);

%%HO3%%
products(15) = r(40) + r(57) - r(73);

%%H2O%%
products(16) = r(1) - r(2) + r(5) - r(6) - r(7) + r(8) + r(11) - r(12) + ...
     r(17) - r(18) - r(21) - 2*r(24) - r(25) - r(27) - r(28) - r(30) + r(35)...
     + r(36) + r(42) + r(44) + r(45) + r(54) - 2*r(58) - r(59) - r(62) - r(64) ...
     + r(67);

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


end

function C0 = icfun(x,arg)

C0=[arg(1,round(x*10) + 1);arg(2,round(x*10) + 1);arg(3,round(x*10) + 1);arg(4,round(x*10) + 1);arg(5,round(x*10) + 1);arg(6,round(x*10) + 1);arg(7,round(x*10) + 1);arg(8,round(x*10) + 1);arg(9,round(x*10) + 1);arg(10,round(x*10) + 1);arg(11,round(x*10) + 1);arg(12,round(x*10) + 1);arg(13,round(x*10) + 1);arg(14,round(x*10) + 1);arg(15,round(x*10) + 1);arg(16,round(x*10) + 1)];

end

function [pL,qL,pR,qR] = bcfun(xL,CL,xR,CR,t)

pR = CR;
pR(16) = CR(16)-4.09e4*0.013;
qR = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

pL = CL;
pL(16) = CL(16)-4.09e4*0.013;
qL = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

end

