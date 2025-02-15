clc
clear

global numberOfSpecies
global D
global generation

global Pressure

Pressure=0.013;

numberOfSpecies = 10;


% Order: [“e-”, “O”, “O+”, “O2”, “O2+”, “O2-”, “O3”, “O3-”, “O4+”, “O4-”]

molarMass = [5.49E-4 15.999 15.999 31.998 31.998 31.998 47.997 47.997 63.996 63.996];
molecularVolume = [0 6.11 6.11 16.3 16.3 16.3 22.4 22.4 28.5 28.5];

for i = 1:numberOfSpecies
    Mab(i) = 2/(1/molarMass(i) + 1/molarMass(4));
    D(i) = 0.00143 * 1e14 * 298.^1.75 / (sqrt(Mab(i)) * ((molecularVolume(i).^1/3) + molecularVolume(4).^1/3).^2 * Pressure); 
end


initialConcentration = [0 0 0 4.09e4*Pressure 0 0 0 0 0 0];  % microMole

name =  ["e-" "O" "O+" "O2" "O2+" "O2-" "O3" "O3-" "O4+" "O4-"];
name2 = ["e^{-}" "O" "O^{+}" "O_{2}" "O_{2}^{+}" "O_{2}^{-}" "O_{3}" "O_{3}^{-}" "O_{4}^{+}" "O_{4}^{-}"];

doseRate = 1e5 * 2.446 * 6.2e-11 / (3.141592 * 2.5e-21) ; % in Gy/s at 200keV 

gValues = [3.27/100 6.10/100 0.1/100 -6.32/100 3.17/100 0 0 0 0 0];  % molecules / eV

density = 1.42e-3 * Pressure; %g/cm3

generation = 1e6 * density * doseRate * gValues / (6.022e23 *1.6e-19); % generation in uM

% Geometry description = 100 nm 1D mesh, index = 1001
% probe size = 0.1 nm, probe interval = 1 nm, 10 point (10nm) scan

mesh = zeros(numberOfSpecies,1001);

for i = 1:numberOfSpecies
    mesh(i,:) = initialConcentration(i);
end

time = linspace(1e-9,1e-4,400001);
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

for i = 1:400001
    if i<=400 || rem(i,100)==0
        u = sol(i,:,:);
        plt=squeeze(u).';
        semilogy(x,[plt(4,:);plt(5,:);plt(2,:);plt(3,:);plt(1,:);plt(9,:);plt(6,:);plt(7,:)],"LineWidth",1);
        line_color = ["#8C8C00";"#FF6300";"#00CECE";"#1A6FDF";"#525252";"#8E8E00";"#CC9900";"#B177DE"];
        colororder(line_color);

        xlabel('X (nm)');
        ylabel('C (uMol)')
        xlim([30 70]);
        ylim([1e-10 1e5]);
        title(sprintf('Time: %.2d sec',time(i)));
        legend(name2(4),name2(5),name2(2),name2(3),name2(1),name2(9),name2(6),name2(7),'Location','southeastoutside');
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
        fprintf(file1,'x (nm)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10));
        for j=1:1000
            fprintf(file1,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',x(j),plt(1,j),plt(2,j),plt(3,j),plt(4,j),plt(5,j),plt(6,j),plt(7,j),plt(8,j),plt(9,j),plt(10,j));
        end
    end

end
close(vout);
fclose(file1);

figure(2);
cplot = squeeze(sol(:,466,:)).';
loglog(time,[cplot(1,:);cplot(2,:);cplot(3,:);cplot(4,:);cplot(5,:);cplot(6,:);cplot(7,:);cplot(8,:);cplot(9,:);cplot(10,:)]);
title('Center');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file2=fopen('center.txt','w');
fprintf(file2,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10));
for j=1:400001
    if rem(j,100)==0 || j<=300
        fprintf(file2,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),cplot(1,j),cplot(2,j),cplot(3,j),cplot(4,j),cplot(5,j),cplot(6,j),cplot(7,j),cplot(8,j),cplot(9,j),cplot(10,j));
    end
end
fclose(file2);

figure(3);
oplot = squeeze(sol(:,510,:)).';
loglog(time,[oplot(1,:);oplot(2,:);oplot(3,:);oplot(4,:);oplot(5,:);oplot(6,:);oplot(7,:);oplot(8,:);oplot(9,:);oplot(10,:)]);
title('Between Beam');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file3=fopen('oversample.txt','w');
fprintf(file3,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10));
for j=1:400001
    if rem(j,100)==0 || j<=300
        fprintf(file3,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),oplot(1,j),oplot(2,j),oplot(3,j),oplot(4,j),oplot(5,j),oplot(6,j),oplot(7,j),oplot(8,j),oplot(9,j),oplot(10,j));
    end
end
fclose(file3);

figure(4);
eplot = squeeze(sol(:,546,:)).';
loglog(time,[eplot(1,:);eplot(2,:);eplot(3,:);eplot(4,:);eplot(5,:);eplot(6,:);eplot(7,:);eplot(8,:);eplot(9,:);eplot(10,:)]);
title('Edge');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file4=fopen('edge.txt','w');
fprintf(file4,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10));
for j=1:400001
    if rem(j,100)==0 || j<=300
        fprintf(file4,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),eplot(1,j),eplot(2,j),eplot(3,j),eplot(4,j),eplot(5,j),eplot(6,j),eplot(7,j),eplot(8,j),eplot(9,j),eplot(10,j));
    end
end
fclose(file4);



figure(5);
dplot = squeeze(sol(:,750,:)).';
loglog(time,[dplot(1,:);dplot(2,:);dplot(3,:);dplot(4,:);dplot(5,:);dplot(6,:);dplot(7,:);dplot(8,:);dplot(9,:);dplot(10,:)]);
title('Distant');
legend(name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10));
xlim([1e-9 1e-4]);
ylim([1e-10 1e5]);

file5=fopen('distant.txt','w');
fprintf(file5,'t (s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',name(1),name(2),name(3),name(4),name(5),name(6),name(7),name(8),name(9),name(10));
for j=1:400001
    if rem(j,100)==0 || j<=300
        fprintf(file5,'%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\t%.2d\n',time(j),dplot(1,j),dplot(2,j),dplot(3,j),dplot(4,j),dplot(5,j),dplot(6,j),dplot(7,j),dplot(8,j),dplot(9,j),dplot(10,j));
    end
end
fclose(file5);


function [c,f,s] = pdefun(x,t,C,dCdx)

global D;
global generation;

c = [1/D(1);1/D(2);1/D(3);1/D(4);1/D(5);1/D(6);1/D(7);1/D(8);1/D(9);1/D(10)];
f = [dCdx(1);dCdx(2);dCdx(3);dCdx(4);dCdx(5);dCdx(6);dCdx(7);dCdx(8);dCdx(9);dCdx(10)];

%%%%%% Probe Generation if the t<=dwell time (continuous probes interlacing) %%%%%%
n = ceil(t/3e-6);
n = rem(n,10);

if x>=(45.45+n) && x<=(45.55+n)
    s = [generation(1)/D(1);generation(2)/D(2);generation(3)/D(3);generation(4)/D(4);generation(5)/D(5);generation(6)/D(6);generation(7)/D(7);generation(8)/D(8);generation(9)/D(9);generation(10)/D(10)];

else
    s = [0;0;0;0;0;0;0;0;0;0];

end


%%%%%% Reaction between the spiecies %%%%%%
k = 1e-6 .*[4.2e10; 1.2e10; 2.0e10; 4.2e10; 2.0e13; 1.2e15; 1.2e15; 1.0e14; 1.2e15; 1.2e15; 1.0e14; 1.2e15; 1.2e15; 1.8e11; 1.2e15; 8.0e6; 4.0e7; 6.0e6]; 

% Reaction Set%

r=zeros(18,1);

r(1) = k(1) * C(5);
r(2) = k(2) * C(3);
r(3) = k(3) * C(1);
r(4) = k(4) * C(6);
r(5) = k(5) * C(3) * C(1);
r(6) = k(6) * C(3) * C(6);
r(7) = k(7) * C(3) * C(10);
r(8) = k(8) * C(5) * C(1);
r(9) = k(9) * C(5) * C(6);
r(10) = k(10) * C(5) * C(10);
r(11) = k(11) * C(9) * C(1);
r(12) = k(12) * C(9) * C(6);
r(13) = k(13) * C(9) * C(10);
r(14) = k(14) * C(6) * C(7);
r(15) = k(15) * C(8) * C(9);
r(16) = k(16) * C(2);
r(17) = k(17) * C(2) * C(2);
r(18) = k(18) * C(2) * C(7);

products=zeros(10,1);

products(1) = -r(3) - r(5) - r(8) - r(11);
products(2) = r(2) + r(5) + r(6) + r(7) + 2 * r(8) + 2 * r(9) + 2 * r(10) + 2 * r(11) + 2 * r(12) + 2 * r(13) - r(16) - 2 * r(17) - r(18);
products(3) = -r(2) - r(5) - r(6) - r(7);
products(4) = r(3) + r(6) + 2 * r(7) + r(9) + 2 * r(10) + r(11) + 2 * r(12) + 3 * r(13) + r(14) + 2 * r(15) + r(16) + r(17) + 2 * r(18) - r(1) - r(2) - 2 * r(3) - r(4) - 2 * r(16);
products(5) = r(2) - r(1) - r(8) - r(9) - r(10);
products(6) = r(3) - r(4) - r(6) - r(9) - r(12) - r(14);
products(7) = r(15) + r(16) - r(14) - r(18);
products(8) = r(14) - r(15);
products(9) = r(1) - r(11) - r(12) - r(13) - r(15);
products(10) = r(4) - r(7) - r(10) - r(13);


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


end

function C0 = icfun(x,arg)

C0=[arg(1,round(x*10) + 1);arg(2,round(x*10) + 1);arg(3,round(x*10) + 1);arg(4,round(x*10) + 1);arg(5,round(x*10) + 1);arg(6,round(x*10) + 1);arg(7,round(x*10) + 1);arg(8,round(x*10) + 1);arg(9,round(x*10) + 1);arg(10,round(x*10) + 1)];

end

function [pL,qL,pR,qR] = bcfun(xL,CL,xR,CR,t)

global Pressure

pR = CR;
pR(4) = CR(4)-4.09e4*Pressure;
qR = [0;0;0;0;0;0;0;0;0;0];

pL = CL;
pL(4) = CL(4)-4.09e4*Pressure;
qL = [0;0;0;0;0;0;0;0;0;0];

end

