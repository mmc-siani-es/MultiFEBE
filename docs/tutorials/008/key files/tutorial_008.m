
% Author: Á.G. Vega-Artiles

% This script is part of the Tutorial 8. The code takes the resulting *.tot
% files from Shh, Svv and Smm and creates copies as txt files so that Matlab
% can open them. Then, the rows corresponding to the contributions of the
% cartesian decomposition for the target resultant forces or moments of the
% "bucket-lid" (Boundary_id == 2) and the "soil-skirt" (Class == 0) are
% summed up. Finally, the absolute values (divided by the minimum absolute
% value) and the argument of these results are plotted together with the
% reference data taken from Liingaard (2007).

close all
clear variables

copyfile 'Shh.dat.tot' 'Shh.dat.tot.txt'
filename1 = 'Shh.dat.tot.txt';
A01 = readmatrix(filename1,'NumHeaderLines',21);
A01 = array2table(A01);
A01 = renamevars(A01,{'A012','A016','A017','A0124','A0125'},{'w','Boundary_id','Class','Re_Fx','Im_Fx'});
rows1 = A01.Class == 0;
rows2 = A01.Boundary_id == 2;
lid1 = A01(rows2,:);
load1 = A01(rows1,:);
w = A01.w(rows1);  
a0 = 5*w/31.6227766;

Re_Fx = zeros(100,1);

for i=1:100
    Re_Fx(i) = lid1.Re_Fx(i)+load1.Re_Fx(i);
end

Im_Fx = zeros(100,1);

for i=1:100
    Im_Fx(i) = lid1.Im_Fx(i)+load1.Im_Fx(i);
end     

Abs_Fx = zeros(100,1);

for i=1:100
    Abs_Fx(i) = sqrt(Re_Fx(i)^2 + Im_Fx(i)^2);
end     

Arg_Fx = zeros(100,1);

for i=1:100
    Arg_Fx(i) = atan2(Im_Fx(i),Re_Fx(i));
end    

copyfile 'Svv.dat.tot' 'Svv.dat.tot.txt'
filename2 = 'Svv.dat.tot.txt';
A02 = readmatrix(filename2,'NumHeaderLines',21);
A02 = array2table(A02);
A02 = renamevars(A02,{'A022','A026','A027','A0228','A0229'},{'w','Boundary_id','Class','Re_Fz','Im_Fz'});
rows3 = A02.Class == 0;
rows4 = A02.Boundary_id == 2;
lid2 = A02(rows4,:);
load2 = A02(rows3,:);

Re_Fz = zeros(100,1);

for i=1:100
    Re_Fz(i) = lid2.Re_Fz(i)+load2.Re_Fz(i);
end

Im_Fz = zeros(100,1);

for i=1:100
    Im_Fz(i) = lid2.Im_Fz(i)+load2.Im_Fz(i);
end     

Abs_Fz = zeros(100,1);

for i=1:100
    Abs_Fz(i) = sqrt(Re_Fz(i)^2 + Im_Fz(i)^2);
end     

Arg_Fz = zeros(100,1);

for i=1:100
    Arg_Fz(i) = atan2(Im_Fz(i),Re_Fz(i));
end    

copyfile 'Smm.dat.tot' 'Smm.dat.tot.txt'
filename3 = 'Smm.dat.tot.txt';
A03 = readmatrix(filename3,'NumHeaderLines',21);
A03 = array2table(A03);
A03 = renamevars(A03,{'A032','A036','A037','A0332','A0333'},{'w','Boundary_id','Class','Re_My','Im_My'});
rows5 = A03.Class == 0;
rows6 = A03.Boundary_id == 2;
lid3 = A03(rows6,:);
load3 = A03(rows5,:);

Re_My = zeros(100,1);

for i=1:100
    Re_My(i) = lid3.Re_My(i)+load3.Re_My(i);
end

Im_My = zeros(100,1);

for i=1:100
    Im_My(i) = lid3.Im_My(i)+load3.Im_My(i);
end     

Abs_My = zeros(100,1);

for i=1:100
    Abs_My(i) = sqrt(Re_My(i)^2 + Im_My(i)^2);
end     

Arg_My = zeros(100,1);

for i=1:100
    Arg_My(i) = atan2(Im_My(i),Re_My(i));
end    

Lshh_abs = importdata('liingaard_fig10_abs_1.dat');
Lshh_arg = importdata('liingaard_fig10_arg_1.dat');
Lsvv_abs = importdata('liingaard_fig5_abs_1.dat');
Lsvv_arg = importdata('liingaard_fig5_arg_1.dat');
Lsmm_abs = importdata('liingaard_fig11_abs_1.dat');
Lsmm_arg = importdata('liingaard_fig11_arg_1.dat');

figure (1);
plot(a0,(Abs_Fx/1e6/5*4)/(Abs_Fx(1)/1e6/5*4),'o-','LineWidth',1.5)
hold on
plot(Lshh_abs(:,1),Lshh_abs(:,2),'LineWidth',1.5)
hold off
xlim([0 10])
xlabel('a_0')
ylabel('|S_{HH}|/K_{HH}')
legend('MultiFEBE','Liingaard et al.')

figure (2);
plot(a0,Arg_Fx,'o-','LineWidth',1.5)
hold on
plot(Lshh_arg(:,1),Lshh_arg(:,2),'LineWidth',1.5)
hold off
xlim([0 10])
ylim([0 pi])
set(gca,'ytick',(0:pi/4:pi))
set(gca,'yticklabels',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
xlabel('a_0')
ylabel('arg(S_{HH})')

figure (3);
plot(a0,(Abs_Fz/1e6/5*4)/(Abs_Fz(1)/1e6/5*4),'o-','LineWidth',1.5)
hold on
plot(Lsvv_abs(:,1),Lsvv_abs(:,2),'LineWidth',1.5)
hold off
xlim([0 10])
xlabel('a_0')
ylabel('|S_{VV}|/K_{VV}')

figure (4);
plot(a0,Arg_Fz,'o-','LineWidth',1.5)
hold on
plot(Lsvv_arg(:,1),Lsvv_arg(:,2),'LineWidth',1.5)
hold off
xlim([0 10])
ylim([0 pi])
set(gca,'ytick',(0:pi/4:pi))
set(gca,'yticklabels',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
xlabel('a_0')
ylabel('arg(S_{VV})')

figure (5);
plot(a0,(Abs_My/1e6/5^3*4)/(Abs_My(1)/1e6/5^3*4),'o-','LineWidth',1.5)
hold on
plot(Lsmm_abs(:,1),Lsmm_abs(:,2),'LineWidth',1.5)
hold off
xlim([0 10])
xlabel('a_0')
ylabel('|S_{MM}|/K_{MM}')

figure (6);
plot(a0,Arg_My,'o-','LineWidth',1.5)
hold on
plot(Lsmm_arg(:,1),Lsmm_arg(:,2),'LineWidth',1.5)
hold off
xlim([0 10])
ylim([0 pi])
set(gca,'ytick',(0:pi/4:pi))
set(gca,'yticklabels',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
xlabel('a_0')
ylabel('arg(S_{MM})')
