

reset

L=3.0;
rho=1.25;
c=343.0;

set terminal pdfcairo color enhanced size 15cm,10cm
set output "Unx0_vs_f.pdf"
set log y
set samples 1000
set xlabel "f [Hz]"
set ylabel "U_n at x=0 m [m]"
plot "<awk '{if ($9==1084) print $0}' room.dat.nso" u 2:(abs($15)) w p t "Numerical solution", \
     abs(1.0/(1.25*2*3.14159265*x*343.0*sin(2*3.14159265*x/343.0*3))) w l t "Analytical solution"
set output
unset log y

omega=2*3.14159265*55.9;

set terminal pdfcairo color enhanced size 15cm,15cm
set output "Un_p_vs_x_fn1.pdf"
set multiplot layout 2,1
set title "f = 55.9 Hz"
set xlabel "x [m]"
set ylabel "U_x [m]"
plot "<awk '{if ($1==10  && $16<0.01 && $17<1.51 && $17>1.49) print $0}' room.dat.eso | LC_ALL=C sort -k 10 -g" u 15:18 w p t "Numerical solution (element nodes)",\
     cos(omega/c*x)/(sin(omega/c*L)*rho*omega*c) w l t "Analytical solution"
set xlabel "x [m]"
set ylabel "p [Pa]"
plot "<awk '{if ($1==10 && $11<0.01 && $12<1.51 && $12>1.49) print $0}' room.dat.nso | LC_ALL=C sort -k 10 -g" u 10:13 w lp t "Numerical solution",\
     sin(omega/c*x)/(sin(omega/c*L)) w l t "Analytical solution"
unset multiplot
set output

omega=2*3.14159265*171.85714286;

set terminal pdfcairo color enhanced size 15cm,15cm
set output "Un_p_vs_x_fn3.pdf"
set multiplot layout 2,1
set title "f = 55.9 Hz"
set xlabel "x [m]"
set ylabel "U_x [m]"
plot "<awk '{if ($1==29  && $16<0.01 && $17<1.51 && $17>1.49) print $0}' room.dat.eso | LC_ALL=C sort -k 10 -g" u 15:18 w p t "Numerical solution (element nodes)",\
     cos(omega/c*x)/(sin(omega/c*L)*rho*omega*c) w l t "Analytical solution"
set xlabel "x [m]"
set ylabel "p [Pa]"
plot "<awk '{if ($1==29 && $11<0.01 && $12<1.51 && $12>1.49) print $0}' room.dat.nso | LC_ALL=C sort -k 10 -g" u 10:13 w lp t "Numerical solution",\
     sin(omega/c*x)/(sin(omega/c*L)) w l t "Analytical solution"
unset multiplot
set output
