reset
set terminal pdfcairo color enhanced size 15cm,10cm
set output "solution1.pdf"
set log y
set samples 1000
set xlabel "f [Hz]"
set ylabel "U_n at x=0 m [m]"
plot "<awk '{if ($9==1084) print $0}' room.dat.nso" u 2:(abs($15)) w p t "Numerical solution", \
     abs(1.0/(1.25*2*3.14159265*x*343.0*sin(2*3.14159265*x/343.0*3))) w l t "Analytical solution"
set output
