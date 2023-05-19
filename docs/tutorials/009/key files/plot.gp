# This is a GNUPlot script file.

set style line 1 dt 1 lw 6 lc rgb "blue"
set style line 2 dt 1 lw 6 lc rgb "red"

set terminal postscript eps  color enhanced font "Times,12" size 16.0cm,8.0cm fontscale 2
set output "utop_shell_hom180.eps"


set log y
set xrange [0:10]
set yrange [0.01:100]

set xtics format "%g" autofreq 1
set ytics format "%g"

set ylabel "|{/Times-Italic u}@_{/Times-Italic x}^{/Times (RNA)}|/|{/Times-Italic u}@_{/Times-Italic x}^{S}(0)|} [-]" offset 0
set xlabel "{/Times-Italic f} [Hz]" offset 0,0

set key opaque box on top Right samplen 1 width -3 opaque box on

plot "<awk '{if($9==50)print $2,$13,$14}' rigid/case.dat.nso"    u 1:(sqrt($2**2+$3**2)) w l ls 2 t "Rigid base",\
     "<awk '{if($9==50)print $2,$13,$14}' flexible/case.dat.nso" u 1:(sqrt($2**2+$3**2)) w l ls 1 t "Flexible base"
     
