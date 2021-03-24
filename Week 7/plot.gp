set xlabel "Time"
set ylabel 'Position'
set title "Average position"
set key left

set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror out scale 0.75
set palette rgb 33,28,10
set palette negative
#set yrange[-0.05:0.6]
#set xrange[-5:5]
fname = "results_mean.txt"

Shadecolor = "#80E0A080"


plot fname u  1:2 w l lw 3 t ""



set terminal pdfcairo color
set output "average_position.pdf"
rep 
set terminal post
