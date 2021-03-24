set xlabel "L"
set ylabel '|Ψ(s)|²'
set title "|Ψ(s)|²"
set key left

set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror out scale 0.75
set palette rgb 33,28,10
set palette negative
set yrange[-0.05:0.6]
set xrange[-5:5]
fname = "results_square.txt"

plot for [i=1:25] fname using 1:(column(5*i+1)) w l lw 3 lt palette frac 1.0*(4*i)/100.0 title ""

#plot fname u  1:4 w l lw 3 


set terminal pdfcairo color
set output "psi_square.pdf"
rep 
set terminal post
