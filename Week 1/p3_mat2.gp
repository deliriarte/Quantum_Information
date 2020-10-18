set xlabel "Size of inner matrices products"
set title "Comparison CPU performance"
#set yrange[0:0.00001]
#set key

plot 'results_p3_Ofast.dat' u 2:6 t "mat(order 2) using Ofast" w lp ls 1.75 linecolor 12, 'results_p3.dat' u 2:6 t "mat(order 2)" w l ls 1.75 linecolor 12

#set title "Comparison CPU performance"

set terminal postscript 
set output "pr3_mat2.ps"
rep 
set terminal wxt
