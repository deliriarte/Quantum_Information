set xlabel "Size of inner matrices products"
set title "Comparison CPU performance"
#set yrange[0:0.00001]
#set key

plot 'results_p3_Ofast.dat' u 2:5 t "mat(order 1) using Ofast" w lp ls 2.75, 'results_p3.dat' u 2:5 t "mat(order 1)" w l ls 2.75

#set title "Comparison CPU performance"

set terminal postscript 
set output "pr3_mat1.ps"
rep 
set terminal wxt
