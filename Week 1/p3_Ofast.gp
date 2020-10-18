set xlabel "Size of inner matrices products"
set title "Comparison CPU performance"

#set key

plot 'results_p3_Ofast.dat' u 2:4 t "Fortran function" w lp ls 3 lt 3 pt 7 linecolor 7, 'results_p3_Ofast.dat' u 2:5 t "mat(order 1)" w lp ls 2.75 , 'results_p3_Ofast.dat' u 2:6 t "mat(order 2)" w lp ls 1.75 linecolor 12

#set title "Comparison CPU performance"

set terminal postscript 
set output "pr3_Ofast.ps"
rep 
set terminal wxt
