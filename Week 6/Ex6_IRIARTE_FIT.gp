set xlabel "Matrix dimension N"
set ylabel 'Cpu time'
set title "Computational time"
set key left
set fit quiet

fname = "cpu_times.txt"

f(x) = a*(x**3)

fit f(x) fname u 1:2 via a

plot fname t fname w p ls 3 pt 7 linecolor 1, f(x) t "Pol Fit" w l ls 3 lt 3 linecolor 6

set print 'param'.fname
print "#f(x)=a+b*x"
print a,a_err

set terminal pdfcairo color
set output "cpu_times.pdf"
rep 
set terminal post
