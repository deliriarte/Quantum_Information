set xlabel "Log(N)"
set ylabel 'Log(Cpu time)'
set title "Fitting"
set key left

fname = system("cat fit_file.txt")
oname = system("cat out_file.txt")

f(x) = a + b * x 

fit f(x) fname u (log($1)):(log($2)) via a, b

plot fname u (log($1)):(log($2)) t fname w lp ls 3 lt 3 pt 7 linecolor 1, f(x) t "Linear Fit" w l ls 3 lt 3 linecolor 6


set print 'logparameters'.fname
print "#f(x)=a+b*x"
print a,a_err
print b,b_err

set terminal pdfcairo color
set output oname
rep 
set terminal post
