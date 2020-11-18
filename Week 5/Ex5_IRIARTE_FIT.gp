set xlabel "Normalized spacing"
set ylabel 'Counts'
set title "Spectral Eigenvalue density"
set key right
set fit quiet

fname = system("cat fit_file.txt")
oname = system("cat out_file.txt")

stats fname u 1
set xrange[STATS_min:STATS_max]
stats fname u 2 
set yrange[STATS_min:STATS_max]

f(x) = a*(x**alpha)*exp(-b*(x**beta))

a = 1
alpha = 1
b = 1
beta = 1

fit f(x) fname u 1:2 via a, alpha, b, beta

plot fname u 1:2 t fname w p ls 3 pt 7 linecolor 1, f(x) t "Pol Fit" w l ls 3 lt 3 linecolor 6


set print 'param'.fname
print "#f(x)=a+b*x"
print a,a_err
print alpha, alpha_err
print b,b_err
print beta, beta_err

set terminal pdfcairo color
set output oname
rep 
set terminal post
