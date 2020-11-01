set xlabel "N"
set ylabel 'Cpu time'
set title "Fitting"
set key left

fname = system("cat fit_file.txt")
oname = system("cat out_file.txt")

f(x) = a + b * x  + c * x**2 + d * x**3

fit f(x) fname u 1:2 via a, b, c, d

plot fname u 1:2 t fname w lp ls 3 lt 3 pt 7 linecolor 1, f(x) t "Fit" w l ls 3 lt 3 linecolor 6

set print 'parameters'.fname
print "#f(x)=a+b*x+c*x**2+d*x**3"
print a,a_err
print b,b_err
print c,c_err
print d,d_err

set terminal pdfcairo color
set output oname
rep 
set terminal post
