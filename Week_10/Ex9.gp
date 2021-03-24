set xrange [0:3.1]


unset key
# border
set style line 101 lc rgb '#808080' lt 1
set border 3 back ls 101
set tics nomirror out scale 0.75
# define grid
set style line 102 lc rgb'#808080' lt 0 lw 1
set grid back ls 102

set ylabel "Energy" 
set xlabel "Lambda" 
set title "Ising model"


# color definitions
set border linewidth 1.5
set style line 1 lt rgb "#F08080" lw 3
set style line 2 lt rgb "blue" lw 2
set style line 3 lt rgb "#66CDAA" lw 3
set style line 4 lt rgb "violet" lw 2

set style line 6 lt rgb "#778899" lw 2
set style line 7 lt rgb "violet" lw 3

N = system("cat N.txt")
fname = system("cat in_file.txt")
oname = system("cat file.txt")

set xtics 
set ytics 

set key left bottom 

plot fname u 1:2 title "1st eigenvalue" w lines ls 1, [0:2] (-1-x**2/4) w l lw 2 lc rgb "#000000" title "Mean field result", [2:4] -x w l lw 2 lc rgb "#000000" title ""





set terminal pdfcairo color
set output oname
rep 
set terminal post
