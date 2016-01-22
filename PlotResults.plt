set size 1,1
set origin 0,0

set multiplot layout 2,2

set xlabel "Position"
set xrange [0:1]

set ylabel "Density"
plot "Results/test4.out" using 1:2 with points notitle,\
    "../HyperbolicPDEs/euler/Results/Exact4.out" using 1:2 with lines notitle

set ylabel "Velocity"
plot "Results/test4.out" using 1:3 with points notitle,\
    "../HyperbolicPDEs/euler/Results/Exact4.out" using 1:3 with lines notitle


set ylabel "Pressure"
plot "Results/test4.out" using 1:4 with points notitle,\
    "../HyperbolicPDEs/euler/Results/Exact4.out" using 1:4 with lines notitle


set ylabel "Internal Energy"
plot "Results/test4.out" using 1:5 with points notitle,\
    "../HyperbolicPDEs/euler/Results/Exact4.out" using 1:5 with lines notitle


unset multiplot

pause(-1)
