set term png

set style line 1  lc rgb '#ff66ff' lt 1 lw 2
set style line 2  lc rgb '#0099ff' lt 1 lw 2 
set style line 3  lc rgb '#ffd633' lt 1 lw 2 
set style line 4  lc rgb '#00ff99' lt 1 lw 2


set xl "t"
set yl "y(t)"

set out "Euler.png"
set title "Euler"
p "Euler0.010000.dat" u 1:2 w lines t "dt = 0.01"   ls 1,\
  "Euler0.100000.dat" u 1:2 w lines t "dt = 0.1"    ls 2,\
  "Euler1.000000.dat" u 1:2 w lines t "dt = 1.0"    ls 3,\
  "Euler0.010000.dat" u 1:3 w lines t "analityczna" ls 4
 
set out "RK2.png"
set title "RK2"
p "RK20.010000.dat" u 1:2 w lines t "dt = 0.01"   ls 1,\
  "RK20.100000.dat" u 1:2 w lines t "dt = 0.1"    ls 2,\
  "RK21.000000.dat" u 1:2 w lines t "dt = 1.0"    ls 3,\
  "RK20.010000.dat" u 1:3 w lines t "analityczna" ls 4
  
set out "RK4.png"
set title "RK4"
p "RK40.010000.dat" u 1:2 w lines t "dt = 0.01"   ls 1,\
  "RK40.100000.dat" u 1:2 w lines t "dt = 0.1"    ls 2,\
  "RK41.000000.dat" u 1:2 w lines t "dt = 1.0"    ls 3,\
  "RK40.010000.dat" u 1:3 w lines t "analityczna" ls 4
  
  
set xl "t"
set yl "y(t)"

set out "Euler_blad.png"
set title "Euler_blad"
p "Euler0.010000.dat" u 1:4 w lines t "dt = 0.01"   ls 1,\
  "Euler0.100000.dat" u 1:4 w lines t "dt = 0.1"    ls 2,\
  "Euler1.000000.dat" u 1:4 w lines t "dt = 1.0"    ls 3,\
  "Euler0.010000.dat" u 1:4 w lines t "analityczna" ls 4
  
set out "RK2_blad.png"
set title "RK2_blad"
p "RK20.010000.dat" u 1:4 w lines t "dt = 0.01"   ls 1,\
  "RK20.100000.dat" u 1:4 w lines t "dt = 0.1"    ls 2,\
  "RK21.000000.dat" u 1:4 w lines t "dt = 1.0"    ls 3,\
  "RK20.010000.dat" u 1:4 w lines t "analityczna" ls 4
  
set out "RK4_blad.png"
set title "RK4_blad"
p "RK40.010000.dat" u 1:4 w lines t "dt = 0.01"   ls 1,\
  "RK40.100000.dat" u 1:4 w lines t "dt = 0.1"    ls 2,\
  "RK41.000000.dat" u 1:4 w lines t "dt = 1.0"    ls 3,\
  "RK40.010000.dat" u 1:4 w lines t "analityczna" ls 4
  
  
  
set title "RLC"
set out "RLC.png"
set yl "RLC"
p "RLC0.500000.dat" u 1:2 w lines t "w_V = 0.5 * w_0" 	ls 1,\
  "RLC0.800000.dat" u 1:2 w lines t "w_V = 0.8 * w_0"	ls 2,\
  "RLC1.000000.dat" u 1:2 w lines t "w_V = 1.0 * w_0"	ls 3,\
  "RLC1.200000.dat" u 1:2 w lines t "w_V = 1.2 * w_0"	ls 4
  
set out "RLC_Q.png"
set yl "RLC1"
p "RLC0.500000.dat" u 1:3 w lines t "w_V = 0.5 * w_0" 	ls 1,\
  "RLC0.800000.dat" u 1:3 w lines t "w_V = 0.8 * w_0"	ls 2,\
  "RLC1.000000.dat" u 1:3 w lines t "w_V = 1.0 * w_0"	ls 3,\
  "RLC1.200000.dat" u 1:3 w lines t "w_V = 1.2 * w_0"	ls 4