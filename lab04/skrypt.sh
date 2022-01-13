set term png size


set style line 1 lc rgb '#ff0000' lt 1 lw 3
set style line 2 lc rgb '#00ff00' lt 1 lw 3 
set style line 3 lc rgb '#0000ff' lt 1 lw 3
set style line 4 lc rgb '#000000' lt 1 lw 3 


set out "Global_S(it).png"
set title "Globalna S(it)"
set xl "it"
set yl "S"
set logscale x
p "Global_S_06.dat" u 1:2 w lines t "ω_G = 0.6"   ls 2,\
  "Global_S_1.dat" u 1:2 w lines t "ω_G = 1.0"   ls 1

set out "Local_S(it).png"
set title "Lokalna S(it)"
set xl "it"
set yl "S"
p "Local_w_1.dat" u 1:2 w lines t "ω_L = 1.0"   ls 1,\
  "Local_w_14.dat" u 1:2 w lines t "ω_L = 1.4"   ls 3,\
  "Local_w_18.dat" u 1:2 w lines t "ω_L = 1.8"   ls 2,\
  "Local_w_19.dat" u 1:2 w lines t "ω_L = 1.9"   ls 4

unset logscale x

set pm3d map
set xl "x"
set yl "y"


set xrange [0:150]
set yrange [0:100]
set out "Global_V(x,y)_06.png"
set title "Globalna - V(x, y), ω_G = 0.6"
splot "Global_V_06.dat" u 1:2:3 
set out "Global_V(x,y)_1.png"
set title "Globalna - V(x, y), ω_G = 1.0"
splot "Global_V_1.dat" u 1:2:3 


set xrange [0:15]
set yrange [0:10]
set out "Global_Er_06.png"
set title "Globalna, blad dla ω_G = 0.6"
splot "Global_Er_06.dat" u 1:2:3 
set out "Global_Er_1.png"
set title "Globalna, blad dla ω_G = 1.0"
splot "Global_Er_1.dat" u 1:2:3 
