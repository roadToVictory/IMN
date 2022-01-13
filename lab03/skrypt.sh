set term png
set grid

set style line 1  lc rgb '#ff0000' lt 1 lw 2
set style line 2  lc rgb '#00ff00' lt 1 lw 2 


######
set out "RK2_dt.png"
set title "dt(t)"
p "rTOL5.dat" u 1:2 w lines t "{/:Italic TOL} = 10^{-5}"   ls 1,\
  "rTOL2.dat" u 1:2 w lines t "{/:Italic TOL} = 10^{-2}"   ls 2

set out "RK2_xt.png"
set title "x(t)"
p "rTOL5.dat" u 1:3 w lines t "{/:Italic TOL} = 10^{-5}"   ls 1,\
  "rTOL2.dat" u 1:3 w lines t "{/:Italic TOL} = 10^{-2}"   ls 2


set out "RK2_vt.png"
set title "v(t)"
p "rTOL5.dat" u 1:4 w lines t "{/:Italic TOL} = 10^{-5}"   ls 1,\
  "rTOL2.dat" u 1:4 w lines t "{/:Italic TOL} = 10^{-2}"   ls 2

  
set out "RK2_vx.png"
set title "v(x)"
p "rTOL5.dat" u 3:4 w lines t "{/:Italic TOL} = 10^{-5}"   ls 1,\
  "rTOL2.dat" u 3:4 w lines t "{/:Italic TOL} = 10^{-2}"   ls 2

######
set out "trapez_xt.png"
set title "x(t)"
p "tTOL5.dat" u 1:3 w lines t "{/:Italic TOL} = 10^{-5}"   ls 1,\
  "tTOL2.dat" u 1:3 w lines t "{/:Italic TOL} = 10^{-2}"   ls 2

set out "trapez_vt.png"
set title "v(t)"
p "tTOL5.dat" u 1:4 w lines t "{/:Italic TOL} = 10^{-5}"   ls 1,\
  "tTOL2.dat" u 1:4 w lines t "{/:Italic TOL} = 10^{-2}"   ls 2

set out "trapez_dt.png"
set title "dt(t)"
p "tTOL5.dat" u 1:2 w lines t "{/:Italic TOL} = 10^{-5}"   ls 1,\
  "tTOL2.dat" u 1:2 w lines t "{/:Italic TOL} = 10^{-2}"   ls 2
  
set out "trapez_vx.png"
set title "v(x)"
p "tTOL5.dat" u 3:4 w lines t "{/:Italic TOL} = 10^{-5}"   ls 1,\
  "tTOL2.dat" u 3:4 w lines t "{/:Italic TOL} = 10^{-2}"   ls 2