set term png
set grid

set style line 1  lc rgb '#ff0000' lw 2
set style line 2  lc rgb '#000000' lw 2 

set xl "czas"
set yl "ilosc osob"

set out "picard.png"
set title "rozprzestrzenianie choroby zakaźnej - metoda Picarda"
plot "picard.dat" u 1:2 w lines t "nosiciele" ls 1,\
  "picard.dat" u 1:3 w lines t "osoby zdrowe" ls 2


set out "Newton.png"
set title "rozprzestrzenianie choroby zakaźnej -metoda Newtona"
plot "Newton.dat" u 1:2 w lines t "nosiciele" ls 1,\
  "Newton.dat" u 1:3 w lines t "osoby zdrowe" ls 2

set out "RK2.png"
set title "rozprzestrzenianie choroby zakaźnej - metoda RK2"
plot "RK2.dat" u 1:2 w lines t "nosiciele" ls 1,\
  "RK2.dat" u 1:3 w lines t "osoby zdrowe" ls 2
