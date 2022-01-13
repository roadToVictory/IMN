set term png enhanced size 800,500 
set size ratio -1

set xl "x"
set yl "y"

set contour
set cntrparam levels 50
unset surface
set view map
unset key



set title "ψ(x,y), Q = -1000"
set out "psi_Q_-1000.png"
splot 'Q_-1000.000000.dat' u 1:2:3 w lines lt -1 palette lw 2 notitle 

set title "ψ(x,y), Q = -4000"
set out "psi_Q_-4000.png"
splot 'Q_-4000.000000.dat' u 1:2:3 w lines lt -1 palette lw 2 notitle 

set title "ψ(x,y), Q = 4000"
set out "psi_Q_4000.png"
splot 'Q_4000.000000.dat' u 1:2:3 w lines lt -1 palette lw 2 notitle 


# dzeta
set title "ζ(x,y), Q = -1000"
set out "dzeta_Q_-1000.png"
splot 'Q_-1000.000000.dat' u 1:2:4 w lines lt -1 palette lw 2 notitle 

set title "ζ(x,y), Q = -4000"
set out "dzeta_Q_-4000.png"
splot 'Q_-4000.000000.dat' u 1:2:4 w lines lt -1 palette lw 2 notitle 

set title "ζ(x,y), Q = 4000"
set out "dzeta_Q_4000.png"
splot 'Q=4000.000000.dat' u 1:2:4 w lines lt -1 palette lw 2 notitle 





set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')
#u
set view map
unset contour
set surface

set title "u(x,y), Q = -1000"
set out "u(x,y)_Q_-1000.png"
splot 'Q_-1000.000000.dat' u 1:2:5 w pm3d notitle

set title "u(x,y), Q = -4000"
set out "u(x,y)_Q_-4000.png"
splot 'Q_-4000.000000.dat' u 1:2:5 w pm3d notitle

set title "u(x,y), Q = 4000"
set out "u(x,y)_Q_4000.png"
splot 'Q=4000.000000.dat' u 1:2:5 w pm3d notitle



# v
set title "v(x,y), Q = -1000"
set out "v(x,y)_Q_-1000.png"
splot 'Q_-1000.000000.dat' u 1:2:6 w pm3d notitle

set title "v(x,y), Q = -4000"
set out "v(x,y)_Q_-4000.png"
splot 'Q_-4000.000000.dat' u 1:2:6 w pm3d notitle

set title "v(x,y), Q = 4000"
set out "v(x,y)_Q_4000.png"
splot 'Q=4000.000000.dat' u 1:2:6 w pm3d notitle
