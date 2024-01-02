#set term postscript eps enhanced color solid

#set log y

set tics font "Helvetica,14"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,14"
set key font "Helvetica,14"
set title font "Helvetica,14"

set dgrid3d 100,100
set hidden3d
splot 'gauss' u 1:5:6 w l



#set yrange [1E-5:0.1]
#set xrange [36:43]


#set output 'A99.eps'
set title 'A=99'
plot\
 'b' i 0 u 3:4 w lp ti 'Wahl' dt 2 lw 2,\
 ''  i 0 u 3:5 w lp ti 'UCD' dt 3 lw 2,\
 ''  i 0 u 3:6 w lp ti 'Minato' dt 4 lw 2,\
 ''  i 0 u 3:7 w lp ti 'Titech' dt 5 lw 2,\
 'hf3d' i 0 u 3:7 w lp ti 'HF3D' lw 2

pause -1

#set output 'A100.eps'
set title 'A=100'
plot\
 'b' i 1 u 3:4 w lp ti 'Wahl' dt 2 lw 2,\
 ''  i 1 u 3:5 w lp ti 'UCD' dt 3 lw 2,\
 ''  i 1 u 3:6 w lp ti 'Minato' dt 4 lw 2,\
 ''  i 1 u 3:7 w lp ti 'Titech' dt 5 lw 2,\
 'hf3d' i 1 u 3:7 w lp ti 'HF3D'

pause -1

#set xrange [49:56]
#set output 'A135.eps'
set title 'A=135'
plot\
 'a' i 0 u 3:4 w lp ti 'Wahl' dt 2 lw 2,\
 ''  i 0 u 3:5 w lp ti 'UCD' dt 3 lw 2,\
 ''  i 0 u 3:6 w lp ti 'Minato' dt 4 lw 2,\
 ''  i 0 u 3:7 w lp ti 'Titech' dt 5 lw 2,\
 'hf3d' i 2 u 5:7 w lp ti 'HF3D'

pause -1
#set output 'A136.eps'
set title 'A=136'
plot\
 'a' i 1 u 3:4 w lp ti 'Wahl' dt 2 lw 2,\
 ''  i 1 u 3:5 w lp ti 'UCD' dt 3 lw 2,\
 ''  i 1 u 3:6 w lp ti 'Minato' dt 4 lw 2,\
 ''  i 1 u 3:7 w lp ti 'Titech' dt 5 lw 2,\
 'hf3d' i 3 u 5:7 w lp ti 'HF3D'

pause -1