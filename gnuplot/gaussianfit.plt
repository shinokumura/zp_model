
set tics font "Helvetica,14"
set xlabel font "Helvetica,14"
set ylabel font "Helvetica,14"
set key font "Helvetica,14"
set title font "Helvetica,14"
set format y "10^{%L}"

sigma = 0.5
mu = 30
amp = 0.001

set log y
# set fit errorvariables

gauss(x)= amp / (sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))
erfs(x) = 0.5 * amp * (erf((x - d + 0.5)/sigma*sqrt(2.)) - erf((x - d - 0.5)/sigma*sqrt(2.)))

# set print "gnuplot_variable.dat"
# set print "gnuplot_variable_yerror.dat"
# set print "gnuplot_variable_erf.dat"
#print "#mass  amp               sigma              mu              ucd_z            deltaZ"

do for [mass = 74.0:109.0]{
    set title 'mass: '.mass
    
    sigma = 0.5
    ucd_z = 92.0/236.0 * mass 
    mu = ucd_z + 0.5
    d = ucd_z + 0.5

    if (mass <= 83) {
        amp = 0.0005
    }
    if (mass >= 84 && mass <= 88) {
        amp = 0.005
    }
    if (mass >= 89  && mass <=100) {
        amp = 0.1
    }
    if (mass >= 101) {
        amp = 0.2
    }
    
    # fit gauss(x) "expU235_average.dat" u ($1 == mass ? $2:NaN):4 via amp, sigma, mu
    # fit gauss(x) "expU235_average.dat" u ($1 == mass ? $2:NaN):4:5 yerror via amp, sigma, mu
    # fit erfs(x) "expU235_average.dat" u ($1 == mass ? $2:NaN):4 via amp, sigma, d
    #fit erfs(x) "expU235_average.dat" u ($1 == mass ? $2:NaN):4:5 yerror via amp, sigma, d
    
    #plot "expU235_average.dat" u ($1 == mass ? $2:1/0):4 ps 2 pt 4, gauss(x) lw 2

    #deltaZ = mu - ucd_z
    #print mass, amp, sigma, mu, ucd_z, deltaZ

    #deltaZ = d - ucd_z
    #print mass, amp, sigma, d, ucd_z, deltaZ
    # pause -1
    
}

set print

###########################
# variables check

set key bottom 
unset log y
unset format y
set yrange [-2:2]
set xlabel "Mass number"
set ylabel "deltaZ"
set xtics 2
set mxtics 2
set grid

set title "deltaZ"
plot 'gnuplot_variable.dat' u 1:6  w l lw 2 dt 2 ti "Gnuplot Gauss no/dy-constrain",\
    'gnuplot_variable_yerror.dat' u 1:6 w l lw 2 dt 3 ti "Gnuplot Gauss w/dy-constrain",\
    'gnuplot_variable_erf.dat' u 1:6  w l lw 2 dt 4 ti "Gnuplot erf",\
    'erf-param' u 1:6  w l lw 2 ti "Python erf",\
    
pause -1

set title "sigma"
set ylabel "sigma"
plot 'gnuplot_variable.dat' u 1:3  w l ti "Gnuplot Gauss - no/dy-constrain",\
    'gnuplot_variable_yerror.dat' u 1:3 w l  ti "Gnuplot Gauss w/dy-constrain",\
    'gnuplot_variable_erf.dat' u 1:3  w l ti "Gnuplot erf",\
    'erf-param' u 1:5  w l ti "Python erf",\
    
pause -1

