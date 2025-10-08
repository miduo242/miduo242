
set terminal pngcairo enhanced size 1000, 800 font "Helvetica,15"


set style line  1 pt 5 lc rgb '#8A2BE2' lw 6 lt 1 # blueviolet
set style line  2 pt 5 lc rgb '#FFC0CB' lw 6 lt 1 # pink
set style line  3 pt 5 lc rgb '#0000FF' lw 6 lt 1 # blue
set style line  4 pt 5 lc rgb '#A52A2A' lw 6 lt 1 # brown
set style line  5 pt 5 lc rgb '#7FFF00' lw 6 lt 1 # chartreuse
set style line  6 pt 5 lc rgb '#00FFFF' lw 6 lt 1 # cyan
set style line  7 pt 5 lc rgb '#006400' lw 6 lt 1 # darkgreen
set style line  8 pt 5 lc rgb '#FF00FF' lw 6 lt 1 # magenta
set style line  9 pt 5 lc rgb '#fde725' lw 6 lt 1 # yellow
set style line  10 pt 5 lc rgb '#cb4679' lw 6 lt 1 # magenta
set style line  11 pt 5 lc rgb '#FF0000' lw 6 lt 1 # red
set style line  12 pt 5 lc rgb '#FF8C00' lw 6 lt 1 # orange
set style line 13 pt 5 lc rgb '#008080' lw 6 lt 1 # teal
set style line 14 pt 5 lc rgb '#4dbeee' lw 6 lt 1 # light-blue

#1         Temperature
set size 1.0, 1.0
set origin 0.0,0.0
set key samplen 2 top left
set output "Temperature.png"
set xlabel "Time [s]"
set ylabel "Temperature [K]"
set title "Martensite (Temperature)"
plot [:] [:] "MartensiteVolumeFraction.dat" u ($1*1e-3):($2) w l ls 1 title "Heat Release"



#2         Temp_VonMises
set size 1.0, 1.0
set origin 0.0,0.0
set key samplen 2 top left
set output "VonMises-Temp.png"
set xlabel "Temperature [K]"
set ylabel "VonMises [MPa]"
set title "Martensite (VonMises-Temp)"
plot [:] [:] "TextData/TemperatureFile.txt" u ($2):($3/1e6) w l ls 1 title "Vm-Temp"


#3         Volume Fraction
set size 1.0, 1.0
set origin 0.0,0.0
set key samplen 2 top left
set output "VolumeFraction.png"
set xlabel "Temperature [K]"
set ylabel "Volume Fraction %"
set title "Martensite (Volume Fraction)"
plot [:] [:] "MartensiteVolumeFraction.dat" u ($2):($3*100) w l ls 1 title "{/Symbol e}_{VF}"


#2         VonMises
set size 1.0, 1.0
set origin 0.0,0.0
set key samplen 2 top left
set output "VonMises.png"
set xlabel "Time [s]"
set ylabel "VonMises [MPa]"
set title "Martensite (VonMises)"
plot [:] [:] "TextData/TemperatureFile.txt" u ($1*1e-3):($3/1e6) w l ls 1 title "Vm-Stress"
