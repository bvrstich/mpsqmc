set output "walkers.svg"
set key right
set pointsize 0.5
set term svg
plot "output/J1J2/4x4/J2=0.3/DT4DW2.txt" u ($1):($2)  w l 
