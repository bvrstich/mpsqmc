set output "energies.svg"
set key left
set pointsize 0.5
set term svg
set key off
plot "energies.txt" u ($1):($3) w l
