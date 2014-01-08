set output "energies.svg"
set key left
set pointsize 0.5
set term svg
set key off
plot "output/Heisenberg1D/ener_L16DT4DW4.txt" u ($1):($3)  w l
