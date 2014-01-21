set output "walkers.svg"
set key right
set pointsize 0.5
set term svg
plot  "output/Heisenberg1D/ener_L16DT32DW2.txt" u ($1):($2)  w l 
