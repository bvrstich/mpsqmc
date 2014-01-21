set output "energies.svg"
set key right
set pointsize 0.5
set term svg
plot  "output/Heisenberg1D/ener_L16DT4DW2.txt" u ($1):($3)  w l 
