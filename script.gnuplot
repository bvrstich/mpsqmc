set output "gauss.svg"
set key left
set pointsize 0.5
set term svg
set key off
plot "test" u ($1):($2) 
