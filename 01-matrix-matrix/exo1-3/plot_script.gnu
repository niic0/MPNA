set terminal pngcairo size 800,600
set output "performance_plot.png"
set title "Temps d\'exécution en fonction de la taille des matrices"
set xlabel "Taille des matrices (N)"
set ylabel "Temps d\'exécution (s)"
set grid
plot 'data.txt' using 1:2 with linespoints title 'Performance'
