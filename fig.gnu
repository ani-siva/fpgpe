set term pdf 
set output "density.pdf"
set notitle
set tics out
unset key
set pm3d map 
set cbtics offset -1.2,0 tc  lt 1
set xtics 3 offset 0,.75  tc lt 1    
set ytics 3 offset 1,0     tc  lt 1  
set mxtics 3    
set mytics 3   
set xrange[-10:10]
set yrange[-10:10] 
set palette defined (  1 "black", 2 "red", 3 "magenta",4 "green", 5 "yellow", 6 "pink",7 "white" ) 
set size square
set xlabel "{/Times-Italic=20 x}" offset 0,1.2  tc  lt 1 
set ylabel "{/Times-Italic=20 y}" offset 1.2,0  tc  lt 1
set colorbox

#splot 'density/den-1.txt'  us ($2):($1):($3)
splot 'density/initial_density.txt'  us ($2):($1):($3)










