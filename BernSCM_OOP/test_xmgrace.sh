#!/bin/sh
xmgrace -printfile ./plots/c4mip_a2_D10QI_BernSCM_t_f_CS25.eps -hdevice EPS -free -noask \
 -pexec 'default font 0' \
-pexec 'title " C4MIP SRES A2";; LEGEND .15000000000000000000, .85000000000000000000 ; LEGEND CHAR SIZE 1 ; yaxis label "Atmospheric CO2 (ppm)" ; yaxis label char size 1.3; xaxis label "Time (yr)"; xaxis label char size 1.3' \
-block ./output/c4mip_a2_D10QI_BernSCM_t_f_CS25.dat -bxy 1:8 \
-pexec 's8 LEGEND "" ; legend box pattern 0 ; legend box fill pattern 0 ; s8 linewidth 2; s8 linestyle 5; s8 color 1' \
-pexec 'with g0; xaxis label ""' \
-pexec 'world 1850, 0, 2100, 1 ; autoscale yaxes; autoticks' \
-pexec 'legend loctype world; legend 1870, 880' \
-graph 1 \
-pexec 'view g0.vx1, g0.vy1, g0.vx2, g0.vy2; world g0.wx1, g0.wy1, g0.wx2, g0.wy2; autoscale onread none; autoticks' \
-pexec 'title " ";; LEGEND .12500000000000000000, .75000000000000000000 ; LEGEND CHAR SIZE 1 ; yaxis label "Global SAT change (K)" ; yaxis label char size 1.3; xaxis label "Time (yr)"; xaxis label char size 1.3' \
-block ./output/c4mip_a2_D10QI_BernSCM_t_f_CS25.dat -bxy 1:2 \
-pexec 's8 LEGEND "" ; legend box pattern 0 ; legend box fill pattern 0 ; s8 linewidth 2; s8 linestyle 5; s8 color 1' \
-pexec 'world 1850, 0, 2100, 1 ; autoscale yaxes; autoticks' \
-graph 2 \
-pexec 'view g0.vx1, g0.vy1, g0.vx2, g0.vy2; world g0.wx1, g0.wy1, g0.wx2, g0.wy2; autoscale onread none; autoticks' \
-pexec 'title " ";; LEGEND .10000000000000000000, .65000000000000000000 ; LEGEND CHAR SIZE 1 ; yaxis label "Land C uptake (GtC/yr)" ; yaxis label char size 1.3; xaxis label "Time (yr)"; xaxis label char size 1.3' \
-block ./output/c4mip_a2_D10QI_BernSCM_t_f_CS25.dat -bxy 1:14 \
-pexec 's8 LEGEND "" ; legend box pattern 0 ; legend box fill pattern 0 ; s8 linewidth 2; s8 linestyle 5; s8 color 1' \
-pexec 'world 1850, 0, 2100, 1 ; autoscale yaxes; autoticks' \
-graph 3 \
-pexec 'view g0.vx1, g0.vy1, g0.vx2, g0.vy2; world g0.wx1, g0.wy1, g0.wx2, g0.wy2; autoscale onread none; autoticks' \
-pexec 'title " ";; LEGEND .07500000000000000000, .55000000000000000000 ; LEGEND CHAR SIZE 1 ; yaxis label "Ocean C uptake (GtC/yr)" ; yaxis label char size 1.3; xaxis label "Time (yr)"; xaxis label char size 1.3' \
-block ./output/c4mip_a2_D10QI_BernSCM_t_f_CS25.dat -bxy 1:13 \
-pexec 's8 LEGEND "" ; legend box pattern 0 ; legend box fill pattern 0 ; s8 linewidth 2; s8 linestyle 5; s8 color 1' \
-pexec 'world 1850, 0, 2100, 1 ; autoscale yaxes; autoticks' \
-pexec 'ARRANGE(2, 2, 0.1, 0.35, 0.2)' \