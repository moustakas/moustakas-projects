; IDL Version 6.0 (linux x86 m32)
; Journal File for ioannis@cerebus.as.arizona.edu
; Working directory: /home/ioannis/catalogs/02garnett
; Date: Tue Sep 27 20:30:10 2005
 
; % No commands containing 'r = rs' exist in the command buffer.
d = rsex('02garnett_table1.dat')
; % No commands containing 'm = ' exist in the command buffer.
ned_webget_diameters, d.galaxy, m
rd = d.rd/d.dist*!radeg*3600.0/1E3
niceprint, m.galaxy, rd, m.rc3_major_axis,  rd/m.rc3_major_axis/2.0
;NGC0224         857.20434        11432.8       0.037488820
;NGC0253         305.27190        1652.50       0.092366686
;NGC0300         206.26479        1312.70       0.078565096
;NGC0598         509.59538        4247.70       0.059984857
;NGC0628         89.680346        628.300       0.071367458
;NGC0925         88.716041        628.300       0.070600066
;NGC1232         60.440382        444.800       0.067941078
;NGC1637         19.187423        238.900       0.040157855
;NGC2403         128.91550        1312.70       0.049103185
;NGC2442         43.424167        329.700       0.065854058
;NGC2805         74.606415        378.600       0.098529337
;NGC2903         62.206843        755.400       0.041174769
;NGC3031         171.88733        1614.90       0.053219186
;NGC3344         60.865021        424.800       0.071639622
;NGC3521         65.890143        657.900       0.050076106
;NGC3621         61.571581        738.200       0.041703860
;NGC4254         40.996729        322.200       0.063620000
;NGC4258         121.00868        1117.30       0.054152275
;NGC4303         52.527060        387.400       0.067794347
;NGC4321         66.619685        444.800       0.074887238
;NGC4395         189.07606        791.000        0.11951711
;NGC5033         72.799339        642.900       0.056617932
;NGC5055         114.59155        755.400       0.075848259
;NGC5194         115.18683        673.200       0.085551717
;NGC5236         160.42817        772.900        0.10378326
;NGC5457         154.69860        1730.40       0.044700241
;NGC6384         91.500924        370.000        0.12364990
;NGC6744         99.165767        1197.20       0.041415708
;NGC6946         101.25726        688.900       0.073491987
;NGC7331         65.567617        628.300       0.052178592
;NGC7793         103.13240        560.000       0.092082498
jj = im_stats(rd/m.rc3_major_axis/2.0,/verbose)
