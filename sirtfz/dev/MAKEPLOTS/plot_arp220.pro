PRO plot_arp220,postscript=postscript,encapsulated=encapsulated
; lam98nov15ucb
; plot of the spectral energy distribution of Arp 220
; Arp220 coordinates (J2000): 15:34:57.3 +23:30:12
; Klaas, U., Haas, M., Heinrichsen, I., \& Schulz, B. 1997, A\&A, 325, L21
; IRAS and ISOPHOT Photometry:

    IF keyword_set(postscript) THEN BEGIN 
      IF keyword_set(encapsulated) THEN BEGIN 
         ps_open,'arp220',/ps_fonts,/encapsulated,/color
         device,xsize=7,ysize=6,/inches,/times
      ENDIF ELSE BEGIN 
         ps_open,'arp220',/ps_fonts,/encapsulated,/color
         device,/times
      ENDELSE 
   ENDIF 

   lam = [3.3,3.6,4.8,7.3,7.7,10.,11.3,12.,12.8,15.,20.,25.,60.,$
          65.,80.,90.,100.,105.,120.,150.,170.,180.,200.]
   fjy = [0.103,0.048,0.053,0.394,0.583,0.127,0.312,0.583,1.54,2.58,3.63,8.28,105.,$
          111.,137.,99.,113.7,114.4,92.5,69.,67.,45.,37.6]

; nb, best fits are with 120K and 47K modified blackbodies (emissivity=1)

    plotsym, 0, 1, /fill
   plot,lam,fjy,xr=[1.,1000.],yr=[0.01,1000.],/xst,/yst,/xlog,/ylog,psym=-8,$
    xtit=textoidl('Wavelength (\mu m)'),ytit=textoidl('Flux F_{\nu} (Jy)')
   xyouts,2.,220,'Arp 220 (Klaas et al. 1997)'
   
   IF keyword_set(postscript) THEN ps_close

END 

   
