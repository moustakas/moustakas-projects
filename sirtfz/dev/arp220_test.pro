pro arp220_test, z=z
; jm01jan17uofa
; test that the flux values measured by my routines in various bands
; match the published flux values for arp 220.

    common cosmology, cosmo_params
    
    sedcube = read_sed_cube()
    arp = sedcube[9]

    cosmo_params = {name: 'Cosmology Parameters', omega_0: 0.3D, omega_lambda: 0.7D, h_100: 0.75}

; plot the SED at arp 220's redshift (0.018126)

    if not keyword_set(z) then z = 0.018126D
    redshift_sed, z, arp.lambda, arp.mlum, wavez, fluxz, /jansky

    lam = [3.3,3.6,4.8,7.3,7.7,10.,11.3,12.,12.8,15.,20.,25.,60.,$
           65.,80.,90.,100.,105.,120.,150.,170.,180.,200.]
    fjy = [0.103,0.048,0.053,0.394,0.583,0.127,0.312,0.583,1.54,2.58,3.63,8.28,105.,$
           111.,137.,99.,113.7,114.4,92.5,69.,67.,45.,37.6]

    plotsym, 0, 1, /fill
    colortable2
    plot, alog10(wavez), alog10(fluxz*1E3), xr=[0,5], yr=[-4,3], xsty=3, ysty=3, color=5
    oplot, alog10(lam), alog10(fjy), ps=-8, color=4

stop
    
return
end
