pro plot_photoz, oflux, obands, zarray, zphot, ztrue, deltaz, bp, $
        sedindx, constant, chi2, pdz, gal_id, user=user
; jm01may13uofa

    common sirtf_simulations

    if keyword_set(user) then stay = ''
    
; redshift the SED to the photometric redshift

;   wavez = *sirtf.sedcube[sedindx].lambda
;   fluxz = *sirtf.sedcube[sedindx].mlum
    
    redshift_sed, zphot, *sirtf.sedcube[sedindx].lambda, $
      *sirtf.sedcube[sedindx].mlum_nu, wavez, fluxz, /mjansky
        
    w1 = 0.1 > min(wavez) ; starting wavelength (mu)
    w2 = 1E4 < max(wavez) ; ending wavelength (mu)

    get_element, wavez, [w1,w2], x1
    nozero = where(oflux[0,*] ne float(0.0))

; SED plot

    plotsym, 0, 1, /fill
    plot, [wavez[x1[0]:x1[1]],wavez[x1[0]:x1[1]]], [min(constant*fluxz[x1[0]:x1[1]]),max(constant*fluxz[x1[0]:x1[1]])], $
      /nodata, xsty=3, ysty=3, color=16, position=[0.22,0.52,0.95,0.93], /xlog, /ylog, xthick=2.0, $
      ythick=2.0, charsize=1.8, charthick=2.0, thick=2.0, $
      xtit=textoidl('\lambda')+' ('+textoidl('\mu')+'m)', ytit=textoidl('f_{\nu}')+' (mJy)', tit='Galaxy '+strn(gal_id)
    oplot, wavez, constant*fluxz, color=3, thick=2.5
    oploterror, sirtf.bandcube[obands].lambda0, oflux[0,nozero], oflux[1,nozero], ps=8, color=4, $
      errcolor=4, thick=2.5
          
; response function plot
           
    plot, [min(zarray),max(zarray)], [min(bp[sedindx,*]),max(bp[sedindx,*])*1.4], /nodata, /noerase, $
      xsty=3, ysty=3, color=16, position=[0.22,0.13,0.95,0.4], xtit='z', ytit='p(z)', xthick=2.0, $
      ythick=2.0, charsize=1.8, charthick=2.0, thick=2.0
    oplot, zarray, bp[sedindx,*], line=0, color=3, thick=2.5
    oplot, [zphot-deltaz,zphot-deltaz], [!y.crange[0],!y.crange[1]], line=2, color=4, thick=2.5
    oplot, [zphot+deltaz,zphot+deltaz], [!y.crange[0],!y.crange[1]], line=2, color=4, thick=2.5
    plotsym, 1, 3, thick=2.5
    plots, [zphot,zphot], 1.3*[max(bp[sedindx,*]),max(bp[sedindx,*])], ps=8, color=6
    
; legend

    text = [textoidl('\chi^{2}')+': '+strn(chi2,format='(F7.2)'), $
            textoidl('p_{\Delta}')+textoidl('_{z}')+' : '+strn(pdz,format='(F6.4)'),$
            textoidl('z_{phot}')+': '+strn(zphot,format='(F4.2)'),$
            textoidl('z_{spec}')+': '+strn(ztrue,format='(F4.2)')]
    if zphot lt max(zarray)/2. then legend, text, box=0, /right, /top, charsize=1.3, charthick=2.0 $
    else legend, text, box=0, /left, /top, charsize=1.3, charthick=2.0
    
    if keyword_set(user) then read, stay
    
return
end
