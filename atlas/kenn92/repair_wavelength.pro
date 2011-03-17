pro repair_wavelength, write=write, doplot=doplot
; jm04apr28uofa
; the wavelength solutions of some (most) of the kennicutt spectra are
; poor

    light = 2.99792458D5        ; speed of light [km/s]
    fwhm2sig = 2.0*sqrt(2.0*alog(2.0))
        
    datapath = atlas_path(/kenn92)+'analysis/'
    
    kenn = read_kenn92(/silent)
    galaxy = kenn.galaxy
    ngalaxy = n_elements(galaxy)

    lpars = read_linepars(linefile='repair_elinelist.dat',linepath=datapath)
    nline = n_elements(lpars)
    
    lineres = replicate(7.0,nline)/(1+kenn.z_obj)
    elineres = light * lineres / lpars.wave / fwhm2sig ; rest-frame sigma width [km/s]

    result = {$
      galaxy: '', $
      line:   strarr(nline), $
      lwave:  fltarr(nline), $
      wshift: fltarr(nline)}
    result = replicate(result,ngalaxy)
    result.galaxy = galaxy
    result.line = lpars.line
    result.lwave = lpars.wave

    pagemaker, nx=3, ny=1, position=pos, xmargin=0.2, xspace=0, /normal

    if keyword_set(doplot) then begin
       window, 0, xs=600, ys=400
       window, 2, xs=600, ys=400
    endif
    
    plotsym, 8, 1.5
    for i = 4L, ngalaxy-1L do begin
;   for i = 0L, ngalaxy-1L do begin

       spec = read_kenn92_specfit(galaxy[i])
       ferr = mrdfits(atlas_path(/kenn92)+'data/'+kenn[i].specfile,1,/silent)
       wave = spec[*,0]
       flux = spec[*,1]
       continuum = spec[*,2]

       eflux = flux - continuum
       norm = max(eflux)
       eflux = eflux / norm
       ferr = ferr / norm
       invvar = 1.0/ferr^2.0

; fit the emission lines without redshift contraints

       fit = ilinefit(eflux,wave,lpars.wave,elineres,zguess=0.0,$
         invvar=invvar,specfit=specfit,windex=lpars.windex,$
         findex=lpars.findex,fvalue=lpars.fvalue,zindex=lpars.zindex,zmaxshift=1000.0)

       wshift = lpars.wave-lpars.wave*(1+fit.linez)
       result[i].wshift = wshift
       
       if keyword_set(doplot) then begin

          wset, 0
          djs_plot, wave, eflux, ps=10, xsty=3, ysty=3, charthick=2.0, $
            charsize=1.0, xthick=2.0, ythick=2.0, ytickname=replicate(' ',10), $
            position=pos[*,0], xrange=[min(wave),4400] ; xrange=[min(wave),min(wave)+1200]
          djs_oplot, wave, specfit, ps=10, color='yellow'
          djs_plot, wave, eflux, ps=10, xsty=3, ysty=3, charthick=2.0, $
            charsize=1.0, xthick=2.0, ythick=2.0, ytickname=replicate(' ',10), $
            position=pos[*,1], xrange=[4750,5100], $ ; xrange=[!x.crange[1],!x.crange[1]+1200], $
            /noerase, title=galaxy[i]
          djs_oplot, wave, specfit, ps=10, color='yellow'
          djs_plot, wave, eflux, ps=10, xsty=3, ysty=3, charthick=2.0, $
            charsize=1.0, xthick=2.0, ythick=2.0, ytickname=replicate(' ',10), $
            position=pos[*,2], xrange=[6450,6850], $ ; xrange=[!x.crange[1],!x.crange[1]+1200], $
            /noerase
          djs_oplot, wave, specfit, ps=10, color='yellow'

          wset, 2
          djs_plot, result.lwave, wshift, xsty=3, ysty=3, charthick=2.0, $
            charsize=2.0, xthick=2.0, ythick=2.0, ps=8, yrange=max(abs(wshift))*[-1,1], $
            title=galaxy[i], xtitle='Wavelength', ytitle='Shift'
          djs_oplot, !x.crange, [0,0], line=0, thick=2.0
          cc = get_kbrd(1)
          
       endif

    endfor

return
end

