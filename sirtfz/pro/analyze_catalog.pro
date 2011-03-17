pro plot_photoz, zinfo, z, zarray, obands, flux, ferror, color

    waitforme = ' '
    
    common sirtf_simulations
    
    wave0 = *sirtf.sedcube[zinfo.tindx].lambda
    lum0 = *sirtf.sedcube[zinfo.tindx].mlum_nu
    redshift_sed, zinfo.zphot_mode, wave0, lum0, wavez, fluxz, /jansky
    
; SED plot

    plot, wavez, fluxz, xrange=1E-4*[min(sirtf.bandcube[obands].lambda_eff)/3,$
;   plot, wavez, zinfo.constant*lum0, xrange=1E-4*[min(sirtf.bandcube[obands].lambda_eff)/3,$
                                     3*max(sirtf.bandcube[obands].lambda_eff)], $
      xsty=3, ysty=3, /ylog, position=[0.22,0.52,0.95,0.93], /xlog, xthick=2.0, $
      ythick=2.0, charsize=1.8, charthick=2.0, thick=2.0, $
      xtit=textoidl('\lambda')+' ('+textoidl('\mu')+'m)', ytit=textoidl('f_{\nu}')+' [Jy]'
    oploterror, sirtf.bandcube[obands].lambda_eff*1E-4, flux, $
      sirtf.bandcube[obands].dlambda*1E-4, ferror, ps=8, color=color.green, $
      errcolor=color.green, thick=2.5
    
; response function plot
           
    plot, zarray, zinfo.zpdf, /noerase, $
      xsty=3, ysty=3, position=[0.22,0.13,0.95,0.4], xtit='Redshift', $
      ytit='Likelihood', xthick=2.0, ythick=2.0, charsize=1.8, charthick=2.0, $
      thick=2.0
    
; legend

    text = [textoidl('\chi^{2}_{\nu}: '+strn(zinfo.chi2_nu,length=6)), $
            textoidl('z_{phot}'+': '+strn(zinfo.zphot_mode,length=5)), $
            textoidl('z_{spec}'+': '+strn(z,length=5))]
    if zinfo.zphot_mode lt max(zarray)/2. then $
      legend, text, box=0, /right, /top, charsize=1.3, charthick=2.0 else $
      legend, text, box=0, /left, /top, charsize=1.3, charthick=2.0

    read, waitforme
    
return
end

pro analyze_catalog, catalog, interactive=interactive
;+
; NAME:
;	ANALYZE_CATALOG
;
; PURPOSE:
;	Analyze the photometric redshift output from ZGALAXY.
;
; CALLING SEQUENCE:
;	catalog  - root name of the photometric catalog [the assumed
;                  extension is .dat for the ZGALAXY output]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;	
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; COMMON BLOCKS:
;	sirtf_simulations
;
; COMMENTS:
;
; PROCEDURES USED:
;	READFAST
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 September, 14, U of A
;-

    common sirtf_simulations
    
    if n_params() eq 0L then begin
       print, 'Syntax - analyze_catalog, catalog'
       return
    endif

    rpath = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='results')
    ppath = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='plots')

; read the photometric redshift results

    cmrestore, rpath+catalog+'.idlsave', pzs, zarray
    readcol, rpath+catalog+'.dat', source, snr, z_mode, z_mean, z_median, sigma_minus, $
      sigma_plus, chi2_nu, constant, type, zspec, format='L,L,F,F,F,F,F,F,F,L,F', /silent

    good = where((zspec gt float(0)) and (finite(z_mode) eq 1B) and (snr ge 2L),ngood)
    if ngood eq 0L then message, 'No spectroscopic redshift information!'
    
    zspec = zspec[good]
    zphot = z_mode[good]
;   zerr = sigma_minus[good]

    read_catalog, catalog+'.cat', oflux, filters, nsources, catinfo
    obands = filter_match(filters,sirtf.bandcube.bandnames)

    plotsym, 0, 1, /fill
    colortable2, color
    plotfaves
    
    if keyword_set(interactive) then for i = 0L, ngood-1L do begin
       
       flux = reform(oflux[0,*,good[i]])
       ferror = reform(oflux[1,*,good[i]])
       
       plot_photoz, pzs[good[i]], catinfo.zspec[good[i]], zarray, obands, flux, ferror, color

    endfor
       
stop
    
;   ploterror, zspec, zphot, zerr, ps=8, color=color.green, errcolor=color.green, $
;     xr=[0,5], yr=[0,5], xsty=3, ysty=3, xtit='Spectroscopic Redshift', ytit='Photometric Redshift'
    plot, zspec, zphot, ps=8, color=color.green, $
      xr=[0,5], yr=[0,5], xsty=3, ysty=3, xtit='Spectroscopic Redshift', ytit='Photometric Redshift'
    oplot, findgen(10), line=2, color=color.red, thick=3.0

    if keyword_set(postscript) then begin
       ps_open, ppath+'dum', /ps_fonts
       device, /inches, /times              
    endif else window, 0, xs=450, ys=450

; plotting preferences

    colortable2, c
    plotsym, 0, 1, /fill
    plotfaves
    
;    binsize = 0.01
;    
;    plothist, zspec, xspec, yspec, binsize=binsize, /noplot
;    plothist, zphot, xphot, yphot, binsize=binsize, /noplot
;    ymax = max(yspec)>max(yphot)
;    
;    plothist, zspec, binsize=binsize, xsty=1, ysty=1, color=16, line=0, $
;      xtit='Redshift', ytit='Number', xthick=2.0, ythick=2.0, $
;      charsize=2.0, charthick=2.0, thick=2.5, yr=[0,1.05*ymax] ;, /peak
;    plothist, zphot, x, y, /overplot, binsize=binsize, line=0, color=4, xsty=1, ysty=1, thick=2.5, $
;      /fill, /fline, fcolor=4, forientation=45, fspacing=0.2 ;, /peak
    
; overplot statistic histograms
       
;   if keyword_set(postscript) then ps_close

; residuals plot    

    resid = zspec-zphot
    
    yminmax = abs(min(resid-zerr))>abs(max(resid+zerr))
    plot, [min(zspec),max(zspec)], [-1.1*yminmax,1.1*yminmax], xsty=3, ysty=3, color=c.white, $
      xtit=textoidl('z_{spec}'), /nodata, position=[0.1,0.15,0.95,0.35], $
      ytit=textoidl('z_{spec}')+'-'+textoidl('z_{phot}')
    oploterror, zspec, resid, zerr, ps=8, color=c.cyan, errcolor=c.cyan
    oplot, [!x.crange[0],!x.crange[1]], [0,0], line=2, thick=2.5, color=c.red
    
;   yminmax = abs(min(zphot-zerr))>abs(max(zphot+zerr))
;   plot, [min(zspec),max(zspec)], [0,1.1*yminmax], xsty=3, ysty=3, color=c.white, $
    plot, [min(zspec),max(zspec)], [min(zspec),max(zspec)], xsty=3, ysty=3, color=c.white, $
      ytit=textoidl('z_{phot}'), /nodata, xtickname=replicate(' ',10), $
      position=[0.1,0.35,0.95,0.95], /noerase
    oploterror, zspec, zphot, zerr, ps=8, color=c.cyan, errcolor=c.cyan
    oplot, findgen(10), line=2, thick=2.5, color=c.red

    plotfaves, /restore
    
    if keyword_set(postscript) then ps_close

stop    


return
end
    
