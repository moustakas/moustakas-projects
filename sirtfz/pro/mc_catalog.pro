function zstatscalc, zspec, zphot

    zstats = fltarr(3)
    
    ntotal = n_elements(zspec)
    deltaz = zspec-zphot

    zstats[0] = total(deltaz)/ntotal                                     ; mean deviation
    zstats[1] = sqrt(total((deltaz-zstats[0])^2)/(ntotal-1))             ; standard deviation
    zstats[2] = sqrt(total((deltaz-zstats[0])^2/(1+zspec)^2)/(ntotal-1)) ; relative deviation

return, zstats
end

pro mc_catalog, filters, fref, zstats_all, zstats_zbins, flimit=flimit, templates=templates, ncatalog=ncatalog, $
      mcnum=mcnum, fmin=fmin, fmax=fmax, zmin=zmin, zmax=zmax, dz=dz, alpha=alpha, $
      mindetect=mindetect, snr=snr, nomonte=nomonte, postscript=postscript, color=color
;+
; NAME:
;	MC_CATALOG
;
; PURPOSE:
;
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;	filters   - bandpass filters in which to compute photometry
;	fref      - filter reference number [zero indexed]
;
; OPTIONAL INPUTS:
;	flimit    - limiting flux corresponding to filters [Jy]
;	templates - SED models
;	ncatalog  - number of galaxies in the simulated catalog
;	mcnum     - number of Monte Carlo realizations of each galaxy
;	fmin      - minimum apparent flux in fref for a galaxy to be
;                   accepted in the catalog [Jy]
;	fmax      - maximum apparent flux [Jy]
;	zmin/zmax - minimum/maximum redshift
;	dz        - redshift bin size
;	alpha     - power of luminosity evolution
;	
; KEYWORD PARAMETERS:
;	nomonte   - do not compute Monte Carlo realizations
;	ps        - generate postscript output
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMON BLOCKS:
;	sirtf_simulations
;
; COMMENTS:
;	Need to add random galaxy types.
;
; EXAMPLE:
;
;
; PROCEDURES USED:
;
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 Sep 16, U of A
;-
    
    common sirtf_simulations
    
    if n_params() lt 2L then begin
       print, 'Syntax - mc_catalog, filters, fref'
       return
    endif

    if not keyword_set(ncatalog) then ncatalog = 50L
    if not keyword_set(templates) then templates = sirtf.templates

    if not keyword_set(mcnum) then mcnum = 100L
    if keyword_set(nomonte) then mcnum = 0L
    
    obands = filter_match(filters,sirtf.bandcube.bandnames)
    nbands = n_elements(obands)
    nseds = n_elements(sirtf.sedcube)
    
    if not keyword_set(snr) then snr = 5.0
    if not keyword_set(fmin) then fmin = snr*sirtf.bandcube[obands[fref]].flimit ; [Jy]
    if not keyword_set(fmax) then fmax = 10.0     ; [10 Jy]
    if not keyword_set(alpha) then alpha = sirtf.evolution.alpha

    if not keyword_set(dz) then dz = 0.01
    if not keyword_set(zmin) then zmin = 0.01
    if not keyword_set(zmax) then zmax = 2.00

    zarray = (findgen((zmax-zmin+dz)/dz))*dz+zmin
    nz = n_elements(zarray)

; plotting preferences
    
; read the model grids.  compute fluxes based on the rest
; monochromatic luminosity and the luminosity distance
    
    read_model_grids, filters, mzarray, mflux, templates=templates
    nmz = n_elements(mzarray)

    zinfo = {zphot_mode  : fltarr(1), $
             zphot_mean  : fltarr(1), $
             zphot_median: fltarr(1), $
             sigmaz      : fltarr(2), $
             chi2_nu     : fltarr(1), $
             zpdf        : dblarr(nmz),$
             constant    : fltarr(1), $
             zindx       : lonarr(1), $
             tindx       : lonarr(1)}

    if mcnum ge 1L then mc_zphot = fltarr(mcnum)
    zspec = fltarr(ncatalog)
    zphot = fltarr(ncatalog)
    zerr = fltarr(ncatalog)
    type = lonarr(ncatalog)
    
    lumrest = fltarr(nbands,nseds)

; segregate galaxies by infrared luminosity
    
    lbolir = sirtf.sedcube.lbolir
    early = where(lbolir lt 1E10,nearly)
    starburst = where((lbolir ge 1E10) and (lbolir lt 1E11),nstarburst)
    ulirg = where(lbolir ge 1E11,nulirg)

    class = lonarr(nseds)
    class[early] = 0L
    class[starburst] = 1L
    class[ulirg] = 2L
    
    nsources = 0L
    while nsources lt ncatalog do begin ; generate the catalog

       z = (zarray[floor(randomu(seed,1)*nz)])[0] ; random redshift
       l_evol = (1.0+z)^alpha                     ; luminosity evolution
       dlum = dluminosity(z)*3.085678D16          ; luminosity distance

       for i = 0L, nseds-1L do $
         lumrest[*,i] = interpol(*sirtf.sedcube[i].mlum_nu,*sirtf.sedcube[i].lambda*(1.0+z),$  ; [W/Hz]
                                 sirtf.bandcube[obands].lambda_eff*1E-4)

       f_nu = (1+z) * l_evol * lumrest / (4.0D * !dpi * dlum * dlum) * 1E26 ; [Jy]
       good = where(f_nu[fref,*] gt fmin,ngood) ; above the minimum flux

       if ngood ne 0L then begin

          print, format='("Source ",i0," of ",i0,".",a1,$)', nsources, ncatalog, string(13b)

          repeat begin

             rantype = (floor(randomu(seed,1)*3))[0] ; random type [early,starburst,ulirg]
             match = where(class eq rantype,nmatch)
             sindx = where(good eq (match[floor(randomu(seed,1)*nmatch)])[0])
;            sindx = floor(randomu(seed,1)*ngood) ; random galaxy          

          endrep until sindx[0] ne -1L
             
          gflux = reform(f_nu[*,good[sindx]])
          ferror = sqrt(sirtf.bandcube[obands].flimit^2 + (0.05*gflux)^2) ; flux error

          zspec[nsources] = z
          type[nsources] = good[sindx]

; convolve the p(z) with the error function via a Monte Carlo analysis

          if keyword_set(nomonte) then begin

             dflux = gflux + randomn(seed,nbands)*ferror
             photoz, dflux, ferror, mflux, mzarray, pdz=pdz, zinfo

             zphot[nsources] = zinfo.zphot_mode
             zerr[nsources] = 0.0
             
          endif else begin

             for k = 0L, mcnum-1L do begin
                
                dflux = gflux + randomn(seed,nbands)*ferror
                photoz, dflux, ferror, mflux, mzarray, pdz=pdz, zinfo

                mc_zphot[k] = zinfo.zphot_mode
                
             endfor

             zphot[nsources] = mean(mc_zphot)
             zerr[nsources] = stddev(mc_zphot)

          endelse

; should we add the individual uncertainties in quadrature

          if z gt 2.0 then stop
          
          nsources = nsources + 1L

            print, nsources, z, zinfo.zphot_median, zinfo.zphot_mode, zinfo.zphot_mean
            print, zinfo.constant*reform(mflux[zinfo.tindx,zinfo.zindx,*])
            
             wave0 = *sirtf.sedcube[zinfo.tindx].lambda
             lum0 = *sirtf.sedcube[zinfo.tindx].mlum_nu
             redshift_sed, zinfo.zphot_mode, wave0, lum0, wavez, fluxz, /jansky

; SED plot

;            plot, wavez, fluxz, xrange=1E-4*[min(sirtf.bandcube[obands].lambda_eff)/3,$
             plot, wave0*(1+zinfo.zphot_mode), lum0*zinfo.constant, xrange=1E-4*[min(sirtf.bandcube[obands].lambda_eff)/3,$
                                              3*max(sirtf.bandcube[obands].lambda_eff)], $
               xsty=3, ysty=3, /ylog, position=[0.22,0.52,0.95,0.93], /xlog, xthick=2.0, $
               ythick=2.0, charsize=1.8, charthick=2.0, thick=2.0, $
               xtit=textoidl('\lambda')+' ('+textoidl('\mu')+'m)', ytit=textoidl('f_{\nu}')+' [Jy]'
             oploterror, sirtf.bandcube[obands].lambda_eff*1E-4, gflux, $
               sirtf.bandcube[obands].dlambda*1E-4, ferror, ps=8, color=color.green, $
               errcolor=color.green, thick=2.5
             
; response function plot
           
             plot, mzarray, zinfo.zpdf, /noerase, $
               xsty=3, ysty=3, position=[0.22,0.13,0.95,0.4], xtit='Redshift', $
               ytit='Likelihood', xthick=2.0, ythick=2.0, charsize=1.8, charthick=2.0, $
               thick=2.0

; legend
;
;             text = [textoidl('\chi^{2}_{\nu}: '+strn(zinfo.chi2_nu,length=6)), $
;                     textoidl('z_{phot}'+': '+strn(zinfo.zphot_mode,length=5)), $
;                     textoidl('z_{spec}'+': '+strn(z,length=5))]
;             if zinfo.zphot_mode lt max(mzarray)/2. then $
;               legend, text, box=0, /right, /top, charsize=1.3, charthick=2.0 else $
;               legend, text, box=0, /left, /top, charsize=1.3, charthick=2.0

          endif

    endwhile

; color-code by type

    irboltype = sirtf.sedcube[type].lbolir
    early = where(irboltype lt 1E10,nearly)
    starburst = where((irboltype ge 1E10) and (irboltype lt 1E11),nstarburst)
    ulirg = where(irboltype ge 1E11,nulirg)

    plotsym, 0, 0.8, /fill
    colortable2, color
    plotcolors = [44,30,78]

    ppath = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='plots/zmonte')
    mini = sirtf.bandcube[obands].mininames
    fname = strjoin(mini,'_')

    if keyword_set(postscript) then begin
       ps_open, ppath+fname, /ps_fonts, /portrait, color=color
       device, /inches, /times, xsize=7, ysize=7
    endif else window, 0, xs=550, ys=550

    plot, [0], [0], xrange=[zmin,zmax], yrange=[zmin,zmax], /nodata, $
      xsty=3, ysty=3, xtit='Model Redshift', ytit='Photometric Redshift'
    if nulirg ne 0L then oploterror, zspec[ulirg], zphot[ulirg], zerr[ulirg], $
      ps=8, color=plotcolors[2], errcolor=plotcolors[2]
    if nstarburst ne 0L then oploterror, zspec[starburst], zphot[starburst], $
      zerr[starburst], ps=8, color=plotcolors[1], errcolor=plotcolors[1]
    if nearly ne 0L then oploterror, zspec[early], zphot[early], zerr[early], $
      ps=8, color=plotcolors[0], errcolor=plotcolors[0]
    oplot, findgen(10), line=0, thick=2.0

    if fmin le 1E-3 then $
      fmintext = strn(fmin*1E6,format='(G0.0)')+' '+textoidl('\mu') else $
      fmintext = strn(fmin*1E3,format='(G0.0)')+' m'

    legend, [textoidl('f_{'+sirtf.bandcube[obands[fref]].mininames+'} > ')+$
                      fmintext+'Jy ('+strn(snr,format='(G0.0)')+textoidl('\sigma')+')',' ',filters,' '], $
      /left, /top, box=0, charsize=1.2, charthick=2.0

    if keyword_set(postscript) then ps_close

; statistics of the fit

    deltaz = zspec-zphot
    badz = where(abs(deltaz) ge 1.0,nbadz,comp=goodz,ncomp=ngoodz)

    percentbad = 100.0*nbadz/ncatalog      ; fraction of catastrophic errors
    
    zstats_all = zstatscalc(zspec[goodz],zphot[goodz])

    zbin1 = where((zspec[goodz] ge 0.0) and (zspec[goodz] lt 1.0),nzbin1)
    zbin2 = where((zspec[goodz] ge 1.0) and (zspec[goodz] lt 2.0),nzbin2)
    zbin3 = where((zspec[goodz] ge 2.0),nzbin3)

    zstats_zbins = fltarr(9)
    
    if nzbin1 ne 0L then zstats_zbins[0:2] = zstatscalc(zspec[goodz[zbin1]],zphot[goodz[zbin1]])    
    if nzbin2 ne 0L then zstats_zbins[3:5] = zstatscalc(zspec[goodz[zbin2]],zphot[goodz[zbin2]])    
    if nzbin3 ne 0L then zstats_zbins[6:8] = zstatscalc(zspec[goodz[zbin3]],zphot[goodz[zbin3]])

;   splog, strn(strjoin(mini,','),length=20), zstats_all, zstats_zbins, /noname, /append, filename = rpath+'zmonte.log'
    
; filter id's

;    id = strn(sirtf.bandcube[obands[fiducial]].lambda0,format='(G0.0)')+'_'+strn(fapar,format='(G0.2)')+'_'
;    
;    filtids = sirtf.bandcube[obands[0]].mininames
;    mininames = sirtf.bandcube[obands].mininames
;    for i = 1L, nbands-1L do filtids = filtids+'_'+mininames[i]
;
;; plot analysis
;
;    resid = ztrue-photoz
;    coef = linfit((1+ztrue),sigmaz)
;
;;   s = sort(sigmaz)
;;   plot, (1+ztrue[s]), sigmaz[s], ps=8, xr=[0,6]
;;   oplot, (1+ztrue[s]), poly((1+ztrue[s]),coef), line=0, thick=2.5
;    
;    plot, [0,6], [0,6], xsty=3, ysty=3, $ ; color=16, $
;      ytit='Photometric Redshift', /nodata, xtickname=replicate(' ',10), $
;      position=[0.1,0.35,0.95,0.95], xthick=3.0, ythick=3.0, $
;      charsize=1.2, charthick=3.0
;    oploterror, ztrue, photoz, sigmaz, ps=8, color=2, errcolor=2, thick=3.0
;    oplot, findgen(10), line=2, thick=3.0;, color=6
;    legend, [textoidl('f_{'+strn(sirtf.bandcube[obands[fiducial]].lambda0,format='(G0.0)')+$
;                      '\mu}'+textoidl('_{m}')+' > '+strn(fapar,format='(G0.2)'))+' mJy ('+strn(fsigma,format='(G0.0)')+$
;             textoidl('\sigma')+')',' ',filters,' ',$
;             textoidl('\sigma')+'(z) = '+strn(coef[1],format='(G0.2)')+$
;             '(1+z)'], $
;      /left, /top, box=0, charsize=1.2, charthick=3.0
;    
;    yminmax = abs(min(resid+sigmaz))>abs(max(resid+sigmaz))
;    plot, [0,6], [-1.1*yminmax,1.1*yminmax], xsty=3, ysty=3, $ ; color=16, $
;      xtit='Redshift', /nodata, position=[0.1,0.15,0.95,0.35], $
;      ytit='Residuals', /noerase, xthick=3.0, ythick=3.0, $
;      charsize=1.2, charthick=3.0
;    oploterror, ztrue, resid, sigmaz, ps=8, color=2, errcolor=2, thick=3.0
;    oplot, [!x.crange[0],!x.crange[1]], [0,0], line=2, thick=3.0
;
;; redshift histograms
;
;    plothist, ztrue, xtrue, ytrue, bin=0.2, /noplot
;    plothist, photoz, xzphot, yzphot, bin=0.2, /noplot
;    ymax = max(ytrue)>max(yzphot)
;    
;    plothist, ztrue, bin=0.2, xsty=1, ysty=1, xr=[0,6], line=0, $
;      xtit='Redshift', ytit='Number', xthick=3.0, ythick=3.0, $
;      charsize=1.2, charthick=3.0, thick=3.0, yr=[0,1.05*ymax], $
;      position=[0.65,0.43,0.9,0.68], /noerase
;    plothist, photoz, x, y, /overplot, bin=0.2, line=0, color=3, xsty=1, ysty=1, thick=3.0, $
;      /fill, /fline, fcolor=3, forientation=45, fspacing=0.2 ;, /peak
;    
;    if keyword_set(ps) then ps_close

;; contour plot
;    
;    bptot = transpose(total(bparray,3))
;
;    window, 2, xs=450, ys=450
;    contour, bptot, mzarray, sirtf.sedcube[sedtype].lum60_sun, xsty=3, ysty=3, $
;      xthick=2.0, ythick=2.0, thick=2.0, charsize=1.2, charthick=2.0, $
;      xtit='Redshift', ytit='

stop
    
return
end
