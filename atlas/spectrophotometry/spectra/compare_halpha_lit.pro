pro compare_halpha_lit, atlas, lit, paper=paper, postscript=postscript
; jm04mar22uofa
; jm04apr20uofa - modified and updated
; jm05aug02uofa - updated

; compare H-alpha + [N II] emission-line fluxes between the spectral
; atlas and the literature; we need to degrade the spectroscopic
; fluxes by re-incorporating the Balmer absorption and re-reddening
; for Galactic extinction

    snrcut = 3.0
    
    pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/'

    if (n_elements(atlas) eq 0L) then atlas = read_integrated()
    if (n_elements(lit) eq 0L) then lit = read_literature_halpha_catalog()
    if keyword_set(paper) then postscript = 1

    lcharsize = 1.2
    pcharsize = 1.5

; match by coordinates

    raref = 15.0*im_hms2dec(atlas.ra)
    deref = im_hms2dec(atlas.dec)

    ra = 15.0*im_hms2dec(lit.ra)
    de = im_hms2dec(lit.dec)
    
    ntot = djs_angle_match(raref,deref,ra,de,dtheta=15.0/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    match = where(mindx ne -1,nmatch)
    srt = sort(mdist[match])
    niceprint, atlas[match[srt]].galaxy, atlas[match[srt]].alt_galaxy, $
      lit[mindx[match[srt]]].galaxy, lit[mindx[match[srt]]].alt_galaxy, mdist[match[srt]]*3600.0

; only keep unique objects

    unique = uniq(lit[mindx[match]].id)
    nmatch = n_elements(unique)

    matchatlas = atlas[match[unique]]
    matchlit = lit[mindx[match[unique]]]
    
    splog, 'There are '+string(nmatch,format='(I0)')+' unique galaxies with literature H-alpha fluxes.'

; initialize some plotting variables

    xmargin = [1.2,0.3]
    ymargin = [0.3,1.2]

    width = 6.0
    height = [3.0,2.5]

    xpage = total(xmargin)+total(width)
    ypage = total(ymargin)+total(height)

    pagemaker, nx=1, ny=2, position=pos, /normal, xspace=0.0, width=width, $
      height=height, xmargin=xmargin, ymargin=ymargin, yspace=0.0, $
      xpage=xpage, ypage=ypage

    fluxrange = [-14,-10]
    ewrange = [0.3,3.1]
    psize = 1.0
    
; ---------------------------------------------------------------------------    
; emission-line flux: Ha+[N II]    
; ---------------------------------------------------------------------------    
    
    if keyword_set(postscript) then begin
       dfpsplot, pspath+'compare_halpha_flux.eps', xsize=xpage, ysize=ypage, /encapsulated
       postthick = 5.0
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
    endelse

    atlasgood = where((matchatlas.h_alpha[0]/matchatlas.h_alpha[1] gt snrcut) and $
      (matchatlas.nii_6584[0]/matchatlas.nii_6548[1] gt snrcut))
;   litgood = where((matchlit.ha_nii gt 0.0) and (matchlit.ha_nii_err gt 0.0))
    litgood = where((matchlit.ha_nii gt 0.0))
    good = cmset_op(atlasgood,'AND',litgood)
    ngood = n_elements(good)

;   good = where((matchatlas.h_alpha[1] gt 0.0) and (matchatlas.nii_6548[1] gt 0.0) and $
;     (matchatlas.nii_6584[1] gt 0.0) and (matchlit.ha_nii gt 0.0),ngood)
    splog, string(ngood,format='(I0)')+' objects satisfy the flux S/N cut.'

; un-correct for Galactic reddening    

    matchatlas[good].h_alpha[0] = matchatlas[good].h_alpha[0]*10^(-0.4*$
      k_lambda(matchatlas[good].h_alpha_wave,/odonnell)*matchatlas[good].ebv_mw)
    matchatlas[good].nii_6548[0] = matchatlas[good].nii_6548[0]*10^(-0.4*$
      k_lambda(matchatlas[good].nii_6548_wave,/odonnell)*matchatlas[good].ebv_mw)
    matchatlas[good].nii_6584[0] = matchatlas[good].nii_6584[0]*10^(-0.4*$
      k_lambda(matchatlas[good].nii_6584_wave,/odonnell)*matchatlas[good].ebv_mw)

    flux = matchatlas[good].h_alpha[0] + matchatlas[good].nii_6548[0] + matchatlas[good].nii_6584[0]
    ferr = sqrt(matchatlas[good].h_alpha[1]^2 + matchatlas[good].nii_6548[1]^2 + matchatlas[good].nii_6584[1]^2)

; take the log of the quantities
    
    ferr = ferr/flux/alog(10.0) 
    flux = alog10(flux)

    photo = where(matchatlas[good].drift_photflag eq 'Y',comp=nophoto)
    yphoto = flux[photo]
    ynophoto = flux[nophoto]

    xphoto = -matchlit[good[photo]].ha_nii
    xnophoto = -matchlit[good[nophoto]].ha_nii
    
    resphoto = yphoto-xphoto
    resnophoto = ynophoto-xnophoto

    im_symbols, 108, psize=psize, /fill
    djs_plot, xphoto, yphoto, ps=8, xsty=3, ysty=3, charsize=pcharsize, $
      charthick=postthick, xthick=postthick, ythick=postthick, xrange=fluxrange, $
      yrange=fluxrange, position=pos[*,0], xtickname=replicate(' ',10), $
      ytitle=textoidl('log (H\alpha+[N II])')+' [This Paper]'
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    im_symbols, 106, psize=psize, fill=0, thick=postthick
    djs_oplot, xnophoto, ynophoto, ps=8

; residuals

    resrange = (max(abs(resphoto))>max(abs(resnophoto)))*[-1.1,1.1]
    
    im_symbols, 108, psize=psize, /fill
    djs_plot, xphoto, resphoto, xsty=3, ysty=3, ps=8, $
      position=pos[*,1], /noerase, yrange=resrange, xrange=fluxrange, $
      charsize=pcharsize, charthick=postthick, xthick=postthick, ythick=postthick, $
      xtitle=textoidl('log (H\alpha+[N II])')+' [Literature]', yminor=3, ytitle='Residuals'
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    im_symbols, 106, psize=psize, fill=0, thick=postthick
    djs_oplot, xnophoto, resnophoto, ps=8

    splog, 'Flux statistics [Photometric, Non-Photometric, Combined]'
    stats = im_stats(resphoto,sigrej=3.0,/verbose)
    junk = im_stats(resnophoto,sigrej=3.0,/verbose,/no_head)
    junk = im_stats([resphoto,resnophoto],sigrej=3.0,/verbose,/no_head)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   legend, textoidl(xstr), /right, /bottom, box=0, charsize=lcharsize, $
;     charthick=postthick, clear=keyword_set(postscript)

    if keyword_set(postscript) then dfpsclose else cc = get_kbrd(1)

; ---------------------------------------------------------------------------    
; emission-line EW: Ha+[N II]
; ---------------------------------------------------------------------------    

    if keyword_set(postscript) then begin
       dfpsplot, pspath+'compare_halpha_ew.eps', xsize=xpage, ysize=ypage, /encapsulated
       postthick = 5.0
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
    endelse

    atlasgood = where((matchatlas.h_alpha_ew[0]/matchatlas.h_alpha_ew[1] gt snrcut) and $
      (matchatlas.nii_6584_ew[0]/matchatlas.nii_6548_ew[1] gt snrcut))
;   litgood = where((matchlit.ha_nii gt 0.0) and (matchlit.ha_nii_err gt 0.0))
    litgood = where((matchlit.ha_nii_ew gt 0.0))
    good = cmset_op(atlasgood,'AND',litgood)
    ngood = n_elements(good)

;   good = where((matchatlas.h_alpha_ew[1] gt 0.0) and (matchatlas.nii_6548_ew[1] gt 0.0) and $
;     (matchatlas.nii_6584_ew[1] gt 0.0) and (matchlit.ha_nii_ew gt 0.0),ngood)

    splog, string(ngood,format='(I0)')+' objects satisfy the EW S/N cut.'

; compare apples to apples    

    ew = matchatlas[good].h_alpha_ew[0] + matchatlas[good].nii_6548_ew[0] + matchatlas[good].nii_6584_ew[0]
    ewerr = sqrt(matchatlas[good].h_alpha_ew[1]^2 + matchatlas[good].nii_6548_ew[1]^2 + matchatlas[good].nii_6584_ew[1]^2)

; take the log of the quantities
    
    ewerr = ewerr/ew/alog(10.0) 
    ew = alog10(ew)

    photo = where(matchatlas[good].drift_photflag eq 'Y',comp=nophoto)
    yphoto = ew[photo]
    ynophoto = ew[nophoto]

    xphoto = alog10(matchlit[good[photo]].ha_nii_ew)
    xnophoto = alog10(matchlit[good[nophoto]].ha_nii_ew)
    
    resphoto = yphoto-xphoto
    resnophoto = ynophoto-xnophoto

    im_symbols, 108, psize=psize, /fill
    djs_plot, xphoto, yphoto, ps=8, xsty=3, ysty=3, charsize=pcharsize, $
      charthick=postthick, xthick=postthick, ythick=postthick, xrange=ewrange, $
      yrange=ewrange, position=pos[*,0], xtickname=replicate(' ',10), $
      ytitle=textoidl('log EW(H\alpha+[N II])')+' [This Paper]'
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    im_symbols, 106, psize=psize, fill=0, thick=postthick
    djs_oplot, xnophoto, ynophoto, ps=8

; residuals

    resrange = (max(abs(resphoto))>max(abs(resnophoto)))*[-1.1,1.1]
    
    im_symbols, 108, psize=psize, /fill
    djs_plot, xphoto, resphoto, xsty=3, ysty=3, ps=8, $
      position=pos[*,1], /noerase, yrange=resrange, xrange=ewrange, $
      charsize=pcharsize, charthick=postthick, xthick=postthick, ythick=postthick, $
      xtitle=textoidl('log EW(H\alpha+[N II])')+' [Literature]', yminor=3, ytitle='Residuals'
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    im_symbols, 106, psize=psize, fill=0, thick=postthick
    djs_oplot, xnophoto, resnophoto, ps=8

    splog, 'EW statistics [Photometric, Non-Photometric, Combined]'
    junk = im_stats(resphoto,sigrej=3.0,/verbose)
    junk = im_stats(resnophoto,sigrej=3.0,/verbose,/no_head)
    stats = im_stats([resphoto,resnophoto],sigrej=3.0,/verbose,/no_head)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   legend, textoidl(xstr), /right, /bottom, box=0, charsize=lcharsize, $
;     charthick=postthick, clear=keyword_set(postscript)
    
    if keyword_set(postscript) then dfpsclose
    
stop    

return
end
