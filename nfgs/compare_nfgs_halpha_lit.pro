pro compare_nfgs_halpha_lit, paper=paper, postscript=postscript
; jm04apr22uofa - written

; compare H-alpha + [N II] emission-line fluxes between the spectral
; atlas and the literature; we need to degrade the spectroscopic
; fluxes by re-incorporating the Balmer absorption and re-reddening
; for Galactic extinction

    pspath = atlas_path(/plots)

    if keyword_set(paper) then begin
       postscript = 1
       pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/'
    endif

    nfgs = read_nfgs(/silent)
    litcat = read_halpha_literature()

; match by coordinates

    raref = 15.0*im_hms2dec(nfgs.ra)
    deref = im_hms2dec(nfgs.dec)

    ra = 15.0*im_hms2dec(litcat.ra)
    de = im_hms2dec(litcat.dec)
    
    ntot = djs_angle_match(raref,deref,ra,de,dtheta=15.0/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    match = where(mindx ne -1,nmatch)
;   srt = sort(mdist[match])
;   niceprint, nfgs[match[srt]].galaxy, litcat[mindx[match[srt]]].name1, $
;     litcat[mindx[match[srt]]].name2, mdist[match[srt]]*3600.0

    splog, 'There are '+strn(nmatch)+' nfgs galaxies with literature H-alpha fluxes.'

    nfgs = nfgs[match]
    litcat = litcat[mindx[match]]

; ---------------------------------------------------------------------------    
; emission-line flux: Ha+[N II]    
; ---------------------------------------------------------------------------    
    
    if keyword_set(postscript) then $
      dfpsplot, pspath+'compare_nfgs_halpha_flux.ps', /square else $
      window, 0, xs=550, ys=550
    
    good = where((nfgs.h_alpha[1] gt 0.0) and (nfgs.nii_6548[1] gt 0.0) and $
      (nfgs.nii_6584[1] gt 0.0) and (litcat.ha_flux gt -900.0),ngood)

    km = nfgs[good]
    cat = litcat[good]

; one object, IRAS08572+3915, has no measured Balmer absorption
; because of its redshift.  replace with the H-beta absorption

    w = where(km.babs_h_alpha_ew[0] eq 0.0,nw)
    if nw ne 0L then km.babs_h_alpha_ew = km.babs_h_beta_ew
    
; un-correct for Balmer absorption     
    
    babscor = km.babs_h_alpha_ew[0]*km.h_alpha_continuum[0]
    babscor_err = sqrt((km.babs_h_alpha_ew[0]*km.h_alpha_continuum[1])^2.0 + $
      (km.babs_h_alpha_ew[1]*km.h_alpha_continuum[0])^2.0)

    splog, 'The median fractional Balmer absorption correction to the flux is '+$
      string(median(100*babscor/km.h_alpha[0]),format='(I0)')+'%.'

    km.h_alpha[0] = km.h_alpha[0] - babscor
    km.h_alpha[1] = sqrt(km.h_alpha[1]^2 + babscor_err^2.0)

; un-correct for Galactic reddening    

    km.h_alpha[0] = km.h_alpha[0]*10^(-0.4*k_lambda(km.h_alpha_wave,/odonnell)*km.ebv_mw)
    km.nii_6548[0] = km.nii_6548[0]*10^(-0.4*k_lambda(km.nii_6548_wave,/odonnell)*km.ebv_mw)
    km.nii_6584[0] = km.nii_6584[0]*10^(-0.4*k_lambda(km.nii_6584_wave,/odonnell)*km.ebv_mw)
    
    litflux = cat.ha_flux
    kmflux = km.h_alpha[0] + km.nii_6548[0] + km.nii_6584[0]
    kmferr = sqrt(km.h_alpha[1]^2 + km.nii_6548[1]^2 + km.nii_6584[1]^2)

; take the log of the quantities
    
    kmferr = kmferr/kmflux/alog(10.0) 
    kmflux = alog10(kmflux)
    
    yrange1 = [min(kmflux)<min(litflux),max(kmflux)>max(litflux)]
    xrange = yrange1

    resid = kmflux - litflux
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    pagemaker, nx=1, ny=2, position=pos, /normal, $
      xspace=0.0, xmargin=[1.8,0.3], ymargin=[0.3,1.3], yspace=0.0

    psize = 1.4
    
    im_symbols, 108, psize=psize
    plot, litflux, kmflux, ps=8, xsty=3, ysty=3, xthick=5.0, $
      ythick=5.0, charthick=5.0, charsize=2.0, position=pos[*,0], $
      yrange=yrange1, xrange=xrange, xtickname=replicate(' ',10), $
      ytitle=textoidl('log H\alpha+[N II]  [NFGS]')
    djs_oplot, !x.crange, !y.crange, line=0, thick=5.0

;   im_legend, ['Photometric','Non-Photometric','Undersampled'], $
;     psym=[108,106,105], fill=[0,0,0], /left, /top, box=0, charsize=1.5, $
;     charthick=5.0, spacing=1.7

; residuals    
    
    im_symbols, 108, psize=psize
    plot, litflux, resid, ps=8, yrange=yrange2, xrange=xrange, $
      position=pos[*,1], /noerase, xthick=5.0, ythick=5.0, xsty=3, ysty=3, $
      charthick=5.0, charsize=2.0, xtitle=textoidl('log H\alpha+[N II]  [Literature]'), $
      ytitle='Residuals'
    djs_oplot, !x.crange, [0,0], line=0, thick=5.0

    splog, 'Statistics:'
    junk = im_stats(resid,sigrej=3.0,/verbose)

    if keyword_set(postscript) then dfpsclose
    
; ---------------------------------------------------------------------------    
; emission-line EW: Ha+[N II]
; ---------------------------------------------------------------------------    
    
    if keyword_set(postscript) then $
      dfpsplot, pspath+'compare_nfgs_halpha_ew.ps', /square else $
      cc = get_kbrd(1)
    
    good = where((nfgs.h_alpha_ew[1] gt 0.0) and (nfgs.nii_6548_ew[1] gt 0.0) and $
      (nfgs.nii_6584_ew[1] gt 0.0) and (litcat.ew gt -900.0),ngood)

    km = nfgs[good]
    cat = litcat[good]
    
; un-correct for Balmer absorption; ignore Galactic reddening

    splog, 'The median fractional Balmer absorption correction to the EW is '+$
      string(median(100*km.babs_h_alpha_ew[0]/km.h_alpha_ew[0]),format='(I0)')+'%.'

    km.h_alpha_ew[0] = km.h_alpha_ew[0] - km.babs_h_alpha_ew[0]
    km.h_alpha[1] = sqrt(km.h_alpha_ew[1]^2 + km.babs_h_alpha_ew[1]^2.0)

; do more stuff here    
    
    litflux = alog10(cat.ew)
    kmflux = km.h_alpha_ew[0] + km.nii_6548_ew[0] + km.nii_6584_ew[0]
    kmferr = sqrt(km.h_alpha_ew[1]^2 + km.nii_6548_ew[1]^2 + km.nii_6584_ew[1]^2)

; take the log    
    
    kmferr = kmferr/kmflux/alog(10.0) 
    kmflux = alog10(kmflux)

    yrange1 = [min(kmflux)<min(litflux),max(kmflux)>max(litflux)]
    xrange = yrange1

    resid = kmflux - litflux
    yrange2 = max(abs(resid))*[-1.1,1.1]
    
    pagemaker, nx=1, ny=2, position=pos, /normal, $
      xspace=0.0, xmargin=[1.6,0.3], ymargin=[0.3,1.3], yspace=0.0

    psize = 1.4
    
    im_symbols, 105, fill=1, psize=psize
    plot, litflux, kmflux, ps=8, xsty=3, ysty=3, xthick=5.0, $
      ythick=5.0, charthick=5.0, charsize=2.0, position=pos[*,0], $
      yrange=yrange1, xrange=xrange, xtickname=replicate(' ',10), $
      ytitle=textoidl('log EW(H\alpha+[N II])  [NFGS]')
    djs_oplot, !x.crange, !y.crange, line=0, thick=5.0

;   im_legend, ['Integrated','Undersampled'], psym=[105,105], fill=[1,0], $
;     /left, /top, box=0, charsize=1.5, charthick=5.0, spacing=1.7

; residuals    

    im_symbols, 105, psize=psize, fill=1
    plot, litflux, resid, ps=8, yrange=yrange2, xrange=xrange, $
      position=pos[*,1], /noerase, xthick=5.0, ythick=5.0, xsty=3, ysty=3, $
      charthick=5.0, charsize=2.0, xtitle=textoidl('log EW(H\alpha+[N II])  [Literature]'), $
      ytitle='Residuals'
    djs_oplot, !x.crange, [0,0], line=0, thick=5.0

    splog, 'Statistics [Combined]'
    junk = im_stats(resid,sigrej=3.0,/verbose)

    if keyword_set(postscript) then dfpsclose

stop
    
return
end
