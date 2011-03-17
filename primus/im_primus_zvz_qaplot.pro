pro primus_zvz_qaplot, extract1, photoinfo1, oned1, template_data1, $
  nophotozlabel, photozlabel, photosuffix, refband
; jm08jul29nyu - make a z versus z QA plot
    
    dzcut = 0.05
    catacut = 0.03
    rmagcut = 30.0 ; 25.0 ; 23.5
    faintcut = 30.0 ; 23.0
    snrcut = 8.0
    
; plotting preliminaries
    
    xpage = 8.5
    ypage = 11.0

;   postthick1 = 2.0
;   postthick2 = 2.0
;   postthick3 = 2.0
    postthick1 = 4.0
    postthick2 = 3.0
    postthick3 = 8.0

    pagemaker, nx=2, ny=4, xpage=xpage, ypage=ypage, $
      xspace=0.0, yspace=[0.6,0.6,0.6], xmargin=[1.1,0.4], $
      ymargin=[0.3,1.1], width=3.5*[1,1], height=[2.7,1.7,1.7,1.7], $
      /normal, position=pos2
    pagemaker, nx=1, ny=1, xpage=xpage, ypage=ypage, $
      xspace=0.0, yspace=0.0, xmargin=[1.4,0.3], $
      ymargin=[0.7,1.1], width=6.8, height=6.8, $
      /normal, position=pos1
    pagemaker, nx=1, ny=3, xpage=xpage, ypage=ypage, $
      xspace=0.0, yspace=[0.0,1.5], xmargin=[1.1,0.4], $
      ymargin=[0.1,1.0], width=7.0, height=[2.8,2.8,2.8], $
      /normal, position=pos3
;   pagemaker, nx=1, ny=2, xpage=xpage, ypage=ypage, $
;     xspace=0.0, yspace=2.0, xmargin=[1.4,0.3], $
;     ymargin=[0.7,1.1], width=6.8, height=3.4*[1,1], $
;     /normal, position=pos3

; define some handy variables    

    gal = where((oned1.zbest_gal gt 0.0) and $
      (oned1.zbest_photoz_gal gt 0.0),ngal)
    extract = extract1[gal]
    photoinfo = photoinfo1[gal]
    oned = oned1[gal]

    tempspec = template_data1[oned.index_gal[0]] ; spectroscopic template
    tempphot = template_data1[oned.index_photoz_gal[0]] ; photometric template
    
    rmag = extract.mag
    snr = total([[extract.sn1],[extract.sn2]],2)/2.0 ; mean S/N

    filters = 'bessell_'+['U','B','V','R','I']+'.par'
    nfilt = n_elements(filters)
    maggies = fltarr(nfilt,ngal)
    for jj = 0L, ngal-1L do begin
       lambda = k_lambda_to_edges(extract[jj].restwave)
       maggies[*,jj] = k_project_filters(lambda,extract[jj].restflux,$
         filterlist=filters,/silent)
    endfor
    ub = reform(-2.5*alog10(maggies[0,*]/maggies[1,*]))
    
    ztrue = extract.z
    zspec = oned.zbest_gal
    zphot = oned.zbest_photoz_gal
    zphot2 = oned.zmin_photoz_gal[1] ; second minimum

    deltazspec = (zspec-ztrue)
    deltazphot = (zphot-ztrue)
    dzspec = deltazspec/(1.0+ztrue)
    dzphot = deltazphot/(1.0+ztrue)

    chi2spec = oned.chi2best_gal
    chi2phot = oned.chi2best_photoz_gal
;   niceprint, ztrue, zspec, zphot, chi2spec, chi2phot
    
    cataspec = where((abs(dzspec) gt catacut) and $
      (rmag lt rmagcut),ncataspec,comp=goodspec)
    cataphot = where((abs(dzphot) gt catacut) and $
      (rmag lt rmagcut),ncataphot,comp=goodphot)
    faint = where(rmag gt faintcut,comp=bright)

    specstats = im_stats(dzspec[goodspec])
    photstats = im_stats(dzphot[goodphot])
    
    labelspec = [$
      '\Delta = '+strtrim(string(specstats.mean,format='(F12.4)'),2)+$
      '\pm'+strtrim(string(specstats.sigma,format='(F12.4)'),2),$
      'F(>'+string(catacut,format='(F4.2)')+') = '+$
      strtrim(string(100.0*ncataspec/float(ngal),$
      format='(F12.1)'),2)+'%']
;     format='(F12.1)'),2)+'% ('+string(ncataspec,format='(I0)')+'/'+$
;     string(ngal,format='(I0)')+')']
    labelphot = [$
      '\Delta = '+strtrim(string(photstats.mean,format='(F12.4)'),2)+$
      '\pm'+strtrim(string(photstats.sigma,format='(F12.4)'),2),$
      'F(>'+string(catacut,format='(F4.2)')+') = '+$
      strtrim(string(100.0*ncataphot/float(ngal),$
      format='(F12.1)'),2)+'%']
;     format='(F12.1)'),2)+'% ('+string(ncataphot,format='(I0)')+'/'+$
;     string(ngal,format='(I0)')+')']
    
; compare the spectro+photoz redshifts against the true redshift 

    charsize1 = 1.2
    charsize2 = 1.0
    zrange = [-0.1,1.3]
    zaxis = findgen((2.0-(-1.0))/0.01)*0.01+(-1.0)
    residrange = 0.35*[-1,1]
    psize1 = 0.55
    psize2 = 0.15
    
    plot, [0], [0], /nodata, charsize=charsize1, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
      xrange=zrange, yrange=zrange, xtickinterval=0.2, $
      xtitle='Redshift [Simulated]', ytitle='Redshift [PRIMUS]', $
      position=pos2[*,0];, title='No PHOTOZ'
    oplot, zaxis, zaxis, line=0, thick=2.0
    oplot, zaxis, zaxis-dzcut, line=5, thick=2.0
    oplot, zaxis, zaxis+dzcut, line=5, thick=2.0
    plotsym, 0, psize1, fill=1, thick=postthick1
    djs_oplot, ztrue, zspec, ps=8, color='dark green'
    legend, nophotozlabel, /left, /top, box=0, $
      charsize=charsize1, charthick=postthick2, /clear
    legend, textoidl(labelspec), /right, /bottom, box=0, $
      charsize=charsize1, charthick=postthick2, /clear

    plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
      xrange=zrange, yrange=zrange, xtickinterval=0.2, $
      xtitle='Redshift [Simulated]', ytitle='', ytickname=replicate(' ',10), $
      position=pos2[*,1];, title='PHOTOZ'
    oplot, zaxis, zaxis, line=0, thick=2.0
    oplot, zaxis, zaxis-dzcut, line=5, thick=2.0
    oplot, zaxis, zaxis+dzcut, line=5, thick=2.0
    plotsym, 8, psize1, fill=1
    djs_oplot, ztrue, zphot, ps=8, color='red'
    legend, photozlabel, /left, /top, box=0, $
      charsize=charsize1, charthick=postthick2, /clear
    legend, textoidl(labelphot), /right, /bottom, box=0, $
      charsize=charsize1, charthick=postthick2, /clear
    
; residuals vs R magnitude

    residtitle = textoidl('\delta'+'z/(1+z)')
    
    rmagrange = [17.5,26.5]
;   rmagrange = [19.5,23.8]
    
    plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
      xrange=rmagrange, yminor=5, yrange=residrange, $
      xtitle='R (mag)', ytitle=residtitle, position=pos2[*,2]
    oplot, !x.crange, [0,0], line=0, thick=2.0
    oplot, !x.crange, catacut*[1,1], line=1, thick=4.0
    oplot, !x.crange, -catacut*[1,1], line=1, thick=4.0
    plotsym, 0, psize1, fill=1
    djs_oplot, rmag, dzspec, ps=8, color='dark green'
    
    plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
      xrange=rmagrange, yminor=5, yrange=residrange, $
      xtitle='R (mag)', ytitle='', ytickname=replicate(' ',10), position=pos2[*,3]
    oplot, !x.crange, [0,0], line=0, thick=2.0
    oplot, !x.crange, catacut*[1,1], line=1, thick=4.0
    oplot, !x.crange, -catacut*[1,1], line=1, thick=4.0
    plotsym, 8, psize1, fill=1
    djs_oplot, rmag, dzphot, ps=8, color='red'
    
; residuals vs S/N

    snrrange = [0,59.0]
    
    plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
      xrange=snrrange, yminor=5, yrange=residrange, $
      xtitle='<S/N>', ytitle=residtitle, position=pos2[*,4]
    oplot, !x.crange, [0,0], line=0, thick=2.0
    oplot, !x.crange, catacut*[1,1], line=1, thick=4.0
    oplot, !x.crange, -catacut*[1,1], line=1, thick=4.0
    plotsym, 0, psize1, fill=1
    djs_oplot, snr, dzspec, ps=8, color='dark green'
    
    plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
      xrange=snrrange, yminor=5, yrange=residrange, $
      xtitle='<S/N>', ytitle='', ytickname=replicate(' ',10), position=pos2[*,5]
    oplot, !x.crange, [0,0], line=0, thick=2.0
    oplot, !x.crange, catacut*[1,1], line=1, thick=4.0
    oplot, !x.crange, -catacut*[1,1], line=1, thick=4.0
    plotsym, 8, psize1, fill=1
    djs_oplot, snr, dzphot, ps=8, color='red'
    
; residuals vs U-B color

    ubrange = [0.21,1.7]
;   ubrange = [-0.5,0.5]
    
    plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
      xrange=ubrange, yminor=5, yrange=residrange, $
      xtitle='U-B (AB mag)', ytitle=residtitle, position=pos2[*,6]
    oplot, !x.crange, [0,0], line=0, thick=2.0
    oplot, !x.crange, catacut*[1,1], line=1, thick=4.0
    oplot, !x.crange, -catacut*[1,1], line=1, thick=4.0
    plotsym, 0, psize1, fill=1
    djs_oplot, ub, dzspec, ps=8, color='dark green'
    
    plot, [0], [0], /nodata, /noerase, charsize=charsize1, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
      xrange=ubrange, yminor=5, yrange=residrange, $
      xtitle='U-B (AB mag)', ytitle='', ytickname=replicate(' ',10), position=pos2[*,7]
    oplot, !x.crange, [0,0], line=0, thick=2.0
    oplot, !x.crange, catacut*[1,1], line=1, thick=4.0
    oplot, !x.crange, -catacut*[1,1], line=1, thick=4.0
    plotsym, 8, psize1, fill=1
    djs_oplot, ub, dzphot, ps=8, color='red'

; now plot individual spectra and chi2 curves

;   check = lindgen(ngal) & ncheck = ngal
    check = where((abs(deltazphot) gt dzcut) and $
      (rmag lt faintcut) and (snr gt snrcut),ncheck)
    srt = reverse(sort(abs(deltazphot[check])))
    check = check[srt]
    
;   jj = 5L & thisindex = 127
;   plot, oned_photoz[these_photoz[check[jj]]].zgrid_gal, $
;     oned_photoz[these_photoz[check[jj]]].chi2_gal, ysty=3

    plotscale = 1E17
    light = 2.99792458D18
    weff = photoinfo[0].weff
    filters = strtrim(photoinfo[0].filterlist,2)
    fnu2flam = 10^(-0.4*48.6)*light/weff^2.0

; pull out the primus spectra and best-fitting model (at the first
; minimum)
    flux = primus2flambda(extract[check],wave=wave,$
      counts2flam=counts2flam)
    modelflux = flux*0.0
    modelflux[*,0,*] = oned[check].fmod_photoz_gal1[*,0]/$
      (counts2flam[*,0,*]+(counts2flam[*,0,*] eq 0.0))*(counts2flam[*,0,*] ne 0.0)
    modelflux[*,1,*] = oned[check].fmod_photoz_gal2[*,0]/$
      (counts2flam[*,1,*]+(counts2flam[*,1,*] eq 0.0))*(counts2flam[*,1,*] ne 0.0)

; true high-resolution model at the true redshift

    npix = n_elements(extract[0].restwave)
    trueflux = fltarr(npix,ncheck)
    truewave = fltarr(npix,ncheck)

    for jj = 0L, ncheck-1L do begin
       trueflux[*,jj] = extract[check[jj]].restflux/(1+ztrue[check[jj]])
       truewave[*,jj] = extract[check[jj]].restwave*(1+ztrue[check[jj]])
    endfor
    
; best-fitting spectroscopic model at the spectroscopic redshift

    npix = n_elements(extract[0].restwave)
    tempspecflux = fltarr(npix,ncheck)
    tempspecwave = fltarr(npix,ncheck)

    for jj = 0L, ncheck-1L do begin
       tempspecflux[*,jj] = tempspec[check[jj]].highflux/(1+zspec[check[jj]])
       tempspecwave[*,jj] = tempspec[check[jj]].highwave*(1+zspec[check[jj]])
    endfor
    
; best-fitting photoz model at the photoz redshift

    npix = n_elements(extract[0].restwave)
    tempphotflux = fltarr(npix,ncheck)
    tempphotwave = fltarr(npix,ncheck)

    for jj = 0L, ncheck-1L do begin
       tempphotflux[*,jj] = tempphot[check[jj]].highflux/(1+zphot[check[jj]])
       tempphotwave[*,jj] = tempphot[check[jj]].highwave*(1+zphot[check[jj]])
    endfor
    
;   dfpsclose & im_window, 0, xratio=0.7

;   for jj = 0L, 1L do begin
    for jj = 0L, ncheck-1L do begin

; calculate some preliminaries and the xyranges
       truenorm = photoinfo[check[jj]].maggies[refband]/$
         (k_project_filters(k_lambda_to_edges(truewave[*,jj]),$
         trueflux[*,jj],filterlist=filters,/silent))[refband]
       tempspecnorm = photoinfo[check[jj]].maggies[refband]/$
         (k_project_filters(k_lambda_to_edges(tempspecwave[*,jj]),$
         tempspecflux[*,jj],filterlist=filters,/silent))[refband]
       tempphotnorm = photoinfo[check[jj]].maggies[refband]/$
         (k_project_filters(k_lambda_to_edges(tempphotwave[*,jj]),$
         tempphotflux[*,jj],filterlist=filters,/silent))[refband]
       norm1 = photoinfo[check[jj]].maggies[refband]/$
         (k_project_filters(k_lambda_to_edges(wave[*,0,jj]),$
         flux[*,0,jj],filterlist=filters,/silent))[refband]
       norm2 = photoinfo[check[jj]].maggies[refband]/$
         (k_project_filters(k_lambda_to_edges(wave[*,1,jj]),$
         flux[*,1,jj],filterlist=filters,/silent))[refband]
       modelnorm1 = photoinfo[check[jj]].maggies[refband]/$
         (k_project_filters(k_lambda_to_edges(wave[*,0,jj]),$
         modelflux[*,0,jj],filterlist=filters,/silent))[refband]
       modelnorm2 = photoinfo[check[jj]].maggies[refband]/$
         (k_project_filters(k_lambda_to_edges(wave[*,1,jj]),$
         modelflux[*,1,jj],filterlist=filters,/silent))[refband]
       yrange = plotscale*[$
         im_min(tempphotnorm*tempphotflux[*,jj],sigrej=3.0,/ignorezero,/nan),$
         im_max(tempphotnorm*tempphotflux[*,jj],sigrej=3.0,/nan)]
;      yrange = plotscale*[$
;        im_min(norm1*flux[*,0,jj],sigrej=4.0,/ignorezero,/nan)<$
;        im_min(norm2*flux[*,1,jj],sigrej=4.0,/ignorezero,/nan),$
;        im_max(norm1*flux[*,0,jj],sigrej=4.0,/nan)>$
;        im_max(norm2*flux[*,1,jj],sigrej=4.0,/nan)]
;      yrange = plotscale*[$
;        im_min(modelnorm1*modelflux[*,0,jj],sigrej=4.0,/ignorezero)<$
;        im_min(modelnorm2*modelflux[*,1,jj],sigrej=4.0,/ignorezero),$
;        im_max(modelnorm1*modelflux[*,0,jj],sigrej=4.0)>$
;        im_max(modelnorm2*modelflux[*,1,jj],sigrej=4.0)]
       xrange = [3300,9900]
; now make the plot
       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
         xthick=postthick1, ythick=postthick1, xsty=1, ysty=3, $
         thick=postthick1, xtitle='', xtickname=replicate(' ',10), $
         ytitle='', $
;        ytitle='Flux (10^{-17} '+flam_units()+')', $
         charsize=1.3, charthick=postthick2, position=pos3[*,0]
;        title='z_{true} = '+strtrim(string(ztrue[check[jj]],format='(F12.4)'),2)+$
;        ', R = '+strtrim(string(rmag[check[jj]],format='(F12.2)'),2)+$
;        ', <S/N> = '+strtrim(string(snr[check[jj]],format='(F12.2)'),2)
       legend, textoidl([$
         'z_{true} = '+strtrim(string(ztrue[check[jj]],format='(F12.4)'),2),$
;        'R = '+strtrim(string(rmag[check[jj]],format='(F12.2)'),2),$
         '<S/N> = '+strtrim(string(snr[check[jj]],format='(F12.2)'),2)]), $
         /right, /bottom, clear=0, box=0, charsize=1.3, charthick=postthick2
; true spectrum
       djs_oplot, truewave[*,jj], plotscale*trueflux[*,jj]*truenorm, $
         thick=postthick2, color=fsc_color('blue',100)
       plotsym, 8, 2.0, fill=1, thick=postthick3
       oploterror, weff, plotscale*photoinfo[check[jj]].maggies*fnu2flam, $
         plotscale*(1.0/sqrt(photoinfo[check[jj]].maggiesivar))*fnu2flam, $
         ps=8, color=djs_icolor('dark green'), thick=postthick1, $
         errcolor=djs_icolor('dark green'), errthick=postthick1
; the data
       djs_plot, [0], [0], /nodata, /noerase, xrange=xrange, yrange=yrange, $
         xthick=postthick1, ythick=postthick1, xsty=1, ysty=3, $
         thick=postthick1, xtitle='Observed Wavelength (\AA)', $
         ytitle='', $
;        ytitle='Flux (10^{-17} '+flam_units()+')', $
         charsize=1.3, charthick=postthick2, position=pos3[*,1]
       djs_oplot, tempphotwave[*,jj], plotscale*tempphotflux[*,jj]*tempphotnorm, $
         thick=postthick2, color=fsc_color('slate grey',101)
; PRIMUS spectra
       djs_oplot, wave[*,0,jj], plotscale*norm1*flux[*,0,jj], line=1, $
         color='purple', thick=postthick1
       djs_oplot, wave[*,1,jj], plotscale*norm2*flux[*,1,jj], line=1, $
         color='red', thick=postthick1
       legend, textoidl(['z_{'+photosuffix+'} = '+strtrim(string(zphot[check[jj]],$
         format='(F12.4)'),2)]), /right, /bottom, clear=0, $
         box=0, charsize=1.3, charthick=postthick2
;      legend, textoidl(['z_{'+photosuffix+'} = '+strtrim(string(zphot[check[jj]],$
;        format='(F12.4)'),2),'z_{PRIMUS} = '+strtrim(string(zspec[check[jj]],$
;        format='(F12.4)'),2)]), /right, /bottom, clear=0, $
;        box=0, charsize=1.3, charthick=postthick2
; best-fitting model spectra
       djs_oplot, wave[*,0,jj], plotscale*modelnorm1*modelflux[*,0,jj], line=0, $
         thick=postthick1
       djs_oplot, wave[*,1,jj], plotscale*modelnorm2*modelflux[*,1,jj], line=0, $
         thick=postthick1
; photometry          
       plotsym, 8, 2.0, fill=1, thick=postthick3
       oploterror, weff, plotscale*photoinfo[check[jj]].maggies*fnu2flam, $
         plotscale*(1.0/sqrt(photoinfo[check[jj]].maggiesivar))*fnu2flam, $
         ps=8, color=djs_icolor('dark green'), thick=postthick1, $
         errcolor=djs_icolor('dark green'), errthick=postthick1
;      djs_oplot, weff, plotscale*photoinfo[these[check[jj]]].maggies*fnu2flam, $
;        ps=8, color=djs_icolor('dark green'), thick=postthick1
; P(z) plots
;         pzspec = exp(-0.5*(oned[these[check[jj]]].chi2_gal-chi2spec[check[jj]])
;         pzspec = pzspec/total(pzspec)
;         yrange = [min(oned[these_photoz[check[jj]]].chi2_gal)<$
;           min(oned[these[check[jj]]].chi2_gal),$
;           max(oned[these_photoz[check[jj]]].chi2_gal)>$
;           max(oned[these[check[jj]]].chi2_gal)]
;         yrange = yrange*[0.9,1.05]
;         djs_plot, [0], [0], /nodata, xrange=[0.0,1.19], yrange=yrange, $
;           xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
;           thick=postthick1, xtitle='Redshift', ytitle='\chi^{2}', $
;           charsize=1.8, charthick=postthick2, position=pos1[*,0]
;         djs_oplot, zphot[check[jj]]*[1,1], !y.crange, line=0, color='red', thick=postthick3
;         djs_oplot, zspec[check[jj]]*[1,1], !y.crange, line=5, color='dark green', thick=postthick3
;         djs_oplot, ztrue[check[jj]]*[1,1], !y.crange, line=3, color='blue', thick=postthick3
;         djs_oplot, oned[these_photoz[check[jj]]].zgrid_gal, $
;           oned[these_photoz[check[jj]]].chi2_gal, $
;           line=0, ps=10, thick=postthick1, color='purple'
;         djs_oplot, oned[these[check[jj]]].zgrid_gal, oned[these[check[jj]]].chi2_gal, $
;           line=2, ps=10, thick=postthick1
; ytitle

       xyouts, pos3[0,0]-0.07, pos3[1,0], textoidl('Flux (10^{-17} ')+flam_units()+')', $
         align=0.5, orientation=90.0, charsize=1.3, charthick=postthick2, /normal
       
; chi2 plots          

       yrange = [min(oned[check[jj]].chi2_photoz_gal)<$
         min(oned[check[jj]].chi2_gal),$
         max(oned[check[jj]].chi2_photoz_gal)>$
         max(oned[check[jj]].chi2_gal)]
       yrange = yrange*[0.9,1.05]

       im_lineid_plot, oned[check[jj]].zgrid_gal, oned[check[jj]].chi2_gal, /noerase, $
;      im_lineid_plot, [0], [0], /nodata, /noerase, $
         [ztrue[check[jj]],zspec[check[jj]],zphot[check[jj]]], $
         textoidl(['z_{true}','z_{PRIMUS}','z_{'+photosuffix+'}']), $
         xrange=[0.0,1.19], yrange=yrange, $
         xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
         thick=postthick1, xtitle='Redshift', ytitle=textoidl('\chi^{2}'), $
         charsize=1.5, charthick=postthick2, position=pos3[*,2], $
         lcharsize=1.5, lcharthick=postthick2, extend_thick=postthick1, $
         line=2, extend_color=fsc_color(['blue','forest green','firebrick'],[90,91,92]), $
         label_color=fsc_color(['blue','forest green','firebrick'],[90,91,92])
       djs_oplot, oned[check[jj]].zgrid_gal, $
         oned[check[jj]].chi2_photoz_gal, $
         line=0, ps=10, thick=postthick1, color='purple'
;      djs_oplot, oned[these[check[jj]]].zgrid_gal, oned[these[check[jj]]].chi2_gal, $
;        line=2, ps=10, thick=postthick1
         
;      djs_plot, [0], [0], /nodata, /noerase, xrange=[0.0,1.19], yrange=yrange, $
;        xthick=postthick1, ythick=postthick1, xsty=1, ysty=1, $
;        thick=postthick1, xtitle='Redshift', ytitle='\chi^{2}', $
;        charsize=1.5, charthick=postthick2, position=pos3[*,1]
;      djs_oplot, zphot[check[jj]]*[1,1], !y.crange*[1.0,0.5], line=0, color='red', thick=postthick3
;      djs_oplot, zspec[check[jj]]*[1,1], !y.crange*[1.0,0.5], line=5, color='dark green', thick=postthick3
;      djs_oplot, ztrue[check[jj]]*[1,1], !y.crange*[1.0,0.5], line=3, color='blue', thick=postthick3

;      legend, textoidl(['z_{sim}','z_{PRIMUS}','z_{'+photosuffix+'}']), /left, /top, $
;        box=0, charsize=1.5, charthick=postthick2, line=[3,5,0], $
;        color=djs_icolor(['blue','dark green','red']), thick=postthick1, $
;        /clear
;      legend, textoidl(['\chi^{2}_{ PRIMUS}','\chi^{2}_{ '+photosuffix+'}']), $
;        /right, /bottom, box=0, charsize=1.3, charthick=postthick2, $
;        line=[2,0], thick=postthick1, clear=1, color=djs_icolor(['default','purple'])

    endfor

return
end
