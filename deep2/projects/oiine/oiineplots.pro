pro oiineplots, keynote=keynote
; jm09apr01nyu - 

    oiinepath = deep2_path(/projects)+'oiine/'
    idlpath = getenv('HOME')+'/idl/projects/deep2/projects/oiine/'
    pspath = deep2_path(/papers)+'oiine/FIG_OIINE/'
    if keyword_set(keynote) then begin
       pspath = deep2_path(/papers)+'oiine/keynote/'
    endif

    encap = 1 ; default is to make EPS files
    if keyword_set(encap) then begin
       ps = 1
       suffix = '.eps' 
    endif else suffix = '.ps'

    deepkcorr = read_oiine_sample(/ancillary)
    deepispec = read_oiine_sample(/ispec)
    sdsskcorr = read_oiine_sample(/ancillary,/sdss)
    sdssispec = read_oiine_sample(/ispec,/sdss)

; ---------------------------------------------------------------------------    
; redshift vs electron density

    if keyword_set(ps) then psfile = pspath+'z_vs_ne'+suffix
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, $
      xmargin=[1.8,0.2], ymargin=[1.0,1.1], width=6.5, $
      height=5.0, thick=8.0, keynote=keynote

; read the output from OIINESTATS
    ssstats = mrdfits(oiinepath+'oiine_stats_sdss.fits.gz',1)
    ddstats = mrdfits(oiinepath+'oiine_stats_deep.fits.gz',1)

    xtitle = 'Redshift'
    ytitle = 'Electron Density (cm^{-3})'
    xrange = [-0.05,2.8]
    yrange = [25,3000]

    ddpsym = 9 & ddcolor = 'dodger blue' & ddline = 5 & ddsymsize = 3.0
    sspsym = 4 & sscolor = 'orange' & ssline = 0 & sssymsize = 4.0
    
    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, $
      xsty=9, ysty=1, xrange=xrange, yrange=yrange, position=pos, /ylog
    timelabel1 = [1.0,3.0,5.0,7.0,9.0,11.0] ; [Gyr]
    axis, /xaxis, xsty=1, xrange=xrange, xtitle='Lookback Time (Gyr)', $
      xtickv=getredshift(getage(0.0)-timelabel1), xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)')
; local starbursts
    polyfill, [0.0,0.,0.05,0.05], [300.0,1200.0,1200.0,300.0], $
      color=fsc_color('firebrick',10)
;   im_oplot_box, 0.01, 900.0, 0.0, xoffset=0.05, yoffset=750.0
    xyouts, 0.15, 600, ' Local!cStarbursts', charsize=1.8, $
      align=0.5, orientation=90
; DEEP2
;   djs_oplot, ddstats.z, ddstats.ne_mean, $
;     psym=symcat(ddpsym), symsize=ddsymsize, color=fsc_color(ddcolor,91)
    oploterror, ddstats.z, ddstats.ne_p50, ddstats.zerr, ddstats.ne_p50-ddstats.ne_p16, $
      psym=symcat(ddpsym,thick=10), symsize=ddsymsize, color=fsc_color(ddcolor,91), $
      errcolor=fsc_color(ddcolor,91), line=ddline, /lobar
    oploterror, ddstats.z, ddstats.ne_p50, ddstats.zerr, ddstats.ne_p84-ddstats.ne_p50, $
      psym=symcat(ddpsym,thick=10), symsize=ddsymsize, color=fsc_color(ddcolor,91), $
      errcolor=fsc_color(ddcolor,91), line=ddline, /hibar
    xyouts, 1.1, 750.0, 'DEEP2', charsize=1.8, align=0.5
; SDSS
    oploterror, ssstats.z, ssstats.ne_p50, ssstats.zerr, ssstats.ne_p50-ssstats.ne_p16, $
      psym=symcat(sspsym,thick=10), symsize=sssymsize, color=fsc_color(sscolor,92), $
      errcolor=fsc_color(sscolor,92), line=ssline, /lobar
    oploterror, ssstats.z, ssstats.ne_p50, ssstats.zerr, ssstats.ne_p84-ssstats.ne_p50, $
      psym=symcat(sspsym,thick=10), symsize=sssymsize, color=fsc_color(sscolor,92), $
      errcolor=fsc_color(sscolor,92), line=ssline, /hibar
    xyouts, 0.15, 80.0, 'SDSS', charsize=1.8, align=0.0
; Lehnert+09
    lz = [2.17,2.21,2.26,2.38,1.51]
    lne = [1200,400,1200,400,260]
    lne_lo = [400,300,400,300,120]
    lne_hi = [700,700,700,700,150]
    lcolor = 'navy' & lpsym = 6
    oploterror, lz, lne, lz*0.0, lne_lo, psym=symcat(lpsym,thick=10), symsize=3.0, $
      color=fsc_color(lcolor,93), errcolor=fsc_color(lcolor,93), /lobar
    oploterror, lz, lne, lz*0.0, lne_lo, psym=symcat(lpsym,thick=10), symsize=3.0, $
      color=fsc_color(lcolor,93), errcolor=fsc_color(lcolor,93), /hibar
    xyouts, 2.1, 600.0, 'LBGs', charsize=1.8, align=1.0
; Nesvadba+07
    nz = [2.57]
    nne = [400]
    nne_lo = [220]
    nne_hi = [500]
    ncolor = 'green' & npsym = 5
    oploterror, nz, nne, nz*0.0, nne_lo, psym=symcat(npsym,thick=10), symsize=3.0, $
      color=fsc_color(ncolor,94), errcolor=fsc_color(ncolor,94), /lobar
    oploterror, nz, nne, nz*0.0, nne_lo, psym=symcat(npsym,thick=10), symsize=3.0, $
      color=fsc_color(ncolor,94), errcolor=fsc_color(ncolor,94), /hibar
    xyouts, 2.57, 130, 'SMG', charsize=1.8, align=0.5

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,suffix,'.pdf'), /sh
       rmfile, psfile
    endif

stop
    
; ---------------------------------------------------------------------------    
; absolute magnitude versus redshift - SDSS

    if keyword_set(ps) then psfile = pspath+'sdss_z_vs_mg'+suffix
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, xmargin=[1.3,0.2]

    xtitle = 'Redshift'
    ytitle = 'M_{0.1g}'
    xrange = [0.03,0.17]
    yrange = [-17.5,-22.9]

    indx = where(sdssispec.dens_edge eq 0)
    xx = sdsskcorr[indx].z
    yy = sdsskcorr[indx].ugriz_absmag[1]
    
    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, position=pos
    djs_oplot, xx, yy, psym=6, symsize=0.12

; overplot the sample definitions    
    par = yanny_readone(idlpath+'sdss_samples.par')
    for ii = 0, n_elements(par)-1 do begin
       djs_oplot, [par[ii].zmin,par[ii].zmax], par[ii].faintcut*[1,1], $
         line=0, color='red', thick=6.0
       djs_oplot, [par[ii].zmin,par[ii].zmax], par[ii].brightcut*[1,1], $
         line=0, color='red', thick=6.0
       djs_oplot, par[ii].zmin*[1,1], [par[ii].faintcut,par[ii].brightcut], $
         line=0, color='red', thick=6.0
       djs_oplot, par[ii].zmax*[1,1], [par[ii].faintcut,par[ii].brightcut], $
         line=0, color='red', thick=6.0
    endfor

    im_plotconfig, /psclose
    
; ---------------------------------------------------------------------------    
; absolute magnitude versus redshift - DEEP2

    if keyword_set(ps) then psfile = pspath+'deep2_z_vs_mg'+suffix
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, xmargin=[1.3,0.2]

    xtitle = 'Redshift'
    ytitle = 'M_{0.1g}'
    xrange = [0.6,1.55]
    yrange = [-18.5,-23.9]

    indx = where(deepispec.dens_edge eq 0)
    xx = deepkcorr[indx].z
    yy = deepkcorr[indx].ugriz_absmag[1]
    
    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, position=pos
    djs_oplot, xx, yy, psym=6, symsize=0.12
    
; overplot the sample definitions    
    par = yanny_readone(idlpath+'deep2_samples.par')
    for ii = 0, n_elements(par)-1 do begin
       djs_oplot, [par[ii].zmin,par[ii].zmax], par[ii].faintcut*[1,1], $
         line=0, color='red', thick=6.0
       djs_oplot, [par[ii].zmin,par[ii].zmax], par[ii].brightcut*[1,1], $
         line=0, color='red', thick=6.0
       djs_oplot, par[ii].zmin*[1,1], [par[ii].faintcut,par[ii].brightcut], $
         line=0, color='red', thick=6.0
       djs_oplot, par[ii].zmax*[1,1], [par[ii].faintcut,par[ii].brightcut], $
         line=0, color='red', thick=6.0
    endfor

    im_plotconfig, /psclose
    
return
end
    
