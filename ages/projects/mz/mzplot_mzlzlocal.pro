pro mzplot_mzlzlocal, ps=ps
; jm09mar19nyu - make various plots pertaining to the local LZ and MZ
; relations
; jm10oct10ucsd - major update    

    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(ps) then begin
       pspath = qapath
       suffix = '.ps'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse
    
; read the data    
    sdssancillary = read_mz_sample(/sdss,/mzhii_ancillary)
    sdssmass = read_mz_sample(/sdss,/mzhii_mass)
    sdssohdust = read_mz_sample(/sdss,/mzhii_log12oh)
    sdssohnodust = read_mz_sample(/sdss,/mzhii_log12oh,/nodust)
;   sdssispec = read_mz_sample(/sdss,/mzhii_ispec)

    agesancillary = read_mz_sample(/mzhii_ancillary)
    agesmass = read_mz_sample(/mzhii_mass)
    agesohdust = read_mz_sample(/mzhii_log12oh)
    agesohnodust = read_mz_sample(/mzhii_log12oh,/nodust)

    massaxis1 = im_array(8.8,11.3,0.01)

    smzbinsize = 0.1
    amzbinsize = 0.1
    slzbinsize = 0.2
    alzbinsize = 0.2

; ---------------------------------------------------------------------------    
; QAplot - compare the local SDSS MZ relations based on the various 
; calibrations; just plot the best-fitting results, not the data
    mzlocal_ews_t04 = mrdfits(mzpath+'mzlocal_sdss_ews_t04.fits.gz',1)
    mzlocal_cor_t04 = mrdfits(mzpath+'mzlocal_sdss_fluxcor_t04.fits.gz',1)

    mzlocal_ews_m91 = mrdfits(mzpath+'mzlocal_sdss_ews_m91.fits.gz',1)
    mzlocal_cor_m91 = mrdfits(mzpath+'mzlocal_sdss_fluxcor_m91.fits.gz',1)

    mzlocal_ews_kk04 = mrdfits(mzpath+'mzlocal_sdss_ews_kk04.fits.gz',1)
    mzlocal_cor_kk04 = mrdfits(mzpath+'mzlocal_sdss_fluxcor_kk04.fits.gz',1)

    mzlocal_mpa = mrdfits(mzpath+'mzlocal_sdss_fluxcor_mpajhu.fits.gz',1)
    
    massrange1 = [8.8,11.6]
    massaxis1 = range(9.0,11.4,500)
    massaxis2 = range(9.0,11.4,500)

    t04color = 'red'    & t04line = 3
    m91color = 'black'  & m91line = 0
    kk04color = 'blue'  & kk04line = 5
    tremonticolor = 'orange' & tremontiline = 4

    psfile = qapath+'mzlocal_compare.ps'
    im_plotconfig, 6, pos, psfile=psfile, height=[4.3,2.8], $
      xmargin=[1.4,0.4], width=6.7

    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xrange=massrange1, yrange=[8.52,9.2], xtickname=replicate(' ',10), $
      ytitle=mzplot_ohtitle()
;   djs_oplot, massaxis2, tremonti_mz(massaxis=massaxis2), $
;     color=fsc_color(tremonticolor,101), line=tremontiline, $
;     thick=8

; EWs    
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal_ews_t04.coeff), $
      color=fsc_color(t04color,101), line=t04line, thick=8
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal_ews_m91.coeff), $
      color=fsc_color(m91color,101), line=m91line, thick=8
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal_ews_kk04.coeff), $
      color=fsc_color(kk04color,101), line=kk04line, thick=8

; cor    
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal_cor_t04.coeff), $
      color=fsc_color(t04color,101), line=t04line, thick=2
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal_cor_m91.coeff), $
      color=fsc_color(m91color,101), line=m91line, thick=2
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal_cor_kk04.coeff), $
      color=fsc_color(kk04color,101), line=kk04line, thick=2
    
    im_legend, ['T04','M91','KK04'], /left, /top, box=0, margin=0, charsize=1.4, $
      line=[t04line,m91line,kk04line], color=[t04color,m91color,kk04color], $
      pspacing=1.9, thick=10
    legend, ['Thick - EWs','Thin - Fluxcor'], /right, /bottom, box=0

; residuals
    djs_plot, [0], [0], /nodata, position=pos[*,1], /noerase, xsty=1, ysty=1, $
      xrange=massrange1, yrange=0.08*[-1,1], xtitle=mzplot_masstitle(), $
      ytitle='Residuals [Fluxcor-EWs]'
    djs_oplot, !x.crange, [0,0], color='grey', line=1
    
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal_cor_t04.coeff)-$
      mz_closedbox(massaxis1,mzlocal_ews_t04.coeff), $
      color=fsc_color(t04color,101), line=t04line, thick=8
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal_cor_m91.coeff)-$
      mz_closedbox(massaxis1,mzlocal_ews_m91.coeff), $
      color=fsc_color(m91color,101), line=m91line, thick=8
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal_cor_kk04.coeff)-$
      mz_closedbox(massaxis1,mzlocal_ews_kk04.coeff), $
      color=fsc_color(kk04color,101), line=kk04line, thick=8
    
    im_plotconfig, /psclose, psfile=psfile, /gzip

stop    
    
; ---------------------------------------------------------------------------    
; paper plot - SDSS LZ relation showing four separate calibrations
    magrange = [-17,-22.7]
    ohrange1 = [8.3,9.3]

    sdsslevels = [0.25,0.5,0.75,0.9,0.95]
    outcolor = 'medium grey'

    ohtitle1 = mzplot_ohtitle(/fluxcor)

    kk04color = 'firebrick'   & kk04psym = 9 & kk04symsize = 0.8 & kk04line = 5
    m91color = 'dodger blue'  & m91psym = 6  & m91symsize  = 0.8 & m91line = 3
    t04color = 'forest green' & t04psym = 5  & t04symsize  = 0.8 & t04line = 4
    mpacolor = 'black'    & mpapsym = 4  & mpasymsize  = 0.8 & mpaline = 0
    tremonticolor = 'purple'    & tremontiline = 2
    
    psfile = pspath+'lzlocal_sdss_calibcompare'+suffix
    im_plotconfig, 5, pos, psfile=psfile, charsize=1.6;, width=[4.5,4.5]

    errcut = 0.2

; reference LZ relation
    refline = mpaline
    refcolor = 'black'
    lzlocalref = mrdfits(mzpath+'lzlocal_sdss_fluxcor_mpajhu.fits.gz',1)

; for the plot, apply an error cut so that the contours get rendered
; correctly     
    
; KK04    
    lzlocal = mrdfits(mzpath+'lzlocal_sdss_fluxcor_kk04.fits.gz',1)

    sinfo = mzlz_grab_info(sdssohnodust,sdssancillary,sdssmass,/kk04,$
      /flux,/nolimit,/errcut)
    mzplot_scatterplot, /sdss, sinfo.mb_ab, sinfo.oh, weight=sinfo.weight, $
      position=pos[*,0], xstyle=1, ystyle=1, xtitle='', xtickname=replicate(' ',10), $
      ytitle=ohtitle1, xrange=magrange, yrange=ohrange1, $
      levels=sdsslevels, outcolor=fsc_color(outcolor,101), xtickinterval=1, $
      ccolor=djs_icolor('grey'), /nogrey
    oplot_lzfit, lzlocal.coeff, linestyle=kk04line, linecolor=kk04color
    oplot_lzfit, lzlocalref.coeff, linestyle=refline, linecolor=refcolor
    im_legend, ['KK04'], /right, /bottom, box=0, charsize=1.3, margin=0, $
      color=kk04color, line=kk04line, pspacing=1.9, thick=10

; M91
    lzlocal = mrdfits(mzpath+'lzlocal_sdss_fluxcor_m91.fits.gz',1)

    sinfo = mzlz_grab_info(sdssohnodust,sdssancillary,sdssmass,/m91,/flux,/nolimit,/errcut)
    mzplot_scatterplot, /sdss, sinfo.mb_ab, sinfo.oh, weight=sinfo.weight, $
      position=pos[*,1], /noerase, xstyle=1, ystyle=1, xtitle='', xtickname=replicate(' ',10), $
      ytitle='', ytickname=replicate(' ',10), xrange=magrange, yrange=ohrange1, $
      levels=sdsslevels, outcolor=fsc_color(outcolor,101), xtickinterval=1, $
      ccolor=djs_icolor('grey'), /nogrey
    oplot_lzfit, lzlocal.coeff, linestyle=m91line, linecolor=m91color
    oplot_lzfit, lzlocalref.coeff, linestyle=refline, linecolor=refcolor
    im_legend, ['M91'], /right, /bottom, box=0, charsize=1.3, margin=0, $
      color=m91color, line=m91line, pspacing=1.9, thick=10

; T04
    lzlocal = mrdfits(mzpath+'lzlocal_sdss_fluxcor_t04.fits.gz',1)

    sinfo = mzlz_grab_info(sdssohnodust,sdssancillary,sdssmass,/t04,/flux,/nolimit,/errcut)
    mzplot_scatterplot, /sdss, sinfo.mb_ab, sinfo.oh, weight=sinfo.weight, $
      position=pos[*,2], /noerase, xstyle=1, ystyle=1, xtitle=mzplot_mbtitle(), $
      ytitle=ohtitle1, xrange=magrange, yrange=ohrange1, $
      levels=sdsslevels, outcolor=fsc_color(outcolor,101), xtickinterval=1, $
      ccolor=djs_icolor('grey'), /nogrey
    oplot_lzfit, lzlocalref.coeff, linestyle=refline, linecolor=refcolor
    oplot_lzfit, lzlocal.coeff, linestyle=t04line, linecolor=t04color
    im_legend, ['T04'], /right, /bottom, box=0, charsize=1.3, margin=0, $
      color=t04color, line=t04line, pspacing=1.9, thick=10

; MPA-JHU
    lzlocal = mrdfits(mzpath+'lzlocal_sdss_fluxcor_mpajhu.fits.gz',1)

    sinfo = mzlz_grab_info(sdssohnodust,sdssancillary,sdssmass,/mpa,/flux,/nolimit,/errcut)
    mzplot_scatterplot, /sdss, sinfo.mb_ab, sinfo.oh, weight=sinfo.weight, $
      position=pos[*,3], /noerase, xstyle=1, ystyle=1, xtitle=mzplot_mbtitle(), $
      ytitle='', ytickname=replicate(' ',10), xrange=magrange, yrange=ohrange1, $
      levels=sdsslevels, outcolor=fsc_color(outcolor,101), xtickinterval=1, $
      ccolor=djs_icolor('grey'), /nogrey
    oplot_lzfit, lzlocal.coeff, linestyle=mpaline, linecolor=mpacolor
    oplot_lzfit, lzlocalref.coeff, linestyle=refline, linecolor=refcolor
;   oplot_lzfit, linestyle=tremontiline, linecolor=tremonticolor, /tremonti
    im_legend, 'MPA-JHU', /right, /bottom, box=0, charsize=1.3, margin=0, $
      color=mpacolor, line=mpaline, pspacing=1.9, thick=10
;   im_legend, ['MPA-JHU','Tremonti+04'], /right, /bottom, box=0, charsize=1.3, margin=0, $
;     color=[mpacolor,tremonticolor], line=[mpaline,tremontiline], pspacing=1.9, thick=10

    im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)

stop    
    
; ---------------------------------------------------------------------------    
; paper plot - SDSS MZ relation showing four separate calibrations;
; see FIT_MZLZLOCAL for the adopted sample cuts
    massrange = [8.4,11.6]
;   ohrange1 = [8.4,9.3]
    ohrange1 = [8.3,9.3]
    verbose = 1
    sdsslevels = [0.25,0.5,0.75,0.9,0.95]
    outcolor = 'medium grey'

    ohtitle1 = mzplot_ohtitle(/fluxcor)

    kk04color = 'firebrick'   & kk04psym = 9 & kk04symsize = 0.8 & kk04line = 5
    m91color = 'dodger blue'  & m91psym = 6  & m91symsize  = 0.8 & m91line = 3
    t04color = 'forest green' & t04psym = 5  & t04symsize  = 0.8 & t04line = 4
    mpacolor = 'chocolate'    & mpapsym = 4  & mpasymsize  = 0.8 & mpaline = 0
    
    psfile = pspath+'mzlocal_sdss_calibcompare'+suffix
    im_plotconfig, 5, pos, psfile=psfile, charsize=1.6;, width=[4.5,4.5]

    errcut = 0.2

; reference MZ relation    
    refline = mpaline
    mzlocalref = mrdfits(mzpath+'mzlocal_sdss_fluxcor_mpajhu.fits.gz',1)
;   maxisref = range(mzlocalref.minmass,massrange[1]-0.2,50)

    maxis = range(8.75,11.4,75)

; for the plot, apply an error cut so that the contours get rendered
; correctly     
    
; KK04    
    mzlocal = mrdfits(mzpath+'mzlocal_sdss_fluxcor_kk04.fits.gz',1)
;   maxis = range(mzlocal.minmass,massrange[1]-0.2,50)

    sinfo = mzlz_grab_info(sdssohnodust,sdssancillary,sdssmass,/kk04,$
      /flux,/nolimit,/errcut)
;   show = where(sinfo.oh_err lt errcut,comp=toss)
    mzplot_scatterplot, /sdss, sinfo.mass, sinfo.oh, weight=sinfo.weight, $
      position=pos[*,0], xstyle=1, ystyle=1, xtitle='', xtickname=replicate(' ',10), $
      ytitle=ohtitle1, xrange=massrange, yrange=ohrange1, $
      levels=sdsslevels, outcolor=fsc_color(outcolor,101), xtickinterval=1, $
      ccolor=djs_icolor('grey'), /nogrey
;   djs_oplot, sinfo.mass[toss], sinfo.oh[toss], psym=6, sym=0.1, color='red'
    djs_oplot, mzlocal.bin_mass, mzlocal.bin_oh, psym=symcat(kk04psym,thick=7), $
      symsize=kk04symsize, thick=6, color=fsc_color(kk04color,101)
    djs_oplot, maxis, mz_closedbox(maxis,mzlocal.coeff), thick=7, $
      line=kk04line;, color=fsc_color(kk04color,101)
;   djs_oplot, mzlocal.bin_mass, mzlocal.bin_oh_mean, psym=6
    djs_oplot, maxis, mz_closedbox(maxis,mzlocalref.coeff), thick=7, $
      line=refline;, color=fsc_color(kk04color,101)
    im_legend, ['KK04'], /right, /bottom, box=0, charsize=1.3, margin=0, $
      psym=-symcat(kk04psym,thick=7), color=kk04color, symsize=kk04symsize*0.4, $
      line=kk04line, pspacing=1.9, thick=10

; M91
    mzlocal = mrdfits(mzpath+'mzlocal_sdss_fluxcor_m91.fits.gz',1)
;   maxis = range(mzlocal.minmass,massrange[1]-0.2,50)

    sinfo = mzlz_grab_info(sdssohnodust,sdssancillary,sdssmass,/m91,/flux,/nolimit,/errcut)
;   gd = where(sinfo.oh_err lt 0.1 and sinfo.ohlimit eq 0)
;   bad = where(sinfo.oh_err gt 0.1 and sinfo.ohlimit eq 0)
;   lim = where(sinfo.ohlimit eq 1)
;   sinfo.oh_err = sqrt(sinfo.oh_err^2+0.05^2)
;   mzplot_scatterplot, /sdss, sinfo.mass[gd], sinfo.oh[gd], weight=sinfo.weight[gd], $;/sinfo.oh_err^2, $
    mzplot_scatterplot, /sdss, sinfo.mass, sinfo.oh, weight=sinfo.weight, $
      position=pos[*,1], /noerase, xstyle=1, ystyle=1, xtitle='', xtickname=replicate(' ',10), $
      ytitle='', ytickname=replicate(' ',10), xrange=massrange, yrange=ohrange1, $
      levels=sdsslevels, outcolor=fsc_color(outcolor,101), xtickinterval=1, $
      ccolor=djs_icolor('grey'), /nogrey
;   djs_oplot, sinfo.mass[bad], sinfo.oh[bad], psym=3, color='blue'
;   djs_oplot, sinfo.mass[lim], sinfo.oh[lim], psym=3, color='dark green'

    djs_oplot, mzlocal.bin_mass, mzlocal.bin_oh, psym=symcat(m91psym,thick=7), $
      symsize=m91symsize, thick=6, color=fsc_color(m91color,101)
    djs_oplot, maxis, mz_closedbox(maxis,mzlocal.coeff), thick=7, $
      line=m91line;, color=fsc_color(m91color,101)
;   djs_oplot, mzlocal.bin_mass, mzlocal.bin_oh_mean, psym=6
    djs_oplot, maxis, mz_closedbox(maxis,mzlocalref.coeff), thick=7, $
      line=refline
    im_legend, ['M91'], /right, /bottom, box=0, charsize=1.3, margin=0, $
      psym=-symcat(m91psym,thick=7), color=m91color, symsize=m91symsize*0.4, $
      line=m91line, pspacing=1.9, thick=10

; T04
    mzlocal = mrdfits(mzpath+'mzlocal_sdss_fluxcor_t04.fits.gz',1)
;   maxis = range(mzlocal.minmass,massrange[1]-0.2,50)

    sinfo = mzlz_grab_info(sdssohnodust,sdssancillary,sdssmass,/t04,/flux,/nolimit,/errcut)
    mzplot_scatterplot, /sdss, sinfo.mass, sinfo.oh, weight=sinfo.weight, $
      position=pos[*,2], /noerase, xstyle=1, ystyle=1, xtitle=mzplot_masstitle(), $
      ytitle=ohtitle1, xrange=massrange, yrange=ohrange1, $
      levels=sdsslevels, outcolor=fsc_color(outcolor,101), xtickinterval=1, $
      ccolor=djs_icolor('grey'), /nogrey
    djs_oplot, mzlocal.bin_mass, mzlocal.bin_oh, psym=symcat(t04psym,thick=7), $
      symsize=t04symsize, thick=6, color=fsc_color(t04color,101)
    djs_oplot, maxis, mz_closedbox(maxis,mzlocal.coeff), thick=7, $
      line=t04line;, color=fsc_color(t04color,101)
;   djs_oplot, mzlocal.bin_mass, mzlocal.bin_oh_mean, psym=6
    djs_oplot, maxis, mz_closedbox(maxis,mzlocalref.coeff), thick=7, $
      line=refline;, color=fsc_color(kk04color,101)
    im_legend, ['T04'], /right, /bottom, box=0, charsize=1.3, margin=0, $
      psym=-symcat(t04psym,thick=7), color=t04color, symsize=t04symsize*0.4, $
      line=t04line, pspacing=1.9, thick=10

; MPA-JHU
    mzlocal = mrdfits(mzpath+'mzlocal_sdss_fluxcor_mpajhu.fits.gz',1)
;   maxis = range(mzlocal.minmass,massrange[1]-0.2,50)

    sinfo = mzlz_grab_info(sdssohnodust,sdssancillary,sdssmass,/mpa,/flux,/nolimit,/errcut)
    mzplot_scatterplot, /sdss, sinfo.mass, sinfo.oh, weight=sinfo.weight, $
      position=pos[*,3], /noerase, xstyle=1, ystyle=1, xtitle=mzplot_masstitle(), $
      ytitle='', ytickname=replicate(' ',10), xrange=massrange, yrange=ohrange1, $
      levels=sdsslevels, outcolor=fsc_color(outcolor,101), xtickinterval=1, $
      ccolor=djs_icolor('grey'), /nogrey
    djs_oplot, mzlocal.bin_mass, mzlocal.bin_oh, psym=symcat(mpapsym,thick=7), $
      symsize=mpasymsize, thick=6, color=fsc_color(mpacolor,101)
;   djs_oplot, mzlocal.bin_mass, mzlocal.bin_oh_mean, psym=6
    djs_oplot, maxis, mz_closedbox(maxis,mzlocal.coeff), thick=7, $
      line=mpaline;, color=fsc_color(mpacolor,101)
;   djs_oplot, maxis, mz_closedbox(maxis,mzlocalref.coeff), thick=7, $
;     line=refline;, color=fsc_color(kk04color,101)
    im_legend, 'MPA-JHU', /right, /bottom, box=0, charsize=1.3, margin=0, $
      psym=-symcat(mpapsym,thick=7), color=mpacolor, symsize=mpasymsize*0.4, $
      line=mpaline, pspacing=1.9, thick=10, position=[11.45,8.43], /data
    im_legend, 'Tremonti+04', /right, /bottom, box=0, charsize=1.3, margin=0, $
      color='purple', line=2, pspacing=1.9, thick=10, position=[11.45,8.36], /data
;   im_legend, ['MPA-JHU','Tremonti+04'], /right, /bottom, box=0, charsize=1.3, margin=0, $
;     psym=-[symcat(mpapsym,thick=7),8], color=[mpacolor,'purple'], symsize=mpasymsize*0.4, $
;     line=[mpaline,1], pspacing=1.9, thick=10

; overplot Tremonti+04
    tmassaxis = range(8.8,11.3,50) ; Kroupa+01 IMF
    tcoeff = [-1.492,1.847,-0.08026] ; T04 published
    djs_oplot, tmassaxis, poly(tmassaxis-0.04,tcoeff), line=2, $
      color=fsc_color('purple',101), thick=10

    im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)

; ---------------------------------------------------------------------------    
; paper plot - SDSS and low-redshift AGES LZ relations at z~0.1 for
; all three calibrations
    magrange1 = [-16.5,-22.7]
    verbose = 1
    magaxis = range(-17,-23,300)
    
    sdsslevels = [0.25,0.5,0.75,0.9,0.95]
    ageslevels = [0.25,0.5,0.75,0.9]

    sline = 0 & scolor = 'firebrick'
    aline = 5 & acolor = 'blue'
    
    for ii = 0, 2 do begin
       case ii of
          0: begin
             t04 = 1 & m91 = 0 & kk04 = 0
             calib = 't04'
             ohrange1 = [8.3,9.35]
          end
          1: begin
             t04 = 0 & m91 = 1 & kk04 = 0
             calib = 'm91'
             ohrange1 = [8.3,9.2]
          end
          2: begin
             t04 = 0 & m91 = 0 & kk04 = 1
             calib = 'kk04'
             ohrange1 = [8.55,9.35]
          end
       endcase
       ohtitle1 = mzplot_ohtitle(t04=t04,m91=m91,kk04=kk04,/fluxcor)

; best-fitting LZ relations
       sdss_lzfit = mrdfits(mzpath+'lzlocal_sdss_fluxcor_'+calib+'.fits.gz',1)
       ages_lzfit = mrdfits(mzpath+'lzlocal_ages_fluxcor_'+calib+'.fits.gz',1)
       splog, ages_lzfit.coeff
       
       sinfo = mzlz_grab_info(sdssohnodust,sdssancillary,sdssmass,$
         t04=t04,m91=m91,kk04=kk04,/nolimit,/flux,/errcut)
       ainfo = mzlz_grab_info(agesohnodust,agesancillary,agesmass,$
         t04=t04,m91=m91,kk04=kk04,zmin=0.05,zmax=0.15,/nolimit,/flux);,/errcut)

;      abin = im_medxbin(ainfo.mb_ab,ainfo.oh,alzbinsize,$
;        weight=ainfo.weight,minpts=15,verbose=verbose,$
;        minx=-23.0)
;      sbin = im_medxbin(sinfo.mb_ab,sinfo.oh,slzbinsize,$
;        weight=sinfo.weight,minpts=500,verbose=verbose,$
;        minx=-23.0)
       
       psfile = pspath+'lzlocal_'+calib+suffix
       im_plotconfig, 12, pos, psfile=psfile, charsize=1.9, $
         height=5.0, width=[4.5,4.5]
; LZ: SDSS
       mzplot_scatterplot, /sdss, sinfo.mb_ab, sinfo.oh, weight=sinfo.weight, $
         position=pos[*,0], xstyle=1, ystyle=1, xtitle=mzplot_mbtitle(), $
         ytitle=ohtitle1, xrange=magrange1, yrange=ohrange1, $
         levels=sdsslevels, $;ccolor=djs_icolor('grey'), $ ; /nogrey, $
         outcolor=fsc_color('medium grey',101)
;       oploterror, abin.xbin, abin.medy, abin.sigymean, psym=symcat(16,thick=7), $
;         symsize=1.1, thick=6, color=fsc_color('blue',101), $
;         errcolor=fsc_color('blue',101)
;       oploterror, sbin.xbin, sbin.medy, sbin.sigymean, psym=symcat(15,thick=7), $
;         symsize=1.1, thick=6, color=fsc_color('firebrick',101), $
;         errcolor=fsc_color('firebrick',101)
       oplot_lzfit, sdss_lzfit.coeff, linestyle=sline, linecolor=scolor
       oplot_lzfit, ages_lzfit.coeff, linestyle=aline, linecolor=acolor

       im_legend, ['SDSS - 0.033<z<0.25'], /left, /top, box=0, thick=8, $
         charsize=1.5, margin=0, line=sline, color=scolor, pspacing=1.9
; LZ: low-redshift AGES
       mzplot_scatterplot, /ages, ainfo.mb_ab, ainfo.oh, weight=ainfo.weight, $
         /noerase, position=pos[*,1], xstyle=1, ystyle=1, xtitle=mzplot_mbtitle(), $
         ytitle='', xrange=magrange1, yrange=ohrange1, ytickname=replicate(' ',10), $
         levels=ageslevels, npix=14;, ccolor=djs_icolor('grey');, /nogrey

;      sixlin, ainfo.mb_ab-lz_pivotmag(), ainfo.oh, a, $
;        siga, b, sigb, weight=ainfo.weight
;      ccoeff = [a[2],b[2]]
;      djs_oplot, magaxis, poly(magaxis-lz_pivotmag(),ccoeff), color='dark green', thick=4

;       oploterror, abin.xbin, abin.medy, abin.sigymean, psym=symcat(16,thick=7), $
;         symsize=1.1, thick=6, color=fsc_color('blue',101), $
;         errcolor=fsc_color('blue',101)
;       oploterror, sbin.xbin, sbin.medy, sbin.sigymean, psym=symcat(15,thick=7), $
;         symsize=1.1, thick=6, color=fsc_color('firebrick',101), $
;         errcolor=fsc_color('firebrick',101)
       oplot_lzfit, sdss_lzfit.coeff, linestyle=sline, linecolor=scolor
       oplot_lzfit, ages_lzfit.coeff, linestyle=aline, linecolor=acolor
       
       im_legend, ['AGES - 0.05<z<0.15'], /left, /top, box=0, thick=8, $
         charsize=1.5, margin=0, line=aline, color=acolor, pspacing=1.9
       im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)
    endfor

stop    
    
; ---------------------------------------------------------------------------    
; QAplot - local MZ relation using the MPA metallicities and compare
; with Tremonti et al. 2004
    mzlocal = mrdfits(mzpath+'mzlocal_sdss_fluxcor_t04.fits.gz',1)

    massrange1 = [7.9,11.7]
    ohrange1 = [8.0,9.45]
    massaxis1 = range(8.5,11.4,500)

    xtremonti = [8.57,8.67,8.76,8.86,8.96,9.06,9.16,9.26,$
      9.36,9.46,9.57,9.66,9.76,9.86,9.96,10.06,10.16,10.26,$
      10.36,10.46,10.56,10.66,10.76,10.86,10.95,11.05,11.15,11.25]
    ytremonti = [8.44,8.48,8.57,8.61,8.63,8.66,8.68,8.71,$
      8.74,8.78,8.82,8.84,8.87,8.90,8.94,8.97,8.99,9.01,$
      9.03,9.05,9.07,9.08,9.09,9.10,9.11,9.11,9.12,9.12]

    good = where(sdssancillary.mass_avg gt 0 and sdssancillary.oh_avg gt 0)
    mass = sdssancillary[good].mass_avg
    oh = sdssancillary[good].oh_avg
    sbin = im_medxbin(mass,oh,0.1,minpts=75,minx=8.8,/verbose)

    psfile = qapath+'mzlocal_mpa.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.8

    mzplot_scatterplot, mass, oh, position=pos, xstyle=1, ystyle=1, $
      xtitle=strjoin([mzplot_masstitle(),' [MPA/JHU]']), $
      ytitle=strjoin([mzplot_ohtitle(),' [MPA/JHU]']), $
      xrange=massrange1, yrange=ohrange1, levels=sdsslevels, /sdss
    djs_oplot, sbin.xbin, sbin.medy, psym=symcat(15,thick=7), $
      symsize=1.1, thick=6, color=fsc_color('dodger blue',101)
    djs_oplot, xtremonti, ytremonti, psym=symcat(9,thick=7), $
      symsize=1.1, thick=6, color=fsc_color('forest green',101)
    djs_oplot, massaxis1, tremonti_mz(massaxis=massaxis1), $
      color=fsc_color('firebrick',101), line=0
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal.coeff), $
      color=fsc_color('red',101), line=5, thick=8

    im_plotconfig, /psclose, psfile=psfile, /gzip
    
; ---------------------------------------------------------------------------    
; paper plot - SDSS and low-redshift AGES MZ relation at z~0.1 for all
; three calibrations
    massrange1 = [7.9,11.7]
;   ohrange1 = [8.25,9.45]
    verbose = 1
    
    sdsslevels = [0.25,0.5,0.75,0.9,0.95]
    ageslevels = [0.25,0.5,0.75,0.9]

    massaxis = range(8.8,11.4,500)
    
;   for ii = 2, 2 do begin
    for ii = 0, 2 do begin
       case ii of
          0: begin
             t04 = 1 & m91 = 0 & kk04 = 0
             ohrange1 = [8.4,9.3]
             calib = 't04'
             ohrange1 = [8.35,9.3]
          end
          1: begin
             t04 = 0 & m91 = 1 & kk04 = 0
             ohrange1 = [8.3,9.15]
             calib = 'm91'
             ohrange1 = [8.4,9.2]
          end
          2: begin
             t04 = 0 & m91 = 0 & kk04 = 1
             ohrange1 = [8.5,9.25]
             calib = 'kk04'
             ohrange1 = [8.4,9.3]
          end
       endcase
       ohtitle1 = mzplot_ohtitle(t04=t04,m91=m91,kk04=kk04,/fluxcor)
       
       sinfo = mzlz_grab_info(sdssohnodust,sdssancillary,sdssmass,$
         t04=t04,m91=m91,kk04=kk04,/nolimit,/flux,/errcut)
       ainfo = mzlz_grab_info(agesohnodust,agesancillary,agesmass,$
         t04=t04,m91=m91,kk04=kk04,zmin=0.05,zmax=0.15,/nolimit,/flux)

       abin = im_medxbin(ainfo.mass,ainfo.oh,amzbinsize,minpts=5,$
         weight=ainfo.weight/ainfo.oh_err^2,/verbose)

       mzlocal = mrdfits(mzpath+'mzlocal_sdss_fluxcor_'+calib+'.fits.gz',1)
       
       psfile = pspath+'mzlocal_'+calib+suffix
       im_plotconfig, 12, pos, psfile=psfile, charsize=1.9, $
         height=5.0, width=[4.5,4.5]
; SDSS
       mzplot_scatterplot, /sdss, sinfo.mass, sinfo.oh, weight=sinfo.weight, $
         position=pos[*,0], xstyle=1, ystyle=1, xtitle=mzplot_masstitle(), $
         ytitle=ohtitle1, xrange=massrange1, yrange=ohrange1, $
         levels=sdsslevels, ccolor=djs_icolor('grey'), /nogrey, $
         outcolor=fsc_color('medium grey',101)

       oploterror, abin.xbin, abin.meany, abin.sigymean, psym=symcat(16,thick=7), $
         symsize=1.1, thick=6, color=fsc_color('blue',101), $
         errcolor=fsc_color('blue',101)
       djs_oplot, mzlocal.bin_mass, mzlocal.bin_oh, $
         psym=symcat(6,thick=7), symsize=0.8, thick=6, color=fsc_color('firebrick',101)

;      maxis = range(8.75,11.4,75)
;      djs_oplot, maxis, mz_closedbox(maxis,mzlocal.coeff), thick=7

;      oploterror, mzlocal.bin_mass, mzlocal.bin_oh, mzlocal.bin_oh_err, $
;        psym=symcat(6,thick=3), symsize=0.5, thick=6, color=fsc_color('firebrick',101), $
;        errcolor=fsc_color('firebrick',101)
       
       im_legend, ['SDSS - 0.033<z<0.25'], /left, /top, box=0, $
         charsize=1.5, margin=0, psym=symcat(6,thick=7), color='firebrick'

; MZ: low-z AGES
       mzplot_scatterplot, /ages, ainfo.mass, ainfo.oh, weight=ainfo.weight, $
         /noerase, position=pos[*,1], xstyle=1, ystyle=1, xtitle=mzplot_masstitle(), $
         ytitle='', xrange=massrange1, yrange=ohrange1, ytickname=replicate(' ',10), $
         levels=ageslevels, npix=16, ccolor=djs_icolor('grey'), /nogrey

;      djs_oplot, ainfo.mass[lim], ainfo.oh[lim], weight=ainfo.weight[lim], psym=6, sym=0.1, color='yellow'
       djs_oplot, mzlocal.bin_mass, mzlocal.bin_oh, $
         psym=symcat(6,thick=7), symsize=0.8, thick=6, color=fsc_color('firebrick',101)
       oploterror, abin.xbin, abin.meany, abin.sigymean, psym=symcat(16,thick=7), $
         symsize=1.1, thick=6, color=fsc_color('blue',101), $
         errcolor=fsc_color('blue',101)
;      oploterror, mzlocal.bin_mass, mzlocal.bin_oh, mzlocal.bin_oh_err, $
;        psym=symcat(6,thick=5), symsize=1.1, thick=6, color=fsc_color('firebrick',101), $
;        errcolor=fsc_color('firebrick',101)

       im_legend, ['AGES - 0.05<z<0.15'], /left, /top, box=0, $
         charsize=1.5, margin=0, psym=symcat(16,thick=7), color='blue'

       im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)
    endfor

return
end


;; ---------------------------------------------------------------------------    
;; SDSS MZ relation (KK04 calibration) with and without the upper
;; limits 
;    all = mzlz_grab_info(sdssohdust,sdssancillary,sdssmass,/kk04)
;    nolim = mzlz_grab_info(sdssohdust,sdssancillary,sdssmass,/kk04,/nolimits)
;
;    allfit = fit_mz_closedbox(all.mass,all.oh,all.weight,oh_err=all.oh_err)
;    nolimfit = fit_mz_closedbox(nolim.mass,nolim.oh,nolim.weight,oh_err=nolim.oh_err)
;
;    ohrange = [8.4,9.3]
;    massrange = [8,12]
;    
;    psfile = qapath+'mzlocal_limits.ps'
;    im_plotconfig, 0, pos, psfile=psfile, height=5.8
;
;; everything    
;    mzplot_scatterplot, /sdss, all.mass, all.oh, weight=all.weight, $
;      position=pos, xstyle=1, ystyle=1, xtitle=mzplot_masstitle(), $
;      ytitle=mzplot_ohtitle(/kk04), xrange=massrange, yrange=ohrange
;    oploterror, allfit.bin_mass, allfit.bin_oh, allfit.bin_oh_err, psym=symcat(15,thick=7), $
;      symsize=1.1, thick=6, color=fsc_color('firebrick',101), $
;      errcolor=fsc_color('firebrick',101)
;    djs_oplot, massaxis1, mz_closedbox(massaxis1,allfit.coeff), line=0, color='orange', thick=4
;    djs_oplot, massaxis1, mz_closedbox(massaxis1,nolimfit.coeff), line=5, color='blue', thick=4
;    im_legend, ['Everything'], /left, /top, box=0, charsize=1.5, margin=0
;    
;; no upper limits
;    mzplot_scatterplot, /sdss, nolim.mass, nolim.oh, weight=nolim.weight, $
;      position=pos, xstyle=1, ystyle=1, xtitle=mzplot_masstitle(), $
;      ytitle=mzplot_ohtitle(/kk04), xrange=massrange, yrange=ohrange
;    oploterror, nolimfit.bin_mass, nolimfit.bin_oh, nolimfit.bin_oh_err, psym=symcat(15,thick=7), $
;      symsize=1.1, thick=6, color=fsc_color('firebrick',101), $
;      errcolor=fsc_color('firebrick',101)
;    djs_oplot, massaxis1, mz_closedbox(massaxis1,allfit.coeff), line=0, color='orange', thick=4
;    djs_oplot, massaxis1, mz_closedbox(massaxis1,nolimfit.coeff), line=5, color='blue', thick=4
;    im_legend, ['No Upper Limits'], /left, /top, box=0, charsize=1.5, margin=0
;    
;    im_plotconfig, /psclose, psfile=psfile, /gzip
;
;stop       
;       


;; ---------------------------------------------------------------------------    
;; compare the various ways of fitting the SDSS MZ relation:
;; polynomial, double power-law, and broken power-law (see
;; FIT_MZLZLOCAL); make the plot just using one calibration and EWs
;; (rather than reddening-corrected fluxes)
;    usefluxcor = 0
;    if usefluxcor then ff = '_fluxcor' else ff = ''
;    closed = mrdfits(mzpath+'mzlocal_sdss'+ff+'_closedbox.fits.gz',1)
;    polyfit = mrdfits(mzpath+'mzlocal_sdss'+ff+'_poly.fits.gz',1)
;    doublepl = mrdfits(mzpath+'mzlocal_sdss'+ff+'_doublepl.fits.gz',1)
;    brokenpl = mrdfits(mzpath+'mzlocal_sdss'+ff+'_brokenpl.fits.gz',1)
;
;;   calib = 'm91' &  t04 = 0 & m91 = 1 & kk04 = 0
;    calib = 'kk04' &  t04 = 0 & m91 = 0 & kk04 = 1
;;   calib = 't04' & t04 = 1 & m91 = 0 & kk04 = 0
;    indx = where(strtrim(polyfit.calib,2) eq calib)
;    if usefluxcor then begin
;       info = mzlz_grab_info(sdssohnodust,/flux,$
;         sdssancillary,sdssmass,t04=t04,m91=m91,kk04=kk04)
;    endif else begin
;       info = mzlz_grab_info(sdssohdust,$ ; EWs
;         sdssancillary,sdssmass,t04=t04,m91=m91,kk04=kk04)
;    endelse
;    ytitle = mzplot_ohtitle(t04=t04,m91=m91,kk04=kk04)
;
;    massrange1 = [8.1,11.9]
;;   ohrange1 = [8.45,9.3] ; KK04
;    ohrange1 = [8.25,9.3] ; T04
;    massaxis = range(8.4,11.6,500)
;
;    pcolor = 'red'   & pline = 3
;    bcolor = 'black' & bline = 4
;    dcolor = 'blue'  & dline = 5
;    ccolor = 'orange'  & cline = 0
;    
;    psfile = pspath+'mz_fitting'+suffix
;    im_plotconfig, 0, pos, psfile=psfile, height=5.8
;
;    mzplot_scatterplot, info.mass, info.oh, /sdss, position=pos[*,0], $
;      xsty=1, ysty=1, xrange=massrange1, yrange=ohrange1, $
;      xtitle=mzplot_masstitle(), ytitle=ytitle, $
;      /nogrey, ccolor=djs_icolor('grey')
;    legend, ['Quadratic Polynomial','Double Power-law',$
;      'Broken Power-law','Closed-Box Model'], $
;      /right, /bottom, box=0, margin=0, charsize=1.4, $
;      line=[pline,dline,bline,cline], color=djs_icolor([pcolor,dcolor,bcolor,ccolor]), $
;      pspacing=1.9, thick=10
;    
;    djs_oplot, polyfit[indx].bin_mass, polyfit[indx].bin_oh, psym=6, $
;      symsize=1.5, thick=6, color=fsc_color('brown',101)
;    djs_oplot, massaxis, mz_brokenpl(massaxis,brokenpl[indx].coeff), $
;      line=bline, color=bcolor, thick=8
;    djs_oplot, massaxis, mz_poly(massaxis,polyfit[indx].coeff), $
;      line=pline, color=pcolor, thick=8
;    djs_oplot, massaxis, mz_doublepl(massaxis,doublepl[indx].coeff), $
;      line=dline, color=dcolor, thick=8
;    djs_oplot, massaxis, mz_closedbox(massaxis,closed[indx].coeff), $
;      line=cline, color=ccolor, thick=8
;    im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)
;
