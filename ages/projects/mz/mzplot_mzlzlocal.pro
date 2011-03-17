pro mzplot_mzlzlocal, ps=ps
; jm09mar19nyu - make various plots pertaining to the local LZ and MZ
; relations
; jm10oct10ucsd - major update    

; read the data    
    sdssancillary = read_mz_sample(/sdss,/mzhii_ancillary)
    sdssmass = read_mz_sample(/sdss,/mzhii_mass)
    sdssohdust = read_mz_sample(/sdss,/mzhii_log12oh)
    sdssohnodust = read_mz_sample(/sdss,/mzhii_log12oh,/nodust)

    agesancillary = read_mz_sample(/mzhii_ancillary)
    agesmass = read_mz_sample(/mzhii_mass)
    agesohdust = read_mz_sample(/mzhii_log12oh)

    mzpath = ages_path(/projects)+'mz/'
    qapath = mzpath+'qaplots/'
    pspath = ages_path(/papers)+'mz/FIG_MZ/'
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

    massaxis1 = im_array(8.8,11.3,0.01)

    smzbinsize = 0.1
    amzbinsize = 0.1
    slzbinsize = 0.2
    alzbinsize = 0.2

; ---------------------------------------------------------------------------    
; SDSS MZ relation (KK04 calibration) with and without the upper
; limits 
    all = mzlz_grab_info(sdssohdust,sdssancillary,sdssmass,/kk04)
    nolim = mzlz_grab_info(sdssohdust,sdssancillary,sdssmass,/kk04,/nolimits)

    allfit = fit_mz_closedbox(all.mass,all.oh,all.weight,oh_err=all.oh_err)
    nolimfit = fit_mz_closedbox(nolim.mass,nolim.oh,nolim.weight,oh_err=nolim.oh_err)

    ohrange = [8.4,9.3]
    massrange = [8,12]
    
    psfile = qapath+'mzlocal_limits.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.8

; everything    
    mzplot_scatterplot, /sdss, all.mass, all.oh, weight=all.weight, $
      position=pos, xstyle=1, ystyle=1, xtitle=mzplot_masstitle(), $
      ytitle=mzplot_ohtitle(/kk04), xrange=massrange, yrange=ohrange
    oploterror, allfit.bin_mass, allfit.bin_oh, allfit.bin_oh_err, psym=symcat(15,thick=7), $
      symsize=1.1, thick=6, color=fsc_color('firebrick',101), $
      errcolor=fsc_color('firebrick',101)
    djs_oplot, massaxis1, mz_closedbox(massaxis1,allfit.coeff), line=0, color='orange', thick=4
    djs_oplot, massaxis1, mz_closedbox(massaxis1,nolimfit.coeff), line=5, color='blue', thick=4
    im_legend, ['Everything'], /left, /top, box=0, charsize=1.5, margin=0
    
; no upper limits
    mzplot_scatterplot, /sdss, nolim.mass, nolim.oh, weight=nolim.weight, $
      position=pos, xstyle=1, ystyle=1, xtitle=mzplot_masstitle(), $
      ytitle=mzplot_ohtitle(/kk04), xrange=massrange, yrange=ohrange
    oploterror, nolimfit.bin_mass, nolimfit.bin_oh, nolimfit.bin_oh_err, psym=symcat(15,thick=7), $
      symsize=1.1, thick=6, color=fsc_color('firebrick',101), $
      errcolor=fsc_color('firebrick',101)
    djs_oplot, massaxis1, mz_closedbox(massaxis1,allfit.coeff), line=0, color='orange', thick=4
    djs_oplot, massaxis1, mz_closedbox(massaxis1,nolimfit.coeff), line=5, color='blue', thick=4
    im_legend, ['No Upper Limits'], /left, /top, box=0, charsize=1.5, margin=0
    
    im_plotconfig, /psclose, psfile=psfile, /gzip

stop       
       
; ---------------------------------------------------------------------------    
; SDSS and low-redshift AGES MZ relation at z~0.1 for all three
; calibrations 
    massrange1 = [7.9,11.7]
    ohrange1 = [8.4,9.3]
;   ohrange1 = [8.25,9.45]
    verbose = 1
    
    sdsslevels = [0.25,0.5,0.75,0.9,0.95]
    ageslevels = [0.25,0.5,0.75,0.9]

    mzlocal = mrdfits(mzpath+'mzlocal_sdss_closedbox.fits.gz',1)
    massaxis = range(8.8,11.4,500)
    
    for ii = 2, 2 do begin
;   for ii = 0, 2 do begin
       t04 = 0 & m91 = 0 & kk04 = 0
       case ii of
          0: t04 = 1
          1: m91 = 1
          2: kk04 = 1
       endcase
       if keyword_set(t04) then calib = '_t04'
       if keyword_set(m91) then calib = '_m91'
       if keyword_set(kk04) then calib = '_kk04'
       ohtitle1 = mzplot_ohtitle(t04=t04,m91=m91,kk04=kk04)
       
       sinfo = mzlz_grab_info(sdssohdust,sdssancillary,sdssmass,$
         t04=t04,m91=m91,kk04=kk04)
       ainfo = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
         t04=t04,m91=m91,kk04=kk04,zmin=0.05,zmax=0.15)
       sbin = im_medxbin(sinfo.mass,sinfo.oh,smzbinsize,$
         weight=sinfo.oh_err/sqrt(sinfo.weight),minpts=500,verbose=verbose,$
;        weight=sinfo.weight,minpts=500,verbose=verbose,$
         minx=8.5)
       abin = im_medxbin(ainfo.mass,ainfo.oh,amzbinsize,$
         weight=ainfo.weight,minpts=15,verbose=verbose,$
         minx=8.5)
;       diff = interpol(sbin.medy,sbin.xbin,abin.xbin)-abin.medy
;       djs_plot, abin.xbin, diff, psym=6, xsty=3, ysty=3, yrange=[-0.1,0.2]
;       ww = where(abin.xbin gt 9.5)
;;      ww = where(abin.xbin gt 5.5)
;       splog, calib, median(diff[ww])
       
       psfile = pspath+'mzlocal'+calib+suffix
       im_plotconfig, 12, pos, psfile=psfile, charsize=1.9, $
         height=5.0, width=[4.5,4.5]
; SDSS
       mzplot_scatterplot, /sdss, sinfo.mass, sinfo.oh, weight=sinfo.oh_err/sqrt(sinfo.weight), $
;      lim = where(sinfo.ohlimit,comp=det)
;      mzplot_scatterplot, /sdss, sinfo.mass[det], sinfo.oh[det], weight=sinfo.weight[det], $
         position=pos[*,0], xstyle=1, ystyle=1, xtitle=mzplot_masstitle(), $
         ytitle=ohtitle1, xrange=massrange1, yrange=ohrange1, $
         levels=sdsslevels;, /annotate
;      djs_oplot, sinfo.mass[lim], sinfo.oh[lim], weight=sinfo.weight[lim], psym=3, color='yellow'
       oploterror, abin.xbin, abin.medy, abin.sigymean, psym=symcat(16,thick=7), $
         symsize=1.1, thick=6, color=fsc_color('blue',101), $
         errcolor=fsc_color('blue',101)
       oploterror, sbin.xbin, sbin.medy, sbin.sigymean, psym=symcat(15,thick=7), $
         symsize=1.1, thick=6, color=fsc_color('firebrick',101), $
         errcolor=fsc_color('firebrick',101)
;      djs_oplot, massaxis, mz_closedbox(massaxis,mzlocal[ii].coeff), $
;        line=0, color='orange', thick=4
;      djs_oplot, abin.xbin, abin.medy, psym=symcat(9,thick=7), $
;        symsize=1.1, thick=6, color=fsc_color('blue',101)
;      djs_oplot, sbin.xbin, sbin.medy, psym=symcat(6,thick=7), $
;        symsize=1.1, thick=6, color=fsc_color('firebrick',101)
       im_legend, ['SDSS - 0.033<z<0.25'], /left, /top, box=0, $
         charsize=1.5, margin=0, psym=symcat(15,thick=7), color='firebrick'
;      im_legend, ['SDSS'], /left, /bottom, box=0, charsize=1.5, $
;        margin=0, psym=symcat(6,thick=7), color='firebrick'
; MZ: low-z AGES
       lim = where(ainfo.ohlimit,comp=det)
       mzplot_scatterplot, /ages, ainfo.mass[det], ainfo.oh[det], weight=ainfo.weight[det], $
         /noerase, position=pos[*,1], xstyle=1, ystyle=1, xtitle=mzplot_masstitle(), $
         ytitle='', xrange=massrange1, yrange=ohrange1, ytickname=replicate(' ',10), $
         levels=ageslevels, npix=16;, /annotate
       djs_oplot, ainfo.mass[lim], ainfo.oh[lim], weight=ainfo.weight[lim], psym=6, sym=0.1, color='yellow'
       oploterror, abin.xbin, abin.medy, abin.sigymean, psym=symcat(16,thick=7), $
         symsize=1.1, thick=6, color=fsc_color('blue',101), $
         errcolor=fsc_color('blue',101)
       oploterror, sbin.xbin, sbin.medy, sbin.sigymean, psym=symcat(15,thick=7), $
         symsize=1.1, thick=6, color=fsc_color('firebrick',101), $
         errcolor=fsc_color('firebrick',101)
;      djs_oplot, massaxis, mz_closedbox(massaxis,mzlocal[ii].coeff), $
;        line=0, color='orange', thick=4
;      djs_oplot, sbin.xbin, sbin.medy, psym=symcat(6,thick=7), $
;        symsize=1.1, thick=6, color=fsc_color('firebrick',101)
;      djs_oplot, abin.xbin, abin.medy, psym=symcat(9,thick=7), $
;        symsize=1.1, thick=6, color=fsc_color('blue',101)
       im_legend, ['AGES - 0.05<z<0.15'], /left, /top, box=0, $
         charsize=1.5, margin=0, psym=symcat(16,thick=7), color='blue'
;      im_legend, ['AGES'], /left, /bottom, box=0, charsize=1.5, $
;        margin=0, psym=symcat(9,thick=7), color='blue'
       im_plotconfig, /psclose
    endfor

stop

; ---------------------------------------------------------------------------    
; compare the local MZ relations based on the various calibrations;
; just plot the best-fitting results, not the data
    mzlocal = mrdfits(mzpath+'mzlocal_sdss_closedbox.fits.gz',1)
    mzlocal_cor = mrdfits(mzpath+'mzlocal_sdss_fluxcor_closedbox.fits.gz',1)
;   usefluxcor = 0
;   if usefluxcor then ff = '_fluxcor' else ff = ''
;   mzlocal = mrdfits(mzpath+'mzlocal_sdss'+ff+'_closedbox.fits.gz',1)
    
    massrange1 = [8.5,11.7]
    ohrange1 = [8.5,9.2]
    massaxis1 = range(9.0,11.4,500)
    massaxis2 = range(9.0,11.4,500)

    t04color = 'red'    & t04line = 3
    m91color = 'black'  & m91line = 0
    kk04color = 'blue'  & kk04line = 5
    tremonticolor = 'orange' & tremontiline = 4

    psfile = pspath+'mzlocal_compare'+suffix
    im_plotconfig, 0, pos, psfile=psfile, height=5.8

    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=massrange1, yrange=ohrange1, xtitle=mzplot_masstitle(), $
      ytitle=mzplot_ohtitle()
    djs_oplot, massaxis2, tremonti_mz(massaxis=massaxis2), $
      color=fsc_color(tremonticolor,101), line=tremontiline, $
      thick=8

; EWs    
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal[0].coeff), $
      color=fsc_color(t04color,101), line=t04line, thick=8
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal[1].coeff), $
      color=fsc_color(m91color,101), line=m91line, thick=8
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal[2].coeff), $
      color=fsc_color(kk04color,101), line=kk04line, thick=8

; cor    
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal_cor[0].coeff), $
      color=fsc_color(t04color,101), line=t04line, thick=2
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal_cor[1].coeff), $
      color=fsc_color(m91color,101), line=m91line, thick=2
    djs_oplot, massaxis1, mz_closedbox(massaxis1,mzlocal_cor[2].coeff), $
      color=fsc_color(kk04color,101), line=kk04line, thick=2
    
    im_legend, strupcase(mzlocal.calib), /left, /top, box=0, margin=0, charsize=1.4, $
      line=[t04line,m91line,kk04line], color=[t04color,m91color,kk04color], $
      pspacing=1.9, thick=10
    im_plotconfig, /psclose

; ---------------------------------------------------------------------------    
; compare the various ways of fitting the SDSS MZ relation:
; polynomial, double power-law, and broken power-law (see
; FIT_MZLZLOCAL); make the plot just using one calibration and EWs
; (rather than reddening-corrected fluxes)
    usefluxcor = 0
    if usefluxcor then ff = '_fluxcor' else ff = ''
    closed = mrdfits(mzpath+'mzlocal_sdss'+ff+'_closedbox.fits.gz',1)
    polyfit = mrdfits(mzpath+'mzlocal_sdss'+ff+'_poly.fits.gz',1)
    doublepl = mrdfits(mzpath+'mzlocal_sdss'+ff+'_doublepl.fits.gz',1)
    brokenpl = mrdfits(mzpath+'mzlocal_sdss'+ff+'_brokenpl.fits.gz',1)

;   calib = 'm91' &  t04 = 0 & m91 = 1 & kk04 = 0
    calib = 'kk04' &  t04 = 0 & m91 = 0 & kk04 = 1
;   calib = 't04' & t04 = 1 & m91 = 0 & kk04 = 0
    indx = where(strtrim(polyfit.calib,2) eq calib)
    if usefluxcor then begin
       info = mzlz_grab_info(sdssohnodust,/flux,$
         sdssancillary,sdssmass,t04=t04,m91=m91,kk04=kk04)
    endif else begin
       info = mzlz_grab_info(sdssohdust,$ ; EWs
         sdssancillary,sdssmass,t04=t04,m91=m91,kk04=kk04)
    endelse
    ytitle = mzplot_ohtitle(t04=t04,m91=m91,kk04=kk04)

    massrange1 = [8.1,11.9]
;   ohrange1 = [8.45,9.3] ; KK04
    ohrange1 = [8.25,9.3] ; T04
    massaxis = range(8.4,11.6,500)

    pcolor = 'red'   & pline = 3
    bcolor = 'black' & bline = 4
    dcolor = 'blue'  & dline = 5
    ccolor = 'orange'  & cline = 0
    
    psfile = pspath+'mz_fitting'+suffix
    im_plotconfig, 0, pos, psfile=psfile, height=5.8

    mzplot_scatterplot, info.mass, info.oh, /sdss, position=pos[*,0], $
      xsty=1, ysty=1, xrange=massrange1, yrange=ohrange1, $
      xtitle=mzplot_masstitle(), ytitle=ytitle, $
      /nogrey, ccolor=djs_icolor('grey')
    legend, ['Quadratic Polynomial','Double Power-law',$
      'Broken Power-law','Closed-Box Model'], $
      /right, /bottom, box=0, margin=0, charsize=1.4, $
      line=[pline,dline,bline,cline], color=djs_icolor([pcolor,dcolor,bcolor,ccolor]), $
      pspacing=1.9, thick=10
    
    djs_oplot, polyfit[indx].bin_mass, polyfit[indx].bin_oh, psym=6, $
      symsize=1.5, thick=6, color=fsc_color('brown',101)
    djs_oplot, massaxis, mz_brokenpl(massaxis,brokenpl[indx].coeff), $
      line=bline, color=bcolor, thick=8
    djs_oplot, massaxis, mz_poly(massaxis,polyfit[indx].coeff), $
      line=pline, color=pcolor, thick=8
    djs_oplot, massaxis, mz_doublepl(massaxis,doublepl[indx].coeff), $
      line=dline, color=dcolor, thick=8
    djs_oplot, massaxis, mz_closedbox(massaxis,closed[indx].coeff), $
      line=cline, color=ccolor, thick=8
    im_plotconfig, /psclose

; ---------------------------------------------------------------------------    
; local MZ relation using the MPA metallicities and compare with
; Tremonti et al. 2004
    usefluxcor = 1
    if usefluxcor then ff = '_fluxcor' else ff = ''
    mzlocal = mrdfits(mzpath+'mzlocal_sdss'+ff+'_brokenpl.fits.gz',1)

    massrange1 = [7.9,11.7]
    ohrange1 = [8.0,9.45]
    massaxis1 = range(8.5,11.4,500)

    xtremonti = [8.57,8.67,8.76,8.86,8.96,9.06,9.16,9.26,$
      9.36,9.46,9.57,9.66,9.76,9.86,9.96,10.06,10.16,10.26,$
      10.36,10.46,10.56,10.66,10.76,10.86,10.95,11.05,11.15,11.25]
    ytremonti = [8.44,8.48,8.57,8.61,8.63,8.66,8.68,8.71,$
      8.74,8.78,8.82,8.84,8.87,8.90,8.94,8.97,8.99,9.01,$
      9.03,9.05,9.07,9.08,9.09,9.10,9.11,9.11,9.12,9.12]
    
;   mass = sdssmass.mass_avg
    mass = sdssancillary.mass_median
    oh = sdssancillary.oh_median+randomn(seed,n_elements(sdssancillary))*0.02
    sbin = im_medxbin(mass,oh,0.1,minpts=75,minx=8.8)

    psfile = pspath+'mzlocal_mpa'+suffix
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
    djs_oplot, massaxis1, mz_brokenpl(massaxis1,mzlocal[0].coeff), $
      color=fsc_color('red',101), line=5, thick=8

    im_plotconfig, /psclose
    
stop    
    
; ---------------------------------------------------------------------------    
; SDSS and low-redshift AGES LZ relations at z~0.1 for all three
; calibrations 
    magrange1 = [-16.3,-22.9]
    ohrange1 = [8.4,9.3]
;   ohrange1 = [8.25,9.45] ; T04
    verbose = 1
    magaxis = range(-17,-23,300)
    
    sdsslevels = [0.25,0.5,0.75,0.9,0.95]
    ageslevels = [0.25,0.5,0.75,0.9]

    sline = 0 & scolor = 'firebrick'
    aline = 5 & acolor = 'blue'
    
    sdss_lzfitall = mrdfits(mzpath+'lzlocal_sdss.fits.gz',1) ; based on EWs
    ages_lzfitall = mrdfits(mzpath+'lzlocal_ages.fits.gz',1) ; based on EWs

    for ii = 2, 0, -1 do begin
;   for ii = 0, 2 do begin
       t04 = 0 & m91 = 0 & kk04 = 0
       case ii of
          0: t04 = 1
          1: m91 = 1
          2: kk04 = 1
       endcase
       if keyword_set(t04) then calib = 't04'
       if keyword_set(m91) then calib = 'm91'
       if keyword_set(kk04) then calib = 'kk04'
       ohtitle1 = mzplot_ohtitle(t04=t04,m91=m91,kk04=kk04)

; best-fitting LZ relations
       indx = where(strtrim(sdss_lzfitall.calib,2) eq calib)
       sdss_lzfit = sdss_lzfitall[indx]
       ages_lzfit = ages_lzfitall[indx]
       
       sinfo = mzlz_grab_info(sdssohdust,sdssancillary,sdssmass,$
         t04=t04,m91=m91,kk04=kk04)
       ainfo = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
         t04=t04,m91=m91,kk04=kk04,zmin=0.05,zmax=0.15)
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
         levels=sdsslevels;, /annotate
;       oploterror, abin.xbin, abin.medy, abin.sigymean, psym=symcat(16,thick=7), $
;         symsize=1.1, thick=6, color=fsc_color('blue',101), $
;         errcolor=fsc_color('blue',101)
;       oploterror, sbin.xbin, sbin.medy, sbin.sigymean, psym=symcat(15,thick=7), $
;         symsize=1.1, thick=6, color=fsc_color('firebrick',101), $
;         errcolor=fsc_color('firebrick',101)
       oplot_lzfit, sdss_lzfit.coeff, band='B', linestyle=sline, linecolor=scolor
       oplot_lzfit, ages_lzfit.coeff, band='B', linestyle=aline, linecolor=acolor

       im_legend, ['SDSS - 0.033<z<0.25'], /left, /top, box=0, thick=8, $
         charsize=1.5, margin=0, line=sline, color=scolor, pspacing=1.9
;      im_legend, ['SDSS - 0.033<z<0.25'], /left, /top, box=0, $
;        charsize=1.5, margin=0, psym=symcat(15,thick=7), color='firebrick'
; LZ: low-redshift AGES
       mzplot_scatterplot, /ages, ainfo.mb_ab, ainfo.oh, weight=ainfo.weight, $
         /noerase, position=pos[*,1], xstyle=1, ystyle=1, xtitle=mzplot_mbtitle(), $
         ytitle='', xrange=magrange1, yrange=ohrange1, ytickname=replicate(' ',10), $
         levels=ageslevels, npix=16;, /annotate
;       oploterror, abin.xbin, abin.medy, abin.sigymean, psym=symcat(16,thick=7), $
;         symsize=1.1, thick=6, color=fsc_color('blue',101), $
;         errcolor=fsc_color('blue',101)
;       oploterror, sbin.xbin, sbin.medy, sbin.sigymean, psym=symcat(15,thick=7), $
;         symsize=1.1, thick=6, color=fsc_color('firebrick',101), $
;         errcolor=fsc_color('firebrick',101)
       oplot_lzfit, sdss_lzfit.coeff, band='B', linestyle=sline, linecolor=scolor
       oplot_lzfit, ages_lzfit.coeff, band='B', linestyle=aline, linecolor=acolor
       
       im_legend, ['AGES - 0.05<z<0.15'], /left, /top, box=0, thick=8, $
         charsize=1.5, margin=0, line=aline, color=acolor, pspacing=1.9
;      im_legend, ['AGES - 0.05<z<0.15'], /left, /top, box=0, $
;        charsize=1.5, margin=0, psym=symcat(16,thick=7), color='blue'
       im_plotconfig, /psclose
    endfor

return
end
