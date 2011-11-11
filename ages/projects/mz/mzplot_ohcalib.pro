pro mzplot_ohcalib, ps=ps
; jm09mar25nyu - (O/H) calibration plots

; read the data    
    sdssancillary = read_mz_sample(/mzhii_ancillary,/sdss)
    sdssmass = read_mz_sample(/mzhii_mass,/sdss)
;   sdssispec = read_mz_sample(/mzhii_ispec,/sdss)
    sdssohdust = read_mz_sample(/mzhii_log12oh,/sdss)
    sdssohnodust = read_mz_sample(/mzhii_log12oh,/nodust,/sdss)

    agesispec = read_mz_sample(/mzhii_ispec)
    agesohdust = read_mz_sample(/mzhii_log12oh)

    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(ps) then begin
       pspath = qapath
       suffix = '.ps'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse

; ------------------------------------------------------------
; Figure A1 - 12+log(O/H) vs R23 for various calibrations

    model_logr23 = range(-0.5,1.0,1500)
    logu = [-3.5,-2.5]
    model_logq = logu+alog10(im_light(/cm))
;   model_logq = alog10([2D7,1D8])
;   logu = model_logq-alog10(im_light(/cm))
    model_logo32 = [-1,-0]
;   model_logo32 = [-0.6,-0.2]

    xrange = [-0.2,1.1]
    yrange = [8.1,9.3]
    
    psfile = pspath+'12oh_vs_r23'+suffix
    im_plotconfig, 0, pos, psfile=psfile, height=5.5

    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle='log (R_{23})', $
      ytitle='12 + log (O/H)'

    djs_oplot, model_logr23, poly(model_logr23,[9.185D,-0.313D,-0.264D,-0.321D]), $
      thick=8, line=0, color='red'
    mzoplot_kk04_models, model_logr23=model_logr23, model_logq=model_logq, $
      linestyle=[2,1], linecolor=['blue','orange'], /nolegend, thick=8
    mzoplot_m91_models, model_logr23=model_logr23, model_logo32=model_logo32, $
      linestyle=[3,5], linecolor=['dark green','magenta'], /nolegend, thick=8

    label = ['T04','KK04: log(U)='+string(logu,format='(F4.1)'),$
     'M91: log(O_{32})='+['','-']+strtrim(string(model_logo32,format='(F4.1)'),2)]
    color = ['red','blue','orange','dark green','magenta']
    line = [0,2,1,3,5]
    
    legend, textoidl(label), /left, /bottom, box=0, charsize=1.4, margin=1, $
      color=djs_icolor(color), textcolor=djs_icolor(color), line=line, $
      thick=8, pspacing=1.8
    
    im_plotconfig, /psclose, psfile=psfile

stop    
    
; ------------------------------------------------------------
; Figure 8 - AGES + SDSS - [NII]/Ha vs R23, illustrating that our
; galaxies belong on the upper branch 

    sdssindx = where((sdssohdust.niiha gt -900.0) and $
      (sdssohdust.r23 gt -900.0),nsdss)
    agesindx = where((agesohdust.niiha gt -900.0) and $
      (agesohdust.r23 gt -900.0),nages)
    splog, nsdss, nages

    levels = [0.5,0.75,0.9,0.975]
;   levels = errorf((findgen(3)+1)/sqrt(2))
    xrange = [-1.65,-0.05]
    yrange = [-0.3,1.2]
    xtitle = textoidl('log ([N II] \lambda6584/H\alpha)')
    ytitle = textoidl('log (R_{23})')
    
    psfile = pspath+'niiha_vs_r23'+suffix
    im_plotconfig, 1, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=[4.3,4.3], height=4.3
; SDSS
    mzplot_scatterplot, /sdss, sdssohdust[sdssindx].niiha, alog10(sdssohdust[sdssindx].r23), $
      position=pos[*,0], xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, levels=levels, /nogrey
    djs_oplot, -1.1*[1,1], [0.2,1.1], line=2, thick=6
    xyouts, -1.3, 0.5, 'Lower!cBranch', align=0.5, charsize=1.2, $
      charthick=2.5
    xyouts, -0.6, 1.05, 'Upper Branch', align=0.5, charsize=1.2, $
      charthick=2.5
    legend, 'SDSS', /left, /bottom, box=0, charsize=1.6, margin=0
; AGES
    mzplot_scatterplot, agesohdust[agesindx].niiha, alog10(agesohdust[agesindx].r23), $
      /noerase, position=pos[*,1], xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle='', ytickname=replicate(' ',10), levels=levels, npix=20, /nogrey

;   plotsym, 6, 0.8, color=djs_icolor('red'), thick=3
;   lim = where(agesispec.bpt_nii_ha_limit gt -900.0)
;   djs_oplot, agesispec[lim].bpt_nii_ha_limit, alog10(agesohdust[lim].r23), $
;     psym=8
    
    djs_oplot, -1.1*[1,1], [0.2,1.1], line=2, thick=6
    xyouts, -1.3, 0.5, 'Lower!cBranch', align=0.5, charsize=1.2, $
      charthick=2.5
    xyouts, -0.6, 1.05, 'Upper Branch', align=0.5, charsize=1.2, $
      charthick=2.5
    legend, 'AGES', /left, /bottom, box=0, charsize=1.6, margin=0
    im_plotconfig, /psclose, psfile=psfile

; ------------------------------------------------------------
; Figure 9 - O/H_cor vs O/H_EW
    for ii = 0, 2 do begin
       case ii of
          0: begin
             t04 = 1 & m91 = 0 & kk04 = 0
             ohrange1 = [8.4,9.3]
             calib = 't04'
          end
          1: begin
             t04 = 0 & m91 = 1 & kk04 = 0
             ohrange1 = [8.3,9.15]
             calib = 'm91'
          end
          2: begin
             t04 = 0 & m91 = 0 & kk04 = 1
             ohrange1 = [8.48,9.21]
             calib = 'kk04'
          end
       endcase
       
       cor = mzlz_grab_info(sdssohnodust,sdssancillary,sdssmass,$
         kk04=kk04,t04=t04,m91=m91,/flux,/nolimit)
       ew = mzlz_grab_info(sdssohdust,sdssancillary,sdssmass,$
         kk04=kk04,t04=t04,m91=m91,/nolimit)
       match, cor.id, ew.id, m1, m2
       ewoh = ew.oh[m2]
       coroh = cor.oh[m1]
       resid = ewoh-coroh
       umb = sdssancillary[ew.id[m2]].k_ubvrijhk_absmag_00[0]-$
         sdssancillary[ew.id[m2]].k_ubvrijhk_absmag_00[1]
;      gmr = sdssancillary[ew.id[m2]].k_ugriz_absmag_01[1]-$
;        sdssancillary[ew.id[m2]].k_ugriz_absmag_01[2]

       splog, calib
       med = im_medxbin(umb,resid,0.05,minpts=60,/verbose)
;      med = im_medxbin(gmr,resid,0.05,minpts=50,/verbose)
;      med = im_medxbin(coroh,resid,0.05,minx=8.5,/verbose)
       splog, calib, mean(resid), median(resid), djsig(resid)
       
       xtitle = textoidl('12 + log (O/H)_{'+strupcase(calib)+', cor}')
       ytitle = textoidl('12 + log (O/H)_{'+strupcase(calib)+', EW}')
       
; make the plot    
       psfile = pspath+'ohcor_vs_ohews_'+calib+suffix
       im_plotconfig, 6, pos, psfile=psfile, yspace=1.0, $
         xmargin=[1.4,0.3], width=6.8, height=[4.0,3.0], charsize=2
; main plot
       mzplot_scatterplot, coroh, ewoh, position=pos[*,0], $
         /sdss, xsty=1, ysty=1, xrange=ohrange1, yrange=ohrange1, $
         xtitle=xtitle, ytitle=ytitle, /nogrey, ccolor=djs_icolor('grey')
       djs_oplot, !x.crange, !y.crange, line=0, thick=7;, color='red'
; residuals
;      mzplot_scatterplot, d4000, resid, /noerase, position=pos[*,1], $
;        /sdss, xsty=1, ysty=1, xrange=[0.8,2.0], yrange=0.4*[-1,1], $
;        xtitle=textoidl('D_{n}(4000)'), ytitle='Residuals (dex)'
       mzplot_scatterplot, umb, resid, /noerase, position=pos[*,1], $
         /sdss, xsty=1, ysty=1, xrange=[0.0,1.6], yrange=0.35*[-1,1], $
         xtitle=textoidl('U - B'), ytitle='Residuals (dex)', $
         /nogrey, ccolor=djs_icolor('grey'), ytickinterval=0.2
       oploterror, med.xbin, med.medy, med.quant75-med.medy, psym=6, $
         color=djs_icolor('navy'), errthick=6, thick=6, $
         errcolor=djs_icolor('navy'), /hibar
       oploterror, med.xbin, med.medy, med.medy-med.quant25, psym=3, $
         color=djs_icolor('navy'), errthick=6, thick=6, $
         errcolor=djs_icolor('navy'), /lobar
       djs_oplot, !x.crange, [0,0], line=0, thick=7;, color='red'
       im_plotconfig, /psclose, psfile=psfile
    endfor

stop    
    
stop
stop
stop    
    
; ------------------------------------------------------------
; 12+log(O/H)_{cor} vs 12+log(O/H)_{EW} measured using different
; assumptions about "alpha" - SDSS

    psfile = pspath+'sdss_12oh_cor_vs_12oh_ews'+suffix
    im_plotconfig, 4, pos, psfile=psfile, width=4.5, height=2.6*[1,1,1], $
      xmargin=[1.2,0.3], ymargin=[0.25,0.85], charsize=1.5

    xrange = [8.52,9.23]
    yrange = xrange

    xtitle = textoidl('12 + log (O/H)_{cor}')
    ytitle = textoidl('12 + log (O/H)_{EW}')

    indx = where((sdssohnodust.zstrong_12oh_kk04 gt -900) and $
      (sdssohdust.zstrong_ew_alpha_unity_12oh_kk04 gt -900) and $
      (sdssohdust.zstrong_ew_alpha_d4000_12oh_kk04 gt -900) and $
      (sdssohdust.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(sdssohnodust.r23branch_kk04,2) eq 'U') and $
      (strtrim(sdssohdust.r23branch_ew_alpha_unity_kk04,2) eq 'U') and $
      (strtrim(sdssohdust.r23branch_ew_alpha_d4000_kk04,2) eq 'U') and $
      (strtrim(sdssohdust.r23branch_ew_alpha_gr_kk04,2) eq 'U'))

    ohcor = sdssohnodust[indx].zstrong_12oh_kk04
    ohcorerr = sdssohnodust[indx].zstrong_12oh_kk04_err
    ohalphaunity = sdssohdust[indx].zstrong_ew_alpha_unity_12oh_kk04
    ohalphaunityerr = sdssohdust[indx].zstrong_ew_alpha_unity_12oh_kk04_err
    ohalphad4000 = sdssohdust[indx].zstrong_ew_alpha_d4000_12oh_kk04
    ohalphad4000err = sdssohdust[indx].zstrong_ew_alpha_d4000_12oh_kk04_err
    ohalphagr = sdssohdust[indx].zstrong_ew_alpha_gr_12oh_kk04
    ohalphagrerr = sdssohdust[indx].zstrong_ew_alpha_gr_12oh_kk04_err

; ##########
; alpha=1.0
; ##########
    
    stats = im_medxbin(ohcor,ohalphaunity-ohcor,0.1,minx=8.6,$
      maxx=9.2,minpts=50L,verbose=1)
    allstats = im_stats(ohalphaunity-ohcor,sigrej=3.0)
    xstr = '\sigma='+strtrim(string(min(stats.sigy),format='(F12.2)'),2)+$
      '-'+strtrim(string(max(stats.sigy),format='(F12.2)'),2)
;   xstr = '\sigma = '+strtrim(string(allstats.sigma_rej,format='(F12.3)'),2)
;   xstr = '\Delta = '+strtrim(string(allstats.median_rej,format='(F12.2)'),2)+$
;     '\pm'+strtrim(string(allstats.sigma_rej,format='(F12.2)'),2)

    mzsdss_hogg_scatterplot, ohcor, ohalphaunity, position=pos[*,0], $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange
    djs_oplot, !x.crange, !y.crange, line=0
    im_legend, '\alpha = 1.0', /left, /top, box=0, margin=0, charsize=1.3
    im_legend, xstr, /right, /bottom, box=0, margin=0

; ##########
; alpha-D(4000)
; ##########
    
    stats = im_medxbin(ohcor,ohalphad4000-ohcor,0.1,minx=8.6,$
      maxx=9.2,minpts=50L,verbose=1)
    allstats = im_stats(ohalphad4000-ohcor,sigrej=3.0)
    xstr = '\sigma='+strtrim(string(min(stats.sigy),format='(F12.2)'),2)+$
      '-'+strtrim(string(max(stats.sigy),format='(F12.2)'),2)
;   xstr = '\sigma = '+strtrim(string(allstats.sigma_rej,format='(F12.3)'),2)
;   xstr = '\Delta = '+strtrim(string(allstats.median_rej,format='(F12.2)'),2)+$
;     '\pm'+strtrim(string(allstats.sigma_rej,format='(F12.2)'),2)

    mzsdss_hogg_scatterplot, ohcor, ohalphad4000, /noerase, position=pos[*,1], $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange
    djs_oplot, !x.crange, !y.crange, line=0
    im_legend, '\alpha-D_{n}(4000)', /left, /top, box=0, margin=0
    im_legend, xstr, /right, /bottom, box=0, margin=0
    
; ##########
; alpha-gr
; ##########
    
    stats = im_medxbin(ohcor,ohalphagr-ohcor,0.1,minx=8.6,$
      maxx=9.2,minpts=50L,verbose=1)
    allstats = im_stats(ohalphagr-ohcor,sigrej=3.0)
    xstr = '\sigma='+strtrim(string(min(stats.sigy),format='(F12.2)'),2)+$
      '-'+strtrim(string(max(stats.sigy),format='(F12.2)'),2)
;   xstr = '\sigma = '+strtrim(string(allstats.sigma_rej,format='(F12.3)'),2)
;   xstr = '\Delta = '+strtrim(string(allstats.median_rej,format='(F12.2)'),2)+$
;     '\pm'+strtrim(string(allstats.sigma_rej,format='(F12.2)'),2)

    mzsdss_hogg_scatterplot, ohcor, ohalphagr, /noerase, position=pos[*,2], $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange
    djs_oplot, !x.crange, !y.crange, line=0, thick=linethick2
    im_legend, '\alpha-^{0.1}(g-r)', /left, /top, box=0, margin=0
    im_legend, xstr, /right, /bottom, box=0, margin=0
    
    im_plotconfig, /psclose, psfile=psfile

; ------------------------------------------------------------
; calibrate the EW-alpha relation; see BUILD_MZ_LOG12OH_SAMPLE

    if keyword_set(ps) then psfile = pspath+'sdss_alpha_vs_d4000_gr'+suffix
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.8

    use_alphafit = mrdfits(mzpath+'sdss_alphafit.fits.gz',1)

    indx = where((sdssohnodust.zstrong_o32 gt -900) and $
      (sdssohnodust.zstrong_ew_alpha_unity_o32 gt -900.0),nindx)
    o32 = alog10(sdssohnodust[indx].zstrong_o32)
    ewo32 = alog10(sdssohnodust[indx].zstrong_ew_alpha_unity_o32)                        ; alpha=1
    alpha = sdssohnodust[indx].zstrong_ew_alpha_unity_o32/sdssohnodust[indx].zstrong_o32 ; linear
    
    alphabest_med = weighted_quantile(alpha,weight,quant=0.5)
    alphabest_sig = (weighted_quantile(alpha,weight,quant=0.84)-$
      weighted_quantile(alpha,weight,quant=0.16))/2.0
    splog, 'Median alpha = ', alphabest_med, '+/-', alphabest_sig

    d4000 = sdssispec[indx].d4000_narrow_model[0]
    d4000err = sdssispec[indx].d4000_narrow_model[1]
    gr = sdssancillary[indx].ugriz_absmag[1]-sdssancillary[indx].ugriz_absmag[2]

    minpts1 = 300
    d4000_running = im_medxbin(d4000,alpha,0.03,minx=1.0,$
      maxx=1.8,minpts=minpts1,verbose=0)
;   alpha_d4000_coeff = poly_fit(d4000_running.binctr,d4000_running.meany,2,$
;     chisq=alpha_d4000_chisq,measure_errors=d4000_running.sigy,sigma=alpha_d4000_coeff_err)
    alpha_d4000_coeff = use_alphafit[0].alpha_coeff

    gr_running = im_medxbin(gr,alpha,0.05,minx=0.15,$
      maxx=1.0,minpts=minpts1,verbose=0)
;   alpha_gr_coeff = poly_fit(gr_running.binctr,gr_running.meany,2,$
;     chisq=alpha_gr_chisq,measure_errors=gr_running.sigy,sigma=alpha_gr_coeff_err)
    alpha_gr_coeff = use_alphafit[1].alpha_coeff

; make the plot    
    ytitle1 = textoidl('\alpha = [EW(O_{32})/(O_{32})_{cor}]')
    xtitle1 = textoidl('D_{n}(4000)')
    xtitle2 = textoidl('^{0.1}(g - r)')
    yrange1 = [0.3,1.9]
    xrange1 = [0.95,1.79]
    xrange2 = [-0.05,1.1]
; D(4000) plot
    mzsdss_hogg_scatterplot, d4000, alpha, position=pos[*,0], $
      xstyle=1, ystyle=1, xrange=xrange1, yrange=yrange1, $
      xtitle=xtitle1, ytitle=ytitle1
    d4000axis = im_array(min(d4000_running.binctr-d4000_running.binsz/2.0),$
      max(d4000_running.binctr+d4000_running.binsz/2.0),0.01)
    oploterror, d4000_running.binctr, d4000_running.meany, $
      d4000_running.sigy84-d4000_running.meany, psym=symcat(15), $
      symsize=1.2, /hibar, color=fsc_color('royal blue',1E8)
    oploterror, d4000_running.binctr, d4000_running.meany, $
      d4000_running.sigy16-d4000_running.meany, $
      psym=symcat(15), symsize=1.2, /lobar, color=fsc_color('royal blue',1E8)
    djs_oplot, d4000axis, poly(d4000axis,alpha_d4000_coeff), line=0, $
      thick=5.0, color='dark red'
; (g-r) plot
    mzsdss_hogg_scatterplot, gr, alpha, /noerase, position=pos[*,1], $
      xstyle=1, ystyle=1, xrange=xrange2, yrange=yrange1, $
      xtitle=xtitle2, ytitle='', ytickname=replicate(' ',10)
    graxis = im_array(min(gr_running.binctr-gr_running.binsz/2.0),$
      max(gr_running.binctr+gr_running.binsz/2.0),0.01)
    oploterror, gr_running.binctr, gr_running.meany, gr_running.sigy84-gr_running.meany, $
      psym=symcat(15), symsize=1.2, /hibar, color=fsc_color('royal blue',1E8)
    oploterror, gr_running.binctr, gr_running.meany, gr_running.sigy16-gr_running.meany, $
      psym=symcat(15), symsize=1.2, /lobar, color=fsc_color('royal blue',1E8)
    djs_oplot, graxis, poly(graxis,alpha_gr_coeff), line=0, $
      thick=5.0, color='dark red'

    im_plotconfig, /psclose, psfile=psfile

; ------------------------------------------------------------
; effect on LZ of EWs vs reddening-corrected fluxes 
    if keyword_set(ps) then psfile = pspath+'glzlocal_ohcalib'+suffix
    im_plotconfig, 5, pos, psfile=psfile, charsize=1.8

; some handy plotting variables    
    massrange1 = [8.1,12.1]
    mrrange1 = [-17.0,-23.9]
    mgrange1 = [-16.1,-23.9]
    ohrange1 = [8.35,9.37] ; [8.3,9.39]

    masstitle1 = textoidl('log (M_{*}h_{70}^2/M'+sunsymbol()+')')
    mgtitle1 = textoidl('M_{0.1g} - 5 log (h_{70})')
    ohtitle1 = textoidl('12 + log (O/H)')

    s0_linecolor = 'black'       & s0_linestyle = 0
    s1_linecolor = 'orange'      & s1_linestyle = 5
    s2_linecolor = 'royal blue'  & s2_linestyle = 1
    s3_linecolor = 'firebrick'   & s3_linestyle = 3

; read the output from FIT_MZLZLOCAL; most of the code below needs to
; match that routine!
    lzfit = mrdfits(mzpath+'lzlocal_sdss_coeffs.fits.gz',1)
    struct_print, lzfit
    pivotmag = lzfit[0].pivotmag

; gLZ: SDSS/cor
    s0 = mzlz_grab_info(sdssohnodust,sdssancillary,/flux)
    mzsdss_hogg_scatterplot, s0.mg_ab, s0.oh, position=pos[*,0], $
      xstyle=1, ystyle=1, xtitle='', ytitle=ohtitle1, $
      xrange=mgrange1, yrange=ohrange1, xtickname=replicate(' ',10)
    oplot_lzfit, lzfit[0].coeff, pivotmag, band='g', $
      linestyle=s0_linestyle, linecolor=s0_linecolor
    im_legend, ['Flux-Cor'], /left, /top, box=0, $
      charsize=1.6, margin=0
    im_legend, '\sigma='+string(lzfit[0].scatter,format='(F4.2)'), $
      /right, /bottom, box=0, charsize=1.4, margin=0
;   im_legend, '\sigma_{(O/H)}='+string(lzfit[0].scatter,format='(F4.2)')+' dex', $
;     /right, /bottom, box=0, charsize=1.4, margin=0
    
; gLZ: SDSS/alpha=1
    s1 = mzlz_grab_info(sdssohdust,sdssancillary,/ewunity)
    mzsdss_hogg_scatterplot, s1.mg_ab, s1.oh, /noerase, position=pos[*,1], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=mgrange1, yrange=ohrange1, ytickname=replicate(' ',10), $
      xtickname=replicate(' ',10)
    oplot_lzfit, lzfit[0].coeff, pivotmag, band='g', $
      linestyle=s0_linestyle, linecolor=s0_linecolor
    oplot_lzfit, lzfit[1].coeff, pivotmag, band='g', $
      linestyle=s1_linestyle, linecolor=s1_linecolor
    im_legend, ['\alpha=1'], /left, /top, box=0, $
      charsize=1.6, margin=0
    im_legend, '\sigma='+string(lzfit[1].scatter,format='(F4.2)'), $
      /right, /bottom, box=0, charsize=1.4, margin=0

; gLZ: SDSS/alpha-D(4000)
    s2 = mzlz_grab_info(sdssohdust,sdssancillary,/ewd4000)
    mzsdss_hogg_scatterplot, s2.mg_ab, s2.oh, /noerase, position=pos[*,2], $
      xstyle=1, ystyle=1, xtitle=mgtitle1, ytitle=ohtitle1, $
      xrange=mgrange1, yrange=ohrange1
    oplot_lzfit, lzfit[0].coeff, pivotmag, band='g', $
      linestyle=s0_linestyle, linecolor=s0_linecolor
    oplot_lzfit, lzfit[2].coeff, pivotmag, band='g', $
      linestyle=s2_linestyle, linecolor=s2_linecolor, thick=10
    im_legend, ['\alpha-D_{n}(4000)'], /left, /top, box=0, $
      charsize=1.6, margin=0
    im_legend, '\sigma='+string(lzfit[2].scatter,format='(F4.2)'), $
      /right, /bottom, box=0, charsize=1.4, margin=0

; gLZ: SDSS/alpha-(g-r)
    s3 = mzlz_grab_info(sdssohdust,sdssancillary,/ewgr)
    mzsdss_hogg_scatterplot, s3.mg_ab, s3.oh, /noerase, position=pos[*,3], $
      xstyle=1, ystyle=1, xtitle=mgtitle1, ytitle='', $
      xrange=mgrange1, yrange=ohrange1, ytickname=replicate(' ',10)
    oplot_lzfit, lzfit[0].coeff, pivotmag, band='g', $
      linestyle=s0_linestyle, linecolor=s0_linecolor
    oplot_lzfit, lzfit[3].coeff, pivotmag, band='g', $
      linestyle=s3_linestyle, linecolor=s3_linecolor
    im_legend, ['\alpha-^{0.1}(g-r)'], /left, /top, box=0, $
      charsize=1.6, margin=0
    im_legend, '\sigma='+string(lzfit[3].scatter,format='(F4.2)'), $
      /right, /bottom, box=0, charsize=1.4, margin=0

    im_plotconfig, /psclose, psfile=psfile

return
end
