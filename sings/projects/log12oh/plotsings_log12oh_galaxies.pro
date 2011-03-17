pro plotsings_log12oh_galaxies, ps=ps
; jm10mar10ucsd - build the galaxy plots

    common plotsings_gal, sdss
    
    metpath = sings_path(/projects)+'log12oh/'
    pspath = sings_path(/papers)+'log12oh/FIG_LOG12OH/'

    version = sings_log12oh_version()
    result = mrdfits(metpath+'sings_log12oh_'+version+'.fits.gz',1,/silent)
    if (n_elements(sdss) eq 0L) then sdss = mrdfits(metpath+$
      'sdss_log12oh_'+version+'.fits.gz',1)
    
    nuclear_all = read_sings_gandalf(/nuclear)
    drift20_all = read_sings_gandalf(/drift20)
    drift56_all = read_sings_gandalf(/drift56)

    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

    earlysym1 = 16 & earlypsize1 = 1.8 & earlycolor1 = 'navy' ; 'dodger blue'
    midsym1   = 15 & midpsize1   = 1.8 & midcolor1   = 'tan' ; 'salmon' ; 'black'
    latesym1  = 14 & latepsize1  = 2.3 & latecolor1  = 'firebrick' ; tan

;   earlysym1 = 16 & earlypsize1 = 1.8 & earlycolor1 = 'dodger blue'
;   midsym1   = 15 & midpsize1   = 1.7 & midcolor1   = 'black'
;   latesym1  = 4 & latepsize1  = 2.0 & latecolor1  = 'firebrick'

    nucsym = 17 & nucpsize = 1.6 & nuccolor = 'dark green'
    d20sym = 16 & d20psize = 1.3 & d20color = 'orange'
    d56sym = 15 & d56psize = 1.1 & d56color = 'navy'

    symthick1 = 6.0
    symthick2 = 8.0
    errthick1 = 5.0
    errthick2 = 7.0

    sfsym    = 16 & sfpsize    = 1.2 & sfcolor = 'dodger blue' ; 'tan'
    agnsym   = 6  & agnpsize   = 1.1 & agncolor = 'black' ; 'navy'
    sfagnsym = 17 & sfagnpsize = 1.5 & sfagncolor = 'firebrick' ; 'salmon'

    levels = [0.5,0.75,0.95]
    cannotation = ['50%','75%','95%']

    snrcut1 = 2.0
    
; ---------------------------------------------------------------------------    
; distribution of light fractions and physical aperture size 
    psfile = pspath+'sings_b_fraction_hist'+suffix
    im_plotconfig, 12, pos, psfile=psfile, width=4.2*[1,1], $
      height=4.2, charsize=1.8, charthick=3.0, xmargin=[1.1,1.1], $
      xspace=0.1

    nucindx = where((result.nuclear_fraction_b[0] gt -900.0),nnucindx)
    d20indx = where((result.drift20_fraction_b[0] gt -900.0),nd20indx)
    d56indx = where((result.drift56_fraction_b[0] gt -900.0),nd56indx)
    indx = where((result.nuclear_fraction_b[0] gt -900.0) or $
      (result.drift20_fraction_b[0] gt -900.0) or $
      (result.drift56_fraction_b[0] gt -900.0),nindx)
    
    xnuc = result[nucindx].nuclear_fraction_b[0]*100
    xd20 = result[d20indx].drift20_fraction_b[0]*100
    xd56 = result[d56indx].drift56_fraction_b[0]*100

    anuc = alog10(result[nucindx].nuclear_aperture_kpc[0])
    ad20 = alog10(result[d20indx].drift20_aperture_kpc[0])
    ad56 = alog10(result[d56indx].drift56_aperture_kpc[0])

; for the radial-strip spectra compute the *rectangular* aperture size
    junk = strtrim(repstr(repstr(result[d56indx].drift56_aperture_kpc_tex,'$',''),'\times',''),2)
    ad56_x = d56indx*0.0 & ad56_y = d56indx*0.0
    for jj = 0, nd56indx-1 do begin
       ad56_x[jj] = (strsplit(junk[jj],' ',/extract))[0]
       ad56_y[jj] = (strsplit(junk[jj],' ',/extract))[1]
    endfor
;   niceprint, ad56_x, ad56_y

    splog, 'Light-fraction statistics (median and full range)'
    print, 'Nuclear: ', median(xnuc), minmax(xnuc), format='(A10,4F8.4)'
    print, 'Drift20: ', median(xd20), minmax(xd20), format='(A10,4F8.4)'
    print, 'Drift56: ', median(xd56), minmax(xd56), format='(A10,4F8.4)'

    splog, 'Aperture statistics (median and full range)'
    print, 'Nuclear: ', median(sqrt(10^anuc)), minmax(sqrt(10^anuc)), format='(A10,4F8.4)'
    print, 'Drift20: ', median(sqrt(10^ad20)), minmax(sqrt(10^ad20)), format='(A10,4F8.4)'
    print, 'Drift56: ', median(sqrt(10^ad56)), minmax(sqrt(10^ad56)), format='(A10,4F8.4)'

;   quant = [0.5,0.1,0.9]
;   splog, 'Light-fraction statistics (median and 90% range)'
;   qq = im_quantile(xnuc,quant=quant) & print, 'Nuclear: ', qq[0], qq[1], qq[2], format='(A10,4F8.4)'
;   qq = im_quantile(xd20,quant=quant) & print, 'Drift20: ', qq[0], qq[1], qq[2], format='(A10,4F8.4)'
;   qq = im_quantile(xd56,quant=quant) & print, 'Drift56: ', qq[0], qq[1], qq[2], format='(A10,4F8.4)'
;
;   splog, 'Aperture statistics (median and 90% range)'
;   qq = im_quantile(sqrt(10^anuc),quant=quant) & print, 'Nuclear: ', qq[0], qq[1], qq[2], format='(A10,4F8.4)'
;   qq = im_quantile(sqrt(10^ad20),quant=quant) & print, 'Drift20: ', qq[0], qq[1], qq[2], format='(A10,4F8.4)'
;   qq = im_quantile(sqrt(10^ad56),quant=quant) & print, 'Drift56: ', qq[0], qq[1], qq[2], format='(A10,4F8.4)'

;   quant = [0.05,0.16,0.25,0.5,0.75,0.84,0.95]    
;   splog, 'Light-fraction statistics (0.05,0.16,0.25,0.5,0.75,0.84,0.95): '
;   print, 'Nuclear: ', im_quantile(xnuc,quant=quant), format='(A10,7F8.4)'
;   print, 'Drift20: ', im_quantile(xd20,quant=quant), format='(A10,7F8.4)'
;   print, 'Drift56: ', im_quantile(xd56,quant=quant), format='(A10,7F8.4)'
;
;   splog, 'Aperture statistics (0.16,0.25,0.5,0.75,0.84): '
;   print, 'Nuclear: ', im_quantile(sqrt(10^anuc),quant=quant), format='(A10,7F8.4)'
;   print, 'Drift20: ', im_quantile(sqrt(10^ad20),quant=quant), format='(A10,7F8.4)'
;   print, 'Drift56: ', im_quantile(sqrt(10^ad56),quant=quant), format='(A10,7F8.4)'

    xstats = im_stats(ad56_x) & ystats = im_stats(ad56_y)
    print, '   Drift56 - mean rectangular size = '+$
      strtrim(string(xstats.mean,format='(F12.2)'),2)+$
      ' x '+strtrim(string(ystats.mean,format='(F12.2)'),2)+' kpc^2'

    xtitle1 = 'B-band Light Fraction (%)'
    xtitle2 = 'log Aperture (kpc^{2})'
    ytitle = 'Number of Galaxies'

    xrange1 = [-5.0,105.0]
    xrange2 = [-3.5,3.0]
    binsize1 = 5.0
    binsize2 = 0.2
    
    im_plothist, xnuc, bin=binsize1, xnucbin, ynucbin, /noplot, histmin=0.0, histmax=100.0
    im_plothist, xd20, bin=binsize1, xd20bin, yd20bin, /noplot, histmin=0.0, histmax=100.0
    im_plothist, xd56, bin=binsize1, xd56bin, yd56bin, /noplot, histmin=0.0, histmax=100.0
    yrange1 = [0.0,(max(ynucbin)>max(yd20bin))>max(yd56bin)]*[1.0,1.03]

    im_plothist, anuc, bin=binsize2, axnucbin, aynucbin, /noplot, histmin=xrange2[0], histmax=xrange2[1]
    im_plothist, ad20, bin=binsize2, axd20bin, ayd20bin, /noplot, histmin=xrange2[0], histmax=xrange2[1]
    im_plothist, ad56, bin=binsize2, axd56bin, ayd56bin, /noplot, histmin=xrange2[0], histmax=xrange2[1]
    yrange2 = [0.0,(max(aynucbin)>max(ayd20bin))>max(ayd56bin)]*[1.0,1.03]

; ####################    
; B-band light fraction    
    djs_plot, [0], [0], /nodata, xsty=5, ysty=5, xrange=xrange1, $
      yrange=yrange1, position=pos[*,0]

    im_plothist, xd56, bin=binsize1, histmin=0.0, histmax=100.0, /overplot, $
      line=3, /fill, /fline, forientation=135, fspacing=0.15, $
      fcolor=fsc_color('dark grey',10), thick=5
    im_plothist, xd20, bin=binsize1, /overplot, histmin=0.0, histmax=100.0, $
      line=5, /fill, /fline, forientation=45, fspacing=0.05, $
      fcolor=djs_icolor('grey'), thick=5;, color=djs_icolor('grey')
    im_plothist, xd56, bin=binsize1, histmin=0.0, histmax=100.0, /overplot, $
      line=3, /fill, /fline, forientation=45, fspacing=0.15, $
      fcolor=fsc_color('dark grey',10), thick=5
    gg = where(ynucbin gt 0,ngg)
    djs_oplot, [xnucbin[gg[0]]-binsize1,xnucbin[gg],xnucbin[gg[ngg-1]]+binsize1], $
      [0,ynucbin[gg],0], psym=10, thick=5
;   im_plothist, xnuc, bin=binsize1, /overplot, histmin=0.0, histmax=100.0, $
;     line=0, thick=7

    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle1, ytitle=ytitle, $
      /xsty, /ysty, xrange=xrange1, yrange=yrange1, position=pos[*,0]

    legend, ['Nuclear','Circumnuclear','Radial Strip'], line=[0,5,3], $
      /right, /top, box=0, charsize=1.5, pspacing=1.8, thick=7
;     color=[djs_icolor('default'),djs_icolor('grey'),fsc_color('dark grey',10)]
    
; ####################    
; projected aperture
    djs_plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, xrange=xrange2, $
      yrange=yrange2, position=pos[*,1]

    im_plothist, ad56, bin=binsize2, histmin=xrange2[0], histmax=xrange2[1], /overplot, $
      line=3, /fill, /fline, forientation=135, fspacing=0.15, $
      fcolor=fsc_color('dark grey',10), thick=5
    im_plothist, ad20, bin=binsize2, /overplot, histmin=xrange2[0], histmax=xrange2[1], $
      line=5, /fill, /fline, forientation=45, fspacing=0.05, $
      fcolor=djs_icolor('grey'), thick=5;, color=djs_icolor('grey')
    im_plothist, ad56, bin=binsize2, histmin=xrange2[0], histmax=xrange2[1], /overplot, $
      line=3, /fill, /fline, forientation=45, fspacing=0.15, $
      fcolor=fsc_color('dark grey',10), thick=5
    gg = where(aynucbin gt 0,ngg)
    djs_oplot, [axnucbin[gg[0]]-binsize2,axnucbin[gg],axnucbin[gg[ngg-1]]+binsize2], $
      [0,aynucbin[gg],0], psym=10, thick=5
;   im_plothist, anuc, bin=binsize2, /overplot, histmin=xrange2[0], histmax=xrange2[1], $
;     line=0, thick=7

    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle2, ytitle='', $
      /xsty, ysty=9, xrange=xrange2, yrange=yrange2, position=pos[*,1], $
      ytickname=replicate(' ',10)
    axis, /yaxis, ysty=1, yrange=yrange, ytitle=textoidl(ytitle)

    im_plotconfig, /psclose

stop    
    
; ---------------------------------------------------------------------------    
; color-magnitude diagram, coded by galaxy type

; SINGS    
    x = result.bvri_absmag[0]
    xerr = result.bvri_absmag_err[0]
    y = result.bv
    yerr = result.bv_err

    early = where((result.t le 0.0),nearly) ; E-->S0
    mid   = where((result.t gt 0.0) and (result.t le 4.0),nmid) ; Sa-->Sbc
    late  = where((result.t gt 5.0),nlate) ; Sc-->Im

    xtitle = textoidl('M_{B}')
    ytitle = textoidl('B - V')
    xrange = [-11.1,-23.3]
    yrange = [-0.4,1.1] ; B-V range

; SDSS
    xsdss = sdss.ubvrijhk_absmag[1]
    ysdss = sdss.ubvrijhk_absmag[1]-sdss.ubvrijhk_absmag[2]

    psfile = pspath+'sings_cmd'+suffix
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.9, $
      xmargin=[1.4,0.3], width=6.7, height=6.7

; plot a random 30% subset of the outliers
    hogg_scatterplot, xsdss, ysdss, position=pos, xsty=1, ysty=1, $
      outliers=1, /internal, label=1, frac_outliers=0.25, $
      outpsym=6, outsymsize=0.01, outcolor='grey', $
      levels=levels, cannotation=cannotation, xtitle=xtitle, ytitle=ytitle, $
      xrange=xrange, yrange=yrange

    oploterror, x[mid], y[mid], xerr[mid], yerr[mid], $
      color=fsc_color(midcolor1,100), errcolor=fsc_color(midcolor1,100), $
      errthick=errthick1, psym=symcat(midsym1,thick=symthick1), symsize=midpsize1
    oploterror, x[late], y[late], xerr[late], yerr[late], $
      color=fsc_color(latecolor1,101), errcolor=fsc_color(latecolor1,101), $
      errthick=errthick1, psym=symcat(latesym1,thick=symthick1), symsize=latepsize1
    oploterror, x[early], y[early], xerr[early], yerr[early], $
      color=fsc_color(earlycolor1,102), errcolor=fsc_color(earlycolor1,102), $
      errthick=errthick1, psym=symcat(earlysym1,thick=symthick1), symsize=earlypsize1

    im_legend, ['E/S0','Sa-Sbc','Sc-Im'], psym=[earlysym1,midsym1,latesym1], $
      /left, /top, box=0, fill=[1,1,1], spacing=1.8, symsize=[1.8,1.5,2.1], $
      color=[earlycolor1,midcolor1,latecolor1], symthick=symthick1

    im_plotconfig, /psclose

; ---------------------------------------------------------------------------    
; light-fraction vs ratio of SF to AGN

    psfile = pspath+'sings_b_fraction_vs_class_fraction'+suffix
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8, $
      height=6.0

    xtitle = 'B-band Light Fraction (%)'
    ytitle = 'Fraction of Galaxies'

    binsize = 20.0
    xrange = [-2.0,102.0]
;   xrange = [0.1,200.0]
    yrange = [-0.05,1.0]
    
    nucindx = where((result.nuclear_fraction_b[0] gt -999.0) and $
      (strtrim(result.nuclear_class,2) ne '?') and $
      (strtrim(result.nuclear_class,2) ne '\nodata'))
    d20indx = where((result.drift20_fraction_b[0] gt -999.0) and $
      (strtrim(result.drift20_class,2) ne '?') and $
      (strtrim(result.drift20_class,2) ne '\nodata'))
    d56indx = where((result.drift56_fraction_b[0] gt -999.0) and $
      (strtrim(result.drift56_class,2) ne '?') and $
      (strtrim(result.drift56_class,2) ne '\nodata'))

    class = strtrim([result[nucindx].nuclear_class,$
      result[d20indx].drift20_class,result[d56indx].drift56_class],2)
    fraction = 100.0*[result[nucindx].nuclear_fraction_b[0],$
      result[d20indx].drift20_fraction_b[0],result[d56indx].drift56_fraction_b[0]]

    sf = where(class eq 'SF',nsf)
    agn = where(class eq 'AGN',nagn)
    sfagn = where(class eq 'SF/AGN',nsfagn)

    im_plothist, fraction[sf], xsf, ysf, bin=binsize, $
      histmin=0.0, histmax=99.0, edge=0L, /noplot
    im_plothist, fraction[agn], xagn, yagn, bin=binsize, $
      histmin=0.0, histmax=99.0, edge=0L, /noplot
    im_plothist, fraction[sfagn], xsfagn, ysfagn, bin=binsize, $
      histmin=0.0, histmax=99.0, edge=0L, /noplot

    norm = total([[yagn],[ysf],[ysfagn]],2) ; normalization constant
    niceprint, xagn, yagn, ysf, ysfagn, norm;, yagnerr, ysferr, ysfagnerr

    ysferr_hi = ysf*0.0 & ysferr_lo = ysf*0.0
    for jj = 0L, n_elements(ysf)-1L do begin
       limits = im_poisson_limits(ysf[jj],0.8413)
       ysferr_lo[jj] = ysf[jj]-limits[0] & ysferr_hi[jj] = limits[1]-ysf[jj]
    endfor
;   niceprint, ysf, ysferr_lo, ysferr_hi, sqrt(ysf)
    ysfagnerr_hi = ysfagn*0.0 & ysfagnerr_lo = ysfagn*0.0
    for jj = 0L, n_elements(ysfagn)-1L do begin
       limits = im_poisson_limits(ysfagn[jj],0.8413)
       ysfagnerr_lo[jj] = ysfagn[jj]-limits[0] & ysfagnerr_hi[jj] = limits[1]-ysfagn[jj]
    endfor
;   niceprint, ysfagn, ysfagnerr_lo, ysfagnerr_hi, sqrt(ysfagn)
    yagnerr_hi = yagn*0.0 & yagnerr_lo = yagn*0.0
    for jj = 0L, n_elements(yagn)-1L do begin
       limits = im_poisson_limits(yagn[jj],0.8413)
       yagnerr_lo[jj] = yagn[jj]-limits[0] & yagnerr_hi[jj] = limits[1]-yagn[jj]
    endfor
;   niceprint, yagn, yagnerr_lo, yagnerr_hi, sqrt(yagn)
    
;   ysferr = sqrt(ysf)/(norm+(norm eq 0.0))
;   yagnerr = sqrt(yagn)/(norm+(norm eq 0.0))
;   ysfagnerr = sqrt(ysfagn)/(norm+(norm eq 0.0))

    ysferr_hi = ysferr_hi/(norm+(norm eq 0.0))
    ysferr_lo = ysferr_lo/(norm+(norm eq 0.0))
    ysfagnerr_hi = ysfagnerr_hi/(norm+(norm eq 0.0))
    ysfagnerr_lo = ysfagnerr_lo/(norm+(norm eq 0.0))
    yagnerr_hi = yagnerr_hi/(norm+(norm eq 0.0))
    yagnerr_lo = yagnerr_lo/(norm+(norm eq 0.0))

    ysf = ysf/(norm+(norm eq 0.0))
    yagn = yagn/(norm+(norm eq 0.0))
    ysfagn = ysfagn/(norm+(norm eq 0.0))

    xsferr = replicate(binsize/2.0,nsf)
    xagnerr = replicate(binsize/2.0,nagn)
    xsfagnerr = replicate(binsize/2.0,nsfagn)
    
    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,0];, /xlog
    djs_oplot, !x.crange, [0,0], line=1
    
    xoff = 0.5 ; 0.5

    with_xerr = 1
    if keyword_set(with_xerr) then xxx = 1.0 else xxx = 0.0
; asymmetric error bars    
    oploterror, xsf, ysf, xsferr*xxx, ysferr_hi, color=fsc_color(sfcolor,100), $
      errcolor=fsc_color(sfcolor,100), errthick=errthick2, /hibar, $
      psym=-symcat(sfsym,thick=symthick2), symsize=2.7*sfpsize
    oploterror, xsf, ysf, xsferr*xxx, ysferr_lo, ps=-8, color=fsc_color(sfcolor,100), $
      errcolor=fsc_color(sfcolor,100), errthick=errthick2, /lobar, $
      psym=-symcat(sfsym,thick=symthick2), symsize=2.7*sfpsize

    oploterror, xsfagn-xoff, ysfagn, xsfagnerr*xxx, ysfagnerr_hi, color=fsc_color(sfagncolor,100), $
      errcolor=fsc_color(sfagncolor,100), errthick=errthick2, /hibar, $
      psym=-symcat(sfagnsym,thick=symthick2), symsize=2.8*sfagnpsize
    oploterror, xsfagn-xoff, ysfagn, xsfagnerr*xxx, ysfagnerr_lo, ps=-8, color=fsc_color(sfagncolor,100), $
      errcolor=fsc_color(sfagncolor,100), errthick=errthick2, /lobar, $
      psym=-symcat(sfagnsym,thick=symthick2), symsize=2.8*sfagnpsize

    oploterror, xagn+xoff, yagn, xagnerr*xxx, yagnerr_hi, color=fsc_color(agncolor,100), $
      errcolor=fsc_color(agncolor,100), errthick=errthick2, /hibar, $
      psym=-symcat(agnsym,thick=symthick2), symsize=3.0*agnpsize
    oploterror, xagn+xoff, yagn, xagnerr*xxx, yagnerr_lo, color=fsc_color(agncolor,100), $
      errcolor=fsc_color(agncolor,100), errthick=errthick2, /lobar, $
      psym=-symcat(agnsym,thick=symthick2), symsize=3.0*agnpsize

; symmetric error bars
;   im_symbols, sfsym, psize=2.7*sfpsize, color=fsc_color(sfcolor,100), fill=1, thick=errthick2
;   oploterror, xsf, ysf, xsferr*xxx, ysferr, ps=-8, color=fsc_color(sfcolor,100), $
;     errcolor=fsc_color(sfcolor,100), errthick=errthick2, thick=errthick2
;   im_symbols, agnsym, psize=3.0*agnpsize, color=fsc_color(agncolor,100), fill=1, thick=errthick2
;   oploterror, xagn+xoff, yagn, xagnerr*xxx, yagnerr, ps=-8, color=fsc_color(agncolor,100), $
;     errcolor=fsc_color(agncolor,100), errthick=errthick2, thick=errthick2
;   im_symbols, sfagnsym, psize=2.1*sfagnpsize, color=fsc_color(sfagncolor,100), fill=1, thick=errthick2
;   oploterror, xsfagn-xoff, ysfagn, xsfagnerr*xxx, ysfagnerr, ps=-8, color=fsc_color(sfagncolor,100), $
;     errcolor=fsc_color(sfagncolor,100), errthick=errthick2, thick=errthick2

;   good = where(ysf gt 0.0,ngood)
;   im_symbols, sfsym, psize=1.4*sfpsize, color=fsc_color(sfcolor,100), fill=1
;   if (ngood ne 0L) then oploterror, xsf[good], ysf[good], replicate(binsize,ngood)/2.0, ysferr[good], ps=-8, $
;     color=fsc_color(sfcolor,100), errcolor=fsc_color(sfcolor,100), errthick=errthick1, thick=postthick1
;   good = where(yagn gt 0.0,ngood)
;   im_symbols, agnsym, psize=1.2*agnpsize, color=fsc_color(agncolor,100), fill=1
;   if (ngood ne 0L) then oploterror, xagn[good], yagn[good], replicate(binsize,ngood)/2.0, yagnerr[good], ps=-8, $
;     color=fsc_color(agncolor,100), errcolor=fsc_color(agncolor,100), errthick=errthick1, thick=postthick1
;   good = where(ysfagn gt 0.0,ngood)
;   im_symbols, sfagnsym, psize=1.2*sfagnpsize, color=fsc_color(sfagncolor,100), fill=1
;   if (ngood ne 0L) then oploterror, xsfagn[good], ysfagn[good], replicate(binsize,ngood)/2.0, ysfagnerr[good], ps=-8, $
;     color=fsc_color(sfagncolor,100), errcolor=fsc_color(sfagncolor,100), errthick=errthick1, thick=postthick1

    im_legend, ['SF','AGN','SF/AGN'], psym=[sfsym,agnsym,sfagnsym], /left, /top, $
      box=0, fill=[1,1,1], spacing=1.8, symsize=1.5, symthick=symthick2, $
      color=[sfcolor,agncolor,sfagncolor]

    im_plotconfig, /psclose

; ---------------------------------------------------------------------------
; BPT diagram       

    psfile = pspath+'sings_niiha_vs_oiiihb'+suffix
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.6,1.1]

    xtitle = textoidl('log ([N II] \lambda6584/H\alpha)')
    ytitle = textoidl('log ([O III] \lambda5007/H\beta)')
    xrange = [-1.4,0.5]     ; [-2.0,0.8]
    yrange = [-1.1,1.1] 
    
; SDSS
    good = where($
      (sdss.nii_6584[0]/sdss.nii_6584[1] gt 5.0) and $
      (sdss.oiii_5007[0]/sdss.oiii_5007[1] gt 5.0) and $
      (sdss.h_beta[0]/sdss.h_beta[1] gt 5.0) and $
      (sdss.h_alpha[0]/sdss.h_alpha[1] gt 5.0),ngood)

    xsdss = alog10(sdss[good].nii_6584[0]/sdss[good].h_alpha[0])
    ysdss = alog10(sdss[good].oiii_5007[0]/sdss[good].h_beta[0])
    
    hogg_scatterplot, xsdss, ysdss, position=pos, xsty=1, ysty=1, $
      outliers=1, internal=1, label=0, $
      outpsym=3, outsymsize=0.05, outcolor=djs_icolor('grey'), $
      levels=levels, cannotation=cannotation, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange

; nuclear
    lineratio, nuclear_all, 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
      nuc_x, nuc_xerr, nuc_y, nuc_yerr, index=nuc_indx, snrcut=snrcut1
    oploterror, nuc_x, nuc_y, nuc_xerr, nuc_yerr, $
      psym=symcat(nucsym), symsize=nucpsize, color=fsc_color(nuccolor,100), $
      errthick=errthick1, errcolor=fsc_color(nuccolor,100)

; circumnuclear    
    lineratio, drift20_all, 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
      d20_x, d20_xerr, d20_y, d20_yerr, index=d20_indx, snrcut=snrcut1
    oploterror, d20_x, d20_y, d20_xerr, d20_yerr, $
      psym=symcat(d20sym), symsize=d20psize, color=fsc_color(d20color,100), $
      errthick=errthick1, errcolor=fsc_color(d20color,100)

; radial-strip
    lineratio, drift56_all, 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
      d56_x, d56_xerr, d56_y, d56_yerr, index=d56_indx, snrcut=snrcut1
    oploterror, d56_x, d56_y, d56_xerr, d56_yerr, $
      psym=symcat(d56sym), symsize=d56psize, color=fsc_color(d56color,100), $
      errthick=errthick1, errcolor=fsc_color(d56color,100)

; overplot the mixing lines
    models = kewley_bpt_lines(/kauffmann,_extra=extra)
    oplot, models.x_nii, models.y_nii, line=0
    models = kewley_bpt_lines(_extra=extra)
    oplot, models.x_nii, models.y_nii, line=2

; overplot the BPT/D coordinate system
;   djs_oplot, [-0.45,0.5], [-0.5,!y.crange[1]], line=1, thick=postthick4
;   djs_oplot, !x.crange, -0.5*[1,1], line=1, thick=postthick1
;   djs_oplot, -0.45*[1,1], !y.crange, line=1, thick=postthick1
    
; label the three regions
    xyouts, -0.9, -0.45, 'SF', align=0.5, charsize=1.6
    xyouts, -0.03, -0.9, 'SF/AGN', align=0.5, charsize=1.6
    xyouts,  -0.4, 0.85, 'AGN',     align=0.5, charsize=1.6
    
    im_legend, ['Nuclear','Circumnuclear','Radial-Strip'], psym=[17,16,15], $
      /left, /bottom, box=0, fill=[1,1,1], charsize=1.6, $
      spacing=1.8, symsize=1.5, color=[nuccolor,d20color,d56color]

    im_plotconfig, /psclose

return
end
    
