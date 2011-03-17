pro talk_08apr_ucsc, postscript=postscript, pdf=pdf
; jm08apr09nyu - plots for my colloquium talk

    common ages_common, ages, ancillary

; default plotting variables - note that POSTSCRIPT over-rides PDF

    postthick1 = 2.0
    postthick2 = 2.0
    postthick3 = 2.0
    postthick4 = 2.0
    charsize1 = 2.0
    charsize2 = 2.0
    charsize3 = 1.5
    charsize4 = 1.8
    textcolor1 = 'white'

    erbcolor = 'forest green'
    liucolor_loz = 'orange'
    liucolor_hiz = 'red'
    shapleycolor_loz = 'orange'
    shapleycolor_hiz = 'red'
    maiercolor = liucolor_hiz   ; same color

    pspath = getenv('PAPERSPATH')+'/literature/'
    datapath2 = pspath+'data/'
    datapath1 = '~/home/research/talks/2008/08apr_ucsc/'

    if keyword_set(pdf) then begin
       postscript = 0L
       pspath = datapath1+'keynote/' ; for keynote presentations
       textcolor1 = 'white'
       charsize1 = 2.2
       charsize2 = 1.6
       charsize3 = 1.2
       charsize4 = 1.4
       postthick1 = 8.0
       postthick2 = 6.0
       postthick3 = 4.0
       postthick4 = 10.0
       erbcolor = 'forest green'
       liucolor_loz = 'orange'
       liucolor_hiz = 'red'
       shapleycolor_loz = 'orange'
       shapleycolor_hiz = 'red'
       maiercolor = liucolor_hiz ; same color
    endif
    if keyword_set(postscript) then begin
       pspath = datapath1
       textcolor1 = 'black'
       postthick1 = 4.0
       postthick2 = 3.0
       postthick3 = 3.0
       postthick4 = 8.0
    endif

    charsize_0 = 1.0
    charsize_1 = 1.1
    charsize_2 = 1.2
    charsize_3 = 1.3
    charsize_4 = 1.4
    charsize_5 = 1.5
    charsize_6 = 1.6
    charsize_7 = 1.7
    charsize_8 = 1.8
    charsize_9 = 1.9
    singlecharsize_0 = 2.0
    singlecharsize_1 = 2.1
    singlecharsize_2 = 2.2
    singlecharsize_3 = 2.3
    singlecharsize_4 = 2.4
    singlecharsize_5 = 2.5
    charsize_30 = 3.0

    ohaxis = findgen((10.0-5.0)/0.01+1)*0.01+5.0
    
    if keyword_set(pdf) then begin
       hoyos05color  = 'forest green' & hoyos05sym = 125 & hoyos05psize = 1.2
       sava05color   = 'blue violet'  & sava05sym  = 115 & sava05psize  = 1.8
       lilly03color  = 'orchid'       & lilly03sym = 109 & lilly03psize = 1.4
       kob04color    = 'orange'       & kob04sym   = 108 & kob04psize   = 0.8
    endif else begin
       hoyos05color = 'forest green' & hoyos05sym = 125 & hoyos05psize = 1.2
       sava05color  = 'blue violet'  & sava05sym  = 115 & sava05psize  = 1.8
       lilly03color = 'orchid'       & lilly03sym = 109 & lilly03psize = 1.4
       kob04color   = 'orange'       & kob04sym   = 108 & kob04psize   = 0.8
    endelse

    if keyword_set(pdf) then begin
       color_normal = 'red'
       color_xray = 'orchid'
       color_nev = 'forest green'
    endif else begin
       color_normal = 'red'
       color_xray = 'dodger blue'
       color_nev = 'forest green'
    endelse
    
    mgtitle = 'M_{0.1g}'
    masstitle = 'log (M_{*}/M'+sunsymbol()+')'
    ohtitle1 = '12 + log (O/H)_{EW}'
    ohtitle2 = '12 + log (O/H)'

    light = 2.99792458D10       ; speed of light [cm/s]

; read some of the data we're going to need

    sdssalldust = read_sdss_main(/ispec)
    sdsskcorr = read_sdss_mz_sample(/mzhii_kcorr)
    sdssohdust = read_sdss_mz_sample(/mzhii_log12oh)
    sdssohnodust = read_sdss_mz_sample(/nodust_mzhii_log12oh)
    ageskcorr = read_ages_mz_sample(/mzhiiplus_ancillary)
    agesidust = read_ages_mz_sample(/mzhiiplus_ispec)
    agesohdust = read_ages_mz_sample(/mzhiiplus_log12oh)
    
; the first mass bin is an upper limit    
    
    erb_all = rsex(datapath2+'06erb.sex')
    ul = where(erb_all.nii lt 0.0,comp=good)
    erb_ul = erb_all[ul]
    erb = erb_all[good]
    
    liu = rsex(datapath2+'08liu_table2.sex')
    liu_n2stack = rsex(datapath2+'08liu_table3.sex')
    liu_o3n2stack = rsex(datapath2+'08liu_table4.sex')
    
    maier = rsex(datapath2+'06maier.sex')

    shapley = rsex(datapath2+'05shapley.sex')

; choose the spectrum to plot; other good candidates: 32, 74, 76, 93,
; 104, 152; choose the spectrum here so that we can plot its
; metallicity on the R23 plot, below
    
    if (n_elements(ages) eq 0L) then $
      ages = read_ages_mz_sample(/mzhiiplus_ispec)
    if (n_elements(ancillary) eq 0L) then $
      ancillary = read_ages_mz_sample(/mzhiiplus_ancillary)

    indx = where(ages.z gt 0.5 and ages.z lt 0.6 and ancillary.main_flag and $
      ages.oii_3727_ew[0] gt 15.0)
    thisages = indx[152]

;   cand = indx[[32,74,76,93,104,152]]
;   niceprint, ages[cand].galaxy, ages[cand].ages_id, ages[cand].h_beta_ew[0], $
;     ages[cand].oii_3727_ew[0], ages[cand].oiii_5007_ew[0]
;   ages_309/127         33593        8.89589        18.2673        3.81152
;   ages_419/082         30279        12.1940        27.7204        7.16050
;   ages_419/255         34837        7.91047        21.1084        5.18883
;   ages_422/270         37914        14.8120        26.7544        11.2311
;   ages_507/224         15748        10.7142        24.9732        6.43218
;   ages_713/099         37454        17.1949        42.0125        17.4742

; read some representative Pegase models

    splog, 'Reading the Pegase models'
    pegpath = getenv('PEGASE_HR_SFHGRID_DIR')+'/MEASURE/'
    peg1 = mrdfits(pegpath+'salp_tau_001.0Gyr.info.fits',1,silent=0)
    peg2 = mrdfits(pegpath+'salp_tau_001.0Gyr.info.fits',1,silent=0)
    peg3 = mrdfits(pegpath+'salp_tau_003.0Gyr.info.fits',1,silent=0)
    pegconst = mrdfits(pegpath+'salp_tau_999.0Gyr.info.fits',1,silent=0)

    rev = reverse(sort(peg1.age)) ; reverse the time array!
    peg1 = peg1[rev]
    peg2 = peg2[rev]
    peg3 = peg3[rev]
    pegconst = pegconst[rev]
    
    peg1_mgalaxy = 2D11 & peg1_zform = 1.0
    peg3_mgalaxy = 5D10 & peg3_zform = 1.2
    pegconst_mgalaxy = 1D10 & pegconst_zform = 2.0

    peg1_zaxis = getredshift(peg1.age/1D3+getage(peg1_zform))
    peg3_zaxis = getredshift(peg3.age/1D3+getage(peg3_zform))
    pegconst_zaxis = getredshift(pegconst.age/1D3+getage(pegconst_zform))

    peg1_good = where(peg1_zaxis gt 0.0)
    peg1 = peg1[peg1_good] & peg1_zaxis = peg1_zaxis[peg1_good]
    peg3_good = where(peg3_zaxis gt 0.0)
    peg3 = peg3[peg3_good] & peg3_zaxis = peg3_zaxis[peg3_good]
    pegconst_good = where(pegconst_zaxis gt 0.0)
    pegconst = pegconst[pegconst_good] & pegconst_zaxis = pegconst_zaxis[pegconst_good]

;   pegzform1 = 1.0
;   peg_zaxis1 = getredshift(peg1.age/1D3+getage(pegzform1))
;   peg_lookback1 = getage(0.0)-getage(peg_zaxis1)
    
;   pegzform2 = 1.0
;   peg_zaxis2 = getredshift(peg1.age/1D3+getage(pegzform2))
;   peg_lookback2 = getage(0.0)-getage(peg_zaxis2)
    
;   peg_good1 = where(peg_zaxis1 gt 0.0)
;   peg_zaxis1 = peg_zaxis1[peg_good1]
;   peg_lookback1 = peg_lookback1[peg_good1]
;   peg1 = peg1[peg_good1]
;   peg2 = peg2[peg_good1]
;   peg3 = peg3[peg_good1]
;   pegconst = pegconst[peg_good1]

; ---------------------------------------------------------------------------    
; compare MZ [R23] and MZ [N2]
; ---------------------------------------------------------------------------    
    
    xsize = 8.5 & ysize = 7.0
    pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, width=7.0, height=5.5, $
      xmargin=[1.2,0.3], ymargin=[0.4,1.1], xpage=xsize, ypage=ysize, $
      position=pos, /normal

    psname = 'mz_r23_vs_mz_n2'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xsize, ysize=ysize, silent=keyword_set(pdf)

    xtitle = 'log (M_{*}/M'+sunsymbol()+')'
    ytitle = '12 + log (O/H)'
    
    xrange = [8.5,11.3]
    yrange = [8.25,9.25]

    sindx = where((sdssohnodust.zstrong_12oh_kk04 gt -900.0) and $
      (strtrim(sdssohnodust.r23branch_kk04,2) eq 'U') and $
      (sdssohdust.zstrong_niiha gt -900))
    smass = sdsskcorr[sindx].mass ; Chabrier IMF!
    soh_n2 = 8.90+0.57*sdssohdust[sindx].zstrong_niiha
    soh_r23 = sdssohnodust[sindx].zstrong_12oh_kk04
    jj = im_stats(soh_r23-soh_n2,/ver)

    im_hogg_scatterplot, smass, soh_n2, outliers=plot_outliers, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor='default', levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,100), $
      contour_color=fsc_color('orange',101), /nogrey, c_linestyle=0, cthick=postthick2

    oplot, [xrange[0],9.0], 8.66*[1,1], line=2, thick=postthick4, color=fsc_color(textcolor1,101)
    xyouts, 8.6, 8.6, textoidl('(O/H)_{\odot}'), align=0.0, charsize=1.5, $
      charthick=postthick2, color=fsc_color(textcolor1,100)

    im_hogg_scatterplot, smass, soh_r23, /noerase, outliers=plot_outliers, xsty=5, ysty=5, $
      label=0, outpsym=8, outsymsize=1.0, outcolor='default', levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle='', ytitle='', $
      xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      contour_color=fsc_color('dodger blue',102), /nogrey, c_linestyle=0, cthick=postthick2

    label = ['Using ([O II]+[O III])/H\beta','Using [N II]/H\alpha']
;   label = ['log {([O II]+[O III])/H\beta}','log ([N II]/H\alpha)']
;   label = ['R_{23}','N2']

    xyouts, 10.78, 8.35, textoidl(label[0]), charsize=1.3, charthick=postthick2, $
      align=0.5, /data, color=fsc_color('dodger blue',101)
    xyouts, 10.78, 8.3, textoidl(label[1]), charsize=1.3, charthick=postthick2, $
      align=0.5, /data, color=fsc_color('orange',102)
    
;   legend, textoidl(label), /right, /bottom, charsize=1.3, $
;     charthick=postthick2, thick=postthick4, box=0, $
;     color=fsc_color(['dodger blue','orange'],[101,102]), $
;     textcolor=fsc_color(['dodger blue','orange'],[101,102]), $
;     spacing=2.2, margin=0

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    
    
; ---------------------------------------------------------------------------    
; mass-metallicity relation    
; ---------------------------------------------------------------------------    
    
    xsize = 8.5 & ysize = 7.0
    pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, width=7.0, height=5.5, $
      xmargin=[1.2,0.3], ymargin=[0.4,1.1], xpage=xsize, ypage=ysize, $
      position=pos, /normal

    psname = 'mz_06erb'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xsize, ysize=ysize, silent=keyword_set(pdf)

    xtitle = 'log (M_{*}/M'+sunsymbol()+')'
    ytitle = '12 + log (O/H)_{N2}'
    
    xrange = [8.6,11.5]
    yrange = [8.1,8.85]

; plot the sdss galaxies

    sindx = where((sdssohdust.zstrong_niiha gt -900))
    smass = sdsskcorr[sindx].mass ; Chabrier IMF!
    soh = 8.90+0.57*sdssohdust[sindx].zstrong_niiha
    
    im_hogg_scatterplot, smass, soh, outliers=plot_outliers, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor='default', levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,100), $
      contour_color=fsc_color(textcolor1,100), /nogrey, c_linestyle=0, cthick=postthick2

;   plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
;     xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
;     charsize=charsize, charthick=postthick2, xsty=1, ysty=1, xrange=xrange, $
;      position=pos[*,0], yrange=yrange, color=fsc_color(textcolor1,150)

    plot_errors = 0

    if keyword_set(plot_errors) then begin
       im_symbols, 105, psize=2.5, thick=postthick3, fill=1, color=fsc_color(erbcolor,150)
       oploterror, alog10(1D10*erb.mass), erb.log12oh, erb.mass_err/erb.mass/alog(10.0), erb.log12oh_err, $
         thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(erbcolor,150), $
         errcolor=fsc_color(erbcolor,150)
       oploterror, alog10(1D10*erb_ul.mass), -erb_ul.log12oh, erb_ul.mass_err/erb_ul.mass/alog(10.0), $
         erb_ul.log12oh_err, thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(erbcolor,150), $
         errcolor=fsc_color(erbcolor,150)
    endif else begin
       im_symbols, 105, psize=3.0, thick=postthick3, fill=1, color=fsc_color(erbcolor,150)
       oplot, alog10(1D10*erb.mass), erb.log12oh, thick=postthick3, psym=8, color=fsc_color(erbcolor,150)
       plots, alog10(1D10*erb_ul.mass), -erb_ul.log12oh, thick=postthick3, psym=8, color=fsc_color(erbcolor,150)
    endelse
    arrow, alog10(1D10*erb_ul.mass), -erb_ul.log12oh, alog10(1D10*erb_ul.mass), $
      -erb_ul.log12oh-0.05, /data, hsize=-0.6, color=fsc_color(erbcolor,150), $
      hthick=postthick3, thick=postthick3
    
; now plot Liu et al. (2008); use the N2 calibrator and differentiate
; between galaxies at z~1 and z~1.4; also overplot the stack

    plot_individual = 1

    if keyword_set(plot_individual) then begin
       
       loz = where((liu.mass gt -900.0) and (liu.log12oh_n2 gt -900.0) and (liu.z lt 1.1))
       hiz = where((liu.mass gt -900.0) and (liu.log12oh_n2 gt -900.0) and (liu.z gt 1.1))
       
       if keyword_set(plot_errors) then begin
          im_symbols, 106, psize=2.0, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
          oploterror, liu[loz].mass, liu[loz].log12oh_n2, liu[loz].mass_err, liu[loz].log12oh_n2_err, $
            thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_loz,150), $
            errcolor=fsc_color(liucolor_loz,150)
          im_symbols, 108, psize=2.0, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
          oploterror, liu[hiz].mass, liu[hiz].log12oh_n2, liu[hiz].mass_err, liu[hiz].log12oh_n2_err, $
            thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_hiz,150), $
            errcolor=fsc_color(liucolor_hiz,150)
       endif else begin
          im_symbols, 106, psize=1.3, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
          oplot, liu[loz].mass, liu[loz].log12oh_n2, thick=postthick3, psym=8, $
            color=fsc_color(liucolor_loz,150)
          im_symbols, 108, psize=1.5, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
          oplot, liu[hiz].mass, liu[hiz].log12oh_n2, thick=postthick3, psym=8, $
            color=fsc_color(liucolor_hiz,150)
       endelse

    endif
       
; stacked spectra    
    
    loz_stack = where((liu_n2stack.mass gt -900.0) and (liu_n2stack.log12oh_n2 gt -900.0) and (liu_n2stack.z lt 1.1))
    hiz_stack = where((liu_n2stack.mass gt -900.0) and (liu_n2stack.log12oh_n2 gt -900.0) and (liu_n2stack.z gt 1.1))
    
    if keyword_set(plot_errors) then begin
       im_symbols, 106, psize=3.0, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
       oploterror, liu_n2stack[loz_stack].mass, liu_n2stack[loz_stack].log12oh_n2, $
         liu_n2stack[loz_stack].mass_err, liu_n2stack[loz_stack].log12oh_n2_err, $
         thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_loz,150), $
         errcolor=fsc_color(liucolor_loz,150)
       im_symbols, 108, psize=3.0, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
       oploterror, liu_n2stack[hiz_stack].mass, liu_n2stack[hiz_stack].log12oh_n2, $
         liu_n2stack[hiz_stack].mass_err, liu_n2stack[hiz_stack].log12oh_n2_err, $
         thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_hiz,150), $
         errcolor=fsc_color(liucolor_hiz,150)
    endif else begin
       im_symbols, 106, psize=3.0, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
       oplot, liu_n2stack[loz_stack].mass, liu_n2stack[loz_stack].log12oh_n2, thick=postthick3, psym=8, $
         color=fsc_color(liucolor_loz,150)
       im_symbols, 108, psize=3.4, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
       oplot, liu_n2stack[hiz_stack].mass, liu_n2stack[hiz_stack].log12oh_n2, thick=postthick3, psym=8, $
         color=fsc_color(liucolor_hiz,150)
    endelse

;   im_symbols, 108, psize=1.5, thick=postthick3, fill=1, color=fsc_color(liucolor_o3n2,150)
;   oploterror, liu_o3n2.mass, liu_o3n2.log12oh_o3n2, liu_o3n2.mass_err, liu_o3n2.log12oh_o3n2_err, $
;     thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_o3n2,150), $
;     errcolor=fsc_color(liucolor_o3n2,150)

; Shapley et al. 2005

    if keyword_set(plot_individual) then begin

       loz = where((shapley.mass gt -900.0) and (shapley.log12oh_n2 gt -900.0) and (shapley.z lt 1.1))
       hiz = where((shapley.mass gt -900.0) and (shapley.log12oh_n2 gt -900.0) and (shapley.z gt 1.1))
       
       if keyword_set(plot_errors) then begin
          im_symbols, 106, psize=2.0, thick=postthick3, fill=1, color=fsc_color(shapleycolor_loz,150)
          oploterror, shapley[loz].mass, shapley[loz].log12oh_n2, shapley[loz].mass_err, shapley[loz].log12oh_n2_err, $
            thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(shapleycolor_loz,150), $
            errcolor=fsc_color(shapleycolor_loz,150)
          im_symbols, 108, psize=2.0, thick=postthick3, fill=1, color=fsc_color(shapleycolor_hiz,150)
          oploterror, shapley[hiz].mass, shapley[hiz].log12oh_n2, shapley[hiz].mass_err, shapley[hiz].log12oh_n2_err, $
            thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(shapleycolor_hiz,150), $
            errcolor=fsc_color(shapleycolor_hiz,150)
       endif else begin
          im_symbols, 106, psize=1.3, thick=postthick3, fill=1, color=fsc_color(shapleycolor_loz,150)
          oplot, shapley[loz].mass, shapley[loz].log12oh_n2, $
            thick=postthick3, psym=8, color=fsc_color(shapleycolor_loz,150)
          im_symbols, 108, psize=1.5, thick=postthick3, fill=1, color=fsc_color(shapleycolor_hiz,150)
          oplot, shapley[hiz].mass, shapley[hiz].log12oh_n2, $
            thick=postthick3, psym=8, color=fsc_color(shapleycolor_hiz,150)
       endelse

    endif    
    
; now overplot the median errors of either the individual points or
; the stacked spectra 

    plot_individual_median_errors = 0
    plot_stack_median_errors = 1

    if keyword_set(plot_individual_median_errors) then begin
       xpos = 10.7 & ypos = 8.3
       im_symbols, 105, psize=3.0, thick=postthick3, fill=1, color=fsc_color(erbcolor,150)
       oploterror, xpos, ypos, median(erb.mass_err/erb.mass/alog(10.0)), median(erb.log12oh_err), $
         thick=postthick3, errthick=postthick4, psym=8, color=fsc_color(erbcolor,150), $
         errcolor=fsc_color(erbcolor,150)
       im_symbols, 108, psize=2.5, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
       oploterror, xpos+0.26, ypos, median([liu[hiz].mass_err,shapley[hiz].mass_err]), $
         median([liu[hiz].log12oh_n2_err,shapley[hiz].log12oh_n2_err]), $
         thick=postthick3, errthick=postthick4, psym=8, color=fsc_color(liucolor_hiz,150), $
         errcolor=fsc_color(liucolor_hiz,150)
       im_symbols, 106, psize=2.5, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
       oploterror, xpos+0.6, ypos, median([liu[loz].mass_err,shapley[loz].mass_err]), $
         median([liu[loz].log12oh_n2_err,shapley[loz].log12oh_n2_err]), $
         thick=postthick3, errthick=postthick4, psym=8, color=fsc_color(liucolor_loz,150), $
         errcolor=fsc_color(liucolor_loz,150)
    endif
    if keyword_set(plot_stack_median_errors) then begin
       xpos = 10.7 & ypos = 8.25
       im_symbols, 105, psize=3.0, thick=postthick3, fill=1, color=fsc_color(erbcolor,150)
       oploterror, xpos+0.6, ypos, median(erb.mass_err/erb.mass/alog(10.0)), median(erb.log12oh_err), $
         thick=postthick3, errthick=postthick4, psym=8, color=fsc_color(erbcolor,150), $
         errcolor=fsc_color(erbcolor,150)
       im_symbols, 108, psize=2.5, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
       oploterror, xpos+0.3, ypos, median([liu_n2stack[hiz_stack].mass_err]), $
         median([liu_n2stack[hiz_stack].log12oh_n2_err]), $
         thick=postthick3, errthick=postthick4, psym=8, color=fsc_color(liucolor_hiz,150), $
         errcolor=fsc_color(liucolor_hiz,150)
       im_symbols, 106, psize=2.5, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
       oploterror, xpos, ypos, median([liu_n2stack[loz_stack].mass_err]), $
         median([liu_n2stack[loz_stack].log12oh_n2_err]), $
         thick=postthick3, errthick=postthick4, psym=8, color=fsc_color(liucolor_loz,150), $
         errcolor=fsc_color(liucolor_loz,150)
    endif
    
; legend    
    
    im_legend, textoidl(['z\sim1','z\sim1.5','z\sim2']), /left, /top, box=0, $
      charsize=singlecharsize_0, $
      charthick=postthick2, psym=[106,108,105], fill=[1,1,1], thick=postthick3, $
      color=fsc_color([liucolor_loz,liucolor_hiz,erbcolor],[101,102,103]), spacing=2.3, $
      textcolor=fsc_color(textcolor1,100), symsize=1.6

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ------------------------------------------------------------
; emission-line versus redshift: when do the various nebular emission
; lines appear in the various optical/near-IR bandpasses?
; ------------------------------------------------------------

    psname = 'emline_vs_redshift'
    xpage = 8.5 & ypage = 8.5
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.1, height=7.1, $
      xmargin=[1.1,0.3], ymargin=[0.3,1.1], xpage=xpage, ypage=ypage, $
      position=pos, /normal

    zmax = 5.3
    zmin = -0.05
    dz = 0.1
    zgrid = findgen((zmax-zmin)/dz+1)*dz+zmin

    ytitle = 'Redshift'

    xrange = [0.0,3.0]
    yrange = minmax(zgrid)
    
    lamgrid = findgen((26E4-3500.0)/1.0+1)*1.0+3500.0
    lamgrid = lamgrid / 1E4

    lamrest = [3727.0,mean([4861.0,5007]),mean([6548.0,6563.0,6584])] / 1E4 ; OII, H-beta, OIII, Ha
    linename = ['[O II]','H\beta+[O III]','H\alpha+[N II]']
    nline = n_elements(lamrest)

    nwindow = 4
    lamopt1 = 0.35 & lamopt2 = 0.95
    lamj1 = 1.1 & lamj2 = 1.4
    lamh1 = 1.5 & lamh2 = 1.8
    lamk1 = 2.0 & lamk2 = 2.3

; when are each of the emission lines within the atmospheric windows?

    zz1 = fltarr(nline,nwindow) & zz2 = fltarr(nline,nwindow)
    for iline = 0L, nline-1L do begin
; optical
       zz1[iline,0] = interpol(zgrid,lamrest[iline]*(1+zgrid),lamopt1)>0.0
       zz2[iline,0] = interpol(zgrid,lamrest[iline]*(1+zgrid),lamopt2)>0.0
; J-band
       zz1[iline,1] = interpol(zgrid,lamrest[iline]*(1+zgrid),lamj1)>0.0
       zz2[iline,1] = interpol(zgrid,lamrest[iline]*(1+zgrid),lamj2)>0.0
; H-band
       zz1[iline,2] = interpol(zgrid,lamrest[iline]*(1+zgrid),lamh1)>0.0
       zz2[iline,2] = interpol(zgrid,lamrest[iline]*(1+zgrid),lamh2)>0.0
; K-band
       zz1[iline,3] = interpol(zgrid,lamrest[iline]*(1+zgrid),lamk1)>0.0
       zz2[iline,3] = interpol(zgrid,lamrest[iline]*(1+zgrid),lamk2)>0.0
    endfor

; initialize the axes, but don't plot anything
    
    plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      xsty=5, ysty=5, position=pos[*,0];, xtitle='', ytitle='', $
;     xtickname=replicate(' ',10), ytickname=replicate(' ',10)

    xw = 0.2
    dx = 0.04
    
; #########################    
; [O II]
    xcen = 0.5
; optical
    polyfill, xcen+[-xw,-xw,+xw,+xw], [zz1[0,0],zz2[0,0],zz2[0,0],zz1[0,0]], $
      color=fsc_color('dodger blue',200), /fill
    xyouts, xcen+dx, mean([zz1[0,0],zz2[0,0]]), 'Optical', orientation=90, $
      align=0.5, charsize=2.0, charthick=postthick1, color=fsc_color(textcolor1,100)
; J
    polyfill, xcen+[-xw,-xw,+xw,+xw], [zz1[0,1],zz2[0,1],zz2[0,1],zz1[0,1]], $
      color=fsc_color('orange',201), /fill
    xyouts, xcen+dx, mean([zz1[0,1],zz2[0,1]]), 'J', orientation=90, $
      align=0.5, charsize=2.2, charthick=postthick1, color=fsc_color(textcolor1,100)
; H
    polyfill, xcen+[-xw,-xw,+xw,+xw], [zz1[0,2],zz2[0,2],zz2[0,2],zz1[0,2]], $
      color=fsc_color('forest green',202), /fill
    xyouts, xcen+dx, mean([zz1[0,2],zz2[0,2]]), 'H', orientation=90, $
      align=0.5, charsize=2.2, charthick=postthick1, color=fsc_color(textcolor1,100)
; K
    polyfill, xcen+[-xw,-xw,+xw,+xw], [zz1[0,3],zz2[0,3],zz2[0,3],zz1[0,3]], $
      color=fsc_color('firebrick',203), /fill
    xyouts, xcen+dx, mean([zz1[0,3],zz2[0,3]]), 'K', orientation=90, $
      align=0.5, charsize=2.2, charthick=postthick1, color=fsc_color(textcolor1,100)
; #########################    
; Hb+[OIII]
    xcen = 1.5
; optical
    polyfill, xcen+[-xw,-xw,+xw,+xw], [zz1[1,0],zz2[1,0],zz2[1,0],zz1[1,0]], $
      color=fsc_color('dodger blue',200), /fill
    xyouts, xcen+dx, mean([zz1[1,0],zz2[1,0]]), 'Optical', orientation=90, $
      align=0.5, charsize=2.0, charthick=postthick1, color=fsc_color(textcolor1,100)
; J
    polyfill, xcen+[-xw,-xw,+xw,+xw], [zz1[1,1],zz2[1,1],zz2[1,1],zz1[1,1]], $
      color=fsc_color('orange',201), /fill
    xyouts, xcen+dx, mean([zz1[1,1],zz2[1,1]]), 'J', orientation=90, $
      align=0.5, charsize=2.2, charthick=postthick1, color=fsc_color(textcolor1,100)
; H
    polyfill, xcen+[-xw,-xw,+xw,+xw], [zz1[1,2],zz2[1,2],zz2[1,2],zz1[1,2]], $
      color=fsc_color('forest green',202), /fill
    xyouts, xcen+dx, mean([zz1[1,2],zz2[1,2]]), 'H', orientation=90, $
      align=0.5, charsize=2.2, charthick=postthick1, color=fsc_color(textcolor1,100)
; K
    polyfill, xcen+[-xw,-xw,+xw,+xw], [zz1[1,3],zz2[1,3],zz2[1,3],zz1[1,3]], $
      color=fsc_color('firebrick',203), /fill
    xyouts, xcen+dx, mean([zz1[1,3],zz2[1,3]]), 'K', orientation=90, $
      align=0.5, charsize=2.2, charthick=postthick1, color=fsc_color(textcolor1,100)
; #########################    
; Ha+[NII]
    xcen = 2.5
; optical
    polyfill, xcen+[-xw,-xw,+xw,+xw], [zz1[2,0],zz2[2,0],zz2[2,0],zz1[2,0]], $
      color=fsc_color('dodger blue',200), /fill
    xyouts, xcen+dx, mean([zz1[2,0],zz2[2,0]]), 'Opt', orientation=90, $
      align=0.5, charsize=2.0, charthick=postthick1, color=fsc_color(textcolor1,100)
; J
    polyfill, xcen+[-xw,-xw,+xw,+xw], [zz1[2,1],zz2[2,1],zz2[2,1],zz1[2,1]], $
      color=fsc_color('orange',201), /fill
    xyouts, xcen+dx, mean([zz1[2,1],zz2[2,1]]), 'J', orientation=90, $
      align=0.5, charsize=2.2, charthick=postthick1, color=fsc_color(textcolor1,100)
; H
    polyfill, xcen+[-xw,-xw,+xw,+xw], [zz1[2,2],zz2[2,2],zz2[2,2],zz1[2,2]], $
      color=fsc_color('forest green',202), /fill
    xyouts, xcen+dx, mean([zz1[2,2],zz2[2,2]]), 'H', orientation=90, $
      align=0.5, charsize=2.2, charthick=postthick1, color=fsc_color(textcolor1,100)
; K
    polyfill, xcen+[-xw,-xw,+xw,+xw], [zz1[2,3],zz2[2,3],zz2[2,3],zz1[2,3]], $
      color=fsc_color('firebrick',203), /fill
    xyouts, xcen+dx, mean([zz1[2,3],zz2[2,3]]), 'K', orientation=90, $
      align=0.5, charsize=2.2, charthick=postthick1, color=fsc_color(textcolor1,100)

; now overplot the axes

    plot, [0], [0], /nodata, /noerase, xrange=xrange, yrange=yrange, $
      charsize=2.3, charthick=postthick2, xthick=postthick1, ythick=postthick1, $
      xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, position=pos[*,0], $
      color=fsc_color(textcolor1,100), xminor=-1, xticklen=1E-6, $
      xtickname=[' ',textoidl(linename[0]),' ',textoidl(linename[1]),' ',$
      textoidl(linename[2]),' ']

;   legend, textoidl(['Opt','J','H','K']), /right, /top, box=0, charsize=2.0, $
;     charthick=postthick2, spacing=2.2, thick=postthick4, $
;     color=fsc_color(['dodger blue','orange','forest green','firebrick'],[50,51,52,53]), $
;     textcolor=fsc_color(['dodger blue','orange','forest green','firebrick'],[50,51,52,53]), $
;     line=[0,0,0,0], fill=[1,1,1,1]

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    
    
; ------------------------------------------------------------
; wavelength versus redshift
; ------------------------------------------------------------

    psname = 'wavelength_vs_redshift'
    xpage = 8.5 & ypage = 8.5
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.1, height=7.1, $
      xmargin=[1.1,0.3], ymargin=[0.3,1.1], xpage=xpage, ypage=ypage, $
      position=pos, /normal

    zmax = 5.1
    zmin = 0.0
    dz = 0.1
    zgrid = findgen((zmax-zmin)/dz+1)*dz+zmin

    lamgrid = findgen((26E4-3500.0)/1.0+1)*1.0+3500.0
    lamgrid = lamgrid / 1E4

;   lamrest = [3727.0,4861.0,6563.0] / 1E4 ; OII, H-beta, OIII, Ha
;   linename = ['[O II]','Hb','Ha']
;   linestyle = [0,3,5]

    lamrest = [3727.0,mean([4861.0,5007]),mean([6563.0,6584])] / 1E4 ; OII, H-beta, OIII, Ha
    linename = ['[O II]','Hb','Ha']
    linestyle = [0,3,5]

;   lamrest = [3727.0,4340.0,4861.0,6563.0] / 1E4 ; OII, H-beta, OIII, Ha
;   linename = ['[O II]','Hg','Hb','Ha']
;   linestyle = [0,1,3,5]

    xtitle = textoidl('log \lambda_{obs} (\mu'+'m)')
    ytitle = 'Redshift'

    xrange = [1600.0,26E3] / 1E4
    yrange = minmax(zgrid)
    
    plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      charsize=2.0, charthick=postthick2, xthick=postthick1, ythick=postthick1, $
      xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, position=pos[*,0], $
      color=fsc_color(textcolor1,100)
    
; fill the optical wavelength range
    
    polyfill, [3500.0,3500.0,9500,9500] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('orchid',200), /fill
    polyfill, [3500.0,3500.0,9500,9500] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('orchid',200), /fill
    xyouts, 6500.0/1E4, 4.2, 'Optical', charsize=2.0, charthick=postthick2, $
      align=0.5, /data, color=fsc_color(textcolor1,100)

    polyfill, [11E3,11E3,14E3,14E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('dodger blue',201), /fill
    polyfill, [11E3,11E3,14E3,14E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('dodger blue',201), /fill
    xyouts, 12.5E3/1E4, 4.2, 'J', charsize=2.0, charthick=postthick2, $
      align=0.5, /data, color=fsc_color(textcolor1,100)

    polyfill, [15E3,15E3,18E3,18E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('forest green',202), /fill
    polyfill, [15E3,15E3,18E3,18E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('forest green',202), /fill
    xyouts, 16.5E3/1E4, 4.2, 'H', charsize=2.0, charthick=postthick2, $
      align=0.5, /data, color=fsc_color(textcolor1,100)

    polyfill, [20E3,20E3,23E3,23E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('orange',203), /fill
    polyfill, [20E3,20E3,23E3,23E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('orange',203), /fill

    xyouts, 21.5E3/1E4, 4.2, 'K', charsize=2.0, charthick=postthick2, $
      align=0.5, /data, color=fsc_color(textcolor1,100)

    for i = 0L, n_elements(lamrest)-1L do djs_oplot, lamrest[i]*(1+zgrid), $
      zgrid, line=linestyle[i], thick=postthick3, color=fsc_color(textcolor1,100)

    xyouts, 1.5, 2.1, textoidl('[O II]'), charsize=1.6, charthick=postthick2, $
      align=0.5, /data, color=fsc_color('firebrick',204), orientation=35
    xyouts, 1.02, 1.4, textoidl('H\beta+[O III]'), charsize=1.6, charthick=postthick2, $
      align=0.5, /data, color=fsc_color(textcolor1,100), orientation=35
    xyouts, 1.02, 0.7, textoidl('H\alpha+[N II]'), charsize=1.6, charthick=postthick2, $
      align=0.5, /data, color=fsc_color(textcolor1,100), orientation=45

; overplot the tick marks    
    
    plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, /noerase, $
      charsize=2.0, charthick=postthick2, xthick=postthick1, ythick=postthick1, $
      xsty=1, ysty=1, xtitle='', ytitle='', position=pos[*,0], color=fsc_color(textcolor1,100)
    
; legend

;   legend, textoidl(['[O II]','H\beta','H\alpha']), /right, $
;     /bottom, charsize=1.8, charthick=postthick, line=linestyle, thick=postthick, $
;     clear=keyword_set(postscript), spacing=1.5, box=0

; when are each of the emission lines within the atmospheric windows?

    for iline = 0L, n_elements(linename)-1L do begin
       print, linename[iline]+':'
       print, ' Optical: <'+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),0.95),format='(F12.1)'),2)
       print, ' J-band : '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.1),format='(F12.1)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.4),format='(F12.1)'),2)
       print, ' H-band : '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.5),format='(F12.1)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.8),format='(F12.1)'),2)
       print, ' K-band : '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),2.0),format='(F12.1)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),2.3),format='(F12.1)'),2)
       print
    endfor
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    
    
; ---------------------------------------------------------------------------    
; BPT diagram
; ---------------------------------------------------------------------------    
    
    xsize = 8.5 & ysize = 7.2
    pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, width=6.7, height=5.7, $
      xmargin=[1.5,0.3], ymargin=[0.4,1.1], xpage=xsize, ypage=ysize, $
      position=pos, /normal

    psname = 'bpt_06erb'        ; bad name!!!
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xsize, ysize=ysize, silent=keyword_set(pdf)

    xtitle = 'log ([N II] \lambda6584/H\alpha)'
    ytitle = 'log ([O III] \lambda5007/H\beta)'
    
    xrange = [-1.6,0.7]
    yrange = [-1.0,1.2]

; plot the sdss galaxies

    sdssindx = where((sdssalldust.nii_6584[0]/sdssalldust.nii_6584[1] gt 5.0) and $
      (sdssalldust.oiii_5007[0]/sdssalldust.oiii_5007[1] gt 5.0) and $
      (sdssalldust.h_alpha[0]/sdssalldust.h_alpha[1] gt 5.0) and $
      (sdssalldust.h_beta[0]/sdssalldust.h_beta[1] gt 5.0))
    
    sniiha = alog10(sdssalldust[sdssindx].nii_6584[0]/sdssalldust[sdssindx].h_alpha[0])
    soiiihb = alog10(sdssalldust[sdssindx].oiii_5007[0]/sdssalldust[sdssindx].h_beta[0])
    
    im_symbols, 108, psize=0.1, fill=1, color=fsc_color(textcolor1,100)
    im_hogg_scatterplot, sniiha, soiiihb, outliers=1, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor=fsc_color(textcolor1,100), $
      levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_1, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,100), $
      contour_color=fsc_color(textcolor1,100), /nogrey, c_linestyle=0, cthick=postthick2

; Kewley curve; label SF and AGN    
    
    models = kewley_bpt_lines(_extra=extra)
;   oplot, models.x_nii, models.y_nii, line=2, thick=postthick3, $
;     color=fsc_color(textcolor1,100)

    xyouts, -1.1, -0.6, 'Star-Forming', align=0.5, charsize=2.0, $
      charthick=postthick2, color=fsc_color(textcolor1,100)
    xyouts, 0.48, 0.7, 'AGN', align=0.5, charsize=2.0, $
      charthick=postthick2, color=fsc_color(textcolor1,100)

; Shapley/z~1.0    
    
    loz = where((shapley.z lt 1.1) and (shapley.oiii gt 0.0))

    niiha_loz = shapley[loz].nii/shapley[loz].ha
    niiha_loz_err = im_compute_error(shapley[loz].nii,shapley[loz].nii_err,$
      shapley[loz].ha,shapley[loz].ha_err,/quotient)

    oiiihb_loz = shapley[loz].oiii/shapley[loz].hb
    oiiihb_loz_err = im_compute_error(shapley[loz].oiii,shapley[loz].oiii_err,$
      shapley[loz].hb,shapley[loz].hb_err,/quotient)
    
    im_symbols, 106, psize=1.5, thick=postthick3, fill=1, color=fsc_color(shapleycolor_loz,150)
    oploterror, alog10(niiha_loz), alog10(oiiihb_loz), niiha_loz_err/niiha_loz/alog(10.0), $
      oiiihb_loz_err/oiiihb_loz/alog(10.0), $
      thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(shapleycolor_loz,150), $
      errcolor=fsc_color(shapleycolor_loz,150)

; Shapley/z~1.4

    hiz = where((shapley.z gt 1.1) and (shapley.oiii gt 0.0))

    niiha_hiz = shapley[hiz].nii/shapley[hiz].ha
    niiha_hiz_err = im_compute_error(shapley[hiz].nii,shapley[hiz].nii_err,$
      shapley[hiz].ha,shapley[hiz].ha_err,/quotient)

    oiiihb_hiz = shapley[hiz].oiii/shapley[hiz].hb
    oiiihb_hiz_err = im_compute_error(shapley[hiz].oiii,shapley[hiz].oiii_err,$
      shapley[hiz].hb,shapley[hiz].hb_err,/quotient)
    
    im_symbols, 108, psize=1.5, thick=postthick3, fill=1, color=fsc_color(shapleycolor_hiz,150)
    oploterror, alog10(niiha_hiz), alog10(oiiihb_hiz), niiha_hiz_err/niiha_hiz/alog(10.0), $
      oiiihb_hiz_err/oiiihb_hiz/alog(10.0), $
      thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(shapleycolor_hiz,150), $
      errcolor=fsc_color(shapleycolor_hiz,150)

; now plot Liu et al. (2008); need to treat the upper limit
; correctly

; Liu/z~1.0    
    
    loz = where((liu.z lt 1.1) and (liu.oiii gt 0.0))

    niiha_loz = liu[loz].nii/liu[loz].ha
    niiha_loz_err = im_compute_error(liu[loz].nii,liu[loz].nii_err,$
      liu[loz].ha,liu[loz].ha_err,/quotient)

    oiiihb_loz = liu[loz].oiii/liu[loz].hb
    oiiihb_loz_err = im_compute_error(liu[loz].oiii,liu[loz].oiii_err,$
      liu[loz].hb,liu[loz].hb_err,/quotient)
    
    im_symbols, 106, psize=1.5, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
    oploterror, alog10(niiha_loz), alog10(oiiihb_loz), niiha_loz_err/niiha_loz/alog(10.0), $
      oiiihb_loz_err/oiiihb_loz/alog(10.0), $
      thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_loz,150), $
      errcolor=fsc_color(liucolor_loz,150)

; Liu/z~1.4

    hiz = where((liu.z gt 1.1) and (liu.oiii gt 0.0))

    niiha_hiz = liu[hiz].nii/liu[hiz].ha
    niiha_hiz_err = im_compute_error(liu[hiz].nii,liu[hiz].nii_err,$
      liu[hiz].ha,liu[hiz].ha_err,/quotient)

    oiiihb_hiz = liu[hiz].oiii/liu[hiz].hb
    oiiihb_hiz_err = im_compute_error(liu[hiz].oiii,liu[hiz].oiii_err,$
      liu[hiz].hb,liu[hiz].hb_err,/quotient)
    
    im_symbols, 108, psize=1.5, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
    oploterror, alog10(niiha_hiz), alog10(oiiihb_hiz), niiha_hiz_err/niiha_hiz/alog(10.0), $
      oiiihb_hiz_err/oiiihb_hiz/alog(10.0), $
      thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_hiz,150), $
      errcolor=fsc_color(liucolor_hiz,150)

; Maier et al. 2006; all the [NII] measurements are 2-sigma upper
; limits (which are negative)

    niiha = -maier.nii/maier.ha

    oiiihb = maier.oiii/maier.hb
    oiiihb_err = im_compute_error(maier.oiii,maier.oiii_err,$
      maier.hb,maier.hb_err,/quotient)
    
    im_symbols, 113, psize=3.5, thick=postthick3, fill=1, color=fsc_color(maiercolor,150)
;   oploterror, alog10(niiha), alog10(oiiihb), oiiihb_err/oiiihb/alog(10.0), $
;     thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(maiercolor,150), $
;     errcolor=fsc_color(maiercolor,150)
    
; legend    
    
    im_legend, textoidl(['z\sim1','z\sim1.5']), /right, /bottom, box=0, $
      charsize=singlecharsize_0, $
      charthick=postthick2, psym=[106,108], fill=[1,1], thick=postthick3, $
      color=fsc_color([liucolor_loz,liucolor_hiz],[101,102]), spacing=2.3, $
      textcolor=fsc_color(textcolor1,100), symsize=1.6

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

    stop    
    
; ---------------------------------------------------------------------------    
; mass-metallicity relation    
; ---------------------------------------------------------------------------    
    
    xsize = 8.5 & ysize = 7.0
    pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, width=7.0, height=5.5, $
      xmargin=[1.2,0.3], ymargin=[0.4,1.1], xpage=xsize, ypage=ysize, $
      position=pos, /normal

    psname = 'mz_06erb'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xsize, ysize=ysize, silent=keyword_set(pdf)

    xtitle = 'log (M_{*}/M'+sunsymbol()+')'
    ytitle = '12 + log (O/H)_{N2}'
    
    xrange = [8.6,11.5]
    yrange = [8.1,8.85]

; plot the sdss galaxies

    sindx = where((sdssohdust.zstrong_niiha gt -900))
    smass = sdsskcorr[sindx].mass ; Chabrier IMF!
    soh = 8.90+0.57*sdssohdust[sindx].zstrong_niiha
    
    im_hogg_scatterplot, smass, soh, outliers=plot_outliers, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor='default', levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,100), $
      contour_color=fsc_color(textcolor1,100), /nogrey, c_linestyle=0, cthick=postthick2

;   plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
;     xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
;     charsize=charsize, charthick=postthick2, xsty=1, ysty=1, xrange=xrange, $
;      position=pos[*,0], yrange=yrange, color=fsc_color(textcolor1,150)

    plot_errors = 0

    if keyword_set(plot_errors) then begin
       im_symbols, 105, psize=2.5, thick=postthick3, fill=1, color=fsc_color(erbcolor,150)
       oploterror, alog10(1D10*erb.mass), erb.log12oh, erb.mass_err/erb.mass/alog(10.0), erb.log12oh_err, $
         thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(erbcolor,150), $
         errcolor=fsc_color(erbcolor,150)
       oploterror, alog10(1D10*erb_ul.mass), -erb_ul.log12oh, erb_ul.mass_err/erb_ul.mass/alog(10.0), $
         erb_ul.log12oh_err, thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(erbcolor,150), $
         errcolor=fsc_color(erbcolor,150)
    endif else begin
       im_symbols, 105, psize=3.0, thick=postthick3, fill=1, color=fsc_color(erbcolor,150)
       oplot, alog10(1D10*erb.mass), erb.log12oh, thick=postthick3, psym=8, color=fsc_color(erbcolor,150)
       plots, alog10(1D10*erb_ul.mass), -erb_ul.log12oh, thick=postthick3, psym=8, color=fsc_color(erbcolor,150)
    endelse
    arrow, alog10(1D10*erb_ul.mass), -erb_ul.log12oh, alog10(1D10*erb_ul.mass), $
      -erb_ul.log12oh-0.05, /data, hsize=-0.6, color=fsc_color(erbcolor,150), $
      hthick=postthick3, thick=postthick3
    
; now plot Liu et al. (2008); use the N2 calibrator and differentiate
; between galaxies at z~1 and z~1.4; also overplot the stack

    plot_individual = 0

    if keyword_set(plot_individual) then begin
       
       loz = where((liu.mass gt -900.0) and (liu.log12oh_n2 gt -900.0) and (liu.z lt 1.1))
       hiz = where((liu.mass gt -900.0) and (liu.log12oh_n2 gt -900.0) and (liu.z gt 1.1))
       
       if keyword_set(plot_errors) then begin
          im_symbols, 106, psize=2.0, thick=postthick3, fill=0, color=fsc_color(liucolor_loz,150)
          oploterror, liu[loz].mass, liu[loz].log12oh_n2, liu[loz].mass_err, liu[loz].log12oh_n2_err, $
            thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_loz,150), $
            errcolor=fsc_color(liucolor_loz,150)
          im_symbols, 108, psize=2.0, thick=postthick3, fill=0, color=fsc_color(liucolor_hiz,150)
          oploterror, liu[hiz].mass, liu[hiz].log12oh_n2, liu[hiz].mass_err, liu[hiz].log12oh_n2_err, $
            thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_hiz,150), $
            errcolor=fsc_color(liucolor_hiz,150)
       endif else begin
          im_symbols, 106, psize=2.5, thick=postthick3, fill=0, color=fsc_color(liucolor_loz,150)
          oplot, liu[loz].mass, liu[loz].log12oh_n2, thick=postthick3, psym=8, $
            color=fsc_color(liucolor_loz,150)
          im_symbols, 108, psize=2.5, thick=postthick3, fill=0, color=fsc_color(liucolor_hiz,150)
          oplot, liu[hiz].mass, liu[hiz].log12oh_n2, thick=postthick3, psym=8, $
            color=fsc_color(liucolor_hiz,150)
       endelse

    endif
       
; stacked spectra    
    
    loz = where((liu_n2stack.mass gt -900.0) and (liu_n2stack.log12oh_n2 gt -900.0) and (liu_n2stack.z lt 1.1))
    hiz = where((liu_n2stack.mass gt -900.0) and (liu_n2stack.log12oh_n2 gt -900.0) and (liu_n2stack.z gt 1.1))
    
    if keyword_set(plot_errors) then begin
       im_symbols, 106, psize=3.0, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
       oploterror, liu_n2stack[loz].mass, liu_n2stack[loz].log12oh_n2, $
         liu_n2stack[loz].mass_err, liu_n2stack[loz].log12oh_n2_err, $
         thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_loz,150), $
         errcolor=fsc_color(liucolor_loz,150)
       im_symbols, 108, psize=3.0, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
       oploterror, liu_n2stack[hiz].mass, liu_n2stack[hiz].log12oh_n2, $
         liu_n2stack[hiz].mass_err, liu_n2stack[hiz].log12oh_n2_err, $
         thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_hiz,150), $
         errcolor=fsc_color(liucolor_hiz,150)
    endif else begin
       im_symbols, 106, psize=3.0, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
       oplot, liu_n2stack[loz].mass, liu_n2stack[loz].log12oh_n2, thick=postthick3, psym=8, $
         color=fsc_color(liucolor_loz,150)
       im_symbols, 108, psize=3.0, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
       oplot, liu_n2stack[hiz].mass, liu_n2stack[hiz].log12oh_n2, thick=postthick3, psym=8, $
         color=fsc_color(liucolor_hiz,150)
    endelse

;   im_symbols, 108, psize=1.5, thick=postthick3, fill=1, color=fsc_color(liucolor_o3n2,150)
;   oploterror, liu_o3n2.mass, liu_o3n2.log12oh_o3n2, liu_o3n2.mass_err, liu_o3n2.log12oh_o3n2_err, $
;     thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(liucolor_o3n2,150), $
;     errcolor=fsc_color(liucolor_o3n2,150)

; Shapley et al. 2005

    if keyword_set(plot_individual) then begin

       loz = where((shapley.mass gt -900.0) and (shapley.log12oh_n2 gt -900.0) and (shapley.z lt 1.1))
       hiz = where((shapley.mass gt -900.0) and (shapley.log12oh_n2 gt -900.0) and (shapley.z gt 1.1))
       
       if keyword_set(plot_errors) then begin
          im_symbols, 106, psize=2.0, thick=postthick3, fill=0, color=fsc_color(shapleycolor_loz,150)
          oploterror, shapley[loz].mass, shapley[loz].log12oh_n2, shapley[loz].mass_err, shapley[loz].log12oh_n2_err, $
            thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(shapleycolor_loz,150), $
            errcolor=fsc_color(shapleycolor_loz,150)
          im_symbols, 108, psize=2.0, thick=postthick3, fill=0, color=fsc_color(shapleycolor_hiz,150)
          oploterror, shapley[hiz].mass, shapley[hiz].log12oh_n2, shapley[hiz].mass_err, shapley[hiz].log12oh_n2_err, $
            thick=postthick3, errthick=postthick3, psym=8, color=fsc_color(shapleycolor_hiz,150), $
            errcolor=fsc_color(shapleycolor_hiz,150)
       endif else begin
          im_symbols, 106, psize=2.5, thick=postthick3, fill=0, color=fsc_color(shapleycolor_loz,150)
          oplot, shapley[loz].mass, shapley[loz].log12oh_n2, $
            thick=postthick3, psym=8, color=fsc_color(shapleycolor_loz,150)
          im_symbols, 108, psize=2.5, thick=postthick3, fill=0, color=fsc_color(shapleycolor_hiz,150)
          oplot, shapley[hiz].mass, shapley[hiz].log12oh_n2, $
            thick=postthick3, psym=8, color=fsc_color(shapleycolor_hiz,150)
       endelse

    endif    
    
; now overplot the median errors

    plot_median_errors = 1
    if keyword_set(plot_median_errors) then begin
       if keyword_set(plot_individual) then begin
          xpos = 10.7 & ypos = 8.3
          im_symbols, 105, psize=3.0, thick=postthick3, fill=1, color=fsc_color(erbcolor,150)
          oploterror, xpos, ypos, median(erb.mass_err/erb.mass/alog(10.0)), median(erb.log12oh_err), $
            thick=postthick3, errthick=postthick4, psym=8, color=fsc_color(erbcolor,150), $
            errcolor=fsc_color(erbcolor,150)
          im_symbols, 108, psize=2.5, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
          oploterror, xpos+0.26, ypos, median([liu[hiz].mass_err,shapley[hiz].mass_err]), $
            median([liu[hiz].log12oh_n2_err,shapley[hiz].log12oh_n2_err]), $
            thick=postthick3, errthick=postthick4, psym=8, color=fsc_color(liucolor_hiz,150), $
            errcolor=fsc_color(liucolor_hiz,150)
          im_symbols, 106, psize=2.5, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
          oploterror, xpos+0.6, ypos, median([liu[loz].mass_err,shapley[loz].mass_err]), $
            median([liu[loz].log12oh_n2_err,shapley[loz].log12oh_n2_err]), $
            thick=postthick3, errthick=postthick4, psym=8, color=fsc_color(liucolor_loz,150), $
            errcolor=fsc_color(liucolor_loz,150)
       endif else begin
          xpos = 10.7 & ypos = 8.25
          im_symbols, 105, psize=3.0, thick=postthick3, fill=1, color=fsc_color(erbcolor,150)
          oploterror, xpos+0.6, ypos, median(erb.mass_err/erb.mass/alog(10.0)), median(erb.log12oh_err), $
            thick=postthick3, errthick=postthick4, psym=8, color=fsc_color(erbcolor,150), $
            errcolor=fsc_color(erbcolor,150)
          im_symbols, 108, psize=2.5, thick=postthick3, fill=1, color=fsc_color(liucolor_hiz,150)
          oploterror, xpos+0.3, ypos, median([liu_n2stack[hiz].mass_err]), $
            median([liu_n2stack[hiz].log12oh_n2_err]), $
            thick=postthick3, errthick=postthick4, psym=8, color=fsc_color(liucolor_hiz,150), $
            errcolor=fsc_color(liucolor_hiz,150)
          im_symbols, 106, psize=2.5, thick=postthick3, fill=1, color=fsc_color(liucolor_loz,150)
          oploterror, xpos, ypos, median([liu_n2stack[loz].mass_err]), $
            median([liu_n2stack[loz].log12oh_n2_err]), $
            thick=postthick3, errthick=postthick4, psym=8, color=fsc_color(liucolor_loz,150), $
            errcolor=fsc_color(liucolor_loz,150)
       endelse
    endif
    
; legend    
    
    im_legend, textoidl(['z\sim1','z\sim1.5','z\sim2']), /left, /top, box=0, $
      charsize=singlecharsize_0, $
      charthick=postthick2, psym=[106,108,105], fill=[1,1,1], thick=postthick3, $
      color=fsc_color([liucolor_loz,liucolor_hiz,erbcolor],[101,102,103]), spacing=2.3, $
      textcolor=fsc_color(textcolor1,100), symsize=1.6

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    
    
; ---------------------------------------------------------------------------    
; (O/H) vs (O/H)
; ---------------------------------------------------------------------------    
    
    xsize = 8.5 & ysize = 7.0
    pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, width=7.0, height=5.5, $
      xmargin=[1.2,0.3], ymargin=[0.4,1.1], xpage=xsize, ypage=ysize, $
      position=pos, /normal

    psname = 'oh_r23_vs_oh_n2'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xsize, ysize=ysize, silent=keyword_set(pdf)

    xtitle = '12 + log (O/H)_{R_{23}}'
    ytitle = '12 + log (O/H)_{N2}'
    
    xrange = [8.3,9.3]
    yrange = xrange
;   xrange = [8.4,9.5]
;   yrange = [8.2,8.9]

; plot the sdss galaxies

    sindx = where((sdssohnodust.zstrong_12oh_kk04 gt -900.0) and $
      (strtrim(sdssohnodust.r23branch_kk04,2) eq 'U') and $
      (sdssohdust.zstrong_niiha gt -900))
    soh_r23 = sdssohnodust[sindx].zstrong_12oh_kk04
    soh_n2 = 8.90+0.57*sdssohdust[sindx].zstrong_niiha
    
    im_hogg_scatterplot, soh_r23, soh_n2, outliers=plot_outliers, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor='default', levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,100), $
      contour_color=fsc_color(textcolor1,100), nogrey=0, c_linestyle=0, cthick=postthick2
    oplot, [7.0,10.0], [7.0,10.0], line=0, thick=postthick1

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    
    
; ------------------------------------------------------------
; integrated spectrum example of an emission-line galaxy
; ------------------------------------------------------------

    galaxy = 'ngc5194'
    psname = 'm51_atlas_example'
    if keyword_set(postscript) or keyword_set(pdf) then begin
       psfile = pspath+psname+'.ps'
       if keyword_set(postscript) then splog, 'Writing '+psfile
    endif else delvarx, psfile

    xmargin = [1.0,0.2] & ymargin = [1.4,1.0]
    xspace = 0.0 & yspace = 0.0
    width = 7.3 & height = 3.7
    xpage = total(height)+total(ymargin)+total(yspace)
    ypage = total(width)+total(xmargin)+total(xspace)

    xspace1 = xspace & yspace1 = yspace
    xmargin1 = xmargin & ymargin1 = ymargin
    width1 = width & height1 = height
    xpage1 = xpage & ypage1 = ypage
    
    arm_plotconfig, landscape=1, nx=1, ny=1, xmargin=xmargin1, $
      ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
      height=height1, coord=pos, xpage=xpage1, ypage=ypage1, $
      psfile=psfile, /writeover, bw=0L ;, /show
    cleanplot, /silent

; version 1
;   linenames = textoidl(['[O II]','H\delta','H\gamma','H\beta','[O III]',$
;     '[N II]','H\alpha','[S II]'])
;   linewaves = [3727,4101,4340,4861,5007,6584,6563,6724]
;   linecolor = ['dodger blue','orange','orange','orange',$
;     'forest green','orange','red','firebrick']
; version 2
    linenames = textoidl(['[O II]','H\delta','H\gamma','H\beta','[O III]','[O III]',$
      '[N II]','[N II]','H\alpha','[S II]','[S II]'])
    linewaves = [3727,4101,4340,4861,4959,5007,6548,6584,6563,6716,6731]
    linecolor = ['dodger blue','orange','orange','orange',$
      'forest green','forest green','red','orange','red','firebrick','firebrick']
; version 3    
;   linenames = textoidl(['[O II]','H\delta','H\gamma','H\beta','[O III]','[O III]',$
;     'He I','[O I]','[N II]','[N II]','H\alpha','[S II]','[S II]'])
;   linewaves = [3727,4101,4340,4861,4959,5007,5876,6300,6548,6584,6563,6716,6731]
    niceprint, linenames, linewaves, linecolor
    
    scale = 1E12
    ytitle = 'Flux (10^{-12} '+flam_units()+')'
    xtitle = 'Rest Wavelength (\AA)'

    specdata = read_atlas_specfit(galaxy,/silent)
    wave = reform(specdata[*,0])
    flux = scale*reform(specdata[*,1])
    yrange = minmax(flux)
    xrange = [min(wave),6830]

    im_lineid_plot, wave, flux, linewaves, linenames, xrange=xrange, yrange=yrange, $
      ps=10, xsty=1, ysty=1, xthick=postthick1, ythick=postthick1, charsize=charsize_8, $
      charthick=postthick2, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      lcharsize=charsize_6, lcharthick=postthick2, extend_thick=postthick2, $
      thick=postthick3, position=pos[*,0], spectrum_color=fsc_color(textcolor1,100), $
      axis_color=fsc_color(textcolor1,100), extend_color=fsc_color('orchid',110), $
      label_color=fsc_color(linecolor,lindgen(n_elements(linecolor))+150)
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

    stop    
    
; ------------------------------------------------------------
; metallicity distribution at z~0.8
; ------------------------------------------------------------

    xpage = 8.5 & ypage = 7.5
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.0, $
      xmargin=[1.6,0.4], ymargin=[0.4,1.1], xpage=xpage, ypage=ypage, $
      position=pos, /normal

    xrange = [8.5,12.3]
    yrange = [8.4,9.35]
;   yrange = [7.71,9.35]

    xtitle = masstitle
    ytitle = ohtitle1

    plot_hoyos = 0
    plot_sava = 0
    plot_lilly = 1
    plot_kob04 = 1

; SDSS - all
    
    sdss_indx = where((sdssohdust.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(sdssohdust.r23branch_ew_alpha_gr_kk04,2) eq 'U'))

    sdss_mg = sdsskcorr[sdss_indx].ugriz_absmag[1]
    sdss_mass = sdsskcorr[sdss_indx].mass + im_convert_imf(/from_chabrier)
    sdss_oh = sdssohdust[sdss_indx].zstrong_ew_alpha_gr_12oh_kk04

; SDSS/AGES - zbin4
    
    sdssages = read_sdssages_mz_sample(/ohdust_zbin1,evolve=0)
;   sdssages = read_sdssages_mz_sample(/ohdust_zbin4,evolve=1)
    sdssages_indx = where((sdssages.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(sdssages.r23branch_ew_alpha_gr_kk04,2) eq 'U'),nsdssages_indx)

    sdssages_mg = sdssages[sdssages_indx].ugriz_absmag[1]
    sdssages_mass = sdssages[sdssages_indx].mass + im_convert_imf(/from_chabrier)
    sdssages_oh = sdssages[sdssages_indx].zstrong_ew_alpha_gr_12oh_kk04
    sdssages_oherr = sdssages[sdssages_indx].zstrong_ew_alpha_gr_12oh_kk04_err
    
; AGES - zbin4    
    
    ages_indx = where((agesohdust.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(agesohdust.r23branch_ew_alpha_gr_kk04,2) eq 'U') and $
      (ageskcorr.z gt 0.60) and (ageskcorr.z lt 0.80) and $
      (agesidust.h_beta_ew[0] gt 3.0),nages_indx)

    mg = ageskcorr[ages_indx].ugriz_absmag[1]
    mass = ageskcorr[ages_indx].mass + im_convert_imf(/from_chabrier)
    oh = agesohdust[ages_indx].zstrong_ew_alpha_gr_12oh_kk04
    oherr = agesohdust[ages_indx].zstrong_ew_alpha_gr_12oh_kk04_err
    weight = ageskcorr[ages_indx].spec_weight

    sfr = agesidust[ages_indx].sfr_oii
;   ml = 10^(mass-(-0.4*(mg-mgsun)))

    ewhb = agesidust[ages_indx].h_beta_ew[0]
    snrhb = agesidust[ages_indx].h_beta[0]/agesidust[ages_indx].h_beta[1]
    ewoii = agesidust[ages_indx].oii_3727_ew[0]
    snroii = agesidust[ages_indx].oii_3727[0]/agesidust[ages_indx].oii_3727[1]

;   ohlo = where(oh lt 8.9)
;   plot, ewhb, snrhb, ps=4, charsize=2, xsty=3, ysty=3
;   plot, snrhb, oh, ps=4, charsize=2, xsty=3, ysty=3
;   plot, snroii, oh, ps=4, charsize=2, xsty=3, ysty=3
;   ages_display_spectrum, agesidust[ages_indx[ohlo]]
    
    normal = where((ageskcorr[ages_indx].x_lum lt -900.0) and $
      (agesidust[ages_indx].nev_3426[0]/agesidust[ages_indx].nev_3426[1] lt 5.0),nnormal)
    xray = where((ageskcorr[ages_indx].x_lum gt -900.0),nxray)
    nev = where((agesidust[ages_indx].nev_3426[0]/agesidust[ages_indx].nev_3426[1] ge 5.0),nnev)
    mgii = where((agesidust[ages_indx].mgii_2800[0]/agesidust[ages_indx].mgii_2800[1] ge 5.0),nmgii)

                                ; now make the plots: (1) just SDSS; (2) all samples with error bars
                                ; and a legend; (3) no error bars; (4) no error bars and the
                                ; high-mass, low-metallicity AGES points on the lower branch; (5)
                                ; same as before, but with Pegase models overlaid
    
                                ; #########################
                                ; (1)    
                                ; #########################
    
    psname = 'mass_vs_12oh_0.8_sdssonly'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)
                                ; SDSS - all    
;    im_hogg_scatterplot, sdssages_mass, sdssages_oh, outliers=0, xsty=1, ysty=1, $
    im_hogg_scatterplot, sdss_mass, sdss_oh, outliers=0, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor='default', levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,101), /nogrey, $
      cthick=postthick1, c_linestyle=0, contour_color=fsc_color(textcolor1,100)
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

                                ; #########################
                                ; (2)    
                                ; #########################
    
    psname = 'mass_vs_12oh_0.8_witherrors'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)
    
                                ; SDSS - all    
;    im_hogg_scatterplot, sdssages_mass, sdssages_oh, outliers=0, xsty=1, ysty=1, $
    im_hogg_scatterplot, sdss_mass, sdss_oh, outliers=0, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor='default', levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,101), /nogrey, $
      cthick=postthick1, c_linestyle=0, contour_color=fsc_color(textcolor1,100)
    
                                ; AGES - 0.6<z<0.8
    plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, xrange=xrange, $
      yrange=yrange, position=pos[*,0], xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), xtitle='', ytitle='', $
      color=fsc_color(textcolor1,100)
                                ;   legend, textoidl('z \sim 0.8'), /left, /top, box=0, $
                                ;     charsize=singlecharsize_0, charthick=postthick2, textcolor=fsc_color(textcolor1,100)
    
    im_symbols, 106, psize=0.8, thick=postthick1, fill=1, color=fsc_color(color_normal,20)
                                ;   djs_oplot, mass, oh, psym=8
    oploterror, mass, oh, oherr, ps=8, thick=postthick1, errthick=postthick1, $
      color=fsc_color(color_normal,20), errcolor=fsc_color(color_normal,20)
    
                                ; literature data
    if keyword_set(plot_hoyos) then begin
       hoyos = read_05hoyos()
       gg = where(hoyos.zstrong_ew_alpha_unity_12oh_kk04 gt -900 and hoyos.lit_log12oh gt -900)
       ohshift = median(hoyos[gg].zstrong_ew_alpha_unity_12oh_kk04-hoyos[gg].lit_log12oh)
       splog, ohshift
       hoyos_good = where((hoyos.m_g gt -900.0) and (hoyos.lit_log12oh gt -900.0) and $
         (hoyos.z gt 0.6) and (hoyos.z lt 1.0),nhoyos_good)
       im_symbols, hoyos05sym, psize=hoyos05psize, fill=1, color=fsc_color(hoyos05color,200)
       oploterror, hoyos[hoyos_good].mass + im_convert_imf(/from_chabrier), $
         hoyos[hoyos_good].lit_log12oh+ohshift, $ ; shift by OHSHIFT
         hoyos[hoyos_good].lit_log12oh_err, ps=8, $ 
         color=fsc_color(hoyos05color,200), errcolor=fsc_color(hoyos05color,200), $
         thick=postthick1, errthick=postthick1
    endif
    
    if keyword_set(plot_sava) then begin
       sava = read_05savaglio()
       sava_good = where((sava.m_g gt -900.0) and (sava.zstrong_12oh_kk04 gt -900.0) and $
         (strtrim(sava.r23branch_kk04,2) eq 'U') and (sava.z gt 0.6) and (sava.z lt 1.0),nsava_good)
       im_symbols, sava05sym, psize=sava05psize, fill=1, color=fsc_color(sava05color,201)
       oploterror, sava[sava_good].mass + im_convert_imf(/from_chabrier), $
         sava[sava_good].zstrong_12oh_kk04, $
         sava[sava_good].zstrong_12oh_kk04_err, ps=8, color=fsc_color(sava05color,201), $
         errcolor=fsc_color(sava05color,201), thick=postthick1, errthick=postthick1
    endif
    
    if keyword_set(plot_lilly) then begin
       lilly = read_03lilly()
       lilly_good = where((lilly.m_g gt -900.0) and (lilly.zstrong_ew_alpha_unity_12oh_kk04 gt -900.0) and $
         (strtrim(lilly.r23branch_ew_alpha_unity_kk04,2) eq 'U') and (lilly.z gt 0.6) and (lilly.z lt 1.0),nlilly_good)
       im_symbols, lilly03sym, psize=lilly03psize, fill=1, color=fsc_color(lilly03color,201)
       oploterror, lilly[lilly_good].mass + im_convert_imf(/from_chabrier), $
         lilly[lilly_good].zstrong_ew_alpha_unity_12oh_kk04, $
         lilly[lilly_good].zstrong_ew_alpha_unity_12oh_kk04_err, ps=8, color=fsc_color(lilly03color,201), $
         errcolor=fsc_color(lilly03color,201), thick=postthick1, errthick=postthick1
    endif
    
    if keyword_set(plot_kob04) then begin
       kob04 = read_04kobulnicky()
       kob04_good = where((kob04.m_g gt -900.0) and (kob04.zstrong_ew_alpha_unity_12oh_kk04 gt -900.0) and $
         (strtrim(kob04.r23branch_ew_alpha_unity_kk04,2) eq 'U') and (kob04.z gt 0.6) and (kob04.z lt 1.0),nkob04_good)
       im_symbols, kob04sym, psize=kob04psize, fill=1, color=fsc_color(kob04color,201)
       oploterror, kob04[kob04_good].mass + im_convert_imf(/from_chabrier), $
         kob04[kob04_good].zstrong_ew_alpha_unity_12oh_kk04, $
         kob04[kob04_good].zstrong_ew_alpha_unity_12oh_kk04_err, ps=8, color=fsc_color(kob04color,201), $
         errcolor=fsc_color(kob04color,201), thick=postthick1, errthick=postthick1
    endif
    
                                ; literature legend
                                ;   legend, textoidl('z \sim 0.8'), /left, /top, box=0, $
                                ;     charsize=singlecharsize_0, charthick=postthick2, textcolor=fsc_color(textcolor1,100)
    
    legendlabel = ['Moustakas et al. 2008','Kobulnicky & Kewley 2004','Lilly et al. 2003']
    legendpsym = [106,kob04sym,lilly03sym]
    legendfill = [1,1,1]
    legendcolor = fsc_color([color_normal,kob04color,lilly03color],[55,56,57])
    
;    legendlabel = ['Moustakas et al. 2008','Kobulnicky et al. 2004','Lilly et al. 2003','Hoyos et al. 2005']
;    legendpsym = [106,kob04sym,lilly03sym,hoyos05sym]
;    legendfill = [1,1,1,1]
;    legendcolor = fsc_color([color_normal,kob04color,lilly03color,hoyos05color],[55,56,57,58])

    postthick_legend = postthick1
    im_legend, textoidl(legendlabel), /left, /top, box=0, charsize=charsize_2, $
      charthick=postthick2, psym=legendpsym, fill=legendfill, symsize=1.2, $
      spacing=1.7, thick=postthick_legend, symthick=postthick_legend, color=legendcolor, $
      /normal, textcolor=fsc_color(textcolor1,100)
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif
    
                                ; #########################
                                ; (3)    
                                ; #########################
    
    psname = 'mass_vs_12oh_0.8_noerrors'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)
    
                                ; SDSS - all    
;    im_hogg_scatterplot, sdssages_mass, sdssages_oh, outliers=0, xsty=1, ysty=1, $
    im_hogg_scatterplot, sdss_mass, sdss_oh, outliers=0, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor='default', levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,101), /nogrey, $
      cthick=postthick1, c_linestyle=0, contour_color=fsc_color(textcolor1,100)
    
                                ; AGES - 0.6<z<0.8
    plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, xrange=xrange, $
      yrange=yrange, position=pos[*,0], xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), xtitle='', ytitle='', $
      color=fsc_color(textcolor1,100)
;    legend, textoidl('z \sim 0.8'), /left, /top, box=0, $
;      charsize=singlecharsize_0, charthick=postthick2, color=fsc_color(textcolor1,100)
    
    pfactor = 1.2
    
    im_symbols, 106, psize=0.8*pfactor, thick=postthick1, fill=1, color=fsc_color(color_normal,20)
    oplot, mass, oh, psym=8, color=fsc_color(color_normal,20)
    
                                ; literature data
    if keyword_set(plot_hoyos) then begin
       hoyos = read_05hoyos()
       gg = where(hoyos.zstrong_ew_alpha_unity_12oh_kk04 gt -900 and hoyos.lit_log12oh gt -900)
       ohshift = median(hoyos[gg].zstrong_ew_alpha_unity_12oh_kk04-hoyos[gg].lit_log12oh)
       splog, ohshift
       hoyos_good = where((hoyos.m_g gt -900.0) and (hoyos.lit_log12oh gt -900.0) and $
         (hoyos.z gt 0.6) and (hoyos.z lt 1.0),nhoyos_good)
       im_symbols, hoyos05sym, psize=hoyos05psize*pfactor, fill=1, color=fsc_color(hoyos05color,200)
       oplot, hoyos[hoyos_good].mass + im_convert_imf(/from_chabrier), $
         hoyos[hoyos_good].lit_log12oh+ohshift, $ ; shift by OHSHIFT
         ps=8, color=fsc_color(hoyos05color,200), thick=postthick1
    endif
    
    if keyword_set(plot_sava) then begin
       sava = read_05savaglio()
       sava_good = where((sava.m_g gt -900.0) and (sava.zstrong_12oh_kk04 gt -900.0) and $
         (strtrim(sava.r23branch_kk04,2) eq 'U') and (sava.z gt 0.6) and (sava.z lt 1.0),nsava_good)
       im_symbols, sava05sym, psize=sava05psize*pfactor, fill=1, color=fsc_color(sava05color,201)
       oplot, sava[sava_good].mass + im_convert_imf(/from_chabrier), $
         sava[sava_good].zstrong_12oh_kk04, $
         ps=8, color=fsc_color(sava05color,201), thick=postthick1
    endif
    
    if keyword_set(plot_lilly) then begin
       lilly = read_03lilly()
       lilly_good = where((lilly.m_g gt -900.0) and (lilly.zstrong_ew_alpha_unity_12oh_kk04 gt -900.0) and $
         (strtrim(lilly.r23branch_ew_alpha_unity_kk04,2) eq 'U') and (lilly.z gt 0.6) and (lilly.z lt 1.0),nlilly_good)
       im_symbols, lilly03sym, psize=lilly03psize*pfactor, fill=1, color=fsc_color(lilly03color,201)
       oplot, lilly[lilly_good].mass + im_convert_imf(/from_chabrier), $
         lilly[lilly_good].zstrong_ew_alpha_unity_12oh_kk04, $
         ps=8, color=fsc_color(lilly03color,201), thick=postthick1
    endif
    
    if keyword_set(plot_kob04) then begin
       kob04 = read_04kobulnicky()
       kob04_good = where((kob04.m_g gt -900.0) and (kob04.zstrong_ew_alpha_unity_12oh_kk04 gt -900.0) and $
         (strtrim(kob04.r23branch_ew_alpha_unity_kk04,2) eq 'U') and (kob04.z gt 0.6) and (kob04.z lt 1.0),nkob04_good)
       im_symbols, kob04sym, psize=kob04psize*pfactor, fill=1, color=fsc_color(kob04color,201)
       oplot, kob04[kob04_good].mass + im_convert_imf(/from_chabrier), $
         kob04[kob04_good].zstrong_ew_alpha_unity_12oh_kk04, $
         ps=8, color=fsc_color(kob04color,201), thick=postthick1
    endif
    
                                ; literature legend
                                ;   legend, textoidl('z \sim 0.8'), /left, /top, box=0, $
                                ;     charsize=singlecharsize_0, charthick=postthick2, textcolor=fsc_color(textcolor1,100)
    
    legendlabel = ['Moustakas et al. 2008','Kobulnicky & Kewley 2004','Lilly et al. 2003']
    legendpsym = [106,kob04sym,lilly03sym]
    legendfill = [1,1,1]
    legendcolor = fsc_color([color_normal,kob04color,lilly03color],[55,56,57])
    
;    legendlabel = ['Moustakas et al. 2008','Kobulnicky et al. 2004','Lilly et al. 2003','Hoyos et al. 2005']
;    legendpsym = [106,kob04sym,lilly03sym,hoyos05sym]
;    legendfill = [1,1,1,1]
;    legendcolor = fsc_color([color_normal,kob04color,lilly03color,hoyos05color],[55,56,57,58])

    postthick_legend = postthick1
    im_legend, textoidl(legendlabel), /left, /top, box=0, charsize=charsize_2, $
      charthick=postthick2, psym=legendpsym, fill=legendfill, symsize=1.2, $
      spacing=1.7, thick=postthick_legend, symthick=postthick_legend, color=legendcolor, $
      /normal, textcolor=fsc_color(textcolor1,100)

; mean error bar

    oploterror, 11.8, 8.55, 0.3, mean(oherr), thick=postthick1, $
      errthick=postthick1, color=fsc_color(textcolor1,150), $
      errcolor=fsc_color(textcolor1,150)
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    
    
                                ; #########################
                                ; (4)    
                                ; #########################
    
    psname = 'mass_vs_12oh_0.8_noerrors_ohlo'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)
    
                                ; SDSS - all    
;    im_hogg_scatterplot, sdssages_mass, sdssages_oh, outliers=0, xsty=1, ysty=1, $
    im_hogg_scatterplot, sdss_mass, sdss_oh, outliers=0, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor='default', levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,101), /nogrey, $
      cthick=postthick1, c_linestyle=0, contour_color=fsc_color(textcolor1,100)
    
                                ; AGES - 0.6<z<0.8
    plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, xrange=xrange, $
      yrange=yrange, position=pos[*,0], xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), xtitle='', ytitle='', $
      color=fsc_color(textcolor1,100)
    legend, textoidl('z \sim 0.8'), /left, /top, box=0, $
      charsize=singlecharsize_0, charthick=postthick2, color=fsc_color(textcolor1,100)
    
    pfactor = 1.2
    
    im_symbols, 106, psize=0.8*pfactor, thick=postthick1, fill=1, color=fsc_color(color_normal,20)
    oplot, mass, oh, psym=8, color=fsc_color(color_normal,20)
    
    lo = where(oh lt 8.9,nlo,comp=hi)
    ages_lo = agesidust[ages_indx[lo]]
    ages_lo.r23branch = 'L'
    mz_log12oh, ages_lo, final_ohdust=ages_ohlo, nmonte=500
    
    im_symbols, 106, psize=0.8*pfactor, thick=postthick1, fill=0, color=fsc_color(color_normal,20)
    oplot, mass[lo], ages_ohlo.zstrong_ew_alpha_gr_12oh_kk04, psym=8, color=fsc_color(color_normal,20)
    for ii = 0L, nlo-1L do oplot, mass[lo[ii]]*[1,1], [oh[lo[ii]],ages_ohlo[ii].zstrong_ew_alpha_gr_12oh_kk04], $
      line=0, thick=1.0, color=fsc_color(textcolor1,100)
    
                                ; literature data
    if keyword_set(plot_hoyos) then begin
       hoyos = read_05hoyos()
       gg = where(hoyos.zstrong_ew_alpha_unity_12oh_kk04 gt -900 and hoyos.lit_log12oh gt -900)
       ohshift = median(hoyos[gg].zstrong_ew_alpha_unity_12oh_kk04-hoyos[gg].lit_log12oh)
       splog, ohshift
       hoyos_good = where((hoyos.m_g gt -900.0) and (hoyos.lit_log12oh gt -900.0) and $
         (hoyos.z gt 0.6) and (hoyos.z lt 1.0),nhoyos_good)
       im_symbols, hoyos05sym, psize=hoyos05psize*pfactor, fill=1, color=fsc_color(hoyos05color,200)
       oplot, hoyos[hoyos_good].mass + im_convert_imf(/from_chabrier), $
         hoyos[hoyos_good].lit_log12oh+ohshift, $ ; shift by OHSHIFT
         ps=8, color=fsc_color(hoyos05color,200), thick=postthick1
    endif
    
    if keyword_set(plot_sava) then begin
       sava = read_05savaglio()
       sava_good = where((sava.m_g gt -900.0) and (sava.zstrong_12oh_kk04 gt -900.0) and $
         (strtrim(sava.r23branch_kk04,2) eq 'U') and (sava.z gt 0.6) and (sava.z lt 1.0),nsava_good)
       im_symbols, sava05sym, psize=sava05psize*pfactor, fill=1, color=fsc_color(sava05color,201)
       oplot, sava[sava_good].mass + im_convert_imf(/from_chabrier), $
         sava[sava_good].zstrong_12oh_kk04, $
         ps=8, color=fsc_color(sava05color,201), thick=postthick1
    endif
    
    if keyword_set(plot_lilly) then begin
       lilly = read_03lilly()
       lilly_good = where((lilly.m_g gt -900.0) and (lilly.zstrong_ew_alpha_unity_12oh_kk04 gt -900.0) and $
         (strtrim(lilly.r23branch_ew_alpha_unity_kk04,2) eq 'U') and (lilly.z gt 0.6) and (lilly.z lt 1.0),nlilly_good)
       im_symbols, lilly03sym, psize=lilly03psize*pfactor, fill=1, color=fsc_color(lilly03color,201)
       oplot, lilly[lilly_good].mass + im_convert_imf(/from_chabrier), $
         lilly[lilly_good].zstrong_ew_alpha_unity_12oh_kk04, $
         ps=8, color=fsc_color(lilly03color,201), thick=postthick1
    endif
    
    if keyword_set(plot_kob04) then begin
       kob04 = read_04kobulnicky()
       kob04_good = where((kob04.m_g gt -900.0) and (kob04.zstrong_ew_alpha_unity_12oh_kk04 gt -900.0) and $
         (strtrim(kob04.r23branch_ew_alpha_unity_kk04,2) eq 'U') and (kob04.z gt 0.6) and (kob04.z lt 1.0),nkob04_good)
       im_symbols, kob04sym, psize=kob04psize*pfactor, fill=1, color=fsc_color(kob04color,201)
       oplot, kob04[kob04_good].mass + im_convert_imf(/from_chabrier), $
         kob04[kob04_good].zstrong_ew_alpha_unity_12oh_kk04, $
         ps=8, color=fsc_color(kob04color,201), thick=postthick1
    endif
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; #########################
; (5)    
; #########################
    
    psname = 'mass_vs_12oh_0.8_pegase'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

; SDSS - all    
;   im_hogg_scatterplot, sdssages_mass, sdssages_oh, outliers=0, xsty=1, ysty=1, $
    im_hogg_scatterplot, sdss_mass, sdss_oh, outliers=0, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor='default', levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,101), /nogrey, $
      cthick=postthick1, c_linestyle=0, contour_color=fsc_color(textcolor1,100)
    
; AGES - 0.6<z<0.8
    plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, xrange=xrange, $
      yrange=yrange, position=pos[*,0], xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), xtitle='', ytitle='', $
      color=fsc_color(textcolor1,100)
;   legend, textoidl('z \sim 0.8'), /left, /top, box=0, $
;     charsize=singlecharsize_0, charthick=postthick2, textcolor=fsc_color(textcolor1,100)

    pfactor = 1.2
    
    im_symbols, 106, psize=0.8*pfactor, thick=postthick1, fill=1, color=fsc_color(color_normal,20)
    oplot, mass, oh, psym=8, color=fsc_color(color_normal,20)

;   im_symbols, 106, psize=0.8*pfactor, thick=postthick1, fill=0, color=fsc_color(color_normal,20)
;   oplot, mass[lo], ages_ohlo.zstrong_ew_alpha_gr_12oh_kk04, psym=8, color=fsc_color(color_normal,20)
;   for ii = 0L, nlo-1L do oplot, mass[lo[ii]]*[1,1], [oh[lo[ii]],ages_ohlo[ii].zstrong_ew_alpha_gr_12oh_kk04], $
;     line=0, thick=1.0, color=fsc_color(textcolor1,100)
    
; literature data
    if keyword_set(plot_hoyos) then begin
       hoyos = read_05hoyos()
       gg = where(hoyos.zstrong_ew_alpha_unity_12oh_kk04 gt -900 and hoyos.lit_log12oh gt -900)
       ohshift = median(hoyos[gg].zstrong_ew_alpha_unity_12oh_kk04-hoyos[gg].lit_log12oh)
       splog, ohshift
       hoyos_good = where((hoyos.m_g gt -900.0) and (hoyos.lit_log12oh gt -900.0) and $
         (hoyos.z gt 0.6) and (hoyos.z lt 1.0),nhoyos_good)
       im_symbols, hoyos05sym, psize=hoyos05psize*pfactor, fill=1, color=fsc_color(hoyos05color,200)
       oplot, hoyos[hoyos_good].mass + im_convert_imf(/from_chabrier), $
         hoyos[hoyos_good].lit_log12oh+ohshift, $ ; shift by OHSHIFT
         ps=8, color=fsc_color(hoyos05color,200), thick=postthick1
    endif
    
    if keyword_set(plot_sava) then begin
       sava = read_05savaglio()
       sava_good = where((sava.m_g gt -900.0) and (sava.zstrong_12oh_kk04 gt -900.0) and $
         (strtrim(sava.r23branch_kk04,2) eq 'U') and (sava.z gt 0.6) and (sava.z lt 1.0),nsava_good)
       im_symbols, sava05sym, psize=sava05psize*pfactor, fill=1, color=fsc_color(sava05color,201)
       oplot, sava[sava_good].mass + im_convert_imf(/from_chabrier), $
         sava[sava_good].zstrong_12oh_kk04, $
         ps=8, color=fsc_color(sava05color,201), thick=postthick1
    endif

    if keyword_set(plot_lilly) then begin
       lilly = read_03lilly()
       lilly_good = where((lilly.m_g gt -900.0) and (lilly.zstrong_ew_alpha_unity_12oh_kk04 gt -900.0) and $
         (strtrim(lilly.r23branch_ew_alpha_unity_kk04,2) eq 'U') and (lilly.z gt 0.6) and (lilly.z lt 1.0),nlilly_good)
       im_symbols, lilly03sym, psize=lilly03psize*pfactor, fill=1, color=fsc_color(lilly03color,201)
       oplot, lilly[lilly_good].mass + im_convert_imf(/from_chabrier), $
         lilly[lilly_good].zstrong_ew_alpha_unity_12oh_kk04, $
         ps=8, color=fsc_color(lilly03color,201), thick=postthick1
    endif
    
    if keyword_set(plot_kob04) then begin
       kob04 = read_04kobulnicky()
       kob04_good = where((kob04.m_g gt -900.0) and (kob04.zstrong_ew_alpha_unity_12oh_kk04 gt -900.0) and $
         (strtrim(kob04.r23branch_ew_alpha_unity_kk04,2) eq 'U') and (kob04.z gt 0.6) and (kob04.z lt 1.0),nkob04_good)
       im_symbols, kob04sym, psize=kob04psize*pfactor, fill=1, color=fsc_color(kob04color,201)
       oplot, kob04[kob04_good].mass + im_convert_imf(/from_chabrier), $
         kob04[kob04_good].zstrong_ew_alpha_unity_12oh_kk04, $
         ps=8, color=fsc_color(kob04color,201), thick=postthick1
    endif
    
; overlay some Pegase models

    labelz = [0.1,0.7,1.0,1.2]
    xlabeloffset = labelz*0.0-0.01
    ylabeloffset = [0.03,0.0,0.0,0.0]
;   labelage = labelz-pegzform
;   labelage = [0.5,1.0,3.0,7.0] ; Gyr
    
;   these = where((peg1.age/1E3 gt min(labelage)) and (peg1.age/1E3 lt max(labelage)))
;   peg1mg = peg1[these].ugriz[1]-2.5*alog10(peg_mgalaxy)
;   peg1oh = peg1[these].log12oh & peg1age = peg1[these].age/1E3

    peg1_mass = alog10(peg1.mstar*peg1_mgalaxy) & peg1_oh = peg1.log12oh
    peg1_mass = peg1_mass/peg1_mass[0]*alog10(peg1_mgalaxy) ; re-normalize (from recycling)
    djs_oplot, peg1_mass, peg1_oh, line=5, thick=postthick3, $
      color=fsc_color('forest green',190)
    xyouts, interpol(peg1_mass,peg1_zaxis,labelz)+xlabeloffset, $
      interpol(peg1_oh,peg1_zaxis,labelz)+ylabeloffset, $
      strtrim(string(labelz,format='(F12.1)'),2), charsize=1.5, align=1.0, $
      charthick=postthick2, color=fsc_color(textcolor1,100)

    peg3_mass = alog10(peg3.mstar*peg3_mgalaxy) & peg3_oh = peg3.log12oh
    peg3_mass = peg3_mass/peg3_mass[0]*alog10(peg3_mgalaxy) ; re-normalize (from recycling)
    djs_oplot, peg3_mass, peg3_oh, line=3, thick=postthick3, $
      color=fsc_color('orange',191)
    xyouts, interpol(peg3_mass,peg3_zaxis,labelz)+xlabeloffset, $
      interpol(peg3_oh,peg3_zaxis,labelz)+ylabeloffset, $
      strtrim(string(labelz,format='(F12.1)'),2), charsize=1.5, align=1.0, $
      charthick=postthick2, color=fsc_color(textcolor1,100)

    pegconst_mass = alog10(pegconst.mstar*pegconst_mgalaxy) & pegconst_oh = pegconst.log12oh
    pegconst_mass = pegconst_mass/pegconst_mass[0]*alog10(pegconst_mgalaxy) ; re-normalize (from recycling)
    djs_oplot, pegconst_mass, pegconst_oh, line=0, thick=postthick3, $
      color=fsc_color('dodger blue',192)
    xyouts, interpol(pegconst_mass,pegconst_zaxis,labelz)+xlabeloffset, $
      interpol(pegconst_oh,pegconst_zaxis,labelz)+ylabeloffset, $
      strtrim(string(labelz,format='(F12.1)'),2), charsize=1.5, align=1.0, $
      charthick=postthick2, color=fsc_color(textcolor1,100)

    legend, textoidl(['\tau=1; z_{f}=1.0','\tau=3; z_{f}=1.2','\tau=\infty; z_{f}=2.0']), $
      /left, /top, box=0, charsize=1.2, $
      charthick=postthick2, spacing=2.2, thick=postthick3, $
      color=fsc_color(['forest green','orange','dodger blue'],[50,51,52]), $
      textcolor=fsc_color(['forest green','orange','dodger blue'],[50,51,52]), $
      line=[5,3,0], fill=[1,1,1], position=[xrange[1],yrange[0]+0.03]

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

    stop    

; ------------------------------------------------------------
; wavelength versus redshift
; ------------------------------------------------------------

    psname = 'wavelength_vs_redshift'
    xpage = 8.5 & ypage = 8.5
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.1, height=7.1, $
      xmargin=[1.1,0.3], ymargin=[0.3,1.1], xpage=xpage, ypage=ypage, $
      position=pos, /normal

    zmax = 5.2
    zmin = 0.0
    dz = 0.1
    zgrid = findgen((zmax-zmin)/dz+1)*dz+zmin

    lamgrid = findgen((26E4-3500.0)/1.0+1)*1.0+3500.0
    lamgrid = lamgrid / 1E4

;   lamrest = [3727.0,4861.0,6563.0] / 1E4 ; OII, H-beta, OIII, Ha
;   linename = ['[O II]','Hb','Ha']
;   linestyle = [0,3,5]

    lamrest = [3727.0,mean([4861.0,5007]),mean([6563.0,6584])] / 1E4 ; OII, H-beta, OIII, Ha
    linename = ['[O II]','Hb','Ha']
    linestyle = [0,3,5]

;   lamrest = [3727.0,4340.0,4861.0,6563.0] / 1E4 ; OII, H-beta, OIII, Ha
;   linename = ['[O II]','Hg','Hb','Ha']
;   linestyle = [0,1,3,5]

    xtitle = textoidl('log \lambda_{obs} (\mu'+'m)')
    ytitle = 'Redshift'

    xrange = [1600.0,26E3] / 1E4
    yrange = minmax(zgrid)
    
    plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      charsize=2.0, charthick=postthick2, xthick=postthick1, ythick=postthick1, $
      xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, position=pos[*,0], $
      color=fsc_color(textcolor1,100)
    
; fill the optical wavelength range
    
    polyfill, [3500.0,3500.0,9500,9500] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('orchid',200), /fill
    polyfill, [3500.0,3500.0,9500,9500] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('orchid',200), /fill
    xyouts, 6500.0/1E4, 4.2, 'Optical', charsize=2.0, charthick=postthick2, $
      align=0.5, /data, color=fsc_color(textcolor1,100)

    polyfill, [11E3,11E3,14E3,14E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('dodger blue',201), /fill
    polyfill, [11E3,11E3,14E3,14E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('dodger blue',201), /fill
    xyouts, 12.5E3/1E4, 4.2, 'J', charsize=2.0, charthick=postthick2, $
      align=0.5, /data, color=fsc_color(textcolor1,100)

    polyfill, [15E3,15E3,18E3,18E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('forest green',202), /fill
    polyfill, [15E3,15E3,18E3,18E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('forest green',202), /fill
    xyouts, 16.5E3/1E4, 4.2, 'H', charsize=2.0, charthick=postthick2, $
      align=0.5, /data, color=fsc_color(textcolor1,100)

    polyfill, [20E3,20E3,23E3,23E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('orange',203), /fill
    polyfill, [20E3,20E3,23E3,23E3] / 1E4, [yrange[0],yrange[1],$
      yrange[1],yrange[0]], color=fsc_color('orange',203), /fill

    xyouts, 21.5E3/1E4, 4.2, 'K', charsize=2.0, charthick=postthick2, $
      align=0.5, /data, color=fsc_color(textcolor1,100)

    for i = 0L, n_elements(lamrest)-1L do djs_oplot, lamrest[i]*(1+zgrid), $
      zgrid, line=linestyle[i], thick=postthick3, color=fsc_color(textcolor1,100)

    xyouts, 1.5, 2.1, textoidl('[O II]'), charsize=1.6, charthick=postthick2, $
      align=0.5, /data, color=fsc_color('firebrick',204), orientation=35
    xyouts, 1.02, 1.4, textoidl('H\beta+[O III]'), charsize=1.6, charthick=postthick2, $
      align=0.5, /data, color=fsc_color(textcolor1,100), orientation=35
    xyouts, 1.02, 0.7, textoidl('H\alpha+[N II]'), charsize=1.6, charthick=postthick2, $
      align=0.5, /data, color=fsc_color(textcolor1,100), orientation=45

; overplot the tick marks    
    
    plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, /noerase, $
      charsize=2.0, charthick=postthick2, xthick=postthick1, ythick=postthick1, $
      xsty=1, ysty=1, xtitle='', ytitle='', position=pos[*,0], color=fsc_color(textcolor1,100)
    
; legend

;   legend, textoidl(['[O II]','H\beta','H\alpha']), /right, $
;     /bottom, charsize=1.8, charthick=postthick, line=linestyle, thick=postthick, $
;     clear=keyword_set(postscript), spacing=1.5, box=0

; when are each of the emission lines within the atmospheric windows?

    for iline = 0L, n_elements(linename)-1L do begin
       print, linename[iline]+':'
       print, ' Optical: <'+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),0.95),format='(F12.1)'),2)
       print, ' J-band : '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.1),format='(F12.1)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.4),format='(F12.1)'),2)
       print, ' H-band : '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.5),format='(F12.1)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.8),format='(F12.1)'),2)
       print, ' K-band : '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),2.0),format='(F12.1)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),2.3),format='(F12.1)'),2)
       print
    endfor
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif
    stop    
    
; ------------------------------------------------------------    
; example AGES spectra
; ------------------------------------------------------------

; plotting preliminaries    
    
    xmargin = [1.3,0.4] & ymargin = [0.5,1.1]
    xspace = 0.1 & yspace = 0.0
    width = 7.0 & height = 4.5
    xpage = total(height)+total(ymargin)+total(yspace)
    ypage = total(width)+total(xmargin)+total(xspace)

    xspace1 = xspace & yspace1 = yspace
    xmargin1 = xmargin & ymargin1 = ymargin
    width1 = width & height1 = height
    xpage1 = xpage & ypage1 = ypage

    nsmooth = 2L
    bestfitcolor = 'red'        ; firebrick
    especcolor = 'forest green'

    lineminwave = 4800
    linemaxwave = 5050

    plotobs = 0
    if keyword_set(plotobs) then begin
       thisminwave = 4400
       thismaxwave = 8300
    endif else begin
       thisminwave = 2915
       thismaxwave = 5495
    endelse

; residual emission-line spectrum with the emission-line fits 
    
    psname = 'ages_example_spectrum_emission'
    if keyword_set(postscript) or keyword_set(pdf) then begin
       psfile = pspath+psname+'.ps'
       if keyword_set(postscript) then splog, 'Writing '+psfile
    endif else delvarx, psfile

    arm_plotconfig, landscape=1, nx=1, ny=1, xmargin=xmargin1, $
      ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
      height=height1, coord=pos, xpage=xpage1, ypage=ypage1, $
      psfile=psfile, /writeover, bw=0L ;, /show
    cleanplot, /silent

    ages_display_spectrum, ages[thisages], /nowindow, plotobs=0, position=pos, $
      labeltype=0L, minwave=lineminwave, maxwave=linemaxwave, color=fsc_color(textcolor1,150), $
      postthick1=postthick1, postthick2=postthick2, speccolor=fsc_color(textcolor1,150), $
      especcolor=fsc_color(especcolor,152), nsmooth=nsmooth, plotthick1=postthick3, $
      plotthick2=postthick1, /noupperaxis, plottype=3L, yrangetype=4L
    xyouts, 4861.0, !y.crange[1]*1.05, textoidl('H\beta'), align=0.5, $
      charsize=1.7, charthick=postthick2, /data, color=fsc_color(textcolor1,150)
    xyouts, 4959.0, !y.crange[1]*1.05, textoidl('[O III]'), align=0.5, $
      charsize=1.7, charthick=postthick2, /data, color=fsc_color(textcolor1,150)
    xyouts, 5007.0, !y.crange[1]*1.05, textoidl('[O III]'), align=0.5, $
      charsize=1.7, charthick=postthick2, /data, color=fsc_color(textcolor1,150)
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; first spectrum - just the data    
    
    psname = 'ages_example_spectrum'
    if keyword_set(postscript) or keyword_set(pdf) then begin
       psfile = pspath+psname+'.ps'
       if keyword_set(postscript) then splog, 'Writing '+psfile
    endif else delvarx, psfile

    arm_plotconfig, landscape=1, nx=1, ny=1, xmargin=xmargin1, $
      ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
      height=height1, coord=pos, xpage=xpage1, ypage=ypage1, $
      psfile=psfile, /writeover, bw=0L ;, /show
    cleanplot, /silent

    ages_display_spectrum, ages[thisages], /nowindow, plotobs=0, position=pos, $
      labeltype=0L, minwave=thisminwave, maxwave=thismaxwave, color=fsc_color(textcolor1,150), $
      postthick1=postthick1, postthick2=postthick2, speccolor=fsc_color(textcolor1,150), $
      nsmooth=nsmooth, plotthick1=postthick3, plotthick2=postthick1, $
      /noupperaxis, plottype=2L
    
;   xyouts, 3727.0, !y.crange[1]*1.05, textoidl('[O II]'), align=0.5, charsize=1.7, charthick=postthick2, /data
;   xyouts, 4861.0, !y.crange[1]*1.05, textoidl('H\beta'), align=0.5, charsize=1.7, charthick=postthick2, /data
;   xyouts, 4959.0, !y.crange[1]*1.05, textoidl('[O III]'), align=0.5, charsize=1.7, charthick=postthick2, /data
;   xyouts, 5007.0, !y.crange[1]*1.05, textoidl('[O III]'), align=0.5, charsize=1.7, charthick=postthick2, /data

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; second spectrum - the data with the best BC03 fit    
    
    psname = 'ages_example_spectrum_with_bestfit'
    if keyword_set(postscript) or keyword_set(pdf) then begin
       psfile = pspath+psname+'.ps'
       if keyword_set(postscript) then splog, 'Writing '+psfile
    endif else delvarx, psfile

    arm_plotconfig, landscape=1, nx=1, ny=1, xmargin=xmargin1, $
      ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
      height=height1, coord=pos, xpage=xpage1, ypage=ypage1, $
      psfile=psfile, /writeover, bw=0L ;, /show
    cleanplot, /silent

    ages_display_spectrum, ages[thisages], /nowindow, plotobs=0, position=pos, $
      labeltype=0L, minwave=thisminwave, maxwave=thismaxwave, color=fsc_color(textcolor1,150), $
      postthick1=postthick1, postthick2=postthick2, speccolor=fsc_color(textcolor1,150), $
      bestfitcolor=fsc_color(bestfitcolor,151), nsmooth=nsmooth, plotthick1=postthick3, $
      plotthick2=postthick1, /noupperaxis, plottype=1L
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

    stop    
    
; ------------------------------------------------------------    
; demonstration of how we measure abundances - upper branch only
; ------------------------------------------------------------

    psname = 'r23_vs_12oh_upper'
    if keyword_set(postscript) or keyword_set(pdf) then begin
       psfile = pspath+psname+'.ps'
       if keyword_set(postscript) then splog, 'Writing '+psfile
    endif else delvarx, psfile

    xmargin = [1.3,0.4] & ymargin = [0.4,1.1]
    xspace = 0.0 & yspace = 0.0
    width = 6.5 & height = 5.5
    xpage = total(height)+total(ymargin)+total(yspace)
    ypage = total(width)+total(xmargin)+total(xspace)

    xspace1 = xspace & yspace1 = yspace
    xmargin1 = xmargin & ymargin1 = ymargin
    width1 = width & height1 = height
    xpage1 = xpage & ypage1 = ypage
    
    arm_plotconfig, landscape=1, nx=1, ny=1, xmargin=xmargin1, $
      ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
      height=height1, coord=pos, xpage=xpage1, ypage=ypage1, $
      psfile=psfile, /writeover, bw=0L ;, /show
    cleanplot, /silent

    r23obs = 4.0 & r23obs_err = 1.5
    o32obs = 0.60 & o32obs_err = 0.03

    logr23obs = alog10(r23obs) & logr23obs_err = r23obs_err/r23obs/alog(10.0)
    logo32obs = alog10(o32obs) & logo32obs_err = o32obs_err/o32obs/alog(10.0)
    
;   logr23obs = 1.0 & logr23obs_err = 0.1
;   logo32obs = -0.5 & logo32obs_err = 0.05
    
    logq = alog10([5D6,4E7,1.5D8]) ; alog10(4E7)
    logu = logq-alog10(light)
;   logr23 = findgen((1.1-(-0.0))/0.001+1)*0.001-0.0
    logr23 = findgen((1.1-(-0.5))/0.001+1)*0.001-0.5
    linestyle = [0,3,5]
    linecolor = ['orange','forest green','dodger blue']
    
; plot the R23-O/H relation and the KK04 theoretical calibration

    ohrange = [8.25,9.27]
    xrange = [-0.1,1.05]
    yrange = ohrange

    xtitle = 'log [([O II] + [O III])/H\beta]'
;   xtitle = 'log (R_{23})'
    ytitle = '12 + log (O/H)'
    
    plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), $
      ytitle=textoidl(ytitle), charsize=charsize1, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, $
      position=pos[*,0], xsty=1, ysty=1, color=fsc_color(textcolor1,150)
    
    for iq = 0L, n_elements(logq)-1L do begin
       
       logoh_upper = 9.72D - 0.777*logr23 - 0.951*logr23^2 - 0.072*logr23^3 - 0.811*logr23^4 - $
         logq[iq]*(0.0737 - 0.0713*logr23 - 0.141*logr23^2 + 0.0373*logr23^3 - 0.058*logr23^4)
       logoh_lower = 9.40D + 4.65D*logr23 - 3.17D*logr23^2 - logq[iq]*(0.272D + 0.547D*logr23 - 0.513D*logr23^2)
       good1 = where((logoh_upper gt logoh_lower))
       good2 = where((logoh_lower[good1] gt 7.5))

       djs_oplot, logr23[good1], logoh_upper[good1], linestyle=linestyle[iq], $
         thick=postthick4, color=fsc_color(linecolor[iq],180+iq)
       djs_oplot, logr23[good1], logoh_lower[good1], linestyle=linestyle[iq], $
         thick=postthick4, color=fsc_color(linecolor[iq],180+iq)
;      djs_oplot, logr23[good1[good2]], logoh_lower[good1[good2]], line=0, thick=postthick

    endfor

    legend, 'log U = '+[string(logu[0],format='(F4.1)'),string(logu[1],format='(F4.1)'),$
      string(logu[2],format='(F4.1)')], /left, /bottom, charthick=postthick2, charsize=1.8, $
      line=linestyle, textcolor=fsc_color(replicate(textcolor1,3),150+findgen(3)), $
      color=fsc_color(linecolor,[180,181,182]), box=0, thick=postthick4
    
; now do stuff

    abund = im_abundance(ages[thisages],/silent,nmonte=500)
    plotsym, 0, 2.0, fill=1, thick=postthick4, color=fsc_color('red',150)
    oploterror, alog10(abund.zstrong_ew_r23), abund.zstrong_ew_12oh_kk04_upper, $
      abund.zstrong_ew_r23_err/abund.zstrong_ew_r23/alog(10.0), $
      abund.zstrong_ew_12oh_kk04_upper_err, ps=8, thick=postthick4, $
      errthick=postthick4, color=fsc_color('red',150), $
      errcolor=fsc_color('red',150)

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ------------------------------------------------------------   
; demonstration of how we measure abundances - full range
; ------------------------------------------------------------

    psname = 'r23_vs_12oh'
    if keyword_set(postscript) or keyword_set(pdf) then begin
       psfile = pspath+psname+'.ps'
       if keyword_set(postscript) then splog, 'Writing '+psfile
    endif else delvarx, psfile

    xmargin = [1.3,0.4] & ymargin = [0.4,1.1]
    xspace = 0.0 & yspace = 0.0
    width = 6.5 & height = 5.5
    xpage = total(height)+total(ymargin)+total(yspace)
    ypage = total(width)+total(xmargin)+total(xspace)

    xspace1 = xspace & yspace1 = yspace
    xmargin1 = xmargin & ymargin1 = ymargin
    width1 = width & height1 = height
    xpage1 = xpage & ypage1 = ypage
    
    arm_plotconfig, landscape=1, nx=1, ny=1, xmargin=xmargin1, $
      ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
      height=height1, coord=pos, xpage=xpage1, ypage=ypage1, $
      psfile=psfile, /writeover, bw=0L ;, /show
    cleanplot, /silent

    r23obs = 4.0 & r23obs_err = 1.5
    o32obs = 0.60 & o32obs_err = 0.03

    logr23obs = alog10(r23obs) & logr23obs_err = r23obs_err/r23obs/alog(10.0)
    logo32obs = alog10(o32obs) & logo32obs_err = o32obs_err/o32obs/alog(10.0)
    
;   logr23obs = 1.0 & logr23obs_err = 0.1
;   logo32obs = -0.5 & logo32obs_err = 0.05
    
    logq = alog10([5D6,4E7,1.5D8]) ; alog10(4E7)
    logu = logq-alog10(light)
;   logr23 = findgen((1.1-(-0.0))/0.001+1)*0.001-0.0
    logr23 = findgen((1.1-(-0.5))/0.001+1)*0.001-0.5
    linestyle = [0,3,5]
    linecolor = ['orange','forest green','dodger blue']
    
; plot the R23-O/H relation and the KK04 theoretical calibration

    ohrange = [7.1,9.4]
    xrange = [-0.1,1.05]
    yrange = ohrange

    xtitle = 'log ([O II] + [O III])/H\beta'
;   xtitle = 'log (R_{23})'
    ytitle = '12 + log (O/H)'
    
    plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), $
      ytitle=textoidl(ytitle), charsize=charsize1, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, $
      position=pos[*,0], xsty=1, ysty=1, color=fsc_color(textcolor1,150)
    
    for iq = 0L, n_elements(logq)-1L do begin
       
       logoh_upper = 9.72D - 0.777*logr23 - 0.951*logr23^2 - 0.072*logr23^3 - 0.811*logr23^4 - $
         logq[iq]*(0.0737 - 0.0713*logr23 - 0.141*logr23^2 + 0.0373*logr23^3 - 0.058*logr23^4)
       logoh_lower = 9.40D + 4.65D*logr23 - 3.17D*logr23^2 - logq[iq]*(0.272D + 0.547D*logr23 - 0.513D*logr23^2)
       good1 = where((logoh_upper gt logoh_lower))
       good2 = where((logoh_lower[good1] gt 7.5))

       djs_oplot, logr23[good1], logoh_upper[good1], linestyle=linestyle[iq], $
         thick=postthick4, color=fsc_color(linecolor[iq],180+iq)
       djs_oplot, logr23[good1], logoh_lower[good1], linestyle=linestyle[iq], $
         thick=postthick4, color=fsc_color(linecolor[iq],180+iq)
;      djs_oplot, logr23[good1[good2]], logoh_lower[good1[good2]], line=0, thick=postthick

    endfor

    xyouts, 0.40, 8.15, 'log U = '+string(logu[0],format='(F4.1)'), align=1.0, $
      charthick=postthick2, charsize=1.8, color=fsc_color(textcolor1,150)
    xyouts, 0.42, 7.85, string(logu[1],format='(F4.1)'), align=0.0, $
      charthick=postthick2, charsize=1.8, color=fsc_color(textcolor1,150)
    xyouts, 0.55, 7.35, string(logu[2],format='(F4.1)'), align=0.0, $
      charthick=postthick2, charsize=1.8, color=fsc_color(textcolor1,150)
    
; now do stuff

    abund = im_abundance(ages[thisages],/silent,nmonte=500)
    plotsym, 0, 2.0, fill=1, thick=postthick4, color=fsc_color('red',150)
    oploterror, alog10(abund.zstrong_ew_r23), abund.zstrong_ew_12oh_kk04_upper, $
      abund.zstrong_ew_r23_err/abund.zstrong_ew_r23/alog(10.0), $
      abund.zstrong_ew_12oh_kk04_upper_err, ps=8, thick=postthick4, $
      errthick=postthick4, color=fsc_color('red',150), $
      errcolor=fsc_color('red',150)

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

    stop    
    
; ---------------------------------------------------------------------------    
; Figure from Tremonti et al. (2007) - ZOOM
; ---------------------------------------------------------------------------    
    
    xsize = 8.0 & ysize = 10.5
    pagemaker, nx=1, ny=3, xspace=0.0, yspace=0.0, width=6.5, height=3.0*[1,1,1], $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=xsize, ypage=ysize, $
      position=pos, /normal

    psname = '07tremonti_myfig1_zoom'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xsize, ysize=ysize, silent=keyword_set(pdf)

    xrange = [2750,2850]
    xtitle = 'Rest Wavelength (\AA)'
;   ytitle = 'Relative Flux'
    ytitle = 'Flux (10^{-17} '+flam_units()+')'

; 0826+4305    

    yrange = [0,20]
    z = 0.603
    ss = rd1dspec('J0826+4305.ms.fits',/silent,datapath=$
      '/Users/ioannis/papers/proposals/HST/2008/')

    plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle='', ytitle='', $
      charsize=charsize1, charthick=postthick2, xsty=1, ysty=1, xrange=xrange, $
      position=pos[*,0], yrange=yrange, color=fsc_color(textcolor1,150), $
      xtickname=replicate(' ',10) ; ytickname=replicate(' ',10), 
    djs_oplot, ss.wave/(1.0+z), 1E17*ss.spec*(1.0+z), ps=10, thick=postthick2, $
      color=fsc_color(textcolor1,150)
    djs_oplot, 2796*[1,1], !y.crange, line=1, thick=postthick1, color=fsc_color('orange',150)
    djs_oplot, 2803*[1,1], !y.crange, line=1, thick=postthick1, color=fsc_color('orange',150)
    legend, 'SDSS J0826+4305', /right, /bottom, charsize=charsize2, $
      charthick=postthick2, box=0, textcolor=fsc_color(textcolor1,150), margin=0

; 0944+0930

    yrange = [0,14]
    z = 0.514
    ss = rd1dspec('J0944+0930.ms.fits',/silent,datapath=$
      '/Users/ioannis/papers/proposals/HST/2008/')

    plot, [0], [0], /nodata, /noerase, xthick=postthick1, ythick=postthick1, $
      xtitle='', ytitle=textoidl(ytitle), ytickinterval=5, $
      charsize=charsize1, charthick=postthick2, xsty=1, ysty=1, xrange=xrange, $
      position=pos[*,1], yrange=yrange, color=fsc_color(textcolor1,150), $
      xtickname=replicate(' ',10) ; ytickname=replicate(' ',10), 
    djs_oplot, ss.wave/(1.0+z), 1E17*ss.spec*(1.0+z), ps=10, thick=postthick2, $
      color=fsc_color(textcolor1,150)
    djs_oplot, 2796*[1,1], !y.crange, line=1, thick=postthick1, color=fsc_color('orange',150)
    djs_oplot, 2803*[1,1], !y.crange, line=1, thick=postthick1, color=fsc_color('orange',150)
    legend, 'SDSS J0944+0930', /right, /bottom, charsize=charsize2, $
      charthick=postthick2, box=0, textcolor=fsc_color(textcolor1,150), margin=0

; 1104+5946
    
    yrange = [0,14]
    z = 0.573
    ss = rd1dspec('J1104+5946.ms.fits',/silent,datapath=$
      '/Users/ioannis/papers/proposals/HST/2008/')

    plot, [0], [0], /nodata, /noerase, xthick=postthick1, ythick=postthick1, $
      xtitle=textoidl(xtitle), ytitle='', ytickinterval=5, $
      charsize=charsize1, charthick=postthick2, xsty=1, ysty=1, xrange=xrange, $
      position=pos[*,2], yrange=yrange, color=fsc_color(textcolor1,150) ;, $
;     ytickname=replicate(' ',10)
    djs_oplot, ss.wave/(1.0+z), 1E17*ss.spec*(1.0+z), ps=10, thick=postthick2, $
      color=fsc_color(textcolor1,150)
    djs_oplot, 2796*[1,1], !y.crange, line=1, thick=postthick1, color=fsc_color('orange',150)
    djs_oplot, 2803*[1,1], !y.crange, line=1, thick=postthick1, color=fsc_color('orange',150)
    legend, 'SDSS J1104+5946', /right, /bottom, charsize=charsize2, $
      charthick=postthick2, box=0, textcolor=fsc_color(textcolor1,150), margin=0

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

    stop    

; ---------------------------------------------------------------------------    
; Figure from Tremonti et al. (2007)
; ---------------------------------------------------------------------------    
    
    psname = '07tremonti_myfig1'
    if keyword_set(postscript) or keyword_set(pdf) then begin
       psfile = pspath+psname+'.ps'
       if keyword_set(postscript) then splog, 'Writing '+psfile
    endif else delvarx, psfile

    xmargin = [1.1,0.4] & ymargin = [0.4,1.1]
    xspace = 0.0 & yspace = 0.0
    width = 7.5 & height = 2.25*[1,1,1]
    xpage = total(height)+total(ymargin)+total(yspace)
    ypage = total(width)+total(xmargin)+total(xspace)

    xspace1 = xspace & yspace1 = yspace
    xmargin1 = xmargin & ymargin1 = ymargin
    width1 = width & height1 = height
    xpage1 = xpage & ypage1 = ypage

    arm_plotconfig, landscape=1, nx=1, ny=3, xmargin=xmargin1, $
      ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
      height=height1, coord=pos, xpage=xpage1, ypage=ypage1, $
      psfile=psfile, /writeover, bw=0L ;, /show
    cleanplot, /silent

    xrange = [2650,4090]
    xtitle = 'Rest Wavelength (\AA)'
;   ytitle = 'Relative Flux'
    ytitle = 'Flux (10^{-17} '+flam_units()+')'

; 0826+4305    

    yrange = [0,20]
    z = 0.603
    ss = rd1dspec('J0826+4305.ms.fits',/silent,datapath=$
      '/Users/ioannis/papers/proposals/HST/2008/')

    plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle='', ytitle='', $
      charsize=charsize1, charthick=postthick2, xsty=1, ysty=1, xrange=xrange, $
      position=pos[*,0], yrange=yrange, color=fsc_color(textcolor1,150), $
      xtickname=replicate(' ',10) ; ytickname=replicate(' ',10), 
    djs_oplot, ss.wave/(1.0+z), 1E17*ss.spec*(1.0+z), ps=10, thick=postthick2, $
      color=fsc_color(textcolor1,150)
    djs_oplot, 2796*[1,1], !y.crange, line=1, thick=postthick1, color=fsc_color('orange',150)
    djs_oplot, 2803*[1,1], !y.crange, line=1, thick=postthick1, color=fsc_color('orange',150)
    legend, 'SDSS J0826+4305', /right, /bottom, charsize=charsize2, $
      charthick=postthick2, box=0, textcolor=fsc_color(textcolor1,150), margin=0

; 0944+0930

    yrange = [0,14]
    z = 0.514
    ss = rd1dspec('J0944+0930.ms.fits',/silent,datapath=$
      '/Users/ioannis/papers/proposals/HST/2008/')

    plot, [0], [0], /nodata, /noerase, xthick=postthick1, ythick=postthick1, $
      xtitle='', ytitle=textoidl(ytitle), ytickinterval=5, $
      charsize=charsize1, charthick=postthick2, xsty=1, ysty=1, xrange=xrange, $
      position=pos[*,1], yrange=yrange, color=fsc_color(textcolor1,150), $
      xtickname=replicate(' ',10) ; ytickname=replicate(' ',10), 
    djs_oplot, ss.wave/(1.0+z), 1E17*ss.spec*(1.0+z), ps=10, thick=postthick2, $
      color=fsc_color(textcolor1,150)
    djs_oplot, 2796*[1,1], !y.crange, line=1, thick=postthick1, color=fsc_color('orange',150)
    djs_oplot, 2803*[1,1], !y.crange, line=1, thick=postthick1, color=fsc_color('orange',150)
    legend, 'SDSS J0944+0930', /right, /bottom, charsize=charsize2, $
      charthick=postthick2, box=0, textcolor=fsc_color(textcolor1,150), margin=0

; 1104+5946
    
    yrange = [0,14]
    z = 0.573
    ss = rd1dspec('J1104+5946.ms.fits',/silent,datapath=$
      '/Users/ioannis/papers/proposals/HST/2008/')

    plot, [0], [0], /nodata, /noerase, xthick=postthick1, ythick=postthick1, $
      xtitle=textoidl(xtitle), ytitle='', ytickinterval=5, $
      charsize=charsize1, charthick=postthick2, xsty=1, ysty=1, xrange=xrange, $
      position=pos[*,2], yrange=yrange, color=fsc_color(textcolor1,150) ;, $
;     ytickname=replicate(' ',10)
    djs_oplot, ss.wave/(1.0+z), 1E17*ss.spec*(1.0+z), ps=10, thick=postthick2, $
      color=fsc_color(textcolor1,150)
    djs_oplot, 2796*[1,1], !y.crange, line=1, thick=postthick1, color=fsc_color('orange',150)
    djs_oplot, 2803*[1,1], !y.crange, line=1, thick=postthick1, color=fsc_color('orange',150)
    legend, 'SDSS J1104+5946', /right, /bottom, charsize=charsize2, $
      charthick=postthick2, box=0, textcolor=fsc_color(textcolor1,150), margin=0

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

    stop    

return
end

