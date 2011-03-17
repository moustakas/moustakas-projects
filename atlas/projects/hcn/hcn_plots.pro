pro hcn_plots, atlas_hcn, sdss, encapsulated=encapsulated, postscript=postscript
; jm06jul20uofa - written
; jm06dec11nyu - update
; jm07sep24nyu - major changes to jive with the changes to
;                WRITE_ATLAS_HCN_SAMPLE  

; read the samples (see WRITE_ATLAS_HCN_SAMPLE)
    
    htmlbase = 'hcn_figures'

    path = '/Users/ioannis/ay/research/projects/atlas/projects/hcn/'
    html_path = '/Users/ioannis/ay/research/projects/atlas/papers/hcn/'
    pspath = html_path+htmlbase+'/'

    if (n_elements(atlas_hcn) eq 0L) then begin
       splog, 'Reading '+path+'atlas_hcn_sample.fits.gz'
       atlas_hcn = mrdfits(path+'atlas_hcn_sample.fits.gz',1,/silent)
    endif
    if (n_elements(sdss) eq 0L) then begin
       splog, 'Reading '+sdss_path()+'sdss_specdata_main_dr4.fits.gz'
       sdss = mrdfits(sdss_path()+'sdss_specdata_main_dr4.fits.gz',1,/silent)
       lineratio, sdss, 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
         x1, xerr1, y1, yerr1, index=indx1, nindex=nindx1, snrcut=5.0
       sdss = (temporary(sdss))[indx1]
    endif

    gao_hcn = rsex(path+'04gao_apj.dat')

; define some plotting variables    
    
    if keyword_set(postscript) then begin
       postthick1 = 4.0 
       postthick2 = 4.0 
       postthick3 = 3.3
       postthick4 = 6.0 
    endif else begin
       postthick1 = 2.0
       postthick2 = 1.0
       postthick3 = 1.0 
       postthick4 = 2.0 
       im_window, 0, xratio=0.5, /square
    endelse

    lsun = 3.826D33       ; bolometric solar luminosity [erg/s]
    kenn_sfr_ha_const = 7.9D-42 ; Salpeter IMF
    kenn_sfr_ir_const = 4.5D-44 ; Salpeter IMF

    lhcnaxis = findgen((10.0-6.0)/0.01+1)*0.01+6.0
    sfraxis = findgen((3.0-(-1.0))/0.01+1)*0.01+(-1.0)

    pivot_lhcn = 8.5
    pivot_sfr = 1.5

    arange1 = [0.0,8.0]
    arange2 = [0.0,6.0]
    sfrrange1 = [-0.2,3.1]
    sfrrange2 = [-0.8,3.0]
    lhcnrange = [6.8,10.0]    
    lirrange = [9.5,13.0]
    l24range = [7.0,11.2] ; [8.5,12.0]    
    l24lhaobsrange = [7.2,11.2] ; [7.0,11.0]    
;   l25range = [7.2,12.0] ; [8.5,12.0]    
;   l24lhaobsrange = [7.2,12.0] ; [7.0,11.0]    

    sfsym = 108 & sfcolor = 'red' & sfpsize = 1.5 & sffill = 1
    sfsym_ul = 114 & sfcolor_ul = 'red' & sfpsize_ul = 4.0 & sffill_ul = 1

    agnsym = 105 & agncolor = 'blue' & agnpsize = 1.8 & agnfill = 1
    agnsym_ul = 114 & agncolor_ul = 'blue' & agnpsize_ul = 4.0 & agnfill_ul = 1

    agnsfsym = 106 & agnsfcolor = 'green' & agnsfpsize = 1.8 & agnsffill = 1
    agnsfsym_ul = 114 & agnsfcolor_ul = 'green' & agnsfpsize_ul = 4.0 & agnsffill_ul = 1

    gaosym = 122 & gaocolor = 'grey' & gaopsize = 1.35 & gaofill = 1

    sdsssym = 108 & sdsscolor = 'grey' & sdsspsize = 0.1 & sdssfill = 1

; ---------------------------------------------------------------------------    
    
    dfpsplot, pspath+'hist_l25_lha_obs.ps', /square
    sf    = where(strtrim(atlas_hcn.class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn.class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn.class,2) eq 'HII/AGN',nagnsf)
    binsize = 0.3
    djs_plot, [0], [0], /nodata, charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, $
      xtitle=textoidl('log [L(25 \mu'+'m)/L(H\alpha)_{obs}]'), ytitle='Number', $
      xrange=[1,5], yrange=[0,7], xsty=1, ysty=1
    im_plothist, atlas_hcn[agnsf].l25-atlas_hcn[agnsf].lha_obs, bin=binsize, /overplot, line=5, thick=5.0
;     /fill, fcolor=fsc_color('grey',10), color=fsc_color('grey',10), line=0, thick=1.0
    im_plothist, atlas_hcn[sf].l25-atlas_hcn[sf].lha_obs, bin=binsize, /overplot, line=0, thick=3.5
    im_plothist, atlas_hcn[agn].l25-atlas_hcn[agn].lha_obs, bin=binsize, /overplot, line=1, thick=6.0
;     /fill, fcolor=fsc_color('dark grey',11), color=fsc_color('dark grey',10), line=0, thick=1.0
    legend, ['SF','SF/AGN','AGN'], line=[0,5,1], thick=4.0, charsize=2.0, charthick=2.0, /right, /top, box=0
    dfpsclose

; ---------------------------------------------------------------------------    
; L(HCN) vs L(24), L(Ha)+a*L(24), and L(IR) (3-panel)
; ---------------------------------------------------------------------------    

    psname = 'lhcn_vs_l24_l24lha_lir'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=10.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=3, xspace=0, yspace=0.0, width=7.0, height=3.0*[1,1,1], $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=10.5, $
      position=pos, /normal

;   im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=10.5, encapsulated=encapsulated
;   pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=7.0, height=4.5*[1,1], $
;     xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=10.5, $
;     position=pos, /normal

    xtitle = 'log L(HCN) (K km s^{-1} pc^{2})'
    xrange = lhcnrange
    yrange1 = l24range
    yrange2 = l24lhaobsrange
    yrange3 = lirrange
    charsize = 1.4
    
; --------------------------------------------------
; L(HCN) vs L(25)

    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.l24 gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.l24 gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn
    xerr = atlas_hcn[indx].gao_lhcn_err
    y = atlas_hcn[indx].l24 + alog10(0.031)
    yerr = atlas_hcn[indx].l24_err
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn
    xerr_ul = atlas_hcn[indx_ul].gao_lhcn_err
    y_ul = atlas_hcn[indx_ul].l24 + alog10(0.031)
    yerr_ul = atlas_hcn[indx_ul].l24_err

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; make the plot    
    
    ytitle = 'log [0.031*L(24) (L_{\odot})]'

    djs_plot, [0], [0], /nodata, xtitle='', xtickname=replicate(' ',10), ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange1, charsize=charsize, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0]
    im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=sffill, color=sfcolor
    oploterror, x[sf], y[sf], xerr[sf], yerr[sf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(sfcolor)
    im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=agnfill, color=agncolor
    oploterror, x[agn], y[agn], xerr[agn], yerr[agn], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agncolor)
    im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=agnsffill, color=agnsfcolor
    oploterror, x[agnsf], y[agnsf], xerr[agnsf], yerr[agnsf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick2

; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

; --------------------------------------------------
; L(HCN) vs L(Ha)_obs+a*L(24)

    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.l24_lha_obs gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.l24_lha_obs gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn
    xerr = atlas_hcn[indx].gao_lhcn_err
    y = atlas_hcn[indx].l24_lha_obs
    yerr = atlas_hcn[indx].l24_lha_obs_err
    photflag = replicate('Y',nindx)
;   photflag = atlas_hcn[indx].photflag
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn
    xerr_ul = atlas_hcn[indx_ul].gao_lhcn_err
    y_ul = atlas_hcn[indx_ul].l24_lha_obs
    yerr_ul = atlas_hcn[indx_ul].l24_lha_obs_err

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; make the plot    
    
    ytitle = 'log [L(H\alpha)_{obs}+0.031*L(24)] (L_{\odot})'
;   ytitle = 'log [L(H\alpha)_{obs}+a*L(24 \mu'+'m)] (L_{\odot})'

    djs_plot, [0], [0], /nodata, /noerase, xtitle='', xtickname=replicate(' ',10), ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange2, charsize=charsize, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,1]
    for ii = 0L, nsf-1L do begin
       im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=(photflag[sf[ii]] eq 'Y'), color=sfcolor
       oploterror, x[sf[ii]], y[sf[ii]], xerr[sf[ii]], yerr[sf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(sfcolor)
    endfor
    for ii = 0L, nagn-1L do begin
       im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=(photflag[agn[ii]] eq 'Y'), color=agncolor
       oploterror, x[agn[ii]], y[agn[ii]], xerr[agn[ii]], yerr[agn[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agncolor)
    endfor
    for ii = 0L, nagnsf-1L do begin
       im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=(photflag[agnsf[ii]] eq 'Y'), color=agnsfcolor
       oploterror, x[agnsf[ii]], y[agnsf[ii]], xerr[agnsf[ii]], yerr[agnsf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    endfor
    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick2

; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

; --------------------------------------------------
; L(HCN) vs L(IR)

    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.lir gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.lir gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn
    xerr = atlas_hcn[indx].gao_lhcn_err
    y = atlas_hcn[indx].lir
    yerr = atlas_hcn[indx].lir_err
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn
    xerr_ul = atlas_hcn[indx_ul].gao_lhcn_err
    y_ul = atlas_hcn[indx_ul].lir
    yerr_ul = atlas_hcn[indx_ul].lir_err

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; make the plot    
    
    ytitle = 'log L(IR) (L_{\odot})]'

    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange3, charsize=charsize, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,2]
    im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=sffill, color=sfcolor
    oploterror, x[sf], y[sf], xerr[sf], yerr[sf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(sfcolor)
    im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=agnfill, color=agncolor
    oploterror, x[agn], y[agn], xerr[agn], yerr[agn], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agncolor)
    im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=agnsffill, color=agnsfcolor
    oploterror, x[agnsf], y[agnsf], xerr[agnsf], yerr[agnsf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick2

; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

    im_openclose, postscript=postscript, /close

stop    
    
; ---------------------------------------------------------------------------    
; SDSS/BPT diagram
; ---------------------------------------------------------------------------    

    psname = 'sdss_niiha_vs_oiiihb'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.3, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.6,1.1], xpage=8.5, ypage=8.3, $
      position=pos, /normal

; ATLAS/HCN

    indx = where((atlas_hcn.oiii_hb gt -900.0) and (atlas_hcn.nii_ha gt -900.0),nindx)
    x = atlas_hcn[indx].nii_ha
    y = atlas_hcn[indx].oiii_hb
    xerr = atlas_hcn[indx].nii_ha_err
    yerr = atlas_hcn[indx].oiii_hb_err
    
    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

; SDSS
    
    xsdss = alog10(sdss.nii_6584[0]/sdss.h_alpha[0])
    ysdss = alog10(sdss.oiii_5007[0]/sdss.h_beta[0])
    
    xtitle = 'log ([N II] \lambda6584/H\alpha)'
    ytitle = 'log ([O III] \lambda5007/H\beta)'

    xrange = [-2.0,1.0]
    yrange = [-1.2,1.2]

;   im_symbols, sdsssym, thick=postthick1, psize=sdsspsize, fill=sdssfill, color=sdsscolor
    hogg_scatterplot, xsdss, ysdss, weight=weight, /outliers, label=0, outpsym=3, outsymsize=1.0, outcolor=sdsscolor, $
      levels=errorf([1.0,2.0,3.0]/2.0),$
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), charsize=2.0, $
      charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0]
    im_symbols, sfsym, thick=postthick4, psize=sfpsize*0.8, fill=sffill, color=sfcolor
    oploterror, x[sf], y[sf], xerr[sf], yerr[sf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(sfcolor)
    im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=agnfill, color=agncolor
    oploterror, x[agn], y[agn], xerr[agn], yerr[agn], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agncolor)
    im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize*0.7, fill=agnsffill, color=agnsfcolor
    oploterror, x[agnsf], y[agnsf], xerr[agnsf], yerr[agnsf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)
    for ii = 0L, nindx-1L do begin
       align = 0.0 & xoff = 0.05 & yoff = 0.0
       if strmatch(galaxy[ii],'*ngc6240*',/fold) then begin
          align = 0.0 & xoff = +0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*ic1623*',/fold) then begin
          align = 1.0 & xoff = -0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*ugc05101*',/fold) then begin
          align = 1.0 & xoff = -0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*iras05189*',/fold) then begin
          align = 1.0 & xoff = -0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*ugc08696*',/fold) then begin
          align = 0.0 & xoff = +0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*iras17208*',/fold) then begin
          align = 1.0 & xoff = -0.07 & yoff = -0.0
       endif
       if strmatch(galaxy[ii],'*ngc3628*',/fold) then begin
          align = 1.0 & xoff = -0.1 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*ngc0660*',/fold) then begin
          align = 0.0 & xoff = +0.05 & yoff = -0.04
       endif
       if strmatch(galaxy[ii],'*ngc3079*',/fold) then begin
          align = 1.0 & xoff = -0.06 & yoff = +0.01
       endif
       if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin
          align = 1.0 & xoff = -0.02 & yoff = -0.07
       endif
       if strmatch(galaxy[ii],'*ngc1144*',/fold) then begin
          align = 1.0 & xoff = -0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*ngc3893*',/fold) then begin
          align = 1.0 & xoff = -0.05 & yoff = 0.0
       endif
       if strmatch(galaxy[ii],'*ngc2903*',/fold) then begin
          align = 1.0 & xoff = -0.05 & yoff = -0.1
       endif
       if strmatch(galaxy[ii],'*ugc04881*',/fold) then begin
          align = 0.0 & xoff = +0.05 & yoff = -0.06
       endif
       if strmatch(galaxy[ii],'*ic5179*',/fold) then begin
          align = 1.0 & xoff = -0.05 & yoff = -0.05
       endif
       if strmatch(galaxy[ii],'*ngc7771*',/fold) then begin
          align = 0.0 & xoff = +0.04 & yoff = +0.025
       endif
       if strmatch(galaxy[ii],'*ic0883*',/fold) then begin
          align = 0.0 & xoff = +0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*mrk0331*',/fold) then begin
          align = 0.0 & xoff = +0.06 & yoff = -0.03
       endif
       if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
          align = 0.0 & xoff = +0.05 & yoff = 0.0
       endif
       if strmatch(galaxy[ii],'*ngc6701*',/fold) then begin
          align = 0.0 & xoff = +0.05 & yoff = +0.01
       endif
       if strmatch(galaxy[ii],'*ngc4414*',/fold) then begin
          align = 1.0 & xoff = +0.0 & yoff = +0.04
       endif
       if strmatch(galaxy[ii],'*ngc0520*',/fold) then begin
          align = 0.0 & xoff = +0.05 & yoff = -0.03
       endif
       if strmatch(galaxy[ii],'*ngc7469*',/fold) then begin
          align = 0.0 & xoff = +0.06 & yoff = -0.02
       endif
       if strmatch(galaxy[ii],'*ngc2146*',/fold) then begin
          align = 1.0 & xoff = -0.02 & yoff = +0.02
       endif
       if strmatch(galaxy[ii],'*ngc3893*',/fold) then begin
          align = 1.0 & xoff = -0.03 & yoff = -0.06
       endif
       if strmatch(galaxy[ii],'*ngc5194*',/fold) then begin
          align = 0.0 & xoff = +0.02 & yoff = -0.09
       endif
       xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
         charsize=1.0, charthick=postthick3, /data, align=align
    endfor

; overplot the mixing lines

    models = kewley_bpt_lines(/kauffmann,_extra=extra)
    oplot, models.x_nii, models.y_nii, line=0, thick=postthick1
    models = kewley_bpt_lines(_extra=extra)
    oplot, models.x_nii, models.y_nii, line=2, thick=postthick1

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
; L(HCN) vs SFR (2-panel)
; ---------------------------------------------------------------------------    

    psname = 'lhcn_vs_sfr_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=5.7, height=4.0*[1,1], $
      xmargin=[1.4,1.4], ymargin=[0.4,1.1], xpage=8.5, ypage=9.5, $
      position=pos, /normal

;   xtitle = 'log L(HCN)'
    xtitle = 'log L(HCN) (K km s^{-1} pc^{2})'
    xrange = lhcnrange
    yrange = sfrrange1
    charsize = 1.6
    
; --------------------------------------------------
; L(HCN) vs SFR[L(IR)]

    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.sfr_lir gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.sfr_lir gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn
    xerr = atlas_hcn[indx].gao_lhcn_err
    y = atlas_hcn[indx].sfr_lir
    yerr = atlas_hcn[indx].sfr_lir_err
    y2 = atlas_hcn[indx].lir
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn
    xerr_ul = atlas_hcn[indx_ul].gao_lhcn_err
    y_ul = atlas_hcn[indx_ul].sfr_lir
    yerr_ul = atlas_hcn[indx_ul].sfr_lir_err

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; Gao & Solomon points    
    
    indx_gao = cmset_op(cmset_op(strtrim(gao_hcn.galaxy,2),'and',/not2,$
      strtrim(atlas_hcn.gao_galaxy,2),/index),'and',where(gao_hcn.lhcn gt 0.0))
    xgao = alog10(gao_hcn[indx_gao].lhcn*1D8)
    ygao = alog10(gao_hcn[indx_gao].lir*1D10) + alog10(kenn_sfr_ir_const) + alog10(lsun)

; fit to the SFR points, *not* the luminosity

    dofit = where((galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
    coeff = robust_linefit(x[dofit]-pivot_lhcn,y[dofit],yfit,sig,coeff_err,/bisect,bisquare_limit=6.0)
    coeff[0] = coeff[0] - coeff[1]*pivot_lhcn

    coeff_slope_one = [-999.0,1.0]
    coeff_slope_one[0] = interpol(poly(lhcnaxis,coeff),lhcnaxis,pivot_lhcn)-coeff_slope_one[1]*pivot_lhcn
    coeff_slope_one[0] = coeff_slope_one[0] - 0.1 ; offset for clarity
    splog, psname, ndofit
    niceprint, coeff, coeff_err;, coeff_slope_one
    
; make the plot    
    
    ytitle = 'log \psi[L(IR)] (M_{\odot} yr^{-1})'
    ytitle2 = 'log L(IR) (L_{\odot})'

    djs_plot, [0], [0], /nodata, xtitle='', xtickname=replicate(' ',10), ytitle=ytitle, xsty=1, ysty=9, $
      xrange=xrange, yrange=yrange, charsize=charsize, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0]
    axis, /yaxis, yrange=interpol(y2,y,!y.crange), ysty=1, ytitle=textoidl(ytitle2), $
      charsize=charsize, charthick=postthick2, ythick=postthick1
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff), line=5, thick=postthick4
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff_slope_one), line=0, thick=postthick1
    im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=sffill, color=sfcolor
    oploterror, x[sf], y[sf], xerr[sf], yerr[sf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(sfcolor)
    im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=agnfill, color=agncolor
    oploterror, x[agn], y[agn], xerr[agn], yerr[agn], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agncolor)
    im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=agnsffill, color=agnsfcolor
    oploterror, x[agnsf], y[agnsf], xerr[agnsf], yerr[agnsf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    im_symbols, gaosym, thick=postthick1, psize=gaopsize, fill=gaofill, color=gaocolor
    djs_oplot, xgao, ygao, psym=8
    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick2

; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

; label some points

    for ii = 0L, nindx-1L do begin
       if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
          align = 1.0 & xoff = -0.14 & yoff = -0.03
          xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
            charsize=1.2, charthick=postthick2, /data, align=align
       endif
       if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin ; = ARP299
          align = 1.0 & xoff = -0.14 & yoff = -0.03
          xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
            charsize=1.2, charthick=postthick2, /data, align=align
       endif
       if strmatch(galaxy[ii],'*iras05189-2524*',/fold) then begin ; = ARP299
          align = 1.0 & xoff = -0.14 & yoff = -0.03
          xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
            charsize=1.1, charthick=postthick2, /data, align=align
       endif
    endfor    

; --------------------------------------------------
; L(HCN) vs SFR[L(Ha)_obs+a*L(24)]

    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.sfr_l24_lha_obs gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.sfr_l24_lha_obs gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn
    xerr = atlas_hcn[indx].gao_lhcn_err
    y = atlas_hcn[indx].sfr_l24_lha_obs
    yerr = atlas_hcn[indx].sfr_l24_lha_obs_err
    photflag = atlas_hcn[indx].photflag
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn
    xerr_ul = atlas_hcn[indx_ul].gao_lhcn_err
    y_ul = atlas_hcn[indx_ul].sfr_l24_lha_obs
    yerr_ul = atlas_hcn[indx_ul].sfr_l24_lha_obs_err

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; fit to the data

    dofit = where((galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
;   dofit = where((photflag eq 'Y') and (galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
    coeff = robust_linefit(x[dofit]-pivot_lhcn,y[dofit],yfit,sig,coeff_err,/bisect,bisquare_limit=6.0)
    coeff[0] = coeff[0] - coeff[1]*pivot_lhcn

    coeff_slope_one = [-999.0,1.0]
    coeff_slope_one[0] = interpol(poly(lhcnaxis,coeff),lhcnaxis,pivot_lhcn)-coeff_slope_one[1]*pivot_lhcn
    coeff_slope_one[0] = coeff_slope_one[0] - 0.1 ; offset for clarity
    splog, psname, ndofit
    niceprint, coeff, coeff_err;, coeff_slope_one
    
; make the plot    
    
    ytitle = 'log \psi[L(H\alpha)_{obs}+a*L(24 \mu'+'m)] (M_{\odot} yr^{-1})'
    ytitle2 = 'log [L(H\alpha)_{obs}+a*L(24 \mu'+'m)] (L_{\odot})'

    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle=ytitle, xsty=1, ysty=9, $
      xrange=xrange, yrange=yrange, charsize=charsize, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,1]
    axis, /yaxis, yrange=interpol(atlas_hcn[indx].l24_lha_obs,y,!y.crange), ysty=1, ytitle=textoidl(ytitle2), $
      charsize=charsize, charthick=postthick2, ythick=postthick1
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff), line=5, thick=postthick4
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff_slope_one), line=0, thick=postthick1
    for ii = 0L, nsf-1L do begin
       im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=(photflag[sf[ii]] eq 'Y'), color=sfcolor
       oploterror, x[sf[ii]], y[sf[ii]], xerr[sf[ii]], yerr[sf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(sfcolor)
    endfor
    for ii = 0L, nagn-1L do begin
       im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=(photflag[agn[ii]] eq 'Y'), color=agncolor
       oploterror, x[agn[ii]], y[agn[ii]], xerr[agn[ii]], yerr[agn[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agncolor)
    endfor
    for ii = 0L, nagnsf-1L do begin
       im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=(photflag[agnsf[ii]] eq 'Y'), color=agnsfcolor
       oploterror, x[agnsf[ii]], y[agnsf[ii]], xerr[agnsf[ii]], yerr[agnsf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    endfor
    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick2

; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

; label some points

    for ii = 0L, nindx-1L do begin
       if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
          align = 1.0 & xoff = -0.14 & yoff = -0.03
          xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
            charsize=1.1, charthick=postthick2, /data, align=align
       endif
       if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin ; = ARP299
          align = 1.0 & xoff = -0.14 & yoff = -0.03
          xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
            charsize=1.1, charthick=postthick2, /data, align=align
       endif
       if strmatch(galaxy[ii],'*iras05189-2524*',/fold) then begin ; = ARP299
          align = 1.0 & xoff = -0.14 & yoff = -0.03
          xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
            charsize=1.1, charthick=postthick2, /data, align=align
       endif
    endfor    

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
; L(HCN)/L(CO) vs L/L(CO) (2-panel)
; ---------------------------------------------------------------------------    

    psname = 'lhcn_lco_vs_l_lco_2panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=10.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=6.5, height=4.5*[1,1], $
      xmargin=[1.4,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=10.5, $
      position=pos, /normal

    xtitle = 'log [L(HCN) / L(CO)]'

    xrange = [-1.9,-0.3]
    charsize = 1.7
    
; --------------------------------------------------
; L(HCN)/L(CO) vs L(IR)/L(CO)

    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.lir gt -900.0) and $
      (atlas_hcn.gao_lco gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.lir gt -900.0) and $
      (atlas_hcn.gao_lco gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn - atlas_hcn[indx].gao_lco
    xerr = sqrt(atlas_hcn[indx].gao_lhcn_err^2.0 + atlas_hcn[indx].gao_lco_err^2.0)
    y = atlas_hcn[indx].lir - atlas_hcn[indx].gao_lco
    yerr = sqrt(atlas_hcn[indx].lir_err^2.0 + atlas_hcn[indx].gao_lco_err^2.0)
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn - atlas_hcn[indx_ul].gao_lco
    xerr_ul = sqrt(atlas_hcn[indx_ul].gao_lhcn_err^2.0 + atlas_hcn[indx_ul].gao_lco_err^2.0)
    y_ul = atlas_hcn[indx_ul].lir - atlas_hcn[indx_ul].gao_lco
    yerr_ul = sqrt(atlas_hcn[indx_ul].lir_err^2.0 + atlas_hcn[indx_ul].gao_lco_err^2.0)

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; Gao & Solomon points    
    
    indx_gao = cmset_op(cmset_op(cmset_op(strtrim(gao_hcn.galaxy,2),'and',/not2,$
      strtrim(atlas_hcn.gao_galaxy,2),/index),'and',where(gao_hcn.lhcn gt 0.0)),$
      'and',where(gao_hcn.lco gt 0.0))
    xgao = alog10(gao_hcn[indx_gao].lhcn*1D8) - alog10(gao_hcn[indx_gao].lco*1D8)
    ygao = alog10(gao_hcn[indx_gao].lir*1D10) - alog10(gao_hcn[indx_gao].lco*1D8)

; make the plot    

    yrange = [0.7,2.7]
    ytitle = 'log [L(IR) / L(CO)]'
;   ytitle = 'log [L(IR) / L(CO)] !c !c (L_{\odot} / K km s^{-1} pc^{2})'

    djs_plot, [0], [0], /nodata, xtitle='', xtickname=replicate(' ',10), ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, charsize=charsize, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0]
    im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=sffill, color=sfcolor
    oploterror, x[sf], y[sf], xerr[sf], yerr[sf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(sfcolor)
    im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=agnfill, color=agncolor
    oploterror, x[agn], y[agn], xerr[agn], yerr[agn], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agncolor)
    im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=agnsffill, color=agnsfcolor
    oploterror, x[agnsf], y[agnsf], xerr[agnsf], yerr[agnsf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    im_symbols, gaosym, thick=postthick1, psize=gaopsize, fill=gaofill, color=gaocolor
    djs_oplot, xgao, ygao, psym=8
    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick2
    
; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

; label some points

;   for ii = 0L, nindx-1L do begin
;      if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
;         align = 1.0 & xoff = -0.05 & yoff = -0.13
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;      if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin ; = ARP299
;         align = 0.0 & xoff = +0.05 & yoff = +0.08
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;   endfor    

; --------------------------------------------------
; L(HCN)/L(CO) vs [L(Ha)_obs+a*L(24)]/L(CO)

    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.sfr_l24_lha_obs gt -900.0) and $
      (atlas_hcn.gao_lco gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.sfr_l24_lha_obs gt -900.0) and $
      (atlas_hcn.gao_lco gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn - atlas_hcn[indx].gao_lco
    xerr = sqrt(atlas_hcn[indx].gao_lhcn_err^2.0 + atlas_hcn[indx].gao_lco_err^2.0)
    y = atlas_hcn[indx].l24_lha_obs - atlas_hcn[indx].gao_lco
    yerr = sqrt(atlas_hcn[indx].l24_lha_obs_err^2.0 + atlas_hcn[indx].gao_lco_err^2.0)
    photflag = atlas_hcn[indx].photflag
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn - atlas_hcn[indx_ul].gao_lco
    xerr_ul = sqrt(atlas_hcn[indx_ul].gao_lhcn_err^2.0 + atlas_hcn[indx_ul].gao_lco_err^2.0)
    y_ul = atlas_hcn[indx_ul].l24_lha_obs - atlas_hcn[indx_ul].gao_lco
    yerr_ul = sqrt(atlas_hcn[indx_ul].l24_lha_obs_err^2.0 + atlas_hcn[indx_ul].gao_lco_err^2.0)

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; make the plot    
    
    yrange = [-1.7,0.7]
    ytitle = 'log [L(H\alpha)_{obs}+a*L(24 \mu'+'m)] / L(CO)]'
;   ytitle = 'log [L(H\alpha)_{obs}+a*L(24 \mu'+'m)] / L(CO)] '+$
;     '!c !c (L_{\odot} / K km s^{-1} pc^{2})'

    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, charsize=charsize, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,1]
    for ii = 0L, nsf-1L do begin
       im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=(photflag[sf[ii]] eq 'Y'), color=sfcolor
       oploterror, x[sf[ii]], y[sf[ii]], xerr[sf[ii]], yerr[sf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(sfcolor)
    endfor
    for ii = 0L, nagn-1L do begin
       im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=(photflag[agn[ii]] eq 'Y'), color=agncolor
       oploterror, x[agn[ii]], y[agn[ii]], xerr[agn[ii]], yerr[agn[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agncolor)
    endfor
    for ii = 0L, nagnsf-1L do begin
       im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=(photflag[agnsf[ii]] eq 'Y'), color=agnsfcolor
       oploterror, x[agnsf[ii]], y[agnsf[ii]], xerr[agnsf[ii]], yerr[agnsf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    endfor
    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick2

; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

; label some points

;   for ii = 0L, nindx-1L do begin
;      if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
;         align = 1.0 & xoff = -0.05 & yoff = -0.13
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;      if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin ; = ARP299
;         align = 0.0 & xoff = +0.05 & yoff = +0.08
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;   endfor    

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
; A(IR) vs A(Ha)
; ---------------------------------------------------------------------------    

    psname = 'a_ir_vs_a_ha'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.2, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.7, height=6.7, $
      xmargin=[1.4,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.2, $
      position=pos, /normal
    
    indx = where((atlas_hcn.a_ir gt -900.0) and (atlas_hcn.a_ha gt -900.0),nindx)
    
    x = atlas_hcn[indx].a_ir
    xerr = atlas_hcn[indx].a_ir_err
    y = atlas_hcn[indx].a_ha
    yerr = atlas_hcn[indx].a_ha_err
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

; make the plot    
    
    xtitle = 'A(IR) (mag)'
    ytitle = 'A(H\alpha) (mag)'

    xrange = arange1
    yrange = arange2

    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, charsize=2.0, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0]
    djs_oplot, [0,10], [0,10], line=0, thick=postthick1
    im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=sffill, color=sfcolor
    oploterror, x[sf], y[sf], xerr[sf], yerr[sf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(sfcolor)
    im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=agnfill, color=agncolor
    oploterror, x[agn], y[agn], xerr[agn], yerr[agn], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agncolor)
    im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=agnsffill, color=agnsfcolor
    oploterror, x[agnsf], y[agnsf], xerr[agnsf], yerr[agnsf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agnsfcolor)

    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
; L(HCN) vs SFR (3-panel)
; ---------------------------------------------------------------------------    

    psname = 'lhcn_vs_sfr_3panel'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=10.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=3, xspace=0, yspace=0.0, width=6.7, height=3.0*[1,1,1], $
      xmargin=[1.4,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=10.5, $
      position=pos, /normal

    xtitle = 'log L(HCN) (K km s^{-1} pc^{2})'
    xrange = lhcnrange
    yrange = sfrrange2
    charsize = 2.0
    
; --------------------------------------------------
; L(HCN) vs SFR[L(IR)]

    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.sfr_lir gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.sfr_lir gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn
    xerr = atlas_hcn[indx].gao_lhcn_err
    y = atlas_hcn[indx].sfr_lir
    yerr = atlas_hcn[indx].sfr_lir_err
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn
    xerr_ul = atlas_hcn[indx_ul].gao_lhcn_err
    y_ul = atlas_hcn[indx_ul].sfr_lir
    yerr_ul = atlas_hcn[indx_ul].sfr_lir_err

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; fit to the data

    dofit = where((galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
    coeff = robust_linefit(x[dofit]-pivot_lhcn,y[dofit],yfit,sig,coeff_err,/bisect,bisquare_limit=6.0)
    coeff[0] = coeff[0] - coeff[1]*pivot_lhcn

    coeff_slope_one = [-999.0,1.0]
    coeff_slope_one[0] = interpol(poly(lhcnaxis,coeff),lhcnaxis,pivot_lhcn)-coeff_slope_one[1]*pivot_lhcn
    splog, psname, ndofit
    niceprint, coeff, coeff_err;, coeff_slope_one
    
; make the plot    
    
    ytitle = 'log \psi[L(IR)]' ;  (M_{\odot} yr^{-1})'

    djs_plot, [0], [0], /nodata, xtitle='', xtickname=replicate(' ',10), ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, charsize=charsize, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0]
;   djs_oplot, lhcnaxis, poly(lhcnaxis,coeff), line=5, thick=postthick4
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff_slope_one), line=0, thick=postthick1
    im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=sffill, color=sfcolor
    oploterror, x[sf], y[sf], xerr[sf], yerr[sf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(sfcolor)
    im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=agnfill, color=agncolor
    oploterror, x[agn], y[agn], xerr[agn], yerr[agn], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agncolor)
    im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=agnsffill, color=agnsfcolor
    oploterror, x[agnsf], y[agnsf], xerr[agnsf], yerr[agnsf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    legend, '(a)', /left, /top, box=0, charsize=charsize, charthick=postthick2
    
; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

; label some points

;   for ii = 0L, nindx-1L do begin
;      if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
;         align = 1.0 & xoff = -0.14 & yoff = -0.03
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;      if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin ; = ARP299
;         align = 1.0 & xoff = -0.14 & yoff = -0.03
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;   endfor    

; --------------------------------------------------
; L(HCN) vs SFR[L(Ha)_obs+a*L(24)]

    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.sfr_l24_lha_obs gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.sfr_l24_lha_obs gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn
    xerr = atlas_hcn[indx].gao_lhcn_err
    y = atlas_hcn[indx].sfr_l24_lha_obs
    yerr = atlas_hcn[indx].sfr_l24_lha_obs_err
    photflag = atlas_hcn[indx].photflag
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn
    xerr_ul = atlas_hcn[indx_ul].gao_lhcn_err
    y_ul = atlas_hcn[indx_ul].sfr_l24_lha_obs
    yerr_ul = atlas_hcn[indx_ul].sfr_l24_lha_obs_err

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; fit to the data

    dofit = where((galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
;   dofit = where((photflag eq 'Y') and (galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
    coeff = robust_linefit(x[dofit]-pivot_lhcn,y[dofit],yfit,sig,coeff_err,/bisect,bisquare_limit=6.0)
    coeff[0] = coeff[0] - coeff[1]*pivot_lhcn

    coeff_slope_one = [-999.0,1.0]
    coeff_slope_one[0] = interpol(poly(lhcnaxis,coeff),lhcnaxis,pivot_lhcn)-coeff_slope_one[1]*pivot_lhcn
    splog, psname, ndofit
    niceprint, coeff, coeff_err;, coeff_slope_one
    
; make the plot    
    
    ytitle = 'log \psi[L(H\alpha)_{obs}+a*L(24)]'; (M_{\odot} yr^{-1})'
;   ytitle = 'log \psi[L(H\alpha)_{obs}+a*L(24 \mu'+'m)]'; (M_{\odot} yr^{-1})'

    djs_plot, [0], [0], /nodata, /noerase, xtitle='', xtickname=replicate(' ',10), ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, charsize=charsize, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,1]
;   djs_oplot, lhcnaxis, poly(lhcnaxis,coeff), line=5, thick=postthick4
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff_slope_one), line=0, thick=postthick1
    for ii = 0L, nsf-1L do begin
       im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=(photflag[sf[ii]] eq 'Y'), color=sfcolor
       oploterror, x[sf[ii]], y[sf[ii]], xerr[sf[ii]], yerr[sf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(sfcolor)
    endfor
    for ii = 0L, nagn-1L do begin
       im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=(photflag[agn[ii]] eq 'Y'), color=agncolor
       oploterror, x[agn[ii]], y[agn[ii]], xerr[agn[ii]], yerr[agn[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agncolor)
    endfor
    for ii = 0L, nagnsf-1L do begin
       im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=(photflag[agnsf[ii]] eq 'Y'), color=agnsfcolor
       oploterror, x[agnsf[ii]], y[agnsf[ii]], xerr[agnsf[ii]], yerr[agnsf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    endfor
    legend, '(b)', /left, /top, box=0, charsize=charsize, charthick=postthick2

; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

; label some points

;   for ii = 0L, nindx-1L do begin
;      if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
;         align = 1.0 & xoff = -0.14 & yoff = -0.03
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;      if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin ; = ARP299
;         align = 1.0 & xoff = -0.14 & yoff = -0.03
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;   endfor    
    
; --------------------------------------------------
; L(HCN) vs SFR[L(Ha)_cor]

    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.sfr_lha_cor gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.sfr_lha_cor gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn
    xerr = atlas_hcn[indx].gao_lhcn_err
    y = atlas_hcn[indx].sfr_lha_cor
    yerr = atlas_hcn[indx].sfr_lha_cor_err
    photflag = atlas_hcn[indx].photflag
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn
    xerr_ul = atlas_hcn[indx_ul].gao_lhcn_err
    y_ul = atlas_hcn[indx_ul].sfr_lha_cor
    yerr_ul = atlas_hcn[indx_ul].sfr_lha_cor_err

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; fit to the data

    dofit = where((galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
;   dofit = where((photflag eq 'Y') and (galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
    coeff = robust_linefit(x[dofit]-pivot_lhcn,y[dofit],yfit,sig,coeff_err,/bisect,bisquare_limit=6.0)
    coeff[0] = coeff[0] - coeff[1]*pivot_lhcn

    coeff_slope_one = [-999.0,1.0]
    coeff_slope_one[0] = interpol(poly(lhcnaxis,coeff),lhcnaxis,pivot_lhcn)-coeff_slope_one[1]*pivot_lhcn
    splog, psname, ndofit
    niceprint, coeff, coeff_err;, coeff_slope_one
    
; make the plot    
    
    ytitle = 'log \psi[L(H\alpha)_{cor}]' ;  (M_{\odot} yr^{-1})'

    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, charsize=charsize, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,2]
;   djs_oplot, lhcnaxis, poly(lhcnaxis,coeff), line=5, thick=postthick4
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff_slope_one), line=0, thick=postthick1
    for ii = 0L, nsf-1L do begin
       im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=(photflag[sf[ii]] eq 'Y'), color=sfcolor
       oploterror, x[sf[ii]], y[sf[ii]], xerr[sf[ii]], yerr[sf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(sfcolor)
    endfor
    for ii = 0L, nagn-1L do begin
       im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=(photflag[agn[ii]] eq 'Y'), color=agncolor
       oploterror, x[agn[ii]], y[agn[ii]], xerr[agn[ii]], yerr[agn[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agncolor)
    endfor
    for ii = 0L, nagnsf-1L do begin
       im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=(photflag[agnsf[ii]] eq 'Y'), color=agnsfcolor
       oploterror, x[agnsf[ii]], y[agnsf[ii]], xerr[agnsf[ii]], yerr[agnsf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    endfor
    legend, '(c)', /left, /top, box=0, charsize=charsize, charthick=postthick2

; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

; label some points

;   for ii = 0L, nindx-1L do begin
;      if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
;         align = 1.0 & xoff = -0.14 & yoff = -0.03
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;      if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin ; = ARP299
;         align = 1.0 & xoff = -0.14 & yoff = -0.03
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;   endfor    
    
    
    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
; SFR[L(IR)] vs SFR[L(Ha)_obs+a*L(24)]
; ---------------------------------------------------------------------------    

    psname = 'sfr_lir_vs_sfr_l24_lha_obs'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.2, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.7, height=6.7, $
      xmargin=[1.4,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.2, $
      position=pos, /normal
    
    indx = where((atlas_hcn.sfr_lir gt -900.0) and (atlas_hcn.sfr_l24_lha_obs gt -900.0),nindx)
    
    x = atlas_hcn[indx].sfr_lir
    xerr = atlas_hcn[indx].sfr_lir_err
    y = atlas_hcn[indx].sfr_l24_lha_obs
    yerr = atlas_hcn[indx].sfr_l24_lha_obs_err
    photflag = atlas_hcn[indx].photflag
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

; fit to the data

    dofit = where((galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
;   dofit = where((photflag eq 'Y') and (galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
    coeff = robust_linefit(x[dofit]-pivot_sfr,y[dofit],yfit,sig,coeff_err,/bisect,bisquare_limit=6.0)
    coeff[0] = coeff[0] - coeff[1]*pivot_sfr

    coeff_slope_one = [-999.0,1.0]
    coeff_slope_one[0] = interpol(poly(sfraxis,coeff),sfraxis,pivot_sfr)-coeff_slope_one[1]*pivot_sfr
    splog, psname, ndofit
    niceprint, coeff, coeff_err;, coeff_slope_one
    
; make the plot    
    
    xtitle = 'log \psi[L(IR)] (M_{\odot} yr^{-1})'
    ytitle = 'log \psi[L(H\alpha)_{obs}+a*L(24 \mu'+'m)] (M_{\odot} yr^{-1})'

    xrange = sfrrange2
    yrange = sfrrange2

    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, charsize=2.0, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0]
    djs_oplot, sfraxis, poly(sfraxis,coeff), line=5, thick=postthick4
    djs_oplot, sfraxis, poly(sfraxis,coeff_slope_one), line=0, thick=postthick1
    for ii = 0L, nsf-1L do begin
       im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=(photflag[sf[ii]] eq 'Y'), color=sfcolor
       oploterror, x[sf[ii]], y[sf[ii]], xerr[sf[ii]], yerr[sf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(sfcolor)
    endfor
    for ii = 0L, nagn-1L do begin
       im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=(photflag[agn[ii]] eq 'Y'), color=agncolor
       oploterror, x[agn[ii]], y[agn[ii]], xerr[agn[ii]], yerr[agn[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agncolor)
    endfor
    for ii = 0L, nagnsf-1L do begin
       im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=(photflag[agnsf[ii]] eq 'Y'), color=agnsfcolor
       oploterror, x[agnsf[ii]], y[agnsf[ii]], xerr[agnsf[ii]], yerr[agnsf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    endfor

; label some points

;   for ii = 0L, nindx-1L do begin
;      if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
;         align = 1.0 & xoff = -0.14 & yoff = -0.03
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;      if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin ; = ARP299
;         align = 1.0 & xoff = -0.14 & yoff = -0.03
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;   endfor    
    
    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
; SFR[L(IR)] vs SFR[L(Ha)_cor]
; ---------------------------------------------------------------------------    

    psname = 'sfr_lir_vs_sfr_lha_cor'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.2, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.7, height=6.7, $
      xmargin=[1.4,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.2, $
      position=pos, /normal
    
    indx = where((atlas_hcn.sfr_lir gt -900.0) and (atlas_hcn.sfr_lha_cor gt -900.0),nindx)
    
    x = atlas_hcn[indx].sfr_lir
    xerr = atlas_hcn[indx].sfr_lir_err
    y = atlas_hcn[indx].sfr_lha_cor
    yerr = atlas_hcn[indx].sfr_lha_cor_err
    photflag = atlas_hcn[indx].photflag
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

; fit to the data

    dofit = where((galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
;   dofit = where((photflag eq 'Y') and (galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
    coeff = robust_linefit(x[dofit]-pivot_sfr,y[dofit],yfit,sig,coeff_err,/bisect,bisquare_limit=6.0)
    coeff[0] = coeff[0] - coeff[1]*pivot_sfr

    coeff_slope_one = [-999.0,1.0]
    coeff_slope_one[0] = interpol(poly(sfraxis,coeff),sfraxis,pivot_sfr)-coeff_slope_one[1]*pivot_sfr
    splog, psname, ndofit
    niceprint, coeff, coeff_err;, coeff_slope_one
    
; make the plot    
    
    xtitle = 'log \psi[L(IR)] (M_{\odot} yr^{-1})'
    ytitle = 'log \psi[L(H\alpha)_{cor}] (M_{\odot} yr^{-1})'

    xrange = sfrrange2
    yrange = sfrrange2

    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, charsize=2.0, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0]
    djs_oplot, sfraxis, poly(sfraxis,coeff), line=5, thick=postthick4
    djs_oplot, sfraxis, poly(sfraxis,coeff_slope_one), line=0, thick=postthick1
    for ii = 0L, nsf-1L do begin
       im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=(photflag[sf[ii]] eq 'Y'), color=sfcolor
       oploterror, x[sf[ii]], y[sf[ii]], xerr[sf[ii]], yerr[sf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(sfcolor)
    endfor
    for ii = 0L, nagn-1L do begin
       im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=(photflag[agn[ii]] eq 'Y'), color=agncolor
       oploterror, x[agn[ii]], y[agn[ii]], xerr[agn[ii]], yerr[agn[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agncolor)
    endfor
    for ii = 0L, nagnsf-1L do begin
       im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=(photflag[agnsf[ii]] eq 'Y'), color=agnsfcolor
       oploterror, x[agnsf[ii]], y[agnsf[ii]], xerr[agnsf[ii]], yerr[agnsf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    endfor

; label some points

;   for ii = 0L, nindx-1L do begin
;      if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
;         align = 1.0 & xoff = -0.14 & yoff = -0.03
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;      if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin ; = ARP299
;         align = 1.0 & xoff = -0.14 & yoff = -0.03
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;   endfor    
    
    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
; L(HCN) vs SFR[L(Ha)_cor]
; ---------------------------------------------------------------------------    

    psname = 'lhcn_vs_sfr_lha_cor'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.2, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.7, height=6.7, $
      xmargin=[1.4,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.2, $
      position=pos, /normal
    
    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.sfr_lha_cor gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.sfr_lha_cor gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn
    xerr = atlas_hcn[indx].gao_lhcn_err
    y = atlas_hcn[indx].sfr_lha_cor
    yerr = atlas_hcn[indx].sfr_lha_cor_err
    photflag = atlas_hcn[indx].photflag
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn
    xerr_ul = atlas_hcn[indx_ul].gao_lhcn_err
    y_ul = atlas_hcn[indx_ul].sfr_lha_cor
    yerr_ul = atlas_hcn[indx_ul].sfr_lha_cor_err

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; fit to the data

    dofit = where((galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
;   dofit = where((photflag eq 'Y') and (galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
    coeff = robust_linefit(x[dofit]-pivot_lhcn,y[dofit],yfit,sig,coeff_err,/bisect,bisquare_limit=6.0)
    coeff[0] = coeff[0] - coeff[1]*pivot_lhcn

    coeff_slope_one = [-999.0,1.0]
    coeff_slope_one[0] = interpol(poly(lhcnaxis,coeff),lhcnaxis,pivot_lhcn)-coeff_slope_one[1]*pivot_lhcn
    splog, psname, ndofit
    niceprint, coeff, coeff_err;, coeff_slope_one
    
; make the plot    
    
    xtitle = 'log L(HCN) (K km s^{-1} pc^{2})'
    ytitle = 'log \psi[L(H\alpha)_{cor}] (M_{\odot} yr^{-1})'

    xrange = lhcnrange
    yrange = sfrrange2

    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, charsize=2.0, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0]
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff), line=5, thick=postthick4
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff_slope_one), line=0, thick=postthick1
    for ii = 0L, nsf-1L do begin
       im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=(photflag[sf[ii]] eq 'Y'), color=sfcolor
       oploterror, x[sf[ii]], y[sf[ii]], xerr[sf[ii]], yerr[sf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(sfcolor)
    endfor
    for ii = 0L, nagn-1L do begin
       im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=(photflag[agn[ii]] eq 'Y'), color=agncolor
       oploterror, x[agn[ii]], y[agn[ii]], xerr[agn[ii]], yerr[agn[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agncolor)
    endfor
    for ii = 0L, nagnsf-1L do begin
       im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=(photflag[agnsf[ii]] eq 'Y'), color=agnsfcolor
       oploterror, x[agnsf[ii]], y[agnsf[ii]], xerr[agnsf[ii]], yerr[agnsf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    endfor

; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

; label some points

;   for ii = 0L, nindx-1L do begin
;      if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
;         align = 1.0 & xoff = -0.14 & yoff = -0.03
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;      if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin ; = ARP299
;         align = 1.0 & xoff = -0.14 & yoff = -0.03
;         xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
;           charsize=1.2, charthick=postthick2, /data, align=align
;      endif
;   endfor    
    
    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
; L(HCN) vs SFR[L(IR)]
; ---------------------------------------------------------------------------    

    psname = 'lhcn_vs_sfr_lir'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.2, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.7, height=6.7, $
      xmargin=[1.4,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.2, $
      position=pos, /normal
    
    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.sfr_lir gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.sfr_lir gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn
    xerr = atlas_hcn[indx].gao_lhcn_err
    y = atlas_hcn[indx].sfr_lir
    yerr = atlas_hcn[indx].sfr_lir_err
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn
    xerr_ul = atlas_hcn[indx_ul].gao_lhcn_err
    y_ul = atlas_hcn[indx_ul].sfr_lir
    yerr_ul = atlas_hcn[indx_ul].sfr_lir_err

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; fit to the data

    dofit = where((galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
    coeff = robust_linefit(x[dofit]-pivot_lhcn,y[dofit],yfit,sig,coeff_err,/bisect,bisquare_limit=6.0)
    coeff[0] = coeff[0] - coeff[1]*pivot_lhcn

    coeff_slope_one = [-999.0,1.0]
    coeff_slope_one[0] = interpol(poly(lhcnaxis,coeff),lhcnaxis,pivot_lhcn)-coeff_slope_one[1]*pivot_lhcn
    splog, psname, ndofit
    niceprint, coeff, coeff_err;, coeff_slope_one
    
; make the plot    
    
    xtitle = 'log L(HCN) (K km s^{-1} pc^{2})'
    ytitle = 'log \psi[L(IR)] (M_{\odot} yr^{-1})'
;   ytitle = 'log \psi[L(8-1000 \mu'+'m)] (M_{\odot} yr^{-1})'

    xrange = lhcnrange
    yrange = sfrrange1

    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, charsize=2.0, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0]
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff), line=5, thick=postthick4
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff_slope_one), line=0, thick=postthick1
    im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=sffill, color=sfcolor
    oploterror, x[sf], y[sf], xerr[sf], yerr[sf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(sfcolor)
    im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=agnfill, color=agncolor
    oploterror, x[agn], y[agn], xerr[agn], yerr[agn], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agncolor)
    im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=agnsffill, color=agnsfcolor
    oploterror, x[agnsf], y[agnsf], xerr[agnsf], yerr[agnsf], psym=8, $
      errthick=postthick4, errcolor=djs_icolor(agnsfcolor)

; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

; label some points

    for ii = 0L, nindx-1L do begin
       if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
          align = 1.0 & xoff = -0.14 & yoff = -0.03
          xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
            charsize=1.2, charthick=postthick2, /data, align=align
       endif
       if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin ; = ARP299
          align = 1.0 & xoff = -0.14 & yoff = -0.03
          xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
            charsize=1.2, charthick=postthick2, /data, align=align
       endif
    endfor    
    
    im_openclose, postscript=postscript, /close

; ---------------------------------------------------------------------------    
; L(HCN) vs SFR[L(Ha)_obs+a*L(24)]
; ---------------------------------------------------------------------------    

    psname = 'lhcn_vs_sfr_l24_lha_obs'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.2, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.7, height=6.7, $
      xmargin=[1.4,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.2, $
      position=pos, /normal
    
    indx = where((atlas_hcn.gao_lhcn gt 0.0) and (atlas_hcn.sfr_l24_lha_obs gt -900.0),nindx)
    indx_ul = where((atlas_hcn.gao_lhcn lt 0.0) and (atlas_hcn.sfr_l24_lha_obs gt -900.0),nindx_ul)
    
    x = atlas_hcn[indx].gao_lhcn
    xerr = atlas_hcn[indx].gao_lhcn_err
    y = atlas_hcn[indx].sfr_l24_lha_obs
    yerr = atlas_hcn[indx].sfr_l24_lha_obs_err
    photflag = atlas_hcn[indx].photflag
    galaxy = strtrim(atlas_hcn[indx].galaxy,2)

    x_ul = -atlas_hcn[indx_ul].gao_lhcn
    xerr_ul = atlas_hcn[indx_ul].gao_lhcn_err
    y_ul = atlas_hcn[indx_ul].sfr_l24_lha_obs
    yerr_ul = atlas_hcn[indx_ul].sfr_l24_lha_obs_err

    sf    = where(strtrim(atlas_hcn[indx].class,2) eq 'HII',nsf)
    agn   = where(strtrim(atlas_hcn[indx].class,2) eq 'AGN',nagn)
    agnsf = where(strtrim(atlas_hcn[indx].class,2) eq 'HII/AGN',nagnsf)

    sf_ul    = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII',nsf_ul)
    agn_ul   = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'AGN',nagn_ul)
    agnsf_ul = where(strtrim(atlas_hcn[indx_ul].class,2) eq 'HII/AGN',nagnsf_ul)

; fit to the data

    dofit = where((galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
;   dofit = where((photflag eq 'Y') and (galaxy ne 'NGC3690') and (galaxy ne 'NGC1614'),ndofit)
    coeff = robust_linefit(x[dofit]-pivot_lhcn,y[dofit],yfit,sig,coeff_err,/bisect,bisquare_limit=6.0)
    coeff[0] = coeff[0] - coeff[1]*pivot_lhcn

    coeff_slope_one = [-999.0,1.0]
    coeff_slope_one[0] = interpol(poly(lhcnaxis,coeff),lhcnaxis,pivot_lhcn)-coeff_slope_one[1]*pivot_lhcn
    splog, psname, ndofit
    niceprint, coeff, coeff_err;, coeff_slope_one
    
; make the plot    
    
    xtitle = 'log L(HCN) (K km s^{-1} pc^{2})'
    ytitle = 'log \psi[L(H\alpha)_{obs}+a*L(24 \mu'+'m)] (M_{\odot} yr^{-1})'

    xrange = lhcnrange
    yrange = sfrrange1

    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, charsize=2.0, charthick=postthick2, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0]
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff), line=5, thick=postthick4
    djs_oplot, lhcnaxis, poly(lhcnaxis,coeff_slope_one), line=0, thick=postthick1
    for ii = 0L, nsf-1L do begin
       im_symbols, sfsym, thick=postthick4, psize=sfpsize, fill=(photflag[sf[ii]] eq 'Y'), color=sfcolor
       oploterror, x[sf[ii]], y[sf[ii]], xerr[sf[ii]], yerr[sf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(sfcolor)
    endfor
    for ii = 0L, nagn-1L do begin
       im_symbols, agnsym, thick=postthick4, psize=agnpsize, fill=(photflag[agn[ii]] eq 'Y'), color=agncolor
       oploterror, x[agn[ii]], y[agn[ii]], xerr[agn[ii]], yerr[agn[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agncolor)
    endfor
    for ii = 0L, nagnsf-1L do begin
       im_symbols, agnsfsym, thick=postthick4, psize=agnsfpsize, fill=(photflag[agnsf[ii]] eq 'Y'), color=agnsfcolor
       oploterror, x[agnsf[ii]], y[agnsf[ii]], xerr[agnsf[ii]], yerr[agnsf[ii]], psym=8, $
         errthick=postthick4, errcolor=djs_icolor(agnsfcolor)
    endfor

; upper limits    
    
    if (nsf_ul ne 0L) then begin
       im_symbols, sfsym_ul, thick=postthick1, psize=sfpsize_ul, fill=sffill_ul, color=sfcolor_ul
       oploterror, x_ul[sf_ul], y_ul[sf_ul], xerr_ul[sf_ul], yerr_ul[sf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(sfcolor_ul)
    endif
    if (nagn_ul ne 0L) then begin
       im_symbols, agnsym_ul, thick=postthick1, psize=agnpsize_ul, fill=agnfill_ul, color=agncolor_ul
       oploterror, x_ul[agn_ul], y_ul[agn_ul], xerr_ul[agn_ul], yerr_ul[agn_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agncolor_ul)
    endif
    if (nagnsf_ul ne 0L) then begin
       im_symbols, agnsfsym_ul, thick=postthick1, psize=agnsfpsize_ul, fill=agnsffill_ul, color=agnsfcolor_ul
       oploterror, x_ul[agnsf_ul], y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], $
         psym=8, errthick=postthick4, errcolor=djs_icolor(agnsfcolor_ul)
    endif

; label some points

    for ii = 0L, nindx-1L do begin
       if strmatch(galaxy[ii],'*ngc1614*',/fold) then begin
          align = 1.0 & xoff = -0.14 & yoff = -0.03
          xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
            charsize=1.2, charthick=postthick2, /data, align=align
       endif
       if strmatch(galaxy[ii],'*ngc3690*',/fold) then begin ; = ARP299
          align = 1.0 & xoff = -0.14 & yoff = -0.03
          xyouts, x[ii]+xoff, y[ii]+yoff, galaxy[ii], $
            charsize=1.2, charthick=postthick2, /data, align=align
       endif
    endfor    
    
    im_openclose, postscript=postscript, /close

; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    

    if keyword_set(postscript) and keyword_set(encapsulated) then begin
       im_ps2html, htmlbase, html_path=html_path, cleanpng=0, npscols=3, _extra=extra
    endif

return
end
