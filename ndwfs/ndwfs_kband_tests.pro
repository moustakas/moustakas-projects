pro ndwfs_kband_tests
; jm09sep14ucsd - test various aspects of the NDWFS K-band data (based
;   on earlier code); compare the NDWFS and FLAMEX photometry against
;   2MASS 

    common ndwfs_kband, ndwfs1, tmass1, flamj1, flamk1
    
    ndwfsdir = getenv('RESEARCHPATH')+'/data/ndwfs/'
    analysis_path = ages_path(/analysis)
    catalogs_path = ages_path(/catalogs)

; read the full set of NDWFS catalogs    
    if (n_elements(ndwfs1) eq 0) then begin
       photofile = file_search(ndwfsdir+'NDWFS_??_??.fits.gz',count=nfile)
       for ii = 0, nfile-1 do begin
          temp = mrdfits(photofile[ii],1)
          keep = where((temp.k_mag_auto gt 0.0) and (temp.k_mag_auto lt 90.0),nkeep)
;         keep = where((temp.bw_flags eq 0) and (temp.r_flags eq 0) and $
;           (temp.i_flags eq 0) and (temp.k_flags eq 0),nkeep)
          splog, 'NKEEP = ', nkeep
          if (nkeep gt 0) then begin
             if (n_elements(ndwfs1) eq 0) then $
               ndwfs1 = temporary(temp[keep]) else $
               ndwfs1 = [temporary(ndwfs1),temporary(temp[keep])]
          endif
       endfor
    endif

; 2MASS    
;   if (n_elements(tmass1) eq 0) then tmass1 = $
;     hogg_mrdfits(vagc_name('object_twomass'),1,nrow=28800L,$
;     columns=['ra','decl','twomass_tag','k_m_ext','k_msig_ext',$
;     'h_m_ext','h_msig_ext','j_m_ext','j_msig_ext'])
    if (n_elements(tmass1) eq 0) then begin
;      % setup twomass_xsc
;      tmassdir = '/global/data/products/NULL/twomass_xsc/v0_13/'
       tmassdir = getenv('TWOMASS_DIR')+'/'
       columns=['ra','decl','twomass_tag','k_m_ext','k_msig_ext',$
         'h_m_ext','h_msig_ext','j_m_ext','j_msig_ext']
       tmassa = hogg_mrdfits(tmassdir+'xsc_aaa.fits.gz',1, $
         nrowchunk=50000,columns=columns)
       tmassb = hogg_mrdfits(tmassdir+'xsc_baa.fits.gz',1, $
         nrowchunk=50000,columns=columns)
       nna = n_elements(tmassa)
       nnb = n_elements(tmassb)
       tmass1 = im_empty_structure(tmassa[0],ncopies=nna+nnb)
       tmass1[0:nna-1] = temporary(tmassa)
       tmass1[nna:nna+nnb-1] = temporary(tmassb)
    endif

; FLAMEX
    if (n_elements(flamj1) eq 0L) then begin
       splog, 'Reading '+catalogs_path+'FLAMEX/BOOTES_j_V1.0.cat'
       flamj1 = rsex(catalogs_path+'FLAMEX/BOOTES_j_V1.0.cat')
    endif
    if (n_elements(flamk1) eq 0L) then begin
       splog, 'Reading '+catalogs_path+'FLAMEX/BOOTES_ks_V1.0.cat'
       flamk1 = rsex(catalogs_path+'FLAMEX/BOOTES_ks_V1.0.cat')
    endif
    
; ###########################################################################
; compare the total K-band constructed in various ways

    ndwfs = parse_ndwfs_phot(ndwfs1,/allndwfs,/nozpoffset)
    keep = where((ndwfs.k_flags eq 0) and (ndwfs.kmag_auto gt 0.0) and $
      (ndwfs.ktot_4 gt 0.0) and (ndwfs.ktot_6 gt 0.0) and $
      (ndwfs.ktot_r_4 gt 0.0) and (ndwfs.ktot_r_6 gt 0.0))
    ndwfs = ndwfs[keep]
    
; make the plot    
    psfile = ndwfsdir+'ndwfs_kband_total.ps'
    im_plotconfig, 6, pos1, psfile=psfile, xmargin=[1.4,0.3], $
      height=[4.5,3.0]

    mrange1 = [11.5,20.5]
    rrange = 0.95*[-1,1]
    
; K_auto vs K_4 from the I-band
    xx = ndwfs.kmag_auto
    yy = ndwfs.ktot_4
    hogg_scatterplot, xx, yy, position=pos1[*,0], /outliers, outcolor='black', $
;   djs_plot, xx, yy, position=pos1[*,0], $
      xsty=1, ysty=1, xrange=mrange1, yrange=mrange1, psym=3, $
      xtitle='', ytitle=textoidl('I_{auto, cor} + (K_{4"}-I_{4"}) (Vega mag)'), $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
    im_legend, '<\Delta'+'K> = '+im_string_stats(yy-xx,sigrej=3.0), $
      /left, /top, box=0
    hogg_scatterplot, xx, yy-xx, position=pos1[*,1], /noerase, /outliers, outcolor='black', $
;   djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, $
      xsty=1, ysty=1, xrange=mrange1, $
      yrange=rrange, psym=3, xtitle=textoidl('K_{auto} (Vega mag)'), $
      ytitle=textoidl('\Delta'+'K (Vega mag)')
    djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
    !p.multi = 0
    
; K_auto vs K_6 from the I-band
    xx = ndwfs.kmag_auto
    yy = ndwfs.ktot_6
    hogg_scatterplot, xx, yy, position=pos1[*,0], /outliers, outcolor='black', $
;   djs_plot, xx, yy, position=pos1[*,0], $
      xsty=1, ysty=1, xrange=mrange1, yrange=mrange1, psym=3, $
      xtitle='', ytitle=textoidl('I_{auto, cor} + (K_{6"}-I_{6"}) (Vega mag)'), $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
    im_legend, '<\Delta'+'K> = '+im_string_stats(yy-xx,sigrej=3.0), $
      /left, /top, box=0

    hogg_scatterplot, xx, yy-xx, position=pos1[*,1], /noerase, /outliers, outcolor='black', $
;   djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, $
      xsty=1, ysty=1, xrange=mrange1, $
      yrange=rrange, psym=3, xtitle=textoidl('K_{auto} (Vega mag)'), $
      ytitle=textoidl('\Delta'+'K (Vega mag)')
    djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
    !p.multi = 0

; K_auto vs K_4 from the R-band
    xx = ndwfs.kmag_auto
    yy = ndwfs.ktot_r_4
    hogg_scatterplot, xx, yy, position=pos1[*,0], /outliers, outcolor='black', $
;   djs_plot, xx, yy, position=pos1[*,0], $
      xsty=1, ysty=1, xrange=mrange1, yrange=mrange1, psym=3, $
      xtitle='', ytitle=textoidl('R_{auto} + (K_{4"}-R_{4"}) (Vega mag)'), $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
    im_legend, '<\Delta'+'K> = '+im_string_stats(yy-xx,sigrej=3.0), $
      /left, /top, box=0
    hogg_scatterplot, xx, yy-xx, position=pos1[*,1], /noerase, /outliers, outcolor='black', $
;   djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, $
      xsty=1, ysty=1, xrange=mrange1, $
      yrange=rrange, psym=3, xtitle=textoidl('K_{auto} (Vega mag)'), $
      ytitle=textoidl('\Delta'+'K (Vega mag)')
    djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
    !p.multi = 0

; K_auto vs K_6 from the R-band
    xx = ndwfs.kmag_auto
    yy = ndwfs.ktot_r_6
    hogg_scatterplot, xx, yy, position=pos1[*,0], /outliers, outcolor='black', $
;   djs_plot, xx, yy, position=pos1[*,0], $
      xsty=1, ysty=1, xrange=mrange1, yrange=mrange1, psym=3, $
      xtitle='', ytitle=textoidl('R_{auto} + (K_{6"}-R_{6"}) (Vega mag)'), $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
    im_legend, '<\Delta'+'K> = '+im_string_stats(yy-xx,sigrej=3.0), $
      /left, /top, box=0
    hogg_scatterplot, xx, yy-xx, position=pos1[*,1], /noerase, /outliers, outcolor='black', $
;   djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, $
      xsty=1, ysty=1, xrange=mrange1, $
      yrange=rrange, psym=3, xtitle=textoidl('K_{auto} (Vega mag)'), $
      ytitle=textoidl('\Delta'+'K (Vega mag)')
    djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'

    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

; ###########################################################################
; compare the NDWFS and FLAMEX zeropoints against 2MASS
    
    psfile = ndwfsdir+'ndwfs_kband_zeropoint.ps'
    im_plotconfig, 6, pos1, psfile=psfile, xmargin=[1.4,0.3], $
      height=[4.5,3.0]

    mrange1 = [11.0,15.5]
    mrange2 = [12.0,16.5]
    mrange3 = [11.0,19.5]
    rrange = 0.95*[-1,1]

    kmks = +0.10
    strkmks = string(kmks,format='(F4.2)')

; FLAMEX vs 2MASS
    spherematch, tmass1.ra, tmass1.decl, flamk1.alphapeak_j2000, $
      flamk1.deltapeak_j2000, 1.5/3600.0, m1, m2
    tmass = tmass1[m1]
    flamk = flamk1[m2]
    flamj = flamj1[m2]

    good = where((tmass.k_m_ext gt 0.0) and $
      (flamj.mag_auto gt 0.0) and (flamj.mag_auto lt 80.0) and $
      (flamj.detprob gt 0.95) and (flamj.sim_magerr_auto lt 1.0) and $
      (flamk.mag_auto gt 0.0) and (flamk.mag_auto lt 80.0) and $
      (flamk.detprob gt 0.95) and (flamk.sim_magerr_auto lt 1.0),ngood)
    splog, 'FLAMEX vs 2MASS ', ngood
    tmass = tmass[good]
    flamk = flamk[good]
    flamj = flamj[good]

; J-band
    xx = tmass.j_m_ext
    yy = flamj.mag_auto
    xxerr = tmass.j_msig_ext
    yyerr = sqrt(tmass.j_msig_ext^2.0+flamj.sim_magerr_auto^2)
    ploterror, xx, yy, xxerr, yyerr, position=pos1[*,0], $
;   djs_plot, xx, yy, position=pos1[*,0], $
      xsty=1, ysty=1, xrange=mrange2, yrange=mrange2, psym=4, $
      xtitle='', ytitle=textoidl('J_{auto} (FLAMEX, Vega mag)'), xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
    im_legend, '<\Delta'+'J> = '+im_string_stats(yy-xx,sigrej=3.0), $
      /left, /top, box=0
    djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, xsty=1, ysty=1, xrange=mrange2, $
      yrange=rrange, psym=4, xtitle='J_{ext} (2MASS, Vega mag)', $
      ytitle='\Delta'+'J (Vega mag)'
    djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
;   djs_oplot, !x.crange, +yyerr*[1,1], line=5, thick=3, color='red'
;   djs_oplot, !x.crange, -yyerr*[1,1], line=5, thick=3, color='red'

; K-band
    xx = tmass.k_m_ext
    yy = flamk.mag_auto
    xxerr = tmass.k_msig_ext
    yyerr = sqrt(tmass.k_msig_ext^2.0+flamk.sim_magerr_auto^2)
    ploterror, xx, yy, xxerr, yyerr, position=pos1[*,0], $
;   djs_plot, xx, yy, position=pos1[*,0], $
      xsty=1, ysty=1, xrange=mrange1, yrange=mrange1, psym=4, $
      xtitle='', ytitle=textoidl('K_{s, auto} (FLAMEX, Vega mag)'), xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
    im_legend, '<\Delta'+'K_{s}> = '+im_string_stats(yy-xx,sigrej=3.0), $
      /left, /top, box=0
    djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, xsty=1, ysty=1, xrange=mrange1, $
      yrange=rrange, psym=4, xtitle='K_{s, ext} (2MASS, Vega mag)', $
      ytitle='\Delta'+'K_{s} (Vega mag)'
    djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
;   djs_oplot, !x.crange, +yyerr*[1,1], line=5, thick=3, color='red'
;   djs_oplot, !x.crange, -yyerr*[1,1], line=5, thick=3, color='red'

; NDWFS vs 2MASS
    spherematch, tmass1.ra, tmass1.decl, ndwfs1.i_alpha_j2000, $
      ndwfs1.i_delta_j2000, 1.5/3600.0, m1, m2
    tmass = tmass1[m1]
    ndwfs = ndwfs1[m2]

    good = where((tmass.k_m_ext gt 0.0) and (ndwfs.k_mag_auto gt 0.0) and $
      (ndwfs.k_mag_auto lt 80.0) and (ndwfs.k_magerr_auto lt 1.0),ngood)
    splog, 'NDWFS vs 2MASS ', ngood
    tmass = tmass[good]
    ndwfs = ndwfs[good]

    xx = tmass.k_m_ext
    yy = ndwfs.k_mag_auto-kmks
    xxerr = tmass.k_msig_ext
    yyerr = sqrt(tmass.k_msig_ext^2.0+ndwfs.k_magerr_auto^2.0)
    ploterror, xx, yy, xxerr, yyerr, position=pos1[*,0], $
;   djs_plot, xx, yy, position=pos1[*,0], $
      xsty=1, ysty=1, xrange=mrange1, yrange=mrange1, psym=4, $
      xtitle='', ytitle=textoidl('K_{auto}-'+strkmks+' (NDWFS, Vega mag)'), xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
    im_legend, '<\Delta'+'K_{s}> = '+im_string_stats(yy-xx,sigrej=3.0), $
      /left, /top, box=0
    im_legend, '<K-K_{s}>_{Vega}=+'+strkmks+' for galaxies', /right, /bottom, box=0, margin=0
    djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, xsty=1, ysty=1, xrange=mrange1, $
      yrange=rrange, psym=4, xtitle='K_{s, ext} (2MASS, Vega mag)', $
      ytitle='\Delta'+'K_{s} (Vega mag)'
    djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
;   djs_oplot, !x.crange, +yyerr*[1,1], line=5, thick=3, color='red'
;   djs_oplot, !x.crange, -yyerr*[1,1], line=5, thick=3, color='red'

; NDWFS vs FLAMEX
    spherematch, flamk1.alphapeak_j2000, flamk1.deltapeak_j2000, $
      ndwfs1.i_alpha_j2000, ndwfs1.i_delta_j2000, 1.5/3600.0, m1, m2
    flamk = flamk1[m1]
    ndwfs = ndwfs1[m2]

    good = where((flamk.mag_auto gt 0.0) and (flamk.mag_auto lt 80.0) and $
      (flamk.detprob gt 0.95) and (flamk.sim_magerr_auto lt 1.0) and $
      (ndwfs.k_mag_auto gt 0.0) and (ndwfs.k_mag_auto lt 80.0) and $
      (ndwfs.k_magerr_auto lt 1.0),ngood)
    splog, 'NDWFS vs FLAMEX ', ngood
    flamk = flamk[good]
    ndwfs = ndwfs[good]

    xx = flamk.mag_auto
    yy = ndwfs.k_mag_auto-kmks
    xxerr = flamk.sim_magerr_auto
    yyerr = sqrt(flamk.sim_magerr_auto^2+ndwfs.k_magerr_auto^2.0)
;   ploterror, xx, yy, xxerr, yyerr, position=pos1[*,0], $
    djs_plot, xx, yy, position=pos1[*,0], symsize=0.12, $
      xsty=1, ysty=1, xrange=mrange3, yrange=mrange3, psym=4, $
      xtitle='', ytitle=textoidl('K_{auto}-'+strkmks+' (NDWFS, Vega mag)'), xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
    range = where(xx lt 16.5) ; NOTE!
    im_legend, '<\Delta'+'K_{s}> = '+im_string_stats(yy[range]-xx[range],sigrej=3.0), $
      /left, /top, box=0
    im_legend, '<K-K_{s}>_{Vega}=+'+strkmks+' for galaxies', /right, /bottom, box=0, margin=0

    djs_plot, xx, yy-xx, position=pos1[*,1], /noerase, xsty=1, ysty=1, xrange=mrange3, $
      yrange=rrange, psym=4, xtitle='K_{s, auto} (FLAMEX, Vega mag)', $
      ytitle='\Delta'+'K_{s} (Vega mag)', symsize=0.12
    djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
;   djs_oplot, !x.crange, +yyerr*[1,1], line=5, thick=3, color='red'
;   djs_oplot, !x.crange, -yyerr*[1,1], line=5, thick=3, color='red'

    im_plotconfig, psfile=psfile, /psclose, /pdf, /pskeep

return
end    
