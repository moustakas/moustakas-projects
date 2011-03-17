pro ndwfs_flamex_kband
; jm09may26nyu - compare the NDWFS and FLAMEX photometry against 2MASS

    common ndwfs_kband, tmass1, flamj1, flamk1, kband1

    catalogs_path = ages_path(/catalogs)
    analysis_path = ages_path(/analysis)
    ndwfsdir = getenv('RESEARCHPATH')+'/data/ndwfs/'

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

; NDWFS/K    
    if (n_elements(kband1) eq 0L) then begin
       photofile = file_search(ndwfsdir+'NDWFS_K_??_??_cat_m.fits.gz',count=nfile)
       for ii = 0, nfile-1 do begin
          splog, 'Reading '+photofile[ii]
          temp = mrdfits(photofile[ii],1)
          keep = where((temp.flags eq 0),nkeep)
          if (nkeep ne 0) then begin
             splog, nkeep
             temp = temp[keep]
             if (n_elements(kband1) eq 0) then $
               kband1 = temporary(temp) else $
               kband1 = [temporary(kband1),temporary(temp)]
          endif
       endfor
    endif
;   if (n_elements(kband1) eq 0L) then begin
;      splog, 'Reading '+analysis_path+'catalog.ndwfsk.fits.gz'
;      kband1 = mrdfits(analysis_path+'catalog.ndwfsk.fits.gz',1,/silent)
;   endif

; make some plots    
    psfile = ndwfsdir+'ndwfs_flamex_kband.ps'
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
    spherematch, tmass1.ra, tmass1.decl, kband1.alpha_j2000, $
      kband1.delta_j2000, 1.5/3600.0, m1, m2
    tmass = tmass1[m1]
    kband = kband1[m2]

    good = where((tmass.k_m_ext gt 0.0) and (kband.mag_auto gt 0.0) and $
      (kband.mag_auto lt 80.0) and (kband.magerr_auto lt 1.0),ngood)
    splog, 'NDWFS vs 2MASS ', ngood
    tmass = tmass[good]
    kband = kband[good]

    xx = tmass.k_m_ext
    yy = kband.mag_auto-kmks
    xxerr = tmass.k_msig_ext
    yyerr = sqrt(tmass.k_msig_ext^2.0+kband.magerr_auto^2.0)
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
      kband1.alpha_j2000, kband1.delta_j2000, 1.5/3600.0, m1, m2
    flamk = flamk1[m1]
    kband = kband1[m2]

    good = where((flamk.mag_auto gt 0.0) and (flamk.mag_auto lt 80.0) and $
      (flamk.detprob gt 0.95) and (flamk.sim_magerr_auto lt 1.0) and $
      (kband.mag_auto gt 0.0) and (kband.mag_auto lt 80.0) and $
      (kband.magerr_auto lt 1.0),ngood)
    splog, 'NDWFS vs FLAMEX ', ngood
    flamk = flamk[good]
    kband = kband[good]

    xx = flamk.mag_auto
    yy = kband.mag_auto-kmks
    xxerr = flamk.sim_magerr_auto
    yyerr = sqrt(flamk.sim_magerr_auto^2+kband.magerr_auto^2.0)
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

    im_plotconfig, psfile=psfile, /psclose, /gzip

return
end
    
