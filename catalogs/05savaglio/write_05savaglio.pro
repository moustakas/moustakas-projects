pro write_05savaglio, outdata, debug=debug
; jm04dec02uofa

    snrcut = 0.0
    Bmsun = k_solar_magnitudes(filterlist='bessell_B.par')
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par',/kurucz))[0]

    root = '05savaglio'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = rsex(path+root+'.dat')
    ngalaxy = n_elements(raw)

    data = init_cat_linefit(ngalaxy=ngalaxy)
    moretags = replicate({m_b_glaze: -999.0, mass_glaze: -999.0, $
      mass_glaze_err: -999.0},ngalaxy)
    data = struct_addtags(data,moretags)

; K-corrections and stellar masses; only trust the results if there is
; photometry in three or more bandpasses; convert to Salpeter

    kcorr = mrdfits(path+'05savaglio_kcorr.fits.gz',1,/silent)
    match, kcorr.galaxy, raw.galaxy, indx1, indx2
    data = struct_addtags(struct_trimtags(data,except='MASS'),$
      struct_trimtags(kcorr[indx1],except=['GALAXY','Z']))

    nband = total((data.abmaggies gt 0.0),1)
    good = where((data.mass gt -900.0) and (nband ge 3.0))
    data[good].mass = data[good].mass + im_convert_imf(/from_chabrier)
    data[good].m_g = data[good].ugriz_absmag[1]
    data[good].m_r = data[good].ugriz_absmag[2]
    data[good].m_b = data[good].ubvri_absmag[1]

; ---------------------------------------------------------------------------    
; store relevant galaxy properties
; ---------------------------------------------------------------------------    

    data.id = lindgen(ngalaxy)+1L
    data.galaxy = raw.galaxy
    data.z = raw.z

    good = where(raw.mass gt -900.0)
    data[good].mass_glaze = raw[good].mass + im_convert_imf(/from_baldry03) ; Salpeter
    data[good].mass_glaze_err = raw[good].mass_err

    good = where(raw.m_b_ab gt -900.0)
    data[good].m_b_glaze = raw[good].m_b_ab - Bvega2ab ; AB-->Vega

; ---------------------------------------------------------------------------
; fluxes    
; ---------------------------------------------------------------------------

    good = where(raw.oii gt -900,ngood)
    if (ngood ne 0L) then data[good].oii_3727 = 1D-18*transpose([ [raw[good].oii], [raw[good].oii_err] ])

    good = where(raw.oiii gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].oiii_5007 = 1D-18*transpose([ [raw[good].oiii], [raw[good].oiii_err] ])
       data[good].oiii_4959 = data[good].oiii_5007 / 3.0
    endif

; predict H-alpha given H-beta and assuming A_V=2.1 [E(B-V)=0.67]
    
    good = where(raw.hb gt -900,ngood) ; corrected for stellar absorption
    if (ngood ne 0L) then begin

       data[good].h_beta = 1D-18*transpose([ [raw[good].hb], [raw[good].hb_err] ])

       ebv = fltarr(ngood)+2.1/k_lambda(5500.0)
       factor = 2.86*10^(0.4*ebv*(k_lambda(4861.0,/odonnell)-k_lambda(6563.0,/odonnell)))
       data[good].h_alpha[0] = data[good].h_beta[0]*factor
       data[good].h_alpha[1] = data[good].h_beta[1]*factor

    endif
    
    good = where(raw.hg gt -900,ngood)
    if (ngood ne 0L) then data[good].h_gamma = 1D-18*transpose([ [raw[good].hg], [raw[good].hg_err] ])

; ---------------------------------------------------------------------------    
; assign R23 branches
; ---------------------------------------------------------------------------    

    splog, 'Assigning R23 branches'
    data.r23branch = 'U'
    
; ---------------------------------------------------------------------------    
; compute the reddening, metallicity, and other miscellaneous
; properties, then write out 
; ---------------------------------------------------------------------------    

    idata = iunred_linedust(data,snrcut=snrcut,/silent)
    ages_mz_log12oh, data, idata, nmonte=500L, snrcut=0.0, $
      logfile=path+'mz_log12oh.05savaglio.log', $
      final_ohdust=log12oh, final_ohnodust=log12oh_nodust

    outdata = struct_addtags(data,log12oh)
    outdata_nodust = struct_addtags(idata,log12oh_nodust)
    
    splog, 'Writing '+path+root+'.fits'
    mwrfits, outdata, path+root+'.fits', /create
    splog, 'Writing '+path+root+'_nodust.fits'
    mwrfits, outdata_nodust, path+root+'_nodust.fits', /create

; ---------------------------------------------------------------------------    
; write out statistics
; ---------------------------------------------------------------------------    

    indx = where((outdata.ubvri_absmag[1] gt -900.0) and (outdata.zstrong_12oh_kk04 gt -900.0),nindx)
    
    print, '##################################################'
    splog, 'Savaglio et al. (2005): '+string(nindx,format='(I0.0)')+' galaxies.'
    stats = im_stats(outdata[indx].z)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<z<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <z>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    stats = im_stats(outdata[indx].m_b)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<M_B<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <M_B>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    print, '##################################################'
    struct_print, struct_trimtags(outdata,select=['GALAXY','M_B','ZSTRONG_12OH_KK04*','R23BRANCH_KK04'])
    print, '##################################################'
       
; ---------------------------------------------------------------------------    
; make some plots
; ---------------------------------------------------------------------------    

; compare the masses and magnitudes; pretty good agreement in the
; masses, but M_B for SA12-5685 disagrees by one magnitude (Savaglio's
; is bigger); not sure why

    if keyword_set(debug) then begin
       
       im_window, 0, xratio=0.4, /square
       niceprint, raw.galaxy, data.galaxy, data.mass_glaze, data.mass
       ploterror, data.mass, data.mass_glaze, data.mass_glaze_err, xrange=[8,11.5], yrange=[8,11.5], $
         ps=4, xsty=3, ysty=3, sym=2
       oplot, !x.crange, !y.crange, thick=2.0
       jj = im_stats(data.mass_glaze-data.mass,/verbose)
       cc = get_kbrd(1)
       
       niceprint, data.galaxy, data.m_b_glaze, data.ubvri_absmag[1] ;, mass.m_b_err
       plot, data.m_b_glaze, data.ubvri_absmag[1], xrange=[-16.5,-22], yrange=[-16.5,-22], $
         ps=4, xsty=3, ysty=3, sym=2
       oplot, !x.crange, !y.crange, thick=2.0
       jj = im_stats(data.m_b_glaze-data.ubvri_absmag[1],/verbose)
       cc = get_kbrd(1)

       niceprint, raw.galaxy, raw.oh, outdata.zstrong_12oh_kk04, outdata.r23branch_kk04
       plot, raw.oh, outdata.zstrong_12oh_kk04, xrange=[8.0,9.2], yrange=[8.0,9.2], $
         ps=4, xsty=3, ysty=3, sym=2
       oplot, !x.crange, !y.crange, thick=2.0
       jj = im_stats(raw.oh-outdata.zstrong_12oh_kk04,/verbose)
       cc = get_kbrd(1)

       niceprint, raw.galaxy, raw.oh, outdata_nodust.zstrong_12oh_kk04, outdata_nodust.r23branch_kk04
       plot, raw.oh, outdata_nodust.zstrong_12oh_kk04, xrange=[8.0,9.2], yrange=[8.0,9.2], $
         ps=4, xsty=3, ysty=3, sym=2
       oplot, !x.crange, !y.crange, thick=2.0
       jj = im_stats(raw.oh-outdata_nodust.zstrong_12oh_kk04,/verbose)
       cc = get_kbrd(1)

       g = where(raw.mass gt -900.0 and raw.oh gt -900.0)
       ploterror, raw[g].mass, raw[g].oh, raw[g].mass_err, raw[g].oh_err, ps=4, $
         xsty=1, ysty=1, xrange=[8,11], yrange=[8.1,9.2]
       cc = get_kbrd(1)
       
       g = where(outdata.mass gt -900.0 and outdata.zstrong_12oh_kk04 gt -900.0)
       ploterror, outdata[g].mass, outdata[g].zstrong_12oh_kk04, raw[g].mass_err, $
         outdata[g].zstrong_12oh_kk04_err, ps=4, xsty=1, ysty=1, xrange=[8,11], yrange=[8.1,9.2]
       cc = get_kbrd(1)
       
; for one object the metallicity differs significantly, SA22-0997
       
       niceprint, lindgen(ngalaxy), outdata.galaxy, raw.oh, outdata.zstrong_12oh_kk04

    endif
       
return
end    
