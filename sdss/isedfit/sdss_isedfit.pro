pro sdss_isedfit, models=models, isedfit=isedfit, qaplot=qaplot, $
  measure=measure, doplots=doplots, clobber=clobber
; jm09feb22nyu - written
; jm09may26nyu - various rewrites    

    common sdss_isedfit, post1, phot1, tmass1
    
    iopath = sdss_path(/isedfit)
    paramfile = iopath+'sdss_isedfit.par'

; ---------------------------------------------------------------------------
; build the model grids    
    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=iopath, /noigm

; ---------------------------------------------------------------------------
; build the sdss sample    
    if (n_elements(phot1) eq 0) then phot1 = $
      hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrow=28800L)

    if (n_elements(tmass1) eq 0) then tmass1 = $
      hogg_mrdfits(vagc_name('object_twomass'),1,nrow=28800L,$
      columns=['ra','decl','twomass_tag','k_m_ext','k_msig_ext',$
      'h_m_ext','h_msig_ext','j_m_ext','j_msig_ext'])

    if (n_elements(post1) eq 0) then post1 = $
      mrdfits(vagc_name('post_catalog',sample='dr72',$
      letter='bsafe',post='25'),1)

    post = post1
    phot = phot1[post.object_position]
    tmass = tmass1[post.object_position]

    keep = where(tmass.k_m_ext gt 0.0)
    post = post[keep]
    phot = phot[keep]
    tmass = tmass[keep]
    
    ngal = 5000
    myseed = 15.5
    keep = long(randomu(myseed,ngal)*n_elements(post))
    post = post[keep]
    phot = phot[keep]
    tmass = tmass[keep]

; ---------------------------------------------------------------------------
; do the SED fitting
    if keyword_set(isedfit) then begin
       sdss_to_maggies, maggies1, ivarmaggies1, calibobj=phot, flux='model'
       twomass_to_maggies, tmass, maggies2, ivarmaggies2
       maggies = fltarr(8,ngal) & ivarmaggies = fltarr(8,ngal)
       maggies[0:4,*] = maggies1 ; ugriz
       maggies[5:7,*] = maggies2 ; JHKs
       ivarmaggies[0:4,*] = ivarmaggies1 ; ugriz
       ivarmaggies[5:7,*] = ivarmaggies2 ; JHKs
; ugriz
       outprefix = 'ugriz'
       use_ivarmaggies = ivarmaggies
       use_ivarmaggies[5:7,*] = 0.0
       isedfit, paramfile, maggies, use_ivarmaggies, post.z, result, $
         iopath=iopath, outprefix=outprefix, nminphot=nminphot, $
         clobber=clobber, debug=0
; ugrizJHKs
       outprefix = 'ugrizJHKs'
       use_ivarmaggies = ivarmaggies
       isedfit, paramfile, maggies, use_ivarmaggies, post.z, result, $
         iopath=iopath, outprefix=outprefix, nminphot=nminphot, $
         clobber=clobber, debug=0
; ugrizKs
       outprefix = 'ugrizKs'
       use_ivarmaggies = ivarmaggies
       use_ivarmaggies[5:6,*] = 0.0
       isedfit, paramfile, maggies, use_ivarmaggies, post.z, result, $
         iopath=iopath, outprefix=outprefix, nminphot=nminphot, $
         clobber=clobber, debug=0
    endif

; build a QAplot
    if keyword_set(qaplot) then begin
; ugriz
       outprefix = 'ugriz'
       isedfit_qaplot, paramfile, iopath=iopath, $
         outprefix=outprefix, clobber=clobber
; ugrizKs
       outprefix = 'ugrizKs'
       isedfit_qaplot, paramfile, iopath=iopath, $
         outprefix=outprefix, clobber=clobber
; ugrizJHKs
       outprefix = 'ugrizJHKs'
       isedfit_qaplot, paramfile, iopath=iopath, $
         outprefix=outprefix, clobber=clobber
    endif

    r1 = mrdfits(iopath+'ugriz_chab_calzetti_sfhgrid01.fits.gz',1)
    r2 = mrdfits(iopath+'ugrizKs_chab_calzetti_sfhgrid01.fits.gz',1)
    r3 = mrdfits(iopath+'ugrizJHKs_chab_calzetti_sfhgrid01.fits.gz',1)

    if keyword_set(doplots) then begin
       
       params = read_isedfit_paramfile(paramfile,iopath=iopath)
;      fp = isedfit_filepaths(params,outprefix=outprefix,iopath=iopath)

       splog, 'Reading '+fp.iopath+fp.isedfit_outfile
       sfh1 = mrdfits(iopath+'ugriz_salp_sfhgrid01_isedfit.fits.gz',1)
       sfh2 = mrdfits(iopath+'ugriz_salp_sfhgrid02_isedfit.fits.gz',1)
       sfh3 = mrdfits(iopath+'ugriz_salp_sfhgrid03_isedfit.fits.gz',1)

       psfile = iopath+'isedfit_vs_kcorrect.ps'
       im_plotconfig, 0, pos, psfile=psfile
       
       massrange = [8.8,12.5]
       agerange = [0.5,13]
       ebvrange = [-0.05,0.7]
       residrange = [-0.6,0.6]
       plotsym, 0, 1, /fill

       kkmass = kcorr[these].mass+im_convert_imf(/from_chabrier)
; ---------------
; these plots show that SFHGRID-2 reproduces the K-correct results
; rather well; the residuals correlate strongly with median AGE
       mymass = sfh2.mass50
       mymass_err = sfh2.mass_err
       myage = sfh2.age50
       thistitle = textoidl('SDSS - Z grid, no dust')

       hogg_scatterplot, kkmass, mymass, xsty=1, ysty=1, xrange=massrange, $
         yrange=massrange, xtitle=textoidl('log (M/M_{\odot}) [K-correct]'), $
         ytitle=textoidl('log (M/M_{\odot}) [isedfit]'), title=thistitle, $
         xnpix=20.0, ynpix=20.0, /outliers, outsymsize=0.4, $
         outcolor=djs_icolor('blue')
;      djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=massrange, $
;        yrange=massrange, xtitle='log (M/M_{\odot}) [K-correct]', $
;        ytitle='log (M/M_{\odot}) [isedfit]', title=thistitle
;      oploterror, kkmass, mymass, mymass_err, psym=8
       djs_oplot, !x.crange, !y.crange, line=0, thick=2, color='red'
       im_legend, im_string_stats(kkmass-mymass,ndecimal=3), $
         /right, /bottom, box=0

       djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=agerange, $
         yrange=residrange, xtitle='Age [Gyr, isedfit]', $
         ytitle='\Delta[log (M/M_{\odot})] [isedfit minus K-correct]', $
         title=thistitle
       djs_oplot, myage, mymass-kkmass, psym=8
       djs_oplot, !x.crange, [0,0], line=0, thick=2, color='red'

; ---------------
; these plots show that...
       mymass = sfh1.mass50
       mymass_err = sfh1.mass_err
       myage = sfh1.age50
       myebv = sfh1.ebv50
       thistitle = textoidl('SDSS - Z=Z_{\odot}, E(B-V) grid')

       hogg_scatterplot, kkmass, mymass, xsty=1, ysty=1, xrange=massrange, $
         yrange=massrange, xtitle=textoidl('log (M/M_{\odot}) [K-correct]'), $
         ytitle=textoidl('log (M/M_{\odot}) [isedfit]'), title=thistitle, $
         xnpix=20.0, ynpix=20.0, /outliers, outsymsize=0.4, $
         outcolor=djs_icolor('blue')
;      djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=massrange, $
;        yrange=massrange, xtitle='log (M/M_{\odot}) [K-correct]', $
;        ytitle='log (M/M_{\odot}) [isedfit]', title=thistitle
;      oploterror, kkmass, mymass, mymass_err, psym=8
       djs_oplot, !x.crange, !y.crange, line=0, thick=2, color='red'
       im_legend, im_string_stats(kkmass-mymass,ndecimal=3), $
         /right, /bottom, box=0

       djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=agerange, $
         yrange=residrange, xtitle='Age [Gyr, isedfit]', $
         ytitle='\Delta[log (M/M_{\odot})] [isedfit minus K-correct]', $
         title=thistitle
       djs_oplot, myage, mymass-kkmass, psym=8
       djs_oplot, !x.crange, [0,0], line=0, thick=2, color='red'

       djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=ebvrange, $
         yrange=residrange, xtitle='E(B-V) [mag, isedfit]', $
         ytitle='\Delta[log (M/M_{\odot})] [isedfit minus K-correct]', $
         title=thistitle
       djs_oplot, myebv, mymass-kkmass, psym=8
       djs_oplot, !x.crange, [0,0], line=0, thick=2, color='red'

; ---------------
; these plots show that...
       mymass = sfh3.mass50
       mymass_err = sfh3.mass_err
       myage = sfh3.age50
       myebv = sfh3.ebv50
       thistitle = textoidl('SDSS - Z grid, E(B-V) grid')

       hogg_scatterplot, kkmass, mymass, xsty=1, ysty=1, xrange=massrange, $
         yrange=massrange, xtitle=textoidl('log (M/M_{\odot}) [K-correct]'), $
         ytitle=textoidl('log (M/M_{\odot}) [isedfit]'), title=thistitle, $
         xnpix=20.0, ynpix=20.0, /outliers, outsymsize=0.4, $
         outcolor=djs_icolor('blue')
;      djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=massrange, $
;        yrange=massrange, xtitle='log (M/M_{\odot}) [K-correct]', $
;        ytitle='log (M/M_{\odot}) [isedfit]', title=thistitle
;      oploterror, kkmass, mymass, mymass_err, psym=8
       djs_oplot, !x.crange, !y.crange, line=0, thick=2, color='red'
       im_legend, im_string_stats(kkmass-mymass,ndecimal=3), $
         /right, /bottom, box=0

       djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=agerange, $
         yrange=residrange, xtitle='Age [Gyr, isedfit]', $
         ytitle='\Delta[log (M/M_{\odot})] [isedfit minus K-correct]', $
         title=thistitle
       djs_oplot, myage, mymass-kkmass, psym=8
       djs_oplot, !x.crange, [0,0], line=0, thick=2, color='red'

       djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=ebvrange, $
         yrange=residrange, xtitle='E(B-V) [mag, isedfit]', $
         ytitle='\Delta[log (M/M_{\odot})] [isedfit minus K-correct]', $
         title=thistitle
       djs_oplot, myebv, mymass-kkmass, psym=8
       djs_oplot, !x.crange, [0,0], line=0, thick=2, color='red'

       im_plotconfig, psfile=psfile, /psclose, /gzip

stop       
       
; ---------------
       mymass = sfh2.mass50
       mymass_err = sfh2.mass_err
       myage = sfh2.age50
       thistitle = 'SDSS - No dust, range of Z'

       djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=massrange, $
         yrange=massrange, xtitle='log (M/M_{\odot}) [K-correct]', $
         ytitle='log (M/M_{\odot}) [isedfit]', title=thistitle
       oploterror, kkmass, mymass, mymass_err, psym=8
       djs_oplot, !x.crange, !y.crange, line=0, thick=2, color='red'
       ss = im_stats(kkmass-mymass,/verbose)
       cc = get_kbrd(1)



stop
; Brinchmann
       djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=massrange, $
         yrange=massrange, xtitle='log (M/M_{\odot}) [Brinchmann+09]', $
         ytitle='log (M/M_{\odot}) [isedfit]'
       oploterror, mm.median, res.mass50-0.25, res.mass_err, psym=8
       djs_oplot, !x.crange, !y.crange, line=0, thick=2, color='red'
       im_legend, 'SFHGRID-'+string(params.sfhgrid,format='(I2.2)'), $
         /left, /top, box=0
       ss = im_stats(mm.median-(res.mass50-0.25),/verbose)
       cc = get_kbrd(1)
       
       
stop       
       
       djs_oplot, kcorr[these].mass, res.mass50-0.25, ps=8, color='red'
       djs_plot, kcorr[these].mass, jj.mass50-0.25, ps=8, xr=[9,12], yr=[9,12], color='red'
       
       
    endif
    
    
    
stop    
    
    if keyword_set(qaplot) then begin
       isedfit_qaplot, paramfile, isedfit, iopath=iopath, $
         galaxy=sdss[these].galaxy
       return
    endif

;   niceprint, isedfit.mass, sdss[these].mass
    
stop    
stop

; select a volume-limited subsample    
;   ss = read_sdss_vagc_mpa(sample='dr7',poststr='32',/postlss)
;   cc = hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrow=28800L)
;   cat = cc[ss.object_position]
;   mwrfits, cat, 'photo.dr7.bsafe.32.fits.gz', /create
    
;   if (n_elements(phot) eq 0L) then phot = mrdfits(iopath+'photo.dr7.bsafe.32.fits.gz',1)
;   cat = read_sdss_vagc_mpa(sample='dr7',poststr='32',/postlss)
;   kcorr = read_sdss_vagc_mpa(sample='dr7',poststr='32',/kcorr)
;   vmax = read_sdss_vagc_mpa(sample='dr7',poststr='32',/vmax)
;   mpa = read_sdss_vagc_mpa(sample='dr7',poststr='32',/mpamass)

;    ngal = 100
;;   maxvmax = max(1.0/vmax.vmax)
;;   chances = (1.0/vmax.vmax)/maxvmax
;;   these = where(randomu(seed,n_elements(vmax)) lt chances,ngal)
;    myseed = 15.5
;    keep = long(randomu(myseed,ngal)*n_elements(post1))
;;   im_plothist, kcorr[keep].z, bin=0.005

return
end
