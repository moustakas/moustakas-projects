pro isedfit_vs_salim, build_parent=build_parent, preliminaries=preliminaries, $
  models=models, isedfit=isedfit, get_colors=get_colors, kcorrect=kcorrect, $
  qaplot=qaplot, clobber=clobber, compare=compare

; echo "build_isedfit_sfhgrid, 5, /cl, /make" | idl > & ./monte5.log &
; echo "isedfit_vs_salim, /ised, /models, /cl" | idl > & doit1.log &
    
    isedfit_dir = getenv('IM_ARCHIVE_DIR')+'/projects/isedfit_vs_salim/'
    
    prefix = 'salim'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    sfhgrid_paramfile = isedfit_dir+prefix+'_sfhgrid.par'
    supergrid_paramfile = isedfit_dir+prefix+'_supergrid.par'
    
    filterlist = [galex_filterlist(),sdss_filterlist()]
    
; --------------------------------------------------
; build a test sample    
    if keyword_set(build_parent) then begin
       sample = 'dr72' & letter = 'bsafe' & poststr = '0'
       post = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/postlss)
       kauff = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/mpamass)
       totsfr = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/totsfr)

       salim = read_07salim()
       spherematch, post.ra, post.dec, salim.ra, salim.dec, 1D/3600.0, m1, m2
       keep = where((post[m1].z gt 0.05) and (post[m1].z lt 0.2) and $
         (kauff[m1].mass_median gt 0.0) and $
         (3D5*abs(post[m1].z-salim[m2].z) lt 100.0),ngal)
       keep = keep[0:9999]
       ngal = n_elements(keep)

       post = post[m1[keep]]
       kauff = kauff[m1[keep]]
       totsfr = totsfr[m1[keep]]

       vagcpath = getenv('VAGC_REDUX')+'/'
       sdssphot = mrdfits(vagcpath+'object_sdss_imaging.fits.gz',$
         row=post.object_position,1)
;      sdssspec = mrdfits(vagcpath+'object_sdss_spectro.fits.gz',$
;        row=post.object_position,1)
       galex = mrdfits(vagcpath+'object_galex_gr6.fits.gz',$
         row=post.object_position,1)
;      if (n_elements(sdsstwomass) eq 0L) then sdsstwomass = $
;        mrdfits(vagcpath+'object_twomass.fits.gz',$
;        row=phot.object_position,1)

       sdss_to_maggies, maggies, ivarmaggies, calib=sdssphot, flux='model' ; flux='cmodel'
       im_galex_to_maggies, galex, gmaggies, givarmaggies
       
       parent = struct_addtags(sdssphot,struct_trimtags(post,$
         select=['z','absm','object_position']))
       parent = struct_addtags(temporary(parent),replicate({maggies: fltarr(7), $
         ivarmaggies: fltarr(7)},ngal))
       parent.maggies = [gmaggies,maggies]
       parent.ivarmaggies = [givarmaggies,ivarmaggies]

       im_mwrfits, parent, 'parent.fits', /clobber
       im_mwrfits, kauff, 'kauffmann.fits', /clobber
       im_mwrfits, totsfr, 'totsfr.fits', /clobber
       im_mwrfits, salim[m2[keep]], 'salim07.fits', /clobber
    endif

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       parent = mrdfits(isedfit_dir+'parent.fits.gz',1)
       kcorr = mf_do_kcorrect(parent.z,parent.maggies,$
         parent.ivarmaggies,/just,maxiter=0,$
         filterlist=[galex_filterlist(),sdss_filterlist()])
;      kk = galex_kcorrect(parent.z,nmgy=parent.maggies*1D9,$
;        ivar=parent.ivarmaggies/1D18,band_shift=0.1,chi2=chi2,$
;        mass=mass,absmag=absmag,rmaggies=rmaggies)
       im_mwrfits, kcorr, 'kcorr.fits', /clobber
    endif
    
; --------------------------------------------------
; do the preliminaries: build the parameter files and the Monte Carlo
; grids
    if keyword_set(preliminaries) then begin
       write_isedfit_paramfile, filterlist, prefix='salim', minz=0.05, $
         maxz=0.2, nzz=30, isedfit_dir=isedfit_dir, /clobber

       nage = 50
       nmonte = 1000
       write_sfhgrid_paramfile, sfhgrid_paramfile, /clobber, $
         nage=nage, nmonte=nmonte
       write_sfhgrid_paramfile, sfhgrid_paramfile, /append, /preset_bursts, $
         nage=nage, nmonte=nmonte
       
       supergrid = [1,2]
       sfhgrid = [1,2]
       write_supergrid_paramfile, supergrid_paramfile, synthmodels='bc03', $
         supergrid=supergrid, sfhgrid=sfhgrid, /clobber

       build_montegrids, sfhgrid_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, /clobber
    endif

; --------------------------------------------------
; build the models
    if keyword_set(models) then begin
       isedfit_models, isedfit_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, clobber=clobber
    endif
    
; --------------------------------------------------
; do the fitting!  
    if keyword_set(isedfit) then begin
;      salim = mrdfits(isedfit_dir+'salim07.fits.gz',1)
;      ised = mrdfits('test1_bc03_chab_charlot_sfhgrid04.fits.gz',1)
;      index = where(ised.tau_avg gt 2 and ised.tau_avg lt 2.1 and $
;        (ised.sfr_avg-salim.sfr_avg) gt 1.3 and (ised.sfr_avg-salim.sfr_avg) lt 1.4)
;      djs_plot, ised.tau_avg, ised.sfr_avg-salim.sfr_avg, ps=3
;      djs_oplot, ised[index].tau_avg, ised[index].sfr_avg-salim[index].sfr_avg, ps=3, color='orange'
;
;      outprefix = 'test2'
       parent = mrdfits(isedfit_dir+'parent.fits.gz',1)

;; use samir's photometry!
;      parent.maggies = mag2maggies(salim.mab,magerr=salim.mab_err,ivarmaggies=iv)
;      parent.ivarmaggies = iv
;      parent.z = salim.z

       isedfit, isedfit_paramfile, parent.maggies, parent.ivarmaggies, parent.z, $
         result, supergrid_paramfile=supergrid_paramfile, sfhgrid_paramfile=sfhgrid_paramfile, $
         isedfit_dir=isedfit_dir, outprefix=outprefix, galchunksize=1000, $
         clobber=clobber, index=index
    endif 

;   mz_kcorrect_qaplot, kk[index], psfile='junk.ps', in_filterlist=sdss_filterlist()

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       jj = mrdfits(isedfit_dir+'test1_bc03_chab_charlot_sfhgrid08.fits.gz',1)
;      kk = mrdfits(isedfit_dir+'kcorr.fits.gz',1)
       index = random_indices(n_elements(jj),50)
;      index = where(abs(jj.mass_50-kk.k_mass) gt 0.3)
       isedfit_qaplot, paramfile, isedfit, isedfit_dir=isedfit_dir, galaxy=galaxy, $
         index=index, clobber=clobber, outprefix=outprefix
    endif

; --------------------------------------------------
    if keyword_set(compare) then begin
       parent = mrdfits(isedfit_dir+'parent.fits.gz',1)
       kcorr = mrdfits(isedfit_dir+'kcorr.fits.gz',1)
       kauff = mrdfits(isedfit_dir+'kauffmann.fits.gz',1)
       salim = mrdfits(isedfit_dir+'salim07.fits.gz',1)
       ngal = n_elements(parent)

       residrange = 0.9*[-1,1]
       sfhgrid = [3,4,5]
       title = ['','','']

;      sfhgrid = [1,2,3,4]
;      title = [$
;        'Dust, Z=0.004-0.05, No bursts',$
;        'Dust, Z=0.02, No bursts',$
;        'No dust, Z=0.004-0.05, No bursts',$
;        'Dust, Z=0.004-0.05, Bursts']
      
;       sfhgrid = [1,2,4,5,6,7]
;       title = [$
;         'Dust, Z=0.004-0.05, No bursts',$
;         'Dust, Z=0.02, No bursts',$
;;        'No dust, Z=0.004-0.05, No bursts',$
;         'Dust, Z=0.004-0.05, Bursts',$
;         'Dust, Z=0.004-0.05, Weaker Bursts',$
;         'Dust, Z=0.004-0.05, Truncated Bursts',$
;         'More Dust, Z=0.004-0.05, No bursts']
      
       psfile = isedfit_dir+'qaplot_ised_vs_salim.ps'
       im_plotconfig, 10, pos, psfile=psfile, charsize=1.6, $
         yspace=[1.0,1.0], xmargin=[1.2,0.4], ymargin=[0.6,1.1];, $
;        height=[3.0,3.0]

       kauffmass = kauff.mass_median-0.07 ; Kroupa-->Chabrier
       salimmass = salim.mass_avg ; Chabrier
       kcorrmass = kcorr.k_mass ; Chabrier

       umr = kcorr.k_fnuvugrizjhk_absmag_01[2]-kcorr.k_fnuvugrizjhk_absmag_01[4]

; ##################################################
; K-correct vs Kauffmann
; mass
       resid = kauffmass-kcorrmass
       hogg_scatterplot, kauffmass, resid, position=pos[*,0], xsty=1, ysty=1, $
         xrange=[9,11.8], yrange=residrange, xtitle=textoidl('log (M/M_{\odot}) [Kauffmann]'), $
         ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
         levels=[0.1,0.25,0.5,0.75,0.9], /internal
       djs_oplot, !x.crange, [0,0], line=0, color='red'
       im_legend, '\Delta = '+im_string_stats(resid,ndec=3), $
         /left, /top, box=0, charsize=1.4
; color
       hogg_scatterplot, umr, resid, $
         /noerase, position=pos[*,1], xsty=1, ysty=1, $
         xrange=[0.9,3.5], yrange=residrange, xtitle=textoidl('^{0.1}(u - r)'), $
         ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
         levels=[0.1,0.25,0.5,0.75,0.9], /internal
       djs_oplot, !x.crange, [0,0], line=0, color='red'
       xyouts, 0.03, 0.5, 'Mass Residuals [Kauffmann-Kcorrect] (dex)', $
         /normal, orientation=90, align=0.5
       !p.multi = 0
       
; K-correct vs Salim
; mass
       resid = salimmass-kcorrmass
       hogg_scatterplot, kauffmass, resid, position=pos[*,0], xsty=1, ysty=1, $
         xrange=[9,11.8], yrange=residrange, xtitle=textoidl('log (M/M_{\odot}) [Kauffmann]'), $
         ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
         levels=[0.1,0.25,0.5,0.75,0.9], /internal
       djs_oplot, !x.crange, [0,0], line=0, color='red'
       im_legend, '\Delta = '+im_string_stats(resid,ndec=3), $
         /left, /top, box=0, charsize=1.4
; color
       hogg_scatterplot, umr, resid, $
         /noerase, position=pos[*,1], xsty=1, ysty=1, $
         xrange=[0.9,3.5], yrange=residrange, xtitle=textoidl('^{0.1}(u - r)'), $
         ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
         levels=[0.1,0.25,0.5,0.75,0.9], /internal
       djs_oplot, !x.crange, [0,0], line=0, color='red'
       xyouts, 0.03, 0.5, 'Mass Residuals [Salim-Kcorrect] (dex)', $
         /normal, orientation=90, align=0.5
       !p.multi = 0
       
; Kauffmann vs Salim
; mass
       resid = kauffmass-salimmass
       hogg_scatterplot, kauffmass, resid, position=pos[*,0], xsty=1, ysty=1, $
         xrange=[9,11.8], yrange=residrange, xtitle=textoidl('log (M/M_{\odot}) [Kauffmann]'), $
         ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
         levels=[0.1,0.25,0.5,0.75,0.9], /internal
       djs_oplot, !x.crange, [0,0], line=0, color='red'
       im_legend, '\Delta = '+im_string_stats(resid,ndec=3), $
         /left, /top, box=0, charsize=1.4
; color
       hogg_scatterplot, umr, resid, $
         /noerase, position=pos[*,1], xsty=1, ysty=1, $
         xrange=[0.9,3.5], yrange=residrange, xtitle=textoidl('^{0.1}(u - r)'), $
         ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
         levels=[0.1,0.25,0.5,0.75,0.9], /internal
       djs_oplot, !x.crange, [0,0], line=0, color='red'
       xyouts, 0.03, 0.5, 'Mass Residuals [Kauffmann-Salim] (dex)', $
         /normal, orientation=90, align=0.5
       !p.multi = 0
       
; now compare isedfit against K-correct, Kauffmann, and Salim
       for ii = 0, n_elements(sfhgrid)-1 do begin
          sfhg = 'sfhgrid'+string(sfhgrid[ii],format='(I2.2)')
          red = '_charlot'
          sfhgridfile = isedfit_dir+'test1_bc03_chab'+red+'_'+sfhg+'.fits.gz'
          splog, 'Reading '+sfhgridfile
          ised = mrdfits(sfhgridfile,1)
          isedmass = ised.mass_50 ; Chabrier
; ##################################################
; isedfit vs K-correct
; mass
          resid = isedmass-kcorrmass
          hogg_scatterplot, isedmass, resid, position=pos[*,0], xsty=1, ysty=1, $
            xrange=[9,11.8], yrange=residrange, xtitle=textoidl('log (M/M_{\odot})'), $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
          im_legend, '\Delta = '+im_string_stats(resid,ndec=3), $
            /left, /top, box=0, charsize=1.4
; color
          hogg_scatterplot, umr, resid, $
            /noerase, position=pos[*,1], xsty=1, ysty=1, $
            xrange=[0.9,3.5], yrange=residrange, xtitle=textoidl('^{0.1}(u - r)'), $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; AGE
          hogg_scatterplot, alog10(ised.age_50), resid, /noerase, position=pos[*,2], xsty=1, ysty=1, $
            xrange=alog10([1,20]), yrange=residrange, xtitle='log (Age) (Gyr)', $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; tau
          hogg_scatterplot, alog10(ised.tau_50>0.4), resid, /noerase, position=pos[*,3], xsty=3, ysty=1, $
            xrange=alog10([1,10]), yrange=residrange, xtitle=textoidl('\tau (Gyr)'), $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; dust
          hogg_scatterplot, alog10(ised.av_50>0.01), resid, /noerase, position=pos[*,4], xsty=3, ysty=1, $
            xrange=alog10([0.1,6]), yrange=residrange, xtitle=textoidl('A_{V} (mag)'), $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; metallicity
          hogg_scatterplot, ised.Z_50, resid, /noerase, position=pos[*,5], xsty=3, ysty=1, $
            xrange=[0.003,0.052], yrange=residrange, xtitle='Metallicity', $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; titles
          xyouts, 0.03, 0.5, 'Mass Residuals [iSEDfit-Kcorrect] (dex)', $
            /normal, orientation=90, align=0.5
          xyouts, 0.5, 0.96, '('+string(sfhgrid[ii],format='(I2.2)')+') '+$
            title[ii], /normal, align=0.5
; ##################################################
; comparison with Kauffmann+
; mass
          resid = isedmass-kauffmass
          hogg_scatterplot, isedmass, resid, position=pos[*,0], xsty=1, ysty=1, $
            xrange=[9,11.8], yrange=residrange, xtitle=textoidl('log (M/M_{\odot})'), $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
          im_legend, '\Delta = '+im_string_stats(resid,ndec=3), $
            /left, /top, box=0, charsize=1.4
; color
          hogg_scatterplot, umr, resid, $
            /noerase, position=pos[*,1], xsty=1, ysty=1, $
            xrange=[0.9,3.5], yrange=residrange, xtitle=textoidl('^{0.1}(u - r)'), $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; AGE
          hogg_scatterplot, alog10(ised.age_50), resid, /noerase, position=pos[*,2], xsty=1, ysty=1, $
            xrange=alog10([1,20]), yrange=residrange, xtitle='log (Age) (Gyr)', $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; tau
          hogg_scatterplot, alog10(ised.tau_50>0.4), resid, /noerase, position=pos[*,3], xsty=3, ysty=1, $
            xrange=alog10([1,10]), yrange=residrange, xtitle=textoidl('\tau (Gyr)'), $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; dust
          hogg_scatterplot, alog10(ised.av_50>0.01), resid, /noerase, position=pos[*,4], xsty=3, ysty=1, $
            xrange=alog10([0.1,6]), yrange=residrange, xtitle=textoidl('A_{V} (mag)'), $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; metallicity
          hogg_scatterplot, ised.Z_50, resid, /noerase, position=pos[*,5], xsty=3, ysty=1, $
            xrange=[0.003,0.052], yrange=residrange, xtitle='Metallicity', $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; titles
          xyouts, 0.03, 0.5, 'Mass Residuals [iSEDfit-Kauffmann] (dex)', $
            /normal, orientation=90, align=0.5
          xyouts, 0.5, 0.96, '('+string(sfhgrid[ii],format='(I2.2)')+') '+$
            title[ii], /normal, align=0.5
; ##################################################
; comparison with Salim+
; mass
          resid = isedmass-salimmass
          hogg_scatterplot, isedmass, resid, position=pos[*,0], xsty=1, ysty=1, $
            xrange=[9,11.8], yrange=residrange, xtitle=textoidl('log (M/M_{\odot})'), $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
          im_legend, '\Delta = '+im_string_stats(resid,ndec=3), $
            /left, /top, box=0, charsize=1.4
; color
          hogg_scatterplot, umr, resid, $
            /noerase, position=pos[*,1], xsty=1, ysty=1, $
            xrange=[0.9,3.5], yrange=residrange, xtitle=textoidl('^{0.1}(u - r)'), $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; AGE
          hogg_scatterplot, alog10(ised.age_50), resid, /noerase, position=pos[*,2], xsty=1, ysty=1, $
            xrange=alog10([1,20]), yrange=residrange, xtitle='log (Age) (Gyr)', $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; tau
          hogg_scatterplot, alog10(ised.tau_50>0.4), resid, /noerase, position=pos[*,3], xsty=3, ysty=1, $
            xrange=alog10([1,10]), yrange=residrange, xtitle=textoidl('\tau (Gyr)'), $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; dust
          hogg_scatterplot, alog10(ised.av_50>0.01), resid, /noerase, position=pos[*,4], xsty=3, ysty=1, $
            xrange=alog10([0.1,6]), yrange=residrange, xtitle=textoidl('A_{V} (mag)'), $
            ytitle='', /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; metallicity
          hogg_scatterplot, ised.Z_50, resid, /noerase, position=pos[*,5], xsty=3, ysty=1, $
            xrange=[0.003,0.052], yrange=residrange, xtitle='Metallicity', $
            ytitle='', ytickname=replicate(' ',10), /outlier, outpsym=6, outcolor=djs_icolor('default'), $
            levels=[0.1,0.25,0.5,0.75,0.9], /internal
          djs_oplot, !x.crange, [0,0], line=0, color='red'
; titles
          xyouts, 0.03, 0.5, 'Mass Residuals [iSEDfit-Salim] (dex)', $
            /normal, orientation=90, align=0.5
          xyouts, 0.5, 0.96, '('+string(sfhgrid[ii],format='(I2.2)')+') '+$
            title[ii], /normal, align=0.5
       endfor
       im_plotconfig, psfile=psfile, /psclose, /gzip ; /pdf
;      spawn, 'rsync -auv '+psfile+'.gz ~/', /sh
    endif
    stop
return
end
