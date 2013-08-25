pro isedfitvsalim, build_parent=build_parent, write_paramfile=write_paramfile, $
  build_grids=build_grids, model_photometry=model_photometry, qaplot_models=qaplot_models, $
  isedfit=isedfit, kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  compare=compare, clobber=clobber
; jm13aug23siena - compare iSEDfit to Salim+07
    
    prefix = 'isedfitvsalim'
    isedfit_dir = getenv('IM_PROJECTS_DIR')+'/isedfit/isedfitvsalim/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    filterlist = [$
      'galex_FUV.par',$
      'galex_NUV.par',$
      'sdss_u0.par',$
      'sdss_g0.par',$
      'sdss_r0.par',$
      'sdss_i0.par',$
      'sdss_z0.par']
    
; ---------------------------------------------------------------------------
; build the parent sample for this test
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
         (kauff[m1].mass_median gt 0.0),ngal)
;      keep = where((post[m1].z gt 0.05) and (post[m1].z lt 0.2) and $
;        (kauff[m1].mass_median gt 0.0) and $
;        (3D5*abs(post[m1].z-salim[m2].z) lt 100.0),ngal)
;      keep = keep[0:14999]
;      ngal = n_elements(keep)
       
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
;      wise = mrdfits(vagcpath+'object_wise.fits.gz',$
;        row=post.object_position,1)
;      if (n_elements(sdsstwomass) eq 0L) then sdsstwomass = $
;        mrdfits(vagcpath+'object_twomass.fits.gz',$
;        row=phot.object_position,1)

       sdss_to_maggies, maggies, ivarmaggies, calib=sdssphot, flux='cmodel'
       im_galex_to_maggies, galex, gmaggies, givarmaggies
;      wise_to_maggies, wise, wmaggies, wivarmaggies, /mpro

       parent = struct_addtags(sdssphot,struct_trimtags(post,$
         select=['z','absm','object_position']))
       parent = struct_addtags(temporary(parent),replicate($
         {maggies: fltarr(7), ivarmaggies: fltarr(7)},ngal))
       parent.maggies = [gmaggies,maggies];,wmaggies[0:1,*]]
       parent.ivarmaggies = [givarmaggies,ivarmaggies];,wivarmaggies[0:1,*]]

       hasuv = where(total(gmaggies gt 0,1) ge 1,ngal)
       
       im_mwrfits, parent[hasuv], isedfit_dir+'isedfitvsalim.fits', /clobber
       im_mwrfits, kauff[hasuv], isedfit_dir+'kauffmann.fits', /clobber
       im_mwrfits, totsfr[hasuv], isedfit_dir+'totsfr.fits', /clobber
       im_mwrfits, salim[m2[keep[hasuv]]], isedfit_dir+'salim07.fits', /clobber
    endif

; --------------------------------------------------
; write the parameter file; these have been chosen to match Salim+07
; as well as possible 
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels='bc03_stelib', $
         imf='chab', redcurve='charlot', /igm, zminmax=[0.05,0.2], zbin=0.005, $
         nmodel=50000L, age=[0.1,13.5], tau=[0.01,1.0], Zmetal=[0.002,0.04], $
         AV=[0.35,2.0], mu=[0.1,4.0], pburst=0.5, interval_pburst=2.0, $
         tburst=[0.1,13.5], /oneovertau, galchunksize=2500L, clobber=clobber
    endif

; --------------------------------------------------
; build the Monte Carlo grids    
    if keyword_set(build_grids) then begin
       isedfit_montegrids, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, clobber=clobber
    endif

; --------------------------------------------------
; calculate the model photometry 
    if keyword_set(model_photometry) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber
    endif

; --------------------------------------------------
; generate the model photometry QAplots
    if keyword_set(qaplot_models) then begin
       cat = mrdfits(isedfit_dir+'isedfitvsalim.fits.gz',1)
       thesefilters = ['galex_NUV','sdss_g0','sdss_r0','sdss_i0']
       isedfit_qaplot_models, isedfit_paramfile, cat.maggies, $
         cat.ivarmaggies, cat.z, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, thesefilters=thesefilters, clobber=clobber
    endif
    
; --------------------------------------------------
; fit!
    if keyword_set(isedfit) then begin
       cat = mrdfits(isedfit_dir+'isedfitvsalim.fits.gz',1)
; use samir's photometry & redshifts!
       salim = mrdfits(isedfit_dir+'salim07.fits.gz',1)
       mm = mag2maggies(salim.mab,magerr=salim.mab_err,ivarmaggies=iv)
       cat.maggies = mm
       cat.ivarmaggies = iv
;      cat.z = salim.z
       isedfit, isedfit_paramfile, cat.maggies, cat.ivarmaggies, $
         cat.z, ra=cat.ra, dec=cat.dec, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       isedfit_qaplot_sed, isedfit_paramfile, nrandom=50, $
         isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber, /xlog
    endif

    splog, 'Compare with Kauffman and Brinchmann!!'
    
; ---------------------------------------------------------------------------
; compare with Salim+07 and MPA+JHU
    if keyword_set(compare) then begin
       ss = mrdfits(isedfit_dir+'salim07.fits.gz',1)
       rr = read_isedfit(isedfit_paramfile,isedfit_dir=isedfit_dir)
       kk = mrdfits(isedfit_dir+'kauffmann.fits.gz',1) ; Kauffmann+03

       psfile = isedfit_dir+'qa_isedfit_vs_salim.ps'
       im_plotconfig, 4, pos, psfile=psfile, yspace=0.0, width=5.0, $
         height=2.5*[1,1,1], charsize=1.5

       im_hogg_scatterplot, rr.mstar_avg, rr.sfr100_avg-ss.sfr_avg, yminor=2, $
         xsty=1, ysty=1, position=pos[*,0], xrange=[9,12], yrange=2.5*[-1,1], $
         /internal, /outlier, /nogrey, outcolor=im_color('dodger blue'), $
         xtickname=replicate(' ',10), $
         ytitle=textoidl('\Delta'+'log (SFR) (M'+sunsymbol()+' yr^{-1})'), $
         title='iSEDfit vs Salim+07'
       im_legend, im_string_stats(rr.sfr100_avg-ss.sfr_avg,type=3,sigrej=3.0), /left, /top, $
         box=0, charsize=1.4, margin=0
       mm = im_medxbin(rr.mstar_avg,rr.sfr100_avg-ss.sfr_avg,0.25,/ver,minpts=50)
       djs_oplot, !x.crange, [0,0], line=0
       djs_oplot, mm.medx, mm.medy, line=0, thick=8
       djs_oplot, mm.medx, mm.quant75, line=5, thick=8
       djs_oplot, mm.medx, mm.quant25, line=5, thick=8

       im_hogg_scatterplot, rr.mstar_avg, rr.mstar_avg-ss.mass_avg, /noerase, yminor=2, $
         xsty=1, ysty=1, position=pos[*,1], xrange=[9,12], yrange=0.5*[-1,1], $
         /internal, /outlier, /nogrey, outcolor=im_color('dodger blue'), $
         xtickname=replicate(' ',10), ytitle=textoidl('\Delta'+'log (M_{*}/M'+sunsymbol()+')')
       im_legend, im_string_stats(rr.mstar_avg-ss.mass_avg,type=3,sigrej=3.0), /left, /top, $
         box=0, charsize=1.4, margin=0
       mm = im_medxbin(rr.mstar_avg,rr.mstar_avg-ss.mass_avg,0.25,/ver,minpts=50)
       djs_oplot, !x.crange, [0,0], line=0
       djs_oplot, mm.medx, mm.medy, line=0, thick=8
       djs_oplot, mm.medx, mm.quant75, line=5, thick=8
       djs_oplot, mm.medx, mm.quant25, line=5, thick=8

       im_hogg_scatterplot, rr.mstar_avg, (rr.sfr100_avg-rr.mstar_avg)-ss.sfrm_avg, /noerase, yminor=2, $
         xsty=1, ysty=1, position=pos[*,2], xrange=[9,12], yrange=3.0*[-1,1], $
         /internal, /outlier, /nogrey, outcolor=im_color('dodger blue'), $
         xtitle=textoidl('log (M_{*}/M'+sunsymbol()+')'), $
         ytitle=textoidl('\Delta'+'log (SFR/M_{*}) (yr^{-1})')
       im_legend, im_string_stats((rr.sfr100_avg-rr.mstar_avg)-ss.sfrm_avg,type=3,sigrej=3.0), /left, /top, $
         box=0, charsize=1.4, margin=0
       mm = im_medxbin(rr.mstar_avg,(rr.sfr100_avg-rr.mstar_avg)-ss.sfrm_avg,0.25,/ver,minpts=50)
       djs_oplot, !x.crange, [0,0], line=0
       djs_oplot, mm.medx, mm.medy, line=0, thick=8
       djs_oplot, mm.medx, mm.quant75, line=5, thick=8
       djs_oplot, mm.medx, mm.quant25, line=5, thick=8

; compare with Brinchmann & Kauffmann
       imfcor = alog10(1.06)
       im_hogg_scatterplot, rr.mstar_avg, rr.sfr100_avg-(kk.sfr_avg-imfcor), yminor=2, $
         xsty=1, ysty=1, position=pos[*,0], xrange=[9,12], yrange=2.5*[-1,1], $
         /internal, /outlier, /nogrey, outcolor=im_color('dodger blue'), $
         xtickname=replicate(' ',10), $
         ytitle=textoidl('\Delta'+'log (SFR) (M'+sunsymbol()+' yr^{-1})'), $
         title='iSEDfit vs MPA/JHU (DR7)'
       im_legend, im_string_stats(rr.sfr100_avg-(kk.sfr_avg-imfcor),type=3,sigrej=3.0), /left, /top, $
         box=0, charsize=1.4, margin=0
       mm = im_medxbin(rr.mstar_avg,rr.sfr100_avg-(kk.sfr_avg-imfcor),0.25,/ver,minpts=50)
       djs_oplot, !x.crange, [0,0], line=0
       djs_oplot, mm.medx, mm.medy, line=0, thick=8
       djs_oplot, mm.medx, mm.quant75, line=5, thick=8
       djs_oplot, mm.medx, mm.quant25, line=5, thick=8
       
       im_hogg_scatterplot, rr.mstar_avg, rr.mstar_avg-(kk.mass_avg-imfcor), /noerase, yminor=2, $
         xsty=1, ysty=1, position=pos[*,1], xrange=[9,12], yrange=0.5*[-1,1], $
         /internal, /outlier, /nogrey, outcolor=im_color('dodger blue'), $
         xtickname=replicate(' ',10), ytitle=textoidl('\Delta'+'log (M_{*}/M'+sunsymbol()+')')
       im_legend, im_string_stats(rr.mstar_avg-(kk.mass_avg-imfcor),type=3,sigrej=3.0), /left, /top, $
         box=0, charsize=1.4, margin=0
       mm = im_medxbin(rr.mstar_avg,rr.mstar_avg-(kk.mass_avg-imfcor),0.25,/ver,minpts=50)
       djs_oplot, !x.crange, [0,0], line=0
       djs_oplot, mm.medx, mm.medy, line=0, thick=8
       djs_oplot, mm.medx, mm.quant75, line=5, thick=8
       djs_oplot, mm.medx, mm.quant25, line=5, thick=8

       im_hogg_scatterplot, rr.mstar_avg, (rr.sfr100_avg-rr.mstar_avg)-(kk.sfr_avg-kk.mass_avg), $
         /noerase, yminor=2, $
         xsty=1, ysty=1, position=pos[*,2], xrange=[9,12], yrange=3.0*[-1,1], $
         /internal, /outlier, /nogrey, outcolor=im_color('dodger blue'), $
         xtitle=textoidl('log (M_{*}/M'+sunsymbol()+')'), $
         ytitle=textoidl('\Delta'+'log (SFR/M_{*}) (yr^{-1})')
       im_legend, im_string_stats((rr.sfr100_avg-rr.mstar_avg)-(kk.sfr_avg-kk.mass_avg),$
         type=3,sigrej=3.0), /left, /top, box=0, charsize=1.4, margin=0
       mm = im_medxbin(rr.mstar_avg,(rr.sfr100_avg-rr.mstar_avg)-(kk.sfr_avg-kk.mass_avg),$
         0.25,/ver,minpts=50)
       djs_oplot, !x.crange, [0,0], line=0
       djs_oplot, mm.medx, mm.medy, line=0, thick=8
       djs_oplot, mm.medx, mm.quant75, line=5, thick=8
       djs_oplot, mm.medx, mm.quant25, line=5, thick=8

       im_plotconfig, psfile=psfile, /psclose, /pdf
    endif

return
end
