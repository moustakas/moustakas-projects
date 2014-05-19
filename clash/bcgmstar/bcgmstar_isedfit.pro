pro bcgmstar_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  qaplot_results=qaplot_results, clobber=clobber
; jm13dec29siena - do SED-fitting

    prefix = 'bcgmstar'
    
    sersicpath = bcgmstar_path(/sersic)
    isedfit_dir = bcgmstar_path(/isedfit)
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; read the sample
    sample = read_bcgmstar_sample(/zsort)
    struct_print, sample
    ncl = n_elements(sample)

    filterlist = bcgmstar_filterlist(short=short)
    nfilt = n_elements(filterlist)

; gather the photometry
;    cat = rsex(isedfit_dir+'flx_iso.dec26')

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       spsmodels = 'fsps_v2.4_miles'
       imf = 'salp'
       nmodel = 5000L
; SFH, age, and metallicity are free
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, use_redshift=sample.z, $
         spsmodels=spsmodels, imf=imf, igm=0, redcurve='none', AV=[0.0,0.0], $
         sfhgrid=1, nmodel=nmodel, age=[1.0,11.2], tau=[0.0,2.0], $
         Zmetal=[0.01,0.03], clobber=clobber
; SFH and age are free; metallicity is fixed at super-solar
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, use_redshift=sample.z, $
         spsmodels=spsmodels, imf=imf, igm=0, redcurve='none', AV=[0.0,0.0], $
         sfhgrid=2, nmodel=nmodel, age=[1.0,11.2], tau=[0.0,2.0], $
         Zmetal=[0.03,0.03], clobber=clobber, /append
; age is free; SFH fixed at an SSP and metallicity fixed at super-solar
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, use_redshift=sample.z, $
         spsmodels=spsmodels, imf=imf, igm=0, redcurve='none', AV=[0.0,0.0], $
         sfhgrid=3, nmodel=nmodel, age=[1.0,11.2], tau=[0.0,0.0], $
         Zmetal=[0.03,0.03], clobber=clobber, /append
; age and metallicity are free; SFH fixed at an SSP
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, use_redshift=sample.z, $
         spsmodels=spsmodels, imf=imf, igm=0, redcurve='none', AV=[0.0,0.0], $
         sfhgrid=4, nmodel=nmodel, age=[1.0,11.2], tau=[0.0,0.0], $
         Zmetal=[0.01,0.03], clobber=clobber, /append
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
; do the fitting!
    if keyword_set(isedfit) then begin
       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          outprefix = prefix+'_'+cluster

          phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
          naper = n_elements(phot[0].photradius_kpc)+1 ; radial bins + integrated
          nband = n_elements(phot)

          maggies = fltarr(nfilt,naper)
          ivarmaggies = fltarr(nfilt,naper)
          for ii = 0, nfilt-1 do begin
             this = where(short[ii] eq strtrim(phot.band,2))
             if this[0] ne -1 then begin
                good = where(phot[this].photradius_kpc_in gt phot[this].rmin_kpc and $
                  phot[this].photradius_kpc_out lt phot[this].rmax_kpc,comp=extrap,ncomp=nextrap)

; take the measure photometry at face-value, but don't give any
; weight to bands that have been extrapolated                
                maggies[ii,*] = [phot[this].maggies_int,phot[this].maggies]
                ivarmaggies[ii,*] = [phot[this].ivarmaggies_int,phot[this].ivarmaggies]
                if nextrap ne 0 then ivarmaggies[ii,1+extrap] = 0.0 ; offset from integrated
;               if cluster eq 'a209' and short[ii] eq 'f390w' then maggies[ii,0] = 0.0 ; crap!
             endif
          endfor
          
          isedfit, isedfit_paramfile, maggies, ivarmaggies, replicate(sample[ic].z,naper), $
            thissfhgrid=thissfhgrid, isedfit_dir=isedfit_dir, isedfit_results=ised, $
            isedfit_post=isedpost, clobber=clobber, outprefix=outprefix
       endfor
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=bessell_filterlist(), band_shift=0.0, $
         clobber=clobber
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
          galaxy = ['Integrated','R = '+strtrim(string(phot[0].photradius_kpc_in,$
            format='(F12.1)'),2)+'-'+strtrim(string(phot[0].photradius_kpc_out,$
            format='(F12.1)'),2)+' kpc']
;         galaxy = ['Integrated','R = '+strtrim(string(phot[0].photradius_kpc,$
;           format='(F12.1)'),2)+' kpc']
          outprefix = prefix+'_'+cluster
          isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
            montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
            clobber=clobber, /xlog, galaxy=galaxy, index=index, $
            outprefix=outprefix, nsigma=1.0
       endfor
    endif
    
; --------------------------------------------------
; generate some QAplots of the results
    if keyword_set(qaplot_results) then begin
       sample = read_bcgmstar_sample()
       ncl = n_elements(sample)

       ncol = 3 ; number of columns
       nrow = ceil(ncl/float(ncol))
       yrange = [1.0,11.0]

       psfile = isedfit_dir+'qa_isedfit.ps'
       im_plotconfig, 0, pos, psfile=psfile, charsize=1.1
       pos = im_getposition(nx=ncol,ny=nrow,yspace=0.0,xspace=0.0,$
         xmargin=[0.9,0.4],width=2.4)

       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          outprefix = prefix+'_'+cluster
          title = repstr(repstr(strupcase(sample[ic].dirname),'_',' '),'ABELL','Abell')

          phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
          rad = phot[0].photradius_kpc

          ised1 = read_isedfit(isedfit_paramfile,outprefix=outprefix,$
            thissfhgrid=1,/silent);,index=lindgen(n_elements(rad))+1)
          ised2 = read_isedfit(isedfit_paramfile,outprefix=outprefix,$
            thissfhgrid=2,/silent);,index=lindgen(n_elements(rad))+1)
          ised3 = read_isedfit(isedfit_paramfile,outprefix=outprefix,$
            thissfhgrid=3,/silent);,index=lindgen(n_elements(rad))+1)
          ised4 = read_isedfit(isedfit_paramfile,outprefix=outprefix,$
            thissfhgrid=4,/silent);,index=lindgen(n_elements(rad))+1)
;         niceprint, rad, ised1.mstar_50, ised1.age_50, ised1.zmetal_50
;if ic eq 1 then stop

          intage1 = ised1[0].age_50
          intage2 = ised2[0].age_50
          intage3 = ised3[0].age_50
          intage4 = ised4[0].age_50

          ised1 = ised1[1:n_elements(rad)]
          ised2 = ised2[1:n_elements(rad)]
          ised3 = ised3[1:n_elements(rad)]
          ised4 = ised4[1:n_elements(rad)]

; minimum of 5 bands fitted          
          good = where(ised1.chi2 lt 1E6 and total(ised1.ivarmaggies gt 0,1) ge 6)
          splog, cluster, max(rad[good])

          if ic ge ncl-3 then begin
             xtitle = 'Equivalent Radius (kpc)'
             delvarx, xtickname
          endif else begin
             xtitle = ''
             xtickname = replicate(' ',10)
          endelse

          if ic mod 3 eq 0 then begin
             ytitle = 'Age (Gyr)'
             delvarx, ytickname
          endif else begin
             ytitle = ''
             ytickname = replicate(' ',10)
          endelse
          
          djs_plot, [0], [0], /nodata, /xlog, noerase=ic gt 0, $
            xrange=[0.4,110], xsty=1, yrange=yrange, position=pos[*,ic], $
            xtickname=xtickname, ytickname=ytickname, xtitle=xtitle, $
            symsize=0.5, ysty=1, ytitle=ytitle
          im_legend, strupcase(cluster), /left, /top, box=0, charsize=1.0, margin=0
;         im_legend, title, /left, /top, box=0, charsize=1.0, margin=0

          djs_oplot, rad[good], ised1[good].age_50, color=cgcolor('dodger blue'), $
            line=0, thick=6, psym=-symcat(15), symsize=0.8
          djs_oplot, rad[good], ised2[good].age_50, color=cgcolor('firebrick'), $
            line=0, thick=6, psym=-symcat(15), symsize=0.8
          djs_oplot, rad[good], ised3[good].age_50, color=cgcolor('orange'), $
            line=0, thick=6, psym=-symcat(15), symsize=0.8

          djs_oplot, rad[good], ised4[good].age_50, color=cgcolor('purple'), $
            line=0, thick=6, psym=-symcat(15), symsize=0.8
          djs_oplot, rad[good], ised4[good].age_50+ised4[good].age_err/2, $
            color=cgcolor('purple'), line=5, thick=6
          djs_oplot, rad[good], ised4[good].age_50-ised4[good].age_err/2, $
            color=cgcolor('purple'), line=5, thick=6

;         djs_oplot, 10^!x.crange, getage(sample[ic].z)*[1,1], line=1

          plots, 1.0, intage1, psym=symcat(6), color=cgcolor('dodger blue')
          plots, 1.0, intage2, psym=symcat(6), color=cgcolor('firebrick')
;         plots, 1.0, intage3, psym=symcat(6), color=cgcolor('orange')
          plots, 1.0, intage4, psym=symcat(6), color=cgcolor('purple')

;         oploterror, rad[good], ised[good].age_50, ised[good].age_err, $
;           psym=symcat(16), color=cgcolor('dodger blue'), $
;           errcolor=cgcolor('dodger blue')
;         if extrap[0] ne -1 then oploterror, rad[extrap], ised[extrap].age_50, $
;           ised[extrap].age_err, psym=symcat(6), color=cgcolor('dodger blue'), $
;           errcolor=cgcolor('dodger blue')
       endfor
       im_plotconfig, psfile=psfile, /psclose, /pdf
    endif
    
return
end
