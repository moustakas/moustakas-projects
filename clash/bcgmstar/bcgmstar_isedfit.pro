pro bcgmstar_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  qaplot_results=qaplot_results, clobber=clobber, parse_massprofiles=parse_massprofiles
; jm13dec29siena - do SED-fitting

    prefix = 'bcgmstar'
    
    sersicpath = bcgmstar_path(/sersic)
    isedfit_dir = bcgmstar_path(/isedfit)
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    massprofpath = bcgmstar_path(/massprofiles)
    lensingpath = bcgmstar_path()+'lensing-profiles-adi/'
    
; read the sample
    sample = read_bcgmstar_sample(/zsort)
;   sample = sample[11]
;   struct_print, sample
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
       nmodel = 10000L
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
          naper = n_elements(phot[0].photradius_kpc)+4 ; integrated,30,50,70 + radial bins
          nband = n_elements(phot)
          
          maggies = fltarr(nfilt,naper)
          ivarmaggies = fltarr(nfilt,naper)
          for ii = 0, nfilt-1 do begin
             this = where(short[ii] eq strtrim(phot.band,2))
             if this[0] ne -1 then begin
                maggies[ii,*] = [phot[this].maggies_int,phot[this].maggies_30,$
                  phot[this].maggies_50,phot[this].maggies_70,phot[this].maggies]
                ivarmaggies[ii,*] = [phot[this].ivarmaggies_int,phot[this].ivarmaggies_30,$
                  phot[this].ivarmaggies_50,phot[this].ivarmaggies_70,phot[this].ivarmaggies]
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
       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          outprefix = prefix+'_'+cluster
          isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
            montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
            absmag_filterlist=bessell_filterlist(), band_shift=0.0, $
            clobber=clobber, outprefix=outprefix
       endfor
    endif 

; --------------------------------------------------
; parse the output
    if keyword_set(parse_massprofiles) then begin
       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
          outfile = massprofpath+cluster+'-massprofile.fits'

          phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
          ised = read_isedfit(isedfit_paramfile,outprefix=prefix+'_'+cluster,$
            isedfit_dir=isedfit_dir,/silent)
          nrad = n_elements(phot[0].photradius_kpc)
          noff = 4 ; offset

          out = {$
            cluster:        cluster,$
            photradius_kpc: phot[0].photradius_kpc, $

            mstar_grid01:     ised[noff:nrad+noff-1,0].mstar_50,$
            mstar_grid02:     ised[noff:nrad+noff-1,1].mstar_50,$
            mstar_grid03:     ised[noff:nrad+noff-1,2].mstar_50,$
            mstar_grid04:     ised[noff:nrad+noff-1,3].mstar_50,$
            mstar_err_grid01: ised[noff:nrad+noff-1,0].mstar_err,$
            mstar_err_grid02: ised[noff:nrad+noff-1,1].mstar_err,$
            mstar_err_grid03: ised[noff:nrad+noff-1,2].mstar_err,$
            mstar_err_grid04: ised[noff:nrad+noff-1,3].mstar_err,$

;           mstar_int_grid01:     ised[0,0].mstar_50,$
;           mstar_int_grid02:     ised[0,1].mstar_50,$
;           mstar_int_grid03:     ised[0,2].mstar_50,$
;           mstar_int_grid04:     ised[0,3].mstar_50,$
;           mstar_int_err_grid01: ised[0,0].mstar_err,$
;           mstar_int_err_grid02: ised[0,1].mstar_err,$
;           mstar_int_err_grid03: ised[0,2].mstar_err,$
;           mstar_int_err_grid04: ised[0,3].mstar_err,$

            mstar:     fltarr(nrad),$ ; average profile
            mstar_err: fltarr(nrad),$

            mstar_30:      0.0,$
            mstar_30_err:  0.0,$
            mstar_50:      0.0,$
            mstar_50_err:  0.0,$
            mstar_70:      0.0,$
            mstar_70_err:  0.0,$
            mstar_int:     0.0,$
            mstar_int_err: 0.0}

          allmstar = [[out.mstar_grid01],[out.mstar_grid02],[out.mstar_grid03],[out.mstar_grid04]]
          allmstar_err = [[out.mstar_err_grid01],[out.mstar_err_grid02],$
            [out.mstar_err_grid03],[out.mstar_err_grid04]]

          out.mstar = mean(allmstar,dim=2)
          out.mstar_err = mean(allmstar_err,dim=2)
;         out.mstar_err = stddev(allmstar,dim=2)

          allmstar_int = [ised[0,0].mstar_50,ised[0,1].mstar_50,ised[0,2].mstar_50,ised[0,3].mstar_50]
          allmstar_30  = [ised[1,0].mstar_50,ised[1,1].mstar_50,ised[1,2].mstar_50,ised[1,3].mstar_50]
          allmstar_50  = [ised[2,0].mstar_50,ised[2,1].mstar_50,ised[2,2].mstar_50,ised[2,3].mstar_50]
          allmstar_70  = [ised[3,0].mstar_50,ised[3,1].mstar_50,ised[3,2].mstar_50,ised[3,3].mstar_50]

          allmstar_int_err = [ised[0,0].mstar_err,ised[0,1].mstar_err,ised[0,2].mstar_err,ised[0,3].mstar_err]
          allmstar_30_err  = [ised[1,0].mstar_err,ised[1,1].mstar_err,ised[1,2].mstar_err,ised[1,3].mstar_err]
          allmstar_50_err  = [ised[2,0].mstar_err,ised[2,1].mstar_err,ised[2,2].mstar_err,ised[2,3].mstar_err]
          allmstar_70_err  = [ised[3,0].mstar_err,ised[3,1].mstar_err,ised[3,2].mstar_err,ised[3,3].mstar_err]

          out.mstar_int = mean(allmstar_int)
          out.mstar_30  = mean(allmstar_30)
          out.mstar_50  = mean(allmstar_50)
          out.mstar_70  = mean(allmstar_70)

          out.mstar_int_err = mean(allmstar_int_err)
          out.mstar_30_err  = mean(allmstar_30_err)
          out.mstar_50_err  = mean(allmstar_50_err)
          out.mstar_70_err  = mean(allmstar_70_err)
;         out.mstar_int_err = stddev(allmstar_int)
;         out.mstar_30_err  = stddev(allmstar_30)
;         out.mstar_50_err  = stddev(allmstar_50)
;         out.mstar_70_err  = stddev(allmstar_70)
          
; pack in the strong-lensing profiles
          lensingsuffix = strupcase(repstr(repstr(cluster,'macs','M'),'clj','cl'))
          lensingfile = lensingpath+'CenterJohnKeiichi_profile_ltm_'+lensingsuffix+'.txt'
          if file_test(lensingfile) eq 0 then begin
             splog, lensingfile+' not found!'
             stop
          endif

; col 1: radius from center, in arcsec
; col 2 =(surface_density);
; col 3: log10(radius);
; col 4: log10(surface_density);
; col 5 irrelevant
; col 6 : projected mass interior to r i.e. M_2D (<r)
; col 7-8 surface density 68.3% range.
; col 9-10 68.3% range on the projected enclosed mass of col 6.
          readcol, lensingfile, radius_arcsec, density, log10radius, $
            log10density, junk, mass_encl, density_lo, density_hi, $
            mass_encl_lo, mass_encl_hi, /silent

          keep = where(radius_arcsec*arcsec2kpc lt 300,nlens)
          out = struct_addtags(out,{$
            lens_radius_kpc: radius_arcsec[keep]*arcsec2kpc,$
            lens_mass_encl:  mass_encl[keep]})

          im_mwrfits, out, outfile, clobber=clobber
       endfor
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

;          maggies = fltarr(nfilt,naper)
;          ivarmaggies = fltarr(nfilt,naper)
;          for ii = 0, nfilt-1 do begin
;             this = where(short[ii] eq strtrim(phot.band,2))
;             if this[0] ne -1 then begin
;                good = where(phot[this].photradius_kpc_in gt phot[this].rmin_kpc and $
;                  phot[this].photradius_kpc_out lt phot[this].rmax_kpc,comp=extrap,ncomp=nextrap)
;
;; take the measured photometry at face-value, but don't give any
;; weight to bands that have been extrapolated                
;                maggies[ii,*] = [phot[this].maggies_int,phot[this].maggies]
;                ivarmaggies[ii,*] = [phot[this].ivarmaggies_int,phot[this].ivarmaggies]
;                if nextrap ne 0 then ivarmaggies[ii,1+extrap] = 0.0 ; offset from integrated
;;               if cluster eq 'a209' and short[ii] eq 'f390w' then maggies[ii,0] = 0.0 ; crap!
;             endif
;          endfor
