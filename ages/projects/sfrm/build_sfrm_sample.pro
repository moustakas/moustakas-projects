function sfrm_kcorrect, in_redshift, in_maggies, $
  in_ivarmaggies, filterlist=filterlist
; compute additional K-corrections for both the AGES and SDSS samples 

    ngal = n_elements(in_redshift)
    nfilt = n_elements(filterlist)
    
    kcorr = {$
      k_maggies:                fltarr(nfilt), $
      k_ivarmaggies:            fltarr(nfilt), $
      k_bestmaggies:            fltarr(nfilt), $
      k_mass:                          -999.0, $
      k_coeffs:                     fltarr(5), $
      k_chi2:                          -999.0, $
      k_cflux:                fltarr(5)-999.0, $
      k_uvflux:               fltarr(2)-999.0, $

      k_galex_absmag_00:         fltarr(2)-999.0,$
      k_galex_absmag_ivar_00:    fltarr(2)-999.0,$
      k_galex_kcorrect_00:       fltarr(2)-999.0,$
      k_galex_absmag_01:         fltarr(2)-999.0,$
      k_galex_absmag_ivar_01:    fltarr(2)-999.0,$
      k_galex_kcorrect_01:       fltarr(2)-999.0,$

      k_ugriz_absmag_00:         fltarr(5)-999.0,$
      k_ugriz_absmag_ivar_00:    fltarr(5)-999.0,$
      k_ugriz_kcorrect_00:       fltarr(5)-999.0,$
      k_ugriz_absmag_01:         fltarr(5)-999.0,$
      k_ugriz_absmag_ivar_01:    fltarr(5)-999.0,$
      k_ugriz_kcorrect_01:       fltarr(5)-999.0,$

      k_ubvrijhk_absmag_00:      fltarr(8)-999.0,$
      k_ubvrijhk_absmag_ivar_00: fltarr(8)-999.0,$
      k_ubvrijhk_kcorrect_00:    fltarr(8)-999.0,$
      k_ubvrijhk_absmag_01:      fltarr(8)-999.0,$
      k_ubvrijhk_absmag_ivar_01: fltarr(8)-999.0,$
      k_ubvrijhk_kcorrect_01:    fltarr(8)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.k_maggies = in_maggies
    kcorr.k_ivarmaggies = in_ivarmaggies
    
; compute k-corrections    
    splog, 'Computing ugriz K-corrections'
    ugriz_kcorrect_00 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,sdss_filterlist(),band_shift=0.0,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=ugriz_absmag_00,$
      ivarabsmag=ugriz_absmag_ivar_00,clineflux=cflux,uvflux=uvflux,$
      /silent,vname=vname) ; AB, band_shift=0.0
    ugriz_kcorrect_01 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,sdss_filterlist(),band_shift=0.1,absmag=ugriz_absmag_01,$
      ivarabsmag=ugriz_absmag_ivar_01,/silent,vname=vname) ; AB, band_shift=0.1

    splog, 'Computing NUV,FUV K-corrections'
    galex_kcorrect_00 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,galex_filterlist(),band_shift=0.0,absmag=galex_absmag_00,$
      ivarabsmag=galex_absmag_ivar_00,/silent,vname=vname) ; AB, band_shift=0.0
    galex_kcorrect_01 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,galex_filterlist(),band_shift=0.1,absmag=galex_absmag_01,$
      ivarabsmag=galex_absmag_ivar_01,/silent,vname=vname) ; AB, band_shift=0.1

    splog, 'Computing UBVRIJHKs K-corrections'
    kcorrect_00 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,[bessell_filterlist(),twomass_filterlist()],band_shift=0.0,$
      absmag=absmag_00,ivarabsmag=absmag_ivar_00,/silent,vname=vname) ; AB, band_shift=0.0
    kcorrect_01 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,[bessell_filterlist(),twomass_filterlist()],band_shift=0.1,$
      absmag=absmag_01,ivarabsmag=absmag_ivar_01,/silent,vname=vname) ; AB, band_shift=0.1

    kcorr.k_bestmaggies = bestmaggies
    kcorr.k_mass = alog10(mass) ; h=0.7, Chabrier
    kcorr.k_coeffs = coeffs
    kcorr.k_chi2 = chi2
    kcorr.k_cflux = cflux
    kcorr.k_uvflux = uvflux

    kcorr.k_galex_absmag_00         = galex_absmag_00
    kcorr.k_galex_absmag_ivar_00    = galex_absmag_ivar_00
    kcorr.k_galex_kcorrect_00       = galex_kcorrect_00
    kcorr.k_galex_absmag_01         = galex_absmag_01
    kcorr.k_galex_absmag_ivar_01    = galex_absmag_ivar_01
    kcorr.k_galex_kcorrect_01       = galex_kcorrect_01

    kcorr.k_ugriz_absmag_00         = ugriz_absmag_00
    kcorr.k_ugriz_absmag_ivar_00    = ugriz_absmag_ivar_00
    kcorr.k_ugriz_kcorrect_00       = ugriz_kcorrect_00
    kcorr.k_ugriz_absmag_01         = ugriz_absmag_01
    kcorr.k_ugriz_absmag_ivar_01    = ugriz_absmag_ivar_01
    kcorr.k_ugriz_kcorrect_01       = ugriz_kcorrect_01

    kcorr.k_ubvrijhk_absmag_00      = absmag_00
    kcorr.k_ubvrijhk_absmag_ivar_00 = absmag_ivar_00
    kcorr.k_ubvrijhk_kcorrect_00    = kcorrect_00
    kcorr.k_ubvrijhk_absmag_01      = absmag_01
    kcorr.k_ubvrijhk_absmag_ivar_01 = absmag_ivar_01
    kcorr.k_ubvrijhk_kcorrect_01    = kcorrect_01

return, kcorr
end

function ages_sfrm_masses, ages_id
; compute a fiducial set of stellar masses and errors for the
;  AGES/SFRM sample

    isedpath = ages_path(/isedfit)

; define the fiducial set of models and stellar masses (SFHGRID02,
; CALZETTI), and the other models
    prefix = 'BwRIJHKsirac_chab_'
    fiducial = prefix+'calzetti_sfhgrid02' ; fiducial models

    moremodels = prefix+[['odonnell','smc']+'_sfhgrid02',$
      ['calzetti','odonnell','smc']+'_sfhgrid01','sfhgrid03']
    nmodels = n_elements(moremodels)

    ngal = n_elements(ages_id)
    result = {$
      mass:              0.0,$  ; fiducial parameters
      mass_err:          0.0,$
      tau:               0.0,$
      tau_err:           0.0,$
      ebv:               0.0,$
      ebv_err:           0.0,$

      model:         moremodels,$
      model_mass:     fltarr(nmodels),$
      model_mass_err: fltarr(nmodels),$
      model_tau:      fltarr(nmodels),$
      model_tau_err:  fltarr(nmodels),$
      model_ebv:      fltarr(nmodels),$
      model_ebv_err:  fltarr(nmodels)}
    result = replicate(result,ngal)

; read the fiducial models
    ised = mrdfits(isedpath+fiducial+'.fits.gz',1,rows=ages_id)
    rest = mrdfits(isedpath+fiducial+'_measure.fits.gz',1,rows=ages_id)
    
    result.mass     = ised.mass50
    result.mass_err = ised.mass_err
    result.tau      = ised.tau50
    result.tau_err  = ised.tau_err
    result.ebv      = ised.ebv50
    result.ebv_err  = ised.ebv_err
    result = struct_addtags(result,rest)
;   result = struct_addtags(result,struct_trimtags(rest,$
;     except=['isedfit_id','zobj','*maggies*']))

; now fill in the additional model results    
    for ii = 0, nmodels-1 do begin
       ised = mrdfits(isedpath+moremodels[ii]+'.fits.gz',1,rows=ages_id)
;      rest = mrdfits(isedpath+moremodels[ii]+'_measure.fits.gz',1,rows=ages_id)
       result.model_mass[ii]     = ised.mass50
       result.model_mass_err[ii] = ised.mass_err
       result.model_tau[ii]      = ised.tau50
       result.model_tau_err[ii]  = ised.tau_err
       result.model_ebv[ii]      = ised.ebv50
       result.model_ebv_err[ii]  = ised.ebv_err
    endfor

return, result
end

pro build_sfrm_sample, sample, sdss=sdss
; jm09nov30ucsd - build the AGES/SFRM parent samples

    common sfrm_sample, sdssphot, sdssgalex, sdsstwomass

; select the parent sample    
    agespath = ages_path(/analysis)
    sfrmpath = ages_path(/projects)+'sfrm/'
    vagcpath = getenv('VAGC_REDUX')+'/'
    h100 = 0.7

    if keyword_set(sdss) then begin
       splog, '#########################'
       splog, 'Building the SDSS comparison sample'

; POSTSTR/25: -17>Mr>-24; 0.01<z<0.25
       sample = 'dr72'
       letter = 'bsafe'
       poststr = '25'
       
       post = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/postlss)
       massoh = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/mpamassoh)
       vmax_noevol = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/vmax_noevol)
       vmax_evol = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/vmax_evol)
       ngal = n_elements(post)

       if (n_elements(sdssphot) eq 0L) then sdssphot = $
         mrdfits(vagcpath+'object_sdss_imaging.fits.gz',$
         row=post.object_position,1)
       if (n_elements(sdssgalex) eq 0L) then sdssgalex = $
         mrdfits(vagcpath+'object_galex_gr45.fits.gz',$
         row=post.object_position,1)
       if (n_elements(sdsstwomass) eq 0L) then sdsstwomass = $
         mrdfits(vagcpath+'object_twomass.fits.gz',$
         row=post.object_position,1)

; build the final catalog by keeping just the tags we want
       sample = struct_trimtags(post,select=['object_position','ra','dec','z'])
       sample = struct_addtags(temporary(sample),struct_trimtags(massoh,$
         except=['ra','dec','z','oh_*']))

       sample = struct_addtags(temporary(sample),struct_trimtags(sdssphot,$
         select=['modelflux*','petroflux*','extinction','petror*']))

       sample = struct_addtags(temporary(sample),im_struct_trimtags(vmax_noevol,$
         select=['zmin_local','zmax_local','vmax'],$
         newtags=['zmin','zmax','vmax']+'_noevol'))
       sample = struct_addtags(temporary(sample),im_struct_trimtags(vmax_evol,$
         select=['zmin_local','zmax_local','vmax'],$
         newtags=['zmin','zmax','vmax']+'_evol'))
       sample.vmax_noevol = sample.vmax_noevol/h100^3.0 ; h=1-->h=0.7
       sample.vmax_evol = sample.vmax_evol/h100^3.0 ; h=1-->h=0.7

       sample = struct_addtags(temporary(sample),$
         struct_trimtags(sdssgalex,select=['nuv_*','fuv_*']))

; compute K-corrections
       splog, 'Computing K-corrections'
       im_galex_to_maggies, sdssgalex, gmaggies, givarmaggies
       sdss_to_maggies, smaggies, sivarmaggies, calib=sdssphot
       twomass_to_maggies, sdsstwomass, tmaggies, tivarmaggies
       maggies = [gmaggies,smaggies,tmaggies]
       ivarmaggies = [givarmaggies,sivarmaggies,tivarmaggies]

       filterlist = [galex_filterlist(),sdss_filterlist(),twomass_filterlist()]
       kcorr = sfrm_kcorrect(sample.z,maggies,ivarmaggies,filterlist=filterlist)
       sample = struct_addtags(temporary(sample),kcorr)

;      bmass = reform(k_sdss_bell(vagc.absmag[0:4])) ; stellar mass
;      bmass = alog10(bmass) - alog10(h100^2)    ; h=1-->0.7

       im_mwrfits, sample, sfrmpath+'sdss_sfrm.fits', /clobber

    endif else begin
; read the parent AGES photometry
       splog, '#########################'
       splog, 'Building the AGES sample'

       agesphot = read_ages(/phot)

; basic parameters of the selection    
       sample_zmin = 0.05
       sample_zmax = 0.75

       vname = 'default.nolines'
       ifilt = 'ndwfs_I.par'
       ibright = 15.0           ; [Vega]
       ifaint = 19.95           ; [Vega]

; define the parent sample
       splog, 'Converting to maggies'
       ages_to_maggies, agesphot, maggies, ivarmaggies, $
         filterlist=filterlist, nozpoffset=nozpoffset
       nfilt = n_elements(filterlist)
       bwindx = where(strmatch(filterlist,'*_Bw*'))
       rindx = where(strmatch(filterlist,'*_R*'))
       iindx = where(strmatch(filterlist,'*_I*'))

       splog, 'Applying the window function'
       window = ages_isin_window(agesphot.ra,agesphot.dec) ; within the 2004 survey region
       istar = (agesphot.psfflux[2]*1E-9 gt 10^(-0.4*19.0)) and $ ; SDSS type=star and r<19
         (agesphot.sdss_star eq 1)

       parent = where($
         (agesphot.z ge sample_zmin) and (agesphot.z le sample_zmax) and $ ; redshift cuts
         (agesphot.gshort gt 0) and (agesphot.gshort ne 2048) and $       ; main galaxy
         (agesphot.field ge 1) and (agesphot.field le 15) and $           ; field assignment cut
;        (agesphot.i_auto le ifaint) and (agesphot.i_auto ge ibright) and $ ; 15<I<19.95
         (agesphot.i_tot le ifaint) and (agesphot.i_tot ge ibright) and $ ; 15<I<19.95
         (istar eq 0) and (window eq 1) and (ivarmaggies[bwindx,*] gt 0.0) and $
         (ivarmaggies[rindx,*] gt 0.0) and (ivarmaggies[iindx,*] gt 0.0),ngal)
       splog, 'Ngal = ', ngal

       sample = agesphot[parent]
       in_maggies = maggies[*,parent]
       in_ivarmaggies = ivarmaggies[*,parent]

       toss = where(strmatch(filterlist,'*_ch[2-4].par*',/fold),ntoss)
       if (ntoss ne 0) then in_ivarmaggies[toss,*] = 0.0
       
;; read and store the isedfit stellar masses
;       splog, 'Adding iSEDfit stellar masses'
;       ised = ages_sfrm_masses(sample.ages_id)
;       sample = struct_addtags(temporary(sample),ised)
       
; compute additional K-corrections; use the same photometry used in
; ISEDFIT_MEASURE 
       splog, 'Computing additional K-corrections'
       kcorr = sfrm_kcorrect(sample.z,in_maggies,$
         in_ivarmaggies,filterlist=filterlist)
       sample = struct_addtags(temporary(sample),kcorr)
       
; compute additional ancillary quantities that we are going to need:
; ZMIN and ZMAX with and without evolution, and the final galaxy
; weight as the product of all the various selection terms
       splog, 'Computing zmin/zmax'
       moretags = replicate({zmin_noevol: 0.0, zmax_noevol: 0.0, $
         zmin_evol: 0.0, zmax_evol: 0.0, final_weight: 0.0},ngal)
       sample = struct_addtags(temporary(sample),moretags)
       sample.final_weight = sample.spec_weight*sample.target_weight*$
         sample.fiber_weight*sample.catalog_weight

       q0 = 0.0 & q1 = 0.0 & qz0 = 0.0 ; no evolution
       zminzmax_noevol = im_zminzmax(sample.z,sample.i_tot,$
         sample.k_coeffs,bright=ibright,faint=ifaint,$
         filter=ifilt,vname=vname,q0=q0,q1=q1,qz0=qz0)

       q0 = 1.6 & q1 = 0.0 & qz0 = 0.1 ; following Eisenstein+09
       zminzmax_evol = im_zminzmax(sample.z,sample.i_tot,$
         sample.k_coeffs,bright=ibright,faint=ifaint,$
         filter=ifilt,vname=vname,q0=q0,q1=q1,qz0=qz0)

       sample.zmin_noevol = zminzmax_noevol.zmin
       sample.zmax_noevol = zminzmax_noevol.zmax
       sample.zmin_evol = zminzmax_evol.zmin
       sample.zmax_evol = zminzmax_evol.zmax

; write out       
       im_mwrfits, sample, sfrmpath+'ages_sfrm.fits', /clobber

    endelse

stop    
    
return
end

;; read the NFRIENDS file (see AGES_FIND_ISOLATED)       
;       friendsfile = agespath+'ages_nfriends.fits.gz'
;       splog, 'Reading '+friendsfile
;       friends = mrdfits(friendsfile,1,rows=sample.ages_id)
;       sample = struct_addtags(temporary(sample),$
;         struct_trimtags(friends,select='nfriends'))

