pro build_mgfe_parent, sdss=sdss
; jm13apr14siena

    mgfepath = getenv('IM_PROJECTS_DIR')+'/ages/mgfe/'
    h100 = 0.7
    
; build the SDSS sample    
    if keyword_set(sdss) then begin
       sample = 'dr72' & letter = 'bsafe' & poststr = '0'
       post = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/postlss)
       
       keep = where((post.z gt 0.05) and (post.z lt 0.2),ngal)
       ngal = n_elements(keep)
       post = post[keep]

       vagcpath = getenv('VAGC_REDUX')+'/'
       sdssphot = mrdfits(vagcpath+'object_sdss_imaging.fits.gz',$
         row=post.object_position,1)
;      sdssspec = mrdfits(vagcpath+'object_sdss_spectro.fits.gz',$
;        row=post.object_position,1)
       galex = mrdfits(vagcpath+'object_galex_gr6.fits.gz',$
         row=post.object_position,1)
       wise = mrdfits(vagcpath+'object_wise.fits.gz',$
         row=post.object_position,1)
;      if (n_elements(sdsstwomass) eq 0L) then sdsstwomass = $
;        mrdfits(vagcpath+'object_twomass.fits.gz',$
;        row=phot.object_position,1)

; use modelmaggies for colors and scale by the cmodelmag in the r-band 
;      sdss_to_maggies, smaggies, sivarmaggies, calib=sdssphot, flux='cmodel'
       sdss_to_maggies, modelmaggies, modelivarmaggies, calib=sdssphot, flux='model'
       sdss_to_maggies, cmodelmaggies, cmodelivarmaggies, calib=sdssphot, flux='cmodel'
       neg = where(modelmaggies[2,*] le 0)
       if (neg[0] ne -1L) then message, 'Bad!'
       factor = rebin(cmodelmaggies[2,*]/modelmaggies[2,*],5,n_elements(sdssphot))
       smaggies = modelmaggies*factor
       sivarmaggies = modelivarmaggies/factor^2
       
       im_galex_to_maggies, galex, gmaggies, givarmaggies, /allow_nuv_nondetect
       wise_to_maggies, wise, wmaggies, wivarmaggies, /mpro

       parent = struct_addtags(sdssphot,struct_trimtags(post,$
         select=['z','absm','object_position']))
       parent = struct_addtags(temporary(parent),replicate($
         {maggies: fltarr(9), ivarmaggies: fltarr(9)},ngal))
       parent.maggies = [gmaggies,smaggies,wmaggies[0:1,*]]
       parent.ivarmaggies = [givarmaggies,sivarmaggies,wivarmaggies[0:1,*]]
    
       im_mwrfits, parent, mgfepath+'sdss_mgfe_parent.fits', /clobber

    endif else begin

;; build the AGES sample       
;       ppxf = read_ages(/ppxf)
;       phot = read_ages(/phot)
;       kcorr = read_ages(/kcorr)
;       
;       phot = phot[ppxf.ages_id]
;       kcorr = kcorr[ppxf.ages_id]
;
;       indx = lindgen(n_elements(phot))
;;      indx = where(ppxf.vdisp ne 165.0 and ppxf.z gt 0.05 and ppxf.z lt 0.75)
;;      indx = where(ppxf.oii_3727_ew[0] lt 1.0 and ppxf.d4000_narrow[0] gt 1.5 and $
;;        ppxf.vdisp ne 165.0 and ppxf.z gt 0.05 and ppxf.z lt 0.7)
;;      qaplot_ages_gandalf_specfit, ppxf[indx[0:10]], ss, psfile='~/junk.ps'
;
;       specfit = read_ages_gandalf_specfit(ppxf[indx])
;
;; write out
;;      im_mwrfits, specfit, mgfepath+'ages_mgfe_specfit.fits', /clobber ; too big!
;       im_mwrfits, specfit, '~/ages_mgfe_specfit.fits', /clobber
;       im_mwrfits, ppxf[indx], mgfepath+'ages_mgfe_ppxf.fits', /clobber
;       im_mwrfits, kcorr[indx], mgfepath+'ages_mgfe_kcorr.fits', /clobber

; area = 8.013 deg^2 (i.e., print, ages_survey_area());
; total volume = jhnvol(0.05,0.75)*8.0135*3600.0 = 1.51329E+07 Mpc^3
; at z=0.05-0.15 the total volume is = jhnvol(0.05,0.15)*8.0135*3600.0
; = 1.86578E+05 Mpc^3

       phot = read_ages(/photo)
       ppxf = read_ages(/ppxf)

; identify main sample galaxies and apply additional sample cuts:
; 0.05<z<0.75; 15<Itot<19.95; reject unfluxed plates and require that
; the spectrum has been fitted by PPXF
       ifaint = mz_ifaint(ibright=ibright,select_filter=select_filter,$ 
         select_vega2ab=iv2ab,/vega)
       ages_to_maggies, phot, allmaggies, allivarmaggies, filterlist=filterlist, bands=bands
       keep = where(strmatch(filterlist,'*ch3*') eq 0 and strmatch(filterlist,'*ch4*') eq 0,nfilt)
       filterlist = filterlist[keep]
       allmaggies = allmaggies[keep,*]
       allivarmaggies = allivarmaggies[keep,*]
       
       sample_zmin = 0.05
       sample_zmax = 0.75
       nminphot = 3
       rejplates = [$
         104,$ ; average fluxing vector
         106,110,209,310,311,$ ; unfluxed
         604]  ; very crummy 
       
       ippxf = intarr(n_elements(phot))
       ippxf[ppxf.ages_id] = 1 ; fitted

       iparent = (phot.imain eq 1) and (phot.i_tot le ifaint) and (phot.i_tot ge ibright) and $
         (ippxf eq 1) and (phot.z ge sample_zmin) and (phot.z le sample_zmax) and $
         (total(allmaggies gt 0,1) ge nminphot) and $
         (phot.pass ne rejplates[0]) and (phot.pass ne rejplates[1]) and $
         (phot.pass ne rejplates[2]) and (phot.pass ne rejplates[3]) and $
         (phot.pass ne rejplates[4]) and (phot.pass ne rejplates[5]) and $
         (phot.pass ne rejplates[6])
       parent = where(iparent,ngal)
       splog, 'Ngal = ', ngal

       ppxf_sample = ppxf[parent]
       sample = struct_trimtags(phot[parent],except=bands+'_*')
       maggies = allmaggies[*,parent]
       ivarmaggies = allivarmaggies[*,parent]

;      ngal = 1001
;      sample = sample[0:1000]  ; TEST!!
;      maggies = maggies[*,0:1000]
;      ivarmaggies = ivarmaggies[*,0:1000]
       
; in addition to the corrections for spectroscopic incompleteness
; (~2.1%), fiber incompleteness (~4.3%), and sparse-sampling, we also
; need to account for the catalog incompleteness (catalog_weight) and
; for the plates that were rejected above (field_weight); from
; Eisenstein, the catalog incompleteness is ~3%-6%, so assume an
; average value of 4%; finally, compute the final galaxy weight as the
; product of all the various selection terms        
       moretags = replicate({catalog_weight: 1.04, field_weight: 1.0, $
         final_weight: 1.0, maggies: fltarr(nfilt), ivarmaggies: fltarr(nfilt), $
         k_chi2: 0.0, k_mass: 0.0, k_coeffs: fltarr(5), k_absmag: fltarr(5), $
         k_ivarabsmag: fltarr(5), k_kcorrect: fltarr(5)},ngal)
       sample = struct_addtags(temporary(sample),moretags)
       sample.maggies = maggies
       sample.ivarmaggies = ivarmaggies
       
       field_weight = ages_upweight_field(rejplates,field=field)
       for ii = 0, n_elements(field)-1 do begin
          these = where(field[ii] eq sample.field,nthese)
          if (nthese ne 0) then sample[these].field_weight = field_weight[ii]
       endfor
       sample.final_weight = sample.spec_weight*sample.target_weight*$
         sample.fiber_weight*sample.catalog_weight*1.0/sample.field_weight ; note 1/FIELD_WEIGHT!

; do a conservative cut on u-r color to select out the candidate
; quiescent galaxies (need to compute K-corrections)
       splog, 'Computing K-corrections'
       kcorr = im_kcorrect(sample.z,sample.maggies,sample.ivarmaggies,$
         filterlist,sdss_filterlist(),band_shift=0.1,chi2=chi2,$
         mass=mass,coeffs=coeffs,bestmaggies=bestmaggies,absmag=absmag,$
         ivarabsmag=ivarabsmag,/silent,h100=h100) ; AB, band_shift=0.1
       sample.k_chi2 = chi2
       sample.k_mass = mass
       sample.k_coeffs = coeffs
       sample.k_absmag = absmag
       sample.k_ivarabsmag = ivarabsmag
       sample.k_kcorrect = kcorr
       
;      qq = where(absmag[0,*]-absmag[1,*] gt 2.0,nqq)
;      splog, 'Nqq ', nqq
;      sample = sample[qq]

;;      qaplot_ages_gandalf_specfit, ppxf_sample[0:10], ss, psfile='~/junk.ps'
;       specfit = read_ages_gandalf_specfit(ppxf_sample)
              
; write out
       im_mwrfits, sample, mgfepath+'ages_mgfe_parent.fits', /clobber
       im_mwrfits, ppxf, mgfepath+'ages_mgfe_ppxf.fits', /clobber
;;     im_mwrfits, specfit, mgfepath+'ages_mgfe_specfit.fits', /clobber ; too big!
;      im_mwrfits, specfit, '~/ages_mgfe_specfit.fits', /clobber

    endelse 

return
end
