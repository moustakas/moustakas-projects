pro get_sdss_bcgs_photometry
; jm10sep28ucsd - crossmatch our local/SDSS sample against the VAGC 

    common sdss, sdss, twomass;, specz
    
    bcgspath = ages_path(/projects)+'bcgs/'
    data = rsex('bcgs_sdss_sample.sex')
;   data = rsex('sdss_sample_v1.sex')
    nobj = n_elements(data)

    if (n_elements(sdss) eq 0L) then sdss = $
      hogg_mrdfits(vagc_name('object_sdss_imaging')+'.gz',1,nrowchunk=75000)
;   if (n_elements(specz) eq 0L) then specz = $
;     hogg_mrdfits(vagc_name('object_sdss_spectro')+'.gz',1,nrowchunk=75000)
    if (n_elements(twomass) eq 0L) then begin
       twomass = hogg_mrdfits(vagc_name('object_twomass')+'.gz',1,$
         nrowchunk=75000, columns=['ra', 'decl','twomass_tag',$
         'k_m_ext','k_msig_ext', 'h_m_ext','h_msig_ext',$
         'j_m_ext','j_msig_ext'])
    endif

; cross-match
    spherematch, sdss.ra, sdss.dec, data.ra, data.dec, $
      5.0/3600.0, m1, m2 & help, m2, data
    missing = lindgen(nobj)
    remove, m2, missing
    struct_print, data[missing]

;   djs_plot, sdss.ra, sdss.dec, psym=3, xsty=3, ysty=3              
;   djs_oplot, data.ra, data.dec, psym=6, color='red'

; use the modelmaggies, but scale the r-band flux to the Petrosian
; flux     
    out = struct_addtags(data[m2],replicate({maggies: fltarr(8), $
      ivarmaggies: fltarr(8)},n_elements(m1)))
    sdss_to_maggies, modelmaggies, imodelmaggies, calib=sdss[m1], flux='model'
    sdss_to_maggies, petromaggies, ipetromaggies, calib=sdss[m1], flux='petro'
    factor = rebin(petromaggies[2,*]/modelmaggies[2,*],5,n_elements(m1))
    maggies = modelmaggies*factor
    imaggies = imodelmaggies/factor^2
    
    out.maggies[0:4] = maggies
    out.ivarmaggies[0:4] = imaggies
    
;    out = struct_addtags(data[m2],im_struct_trimtags(sdss[m1],$
;      except=['ra','dec',flux_type+'flux',flux_type+'flux_ivar'],$
;      newtags=
;    out = struct_addtags(out,replicate({object_position: 0L},n_elements(m1)))
;    out.object_position = m1

    twomass_to_maggies, twomass[m1], mm, ii
    out.maggies[5:7] = mm
    out.ivarmaggies[5:7] = ii

    im_mwrfits, out, bcgspath+'bcgs_sdss_phot.fits', /clobber

return
end
    
