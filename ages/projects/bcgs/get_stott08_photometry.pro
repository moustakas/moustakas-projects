pro get_stott08_photometry
; crossmatch Stott+08 against the VAGC

    common sdss, sdss, twomass
    
    bcgspath = ages_path(/projects)+'bcgs/'
    stott = mrdfits(bcgspath+'08stott.fits.gz',1)
    nobj = n_elements(stott)

    if (n_elements(flux_type) eq 0) then flux_type = 'model' ; = 'petro'
    if (n_elements(sdss) eq 0L) then begin
       object_sdss_imaging_file = vagc_name('object_sdss_imaging')+'.gz'
       sdss = hogg_mrdfits(object_sdss_imaging_file,1,$
         nrowchunk=75000,columns=['object_position','extinction',$
         'ra','dec',flux_type+'flux',flux_type+'flux_ivar'])
    endif

; cross-match
    spherematch, sdss.ra, sdss.dec, stott._raj2000, stott._dej2000, $
      3.0/3600.0, m1, m2
    djs_plot, sdss.ra, sdss.dec, psym=3, xsty=3, ysty=3              
    djs_oplot, stott._raj2000, stott._dej2000, psym=6, color='red'

    out = struct_addtags(stott[m2],replicate({maggies: fltarr(8), $
      ivarmaggies: fltarr(8)},n_elements(m1)))
    sdss_to_maggies, mm, ii, calib=sdss[m1], flux=flux_type
    out.maggies[0:4] = mm
    out.ivarmaggies[0:4] = ii
    
;    out = struct_addtags(stott[m2],im_struct_trimtags(sdss[m1],$
;      except=['ra','dec',flux_type+'flux',flux_type+'flux_ivar'],$
;      newtags=
;    out = struct_addtags(out,replicate({object_position: 0L},n_elements(m1)))
;    out.object_position = m1

; now read and add 2MASS photometry    
    if (n_elements(twomass) eq 0L) then begin
       if (file_test(vagc_name('object_twomass')+'.gz')) then begin
          twomass = hogg_mrdfits(vagc_name('object_twomass')+'.gz',1,$
            nrowchunk=75000, columns=['ra', $
            'decl','twomass_tag','k_m_ext','k_msig_ext', $
            'h_m_ext','h_msig_ext','j_m_ext','j_msig_ext'])
       endif
    endif

    twomass_to_maggies, twomass[m1], mm, ii
    out.maggies[5:7] = mm
    out.ivarmaggies[5:7] = ii

    im_mwrfits, out, bcgspath+'sdss_08stott.fits', /clobber
    
stop    

return
end
    
