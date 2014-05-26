function deep2_get_oiiflux, ppxf, cflux_3727_rest=cflux_3727_rest
; jm14mar11siena - given a pPXF-style emission-line catalog for DEEP2
; and the observed continuum flux at 3727 A, return a data structure
; which contains high-level [OII] stuff where, in particular, the
; fluxes have been converted to physical units
;
; recall that all the measured quantities are in the rest-frame,
; including OII_3727_EW [A], OII_3727_AMP [ADU], OII_3727_CONTINUUM
; [ADU]
;
; the integrated flux is redshift-independent (other than the distance
; modulus dependence) OII_3727 [ADU]

    ngal = n_elements(ppxf)
    if n_elements(snrcut) eq 0 then snrcut = 1.5

    oii = struct_trimtags(ppxf,select=['oii_3727*'],except=$
      ['*wave*','*continuum*','*linez*','*sigma*'])
    oii = struct_addtags(im_struct_trimtags(ppxf,select='sigma_forbidden',$
      newtags='sigma_kms'),temporary(oii))

    zobj = ppxf.z
    cflux_3727_obs = cflux_3727_rest/(1.0+zobj) ; [erg/s/cm2/A]
    
; [OII] 3726    
    good = where(oii.oii_3727_1[1] gt 0,ngood,comp=crap,ncomp=ncrap)
    if ngood ne 0L then begin
       Fc_3727_rest = cflux_3727_rest[good] ; *rest-frame* [erg/s/cm2/A]
       Fc_3727_obs = cflux_3727_obs[good]   ; *observed-frame* [erg/s/cm2/A]

       oii[good].oii_3727_1 = oii[good].oii_3727_1_ew*rebin(reform(Fc_3727_rest,1,ngood),2,ngood) ; [erg/s/cm2]
       oii[good].oii_3727_1_limit = oii[good].oii_3727_1_ew_limit*Fc_3727_rest
       
       oii[good].oii_3727_1_amp = oii[good].oii_3727_1_amp*$ ; *observed* [erg/s/cm2/A]
         rebin(reform(Fc_3727_obs,1,ngood),2,ngood)/$ 
         rebin(reform(ppxf[good].oii_3727_1_continuum[0],1,ngood),2,ngood)
    endif

; [OII] 3729
    good = where(oii.oii_3727_2[1] gt 0,ngood,comp=crap,ncomp=ncrap)
    if ngood ne 0L then begin
       Fc_3727_rest = cflux_3727_rest[good] ; *rest-frame* [erg/s/cm2/A]
       Fc_3727_obs = cflux_3727_obs[good]   ; *observed-frame* [erg/s/cm2/A]

       oii[good].oii_3727_2 = oii[good].oii_3727_2_ew*rebin(reform(Fc_3727_rest,1,ngood),2,ngood) ; [erg/s/cm2]
       oii[good].oii_3727_2_limit = oii[good].oii_3727_2_ew_limit*Fc_3727_rest
       
       oii[good].oii_3727_2_amp = oii[good].oii_3727_2_amp*$ ; *observed* [erg/s/cm2/A]
         rebin(reform(Fc_3727_obs,1,ngood),2,ngood)/$ 
         rebin(reform(ppxf[good].oii_3727_2_continuum[0],1,ngood),2,ngood)
    endif

; [OII] 3726,3729
    good = where(oii.oii_3727[1] gt 0,ngood,comp=crap,ncomp=ncrap)
    if ngood ne 0L then begin
       Fc_3727_rest = cflux_3727_rest[good] ; *rest-frame* [erg/s/cm2/A]
       Fc_3727_obs = cflux_3727_obs[good]   ; *observed-frame* [erg/s/cm2/A]

       oii[good].oii_3727 = oii[good].oii_3727_ew*rebin(reform(Fc_3727_rest,1,ngood),2,ngood) ; [erg/s/cm2]
       oii[good].oii_3727_limit = oii[good].oii_3727_ew_limit*Fc_3727_rest
       
       oii[good].oii_3727_amp = oii[good].oii_3727_amp*$ ; *observed* [erg/s/cm2/A]
         rebin(reform(Fc_3727_obs,1,ngood),2,ngood)/$ 
         rebin(reform(ppxf[good].oii_3727_continuum[0],1,ngood),2,ngood)
    endif

return, oii
end
