pro write_ediscs_kplusa_sample, ancillary, ispec, write=write
; jm08may14nyu - written, based on WRITE_EDISCS_MZ_SAMPLE

    red, h100=0.7, omega0=0.3, omega_lambda=0.7
    h100 = redh100()

    lsun = 3.826D33
    hahb = 2.86 ; case B
    haoii = 1.0 ; assumed intrinsic ratio
    haconst = 7.9D-42 ; K98 conversion L(Ha) --> SFR(Ha) (Salpeter 0.1-100)

    kplusapath = ediscs_path(/projects)+'kplusa/'

    stime0 = systime(1)
    if keyword_set(write) then begin
       splogfile = kplusapath+'write_ediscs_kplusa_sample.log'
       splog, filename=splogfile
       splog, 'Log file '+splogfile+' opened '+systime()
    endif
    splog, 'IDL version: ' + string(!version,format='(99(A," "))')
    spawn, 'uname -a', uname
    splog, 'UNAME: '+uname[0]

    if keyword_set(write) then begin
       postscript = 1L
       postthick1 = 4.0
       postthick2 = 8.0
    endif else begin
       postthick1 = 1.8
       postthick2 = 2.0
       im_window, 0, xratio=0.45, /square
    endelse

; ---------------------------------------------------------------------------    
; read the data
; ---------------------------------------------------------------------------    

    splog, 'Reading the data...'
    if (n_elements(ancillary) eq 0L) then ancillary = read_ediscs(/ancillary)
    if (n_elements(ispec) eq 0L) then begin
       ispec1 = read_ediscs(/specfit)
       ispecmore = {$
         oii_lum:          -999.0, $
         hb_lum:           -999.0, $
         sfr_oii:          -999.0, $
         sfr_hb:           -999.0}
       ispecmore = replicate(ispecmore,n_elements(ispec1))
       ispec = struct_addtags(ispec1,ispecmore)
    endif
    ngalaxy = n_elements(ancillary)

; ---------------------------------------------------------------------------
; compute emission-line star-formation rates
; ---------------------------------------------------------------------------

    splog, 'Computing emission-line star-formation rates'
    dlum = dluminosity(ancillary.z,/cm)

    hbgood = where((ancillary.cflux_4861 gt -900.0) and $
      (ispec.h_beta_ew[0] gt 0.0) and (ispec.h_beta_ew[1] gt 0.0),nhbgood)
    if (nhbgood ne 0L) then begin
       hb_lum_factor = ancillary[hbgood].cflux_4861*(4.0*!dpi*dlum[hbgood]^2.0)/lsun
       ispec[hbgood].hb_lum = alog10(ispec[hbgood].h_beta_ew[0]*hb_lum_factor)
       ispec[hbgood].sfr_hb = ispec[hbgood].hb_lum+alog10(lsun)+alog10(hahb)+alog10(haconst)
    endif
       
    oiigood = where((ancillary.cflux_3727 gt -900.0) and $
      (ispec.oii_3727_ew[0] gt 0.0) and (ispec.oii_3727_ew[1] gt 0.0),noiigood)
    if (noiigood ne 0L) then begin
       oii_lum_factor = ancillary[oiigood].cflux_3727*(4.0*!dpi*dlum[oiigood]^2.0)/lsun
       ispec[oiigood].oii_lum = alog10(ispec[oiigood].oii_3727_ew[0]*oii_lum_factor)
       ispec[oiigood].sfr_oii = ispec[oiigood].oii_lum+alog10(lsun)+alog10(haoii)+alog10(haconst)
    endif
       
; ---------------------------------------------------------------------------    
; define our samples
; ---------------------------------------------------------------------------    

    help, ancillary, ispec
    help, where(ancillary.mass lt -900)
    help, where((ispec.oii_3727_ew[1] le 0.0) and (ispec.oii_3727_ew[0]/ispec.oii_3727_ew[1] le 1.0))
    help, where((ispec.lick_hd_a_model[1] le 0.0) and (ispec.lick_hd_a_model[0]/ispec.lick_hd_a_model[1] le 1.0))
    help, where((ispec.d4000_narrow_model[1] le 0.0) and (ispec.d4000_narrow_model[0]/ispec.d4000_narrow_model[1] le 1.0))

; cluster sample - all galaxies    
    cluster1 = where((ancillary.memberflag eq 1) and (ancillary.mass gt -900.0) and $
;     (ispec.oii_3727_ew[1] gt 0.0) and (ispec.oii_3727_ew[0]/ispec.oii_3727_ew[1] gt 1.0) and $
      (ispec.lick_hd_a_model[1] gt 0.0) and (ispec.lick_hd_a_model[0]/ispec.lick_hd_a_model[1] gt 1.0) and $
      (ispec.d4000_narrow_model[1] gt 0.0) and (ispec.d4000_narrow_model[0]/ispec.d4000_narrow_model[1] gt 1.0),ncluster1)
    splog, 'Cluster Sample-1: '+string(ncluster1,format='(I0)')+' galaxies'
; field sample
    field1 = where((ancillary.memberflag eq 0) and (ancillary.mass gt -900.0) and $
;     (ispec.oii_3727_ew[1] gt 0.0) and (ispec.oii_3727_ew[0]/ispec.oii_3727_ew[1] gt 1.0) and $
      (ispec.lick_hd_a_model[1] gt 0.0) and (ispec.lick_hd_a_model[0]/ispec.lick_hd_a_model[1] gt 1.0) and $
      (ispec.d4000_narrow_model[1] gt 0.0) and (ispec.d4000_narrow_model[0]/ispec.d4000_narrow_model[1] gt 1.0),nfield1)
    splog, 'Field  Sample-1: '+string(nfield1,format='(I0)')+' galaxies'

;   field1 = where((ancillary.memberflag eq 0) and (ancillary.mass gt -900.0) and $
;     (ispec.oii_3727_ew[1] gt 0.0) and (ispec.lick_hd_a[1] gt 0.0) and $
;     (ispec.d4000_narrow[1] gt 0.0) and (ancillary.ubvrijhk_absmag[1] lt -17.5),nfield1)

;   plot, ancillary[sample1].z, ancillary[sample1].ubvrijhk_absmag[1], ps=4, $
;     xr=[0.35,1.05], yr=[-15,-24], xsty=3, ysty=3, charsize=2

; ---------------------------------------------------------------------------    
; write out the various samples
; ---------------------------------------------------------------------------    
    
    if keyword_set(write) then begin

       outfile_cluster1 = kplusapath+'ediscs_kplusa_cluster1.fits'
       splog, 'Writing '+outfile_cluster1
       mwrfits, struct_addtags(ancillary[cluster1],$
         struct_trimtags(ispec[cluster1],except=['SPECFILE','GALAXY'])), $
         outfile_cluster1, /create

       outfile_field1 = kplusapath+'ediscs_kplusa_field1.fits'
       splog, 'Writing '+outfile_field1
       mwrfits, struct_addtags(ancillary[field1],$
         struct_trimtags(ispec[field1],except=['SPECFILE','GALAXY'])), $
         outfile_field1, /create

       splog, /close
       
    endif

stop

return
end
