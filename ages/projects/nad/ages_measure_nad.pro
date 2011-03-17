pro ages_measure_nad, write=write, debug=debug
; jm06mar06uofa - measure the NaD absorption-line velocity in the AGES
;                 spectra

    rootpath = ages_path(/rawdata)
    analysis_path = ages_path(/analysis)
    nadpath = ages_path(/projects)+'nad/'
    
    vredux_11 = '0300'
    vredux_20 = '0051'
    vredux = [vredux_11,vredux_20]

    datapath_11 = rootpath+vredux_11+'/'
    datapath_20 = rootpath+vredux_20+'/'
    datapath = [datapath_11,datapath_20]

    nfiber = 300L
    fwhmres = 6.0 ; instrumental resolution [Angstrom]

    snrcut1 = 5.0
    snrcut2 = 3.0

    nadwave = 5890.0 ; [Angstroms]
    nadwidth = 35.0  ; [Angstroms]
    
; initialize the output structure

    result_template = {$
      zmerge_catid:              0L,$
      zmerge_ra:               0.0D,$
      zmerge_dec:              0.0D,$
      pass:                       0,$
      aper:                       0,$
      z:                        0.0,$
      nad_ew:             [0.0,0.0],$ ; EW(NaD)
      nad_center_gauss: [0.0D,0.0D],$ ; Gaussian-weighted center
      nad_center_flux:  [0.0D,0.0D]}  ; flux-weighted center
    result1 = result_template

; read the ZMERGE redshift catalog
    
    zmerge = mrdfits(analysis_path+'catalog.zmerge.fits.gz',1,/silent)

    t0 = systime(1)
    for ipath = 0L, n_elements(datapath)-1L do begin
    
       sphectfits = file_basename(file_search(datapath[ipath]+'spHect-???-'+vredux[ipath]+'.fits',count=npassfield))
       spzbestfits = file_basename(file_search(datapath[ipath]+'spZbest-???-'+vredux[ipath]+'.fits'))
       passfield = strmid(sphectfits,7,3) ; NOT GENERAL

       for jpass = 0L, npassfield-1L do begin

          splog, 'Unpacking PASSFIELD '+passfield[jpass]+'.'
          pass = long(passfield[jpass])
          
; read the data          
          
          bigwave   = mrdfits(datapath[ipath]+sphectfits[jpass],0,/silent)
          bigflux   = mrdfits(datapath[ipath]+sphectfits[jpass],1,/silent)
          biginvvar = mrdfits(datapath[ipath]+sphectfits[jpass],2,/silent)
    
; loop on each fiber

          for ifiber = 0L, nfiber-1L do begin

             print, format='("Fiber ",I0,"/",I0,".",A1,$)', ifiber+1, nfiber, string(13b)
             
; is this object in the ZMERGE catalog?

             zmatch = where((pass eq zmerge.zmerge_pass) and (ifiber+1L eq zmerge.zmerge_aper),nzmatch)
             if (nzmatch ne 0L) then begin

                if (nzmatch gt 1L) then message, 'This should not happen.'

                wave = reform(bigwave[*,ifiber])/(1+zmerge[zmatch].zmerge_z)
                flux = reform(bigflux[*,ifiber])
                invvar = reform(biginvvar[*,ifiber])

                if (min(wave) lt (nadwave-2*nadwidth)) and (max(wave) gt (nadwave+2*nadwidth)) then begin

                   flux = 1D-17*djs_maskinterp(flux,(invvar le 0.0),wave,/const)
                   invvar = djs_maskinterp(invvar,(invvar le 0.0),wave,/const)
                   ferr = 1D-17*1.0/sqrt(invvar)

;                  nadindx = where((wave gt nadwave-3.0*fwhmres) and (wave lt nadwave+3.0*fwhmres),nnadindx)
;                  nad_signal = total(flux[nadindx]/ferr[nadindx]^2)/total(1.0/ferr[nadindx]^2)
;                  nad_noise = 1.0/sqrt(total(1.0/ferr[nadindx]^2))
;                  nad_snr = nad_signal/nad_noise
;                  nad_snr = total(flux[nadindx])/sqrt(total(ferr[nadindx]^2))
                   
                   indices = spectral_indices(wave,flux,ferr=ferr,indexpath=nadpath,indexfile='nadindex.txt',/silent)

                   if (abs(indices.lick_nad[0]/indices.lick_nad[1]) gt snrcut1) then begin
                   
                      nad = im_splotew(wave,flux,ferr,nadwave,boxwidth=nadwidth,/absorption,doplot=debug,/silent)
                      if keyword_set(debug) then cc = get_kbrd(1)

                      if (abs(nad.ew[0]/nad.ew[1]) gt snrcut2) then begin
                         result1.nad_ew           = nad.ew
                         result1.nad_center_gauss = nad.linecenter
                         result1.nad_center_flux  = nad.linecenter_flux
                         result1.zmerge_catid     = zmerge[zmatch].zmerge_catid
                         result1.zmerge_ra        = zmerge[zmatch].zmerge_ra
                         result1.zmerge_dec       = zmerge[zmatch].zmerge_dec
                         result1.z                = zmerge[zmatch].zmerge_z
                         if (n_elements(result) eq 0L) then result = result1 else result = struct_append(result,result1)
                         result1 = result_template ; zero out
                      endif

                   endif
                endif
             endif 
          endfor                ; close FIBER loop
       endfor                   ; close PASSFIELD loop
    endfor                      ; close DATAPATH loop
    splog, format='("Total time to measure NaD = ",G0," minutes.")', (systime(1)-t0)/60.0

; write out    

    if keyword_set(write) then begin
       nadfile = 'ages_nad.fits'
       splog, 'Writing '+nadpath+nadfile+'.'
       mwrfits, result, nadpath+nadfile, /create
       spawn, ['gzip -f '+nadpath+nadfile], /sh
    endif

stop    
    
return
end
    
