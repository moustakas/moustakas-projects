pro ages_reclassify_bad_redshifts, info, write=write
; jm06feb05uofa - parse the output from AGES_RERUN_SPZBEST

    dzthresh = 0.02
    
    analysis_path = ages_path(/analysis)
    datapath = ages_path()+'spzbest_reclassify/'
    if (n_elements(info) eq 0L) then info = mrdfits(analysis_path+'ages_flag_bad_redshifts.fits.gz',1,/silent)

    allplatefiles = info[uniq(info.platefile,sort(info.platefile))].platefile
    nplate = n_elements(allplatefiles)

;   for i = 0L, 4L do begin
    for i = 0L, nplate-1L do begin

       match = where(allplatefiles[i] eq info.platefile,nmatch)
       info1 = struct_addtags(info[match],replicate({reclassified: 0L, zans_indx: 0L},nmatch))
       
       if (info1[0].pass lt 401) then rerun = '0300' else rerun = '0051'
       
       platefile = info1[0].platefile
       pass = strmid(platefile,7,1)
       field = strmid(platefile,8,2)
       fiberid = info1.aper

; see HS_REDUCE1D       
       
       platemjd = string(pass, format='(i1.1)') + string(field, format='(i2.2)')
       outkey = platemjd + '-' + rerun
       zallfile = 'spZall-'+outkey+'.fits'
       zbestfile = 'spZbest-'+outkey+'.fits'

       res_all = reform(mrdfits(datapath+zallfile,1,/silent),51L,nmatch)
       zans = mrdfits(datapath+zbestfile,1,/silent)

; loop on each fiber

       for ifiber = 0L, nmatch-1L do begin

          dzall = abs(info1[ifiber].z-res_all[*,ifiber].z)
          minz = min(dzall,mindx)
          dz = dzall[mindx]

          print, mindx, info1[ifiber].z, res_all[mindx,ifiber].z, dz, dz/(1+info1[ifiber].z) gt dzthresh

; only retain re-classified galaxies if dz/(1+z)<DZTHRESH; remove the
; remaining objects from further analysis
          
          if (dz/(1+info1[ifiber].z) lt dzthresh) then begin
             zans_reclass1 = res_all[mindx,ifiber]
             info1[ifiber].reclassified = 1L
             info1[ifiber].zans_indx = mindx
             info1[ifiber].newclass = zans_reclass1.class
             info1[ifiber].z_spzbest = zans_reclass1.z
             if (n_elements(outinfo) eq 0L) then outinfo = info1[ifiber] else $
               outinfo = struct_append(outinfo,info1[ifiber])
;            zans_reclass1 = struct_addtags(info1[ifiber],res_all[mindx,ifiber])
;            if (n_elements(zans_reclass) eq 0L) then zans_reclass = zans_reclass1 else $
;              zans_reclass = struct_append(zans_reclass,zans_reclass1)
          endif
          
       endfor 

    endfor

    nall = n_elements(info)        ; number of objects needing reclassification
    nreclass = n_elements(outinfo) ; number of objects successfully reclassified

    splog, 'Successfully re-classified '+string(nreclass,format='(I0)')+'/'+string(nall,format='(I0)')+' objects.'
    struct_print, struct_trimtags(outinfo,except=['WAVE','FLUX'])
    
    if keyword_set(write) then begin
       outfile = 'ages_reclassify_bad_redshifts.fits'
       splog, 'Writing '+analysis_path+outfile+'.gz.'
       mwrfits, outinfo, analysis_path+outfile, /create
       spawn, ['gzip -f '+analysis_path+outfile], /sh
    endif

stop
    
return
end
    
