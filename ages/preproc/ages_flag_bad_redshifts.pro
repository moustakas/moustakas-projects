pro ages_flag_bad_redshifts, write=write
; jm06jan19uofa - flag objects for which the SPZBEST redshift differs
;                 from the "official" AGES redshift by more than
;                 DZTHRESH percent; these objects are important
;                 because they may be misclassified ("QSO" into
;                 "GALAXY") by the Schlegel classification scheme; see
;                 AGES_RERUN_SPZBEST 

    dzthresh = 0.02 ; photometric redshifts having dz[=(z1-z2)/(1+z)]<0.02 are considered "accurate"
    
    rootpath = ages_path(/rawdata)
    analysis_path = ages_path(/analysis)

    vredux_11 = '0301'
    vredux_20 = '0051'
    vredux = [vredux_11,vredux_20]

    datapath_11 = rootpath+vredux_11+'/'
    datapath_20 = rootpath+vredux_20+'/'
    datapath = [datapath_11,datapath_20]

    badplates = [$
      '106',$ ; not fluxed
      '110',$ ; not fluxed
      '209',$ ; not fluxed
      '310',$ ; not fluxed
      '311' $ ; not fluxed
      ]
    badplates = '' ; note!
    
    info = mrdfits(analysis_path+'catalog.zmerge.fits.gz',1,/silent)

    outinfo_template = {$
      platefile:  '', $
      pass:       0L, $
      aper:       0L, $
      class:      '', $
      z_spzbest: 0.0, $
      z:         0.0, $
      newclass:   '', $
      wave:      fltarr(4608), $
      flux:      fltarr(4608)}

    t0 = systime(1)
    for ipath = 0L, n_elements(datapath)-1L do begin
    
       sphectfits = file_basename(file_search(datapath[ipath]+'spHect-???-'+vredux[ipath]+'.fits.gz',count=npassfield))
       spzbestfits = file_basename(file_search(datapath[ipath]+'spZbest-???-'+vredux[ipath]+'.fits.gz'))
       passfield = strmid(sphectfits,7,3) ; NOT GENERAL

;      for jpass = 10L, npassfield-1L do begin
       for jpass = 0L, npassfield-1L do begin

;         splog, 'Unpacking PASSFIELD '+passfield[jpass]+'.'

; skip unfluxed plates          
          
          badmatch = where(strmatch(badplates,passfield[jpass]) eq 1B,nbadmatch)
          if (nbadmatch eq 1L) then splog, 'Skipping bad plate '+badplates[badmatch]+'.' else begin

; read the data          
             
             bigwave   = mrdfits(datapath[ipath]+sphectfits[jpass],0,/silent)
             bigflux   = mrdfits(datapath[ipath]+sphectfits[jpass],1,/silent)
             biginvvar = mrdfits(datapath[ipath]+sphectfits[jpass],2,/silent)
             plugmap   = mrdfits(datapath[ipath]+sphectfits[jpass],5,/silent)
             spzbest   = mrdfits(datapath[ipath]+spzbestfits[jpass],1,/silent)

; identify every fiber with a well-measured ZMERGE redshift             
             
             fiberindx = where(info.zmerge_pass eq passfield[jpass],nfiberindx)
             if (nfiberindx eq 0L) then message, 'This should not happen!'
             thesefibers = info[fiberindx].zmerge_aper-1L ; zeroindex

;            srt = sort(info[fiberindx].zmerge_aper)
;            struct_print, info[fiberindx[srt]]
;            struct_print, info[fiberindx[srt[229:233]]]

; this happens sometimes because there are duplicate entries in the
; zmerge catalog; but we don't care             
;            nunique = n_elements(uniq(info[fiberindx].zmerge_aper,$
;              sort(info[fiberindx].zmerge_aper)))
;            if (nfiberindx ne nunique) then message, 'This should not happen!'
             
;            zbadflag = where((abs(spzbest[thesefibers].z-info[fiberindx].zmerge_z) ge 0.01),nzbadflag)
;            zbadflag = where((strtrim(spzbest[thesefibers].class,2) eq 'GALAXY') and $
;              (abs(spzbest[thesefibers].z-info[fiberindx].zmerge_z) ge 0.01),nzbadflag)
;            zbadflag = where(((strtrim(spzbest[thesefibers].class,2) eq 'GALAXY') or $
;              (strtrim(spzbest[thesefibers].class,2) eq 'QSO')) and $
;              (info[fiberindx].zmerge_z gt 0.0) and $
;              (abs(spzbest[thesefibers].z-info[fiberindx].zmerge_z) ge 0.01),nzbadflag)

             deltaz = abs(spzbest[thesefibers].z-info[fiberindx].zmerge_z)
             dz = deltaz/(1.0+mean([info[fiberindx].zmerge_z,spzbest[thesefibers].z]))
;            plot, deltaz, dz, ps=4, /xlog, /ylog

             zbadflag = where((info[fiberindx].zmerge_z gt 0.0) and (dz ge dzthresh),nzbadflag)

             if (nzbadflag ne 0L) then begin

;               splog, 'The following '+string(nzbadflag,format='(I0)')+$
;                 ' fiber(s) may be mis-classified:'
;               niceprint, info[fiberindx[zbadflag]].zmerge_pass, info[fiberindx[zbadflag]].zmerge_aper, $
;                 spzbest[thesefibers[zbadflag]].class, spzbest[thesefibers[zbadflag]].z, $
;                 info[fiberindx[zbadflag]].zmerge_z

; store the results                
                
                outinfo1 = replicate(outinfo_template,nzbadflag)
                outinfo1.platefile = sphectfits[jpass]
                
                for iflag = 0L, nzbadflag-1L do begin

                   outinfo1[iflag].pass = info[fiberindx[zbadflag[iflag]]].zmerge_pass
                   outinfo1[iflag].aper = info[fiberindx[zbadflag[iflag]]].zmerge_aper
                   outinfo1[iflag].class = spzbest[thesefibers[zbadflag[iflag]]].class
                   outinfo1[iflag].newclass = outinfo1[iflag].class ; for classification, later
                   outinfo1[iflag].z_spzbest = spzbest[thesefibers[zbadflag[iflag]]].z
                   outinfo1[iflag].z = info[fiberindx[zbadflag[iflag]]].zmerge_z
                   outinfo1[iflag].wave = bigwave[*,thesefibers[zbadflag[iflag]]]
                   outinfo1[iflag].flux = bigflux[*,thesefibers[zbadflag[iflag]]]

                   struct_print, struct_trimtags(outinfo1[iflag],except=['WAVE','FLUX']), /no_head
;                  plot, outinfo1[iflag].wave, outinfo1[iflag].flux, ps=10, xsty=3, ysty=3
;                  cc = get_kbrd(1)
                   
                endfor

                if (n_elements(outinfo) eq 0L) then outinfo = outinfo1 else outinfo = struct_append(outinfo,outinfo1)

;               for irej = 0L, nzbadflag-1L do printf, lun1, info[fiberindx[zbadflag[irej]]].galaxy, $
;                 passfield[jpass], thesefibers[zbadflag[irej]], spzbest[thesefibers[zbadflag[irej]]].z, $
;                 info[fiberindx[zbadflag[irej]]].zmerge_z, format='(A14,2x,I3.3,2x,I3.3,2x,F8.5,2x,F8.5)'
             endif

          endelse

       endfor 

    endfor

    if keyword_set(write) then begin

; NOTE!  I am only interested in galaxies that may have been
; classified as QSO's, so only write out the objects at z<1 since
; everything else will be a QSO regardless of whether the redshift was
; right or not       

       flag = lindgen(n_elements(outinfo))
;      flag = where(outinfo.z gt 1.2,nflag)
;      flag = where(outinfo.z le 1.2,nflag)

       outfile = 'ages_flag_bad_redshifts.fits'
       splog, 'Writing '+analysis_path+outfile+'.gz.'
       mwrfits, outinfo[flag], analysis_path+outfile, /create
       spawn, ['gzip -f '+analysis_path+outfile], /sh

       openw, lun, analysis_path+repstr(outfile,'.fits','.txt'), /get_lun
       struct_print, struct_trimtags(outinfo[flag],except=['WAVE','FLUX']), lun=lun
       free_lun, lun

    endif

;   z = outinfo.z
;   deltaz = abs(outinfo.z_spzbest-outinfo.z)
;   dz = deltaz/(1.0+z)
;   plot, deltaz, dz, ps=4, /xlog, /ylog
;   plot, z, dz, ps=4, /xlog, /ylog

stop
    
return
end
    
