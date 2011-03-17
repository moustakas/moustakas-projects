pro unpack_vimos_single, outpath=outpath, datapath=datapath, calibpath=calibpath, $
  recentdata=recentdata, unpack=unpack
; jm07jan15nyu - written

    imageinfo = vimos_forage(file_search(datapath+'VIMOS*.fits',count=nimage))
    if keyword_set(recentdata) then begin
       caliblist = [file_search(calibpath+'DET/VIMOS*.fits'),$
         file_search(calibpath+'IMG_FLAT/VIMOS*.fits'),$
         file_search(calibpath+'IMG_STD/VIMOS*.fits')]
       str = strsplit(datapath,'/',/extract)
       datapath_reduced = '/'+strjoin(str[0L:n_elements(str)-2L],'/')+'/sci_proc/'
       allfiles = file_basename(file_search(datapath_reduced+'/*.fits'))
       for ii = 0L, nimage-1L do imageinfo[ii].file_reduced = $
         allfiles[where(strmatch(allfiles,'*'+repstr(repstr(imageinfo[ii].file,'.fits',''),'VIMOS.','')+'*'))]
    endif else begin
       caliblist = file_search(calibpath+'VIMOS*.fits')
       str = strsplit(datapath,'/',/extract)
       datapath_reduced = '/'+strjoin(str[0L:n_elements(str)-2L],'/')+'/reduced/'
       allfiles = file_basename(file_search(datapath_reduced+'/*.fits'))
       for ii = 0L, nimage-1L do imageinfo[ii].file_reduced = $
         allfiles[where(strmatch(allfiles,'*'+repstr(repstr(imageinfo[ii].file,'.fits',''),'VIMOS.','')+'*'))]
    endelse
    calibinfo = vimos_forage(caliblist)
    ncalib = n_elements(calibinfo)
    if (nimage eq 0L) or (ncalib eq 0L) then message, 'Problem here!'
;   struct_print, calibinfo

; object data; ; remove R-band frames from 2003-12-24 with crap seeing
; (see 072.A-0367A_01/GEN_INFO/qc0_report.txt)
    
;   ignore = ['junk'] 

    ignore = [$
      'VIMOS.2003-12-24T07:50:43.548.fits',$ ; the first four exposures
      'VIMOS.2003-12-24T07:50:43.549.fits',$ ; *could* be used (seeing 
      'VIMOS.2003-12-24T07:50:43.550.fits',$ ; was 0.65")
      'VIMOS.2003-12-24T07:50:43.551.fits',$ ; 

      'VIMOS.2003-12-24T08:11:31.819.fits',$
      'VIMOS.2003-12-24T08:11:31.820.fits',$
      'VIMOS.2003-12-24T08:06:16.616.fits',$
      'VIMOS.2003-12-24T08:06:16.617.fits',$
      'VIMOS.2003-12-24T08:06:16.698.fits',$
      'VIMOS.2003-12-24T08:06:16.699.fits',$
      'VIMOS.2003-12-24T08:01:05.708.fits',$
      'VIMOS.2003-12-24T08:01:05.709.fits',$
      'VIMOS.2003-12-24T08:01:05.710.fits',$
      'VIMOS.2003-12-24T08:01:05.711.fits',$
      'VIMOS.2003-12-24T07:55:55.923.fits',$
      'VIMOS.2003-12-24T08:11:31.901.fits',$
      'VIMOS.2003-12-24T07:55:55.924.fits',$
      'VIMOS.2003-12-24T08:11:31.902.fits',$
      'VIMOS.2003-12-24T07:55:55.927.fits',$
      'VIMOS.2003-12-24T07:55:55.928.fits']

    for i = 0L, nimage-1L do begin
       file = strcompress(imageinfo[i].file)
       skip = 0
       for j = 0L, n_elements(ignore)-1L do if strmatch(file,'*'+ignore[j]+'*') then skip = 1
       if (skip eq 0L) then begin
          root = 'sg1120'
          outfile = 'a.'+root+'_'+strtrim(imageinfo[i].date,2)+'_'+$
            strtrim(imageinfo[i].filter,2)+'_Q'+strtrim(imageinfo[i].quadrant,2)+'.fits'
          if keyword_set(unpack) then begin
             spawn, '/bin/cp -vf '+datapath+file+' '+outpath+outfile
             im = mrdfits(outpath+outfile,0,hdr,/silent)
;            hdr = headfits(outpath+outfile)
;            hdr = hdr[where(strcompress(hdr,/remove) ne '')]
             hdr = hdr[where((strcompress(hdr,/remove) ne '') and (strmatch(hdr,'*HIERARCH*',/fold) eq 0B) and $
               (strmatch(hdr,'*COMMENT*',/fold) eq 0B) and (strmatch(hdr,'*HISTORY*',/fold) eq 0B))]
;            sxdelpar, hdr, 'BSCALE' ; DO NOT DELETE THESE KEYWORDS
;            sxdelpar, hdr, 'BZERO'
             sxdelpar, hdr, 'FILTER1' & sxdelpar, hdr, 'FILTER2' & sxdelpar, hdr, 'FILTER3' & sxdelpar, hdr, 'FILTER4'
             sxaddpar, hdr, 'FILTER', strtrim(imageinfo[i].filter,2), ' Response Function', after='AIRMASS'
             sxaddpar, hdr, 'QUADRANT', 'Q'+strtrim(imageinfo[i].quadrant,2), ' Quadrant Number', after='FILTER'
             sxaddpar, hdr, 'PIXSCALE', float(imageinfo[i].pixscale), ' Pixel scale (arcsec/pixel)', after='QUADRANT'
             sxaddpar, hdr, 'SEEING', float(imageinfo[i].seeing), ' FWHM seeing (arcsec)', after='PIXSCALE'
             sxaddpar, hdr, 'GAIN', float(imageinfo[i].gain), ' Gain (e-/ADU)', after='SEEING'
             sxaddpar, hdr, 'RDNOISE', float(imageinfo[i].rdnoise), ' Read-noise (e-)', after='GAIN'
; grab the magnitude zero-point
             moreinfo = vimos_forage(datapath_reduced+imageinfo[i].file_reduced)
;            help, moreinfo, /str
             sxaddpar, hdr, 'MAG0', float(moreinfo.mag_zero), ' VIMOS pipeline magnitude zero-point', after='RDNOISE'
             sxaddpar, hdr, 'EMAG0', float(moreinfo.magerr_zero), ' error in MAG0', after='MAG0'
; write out
             mwrfits, im, outpath+outfile, hdr, /create
;            modfits, outpath+outfile, 0, hdr ; jm08jul18nyu - modfits is broken
          endif else splog, file+' --> '+outpath+outfile
       endif else splog, 'Skipping '+file
    endfor

; calibration data    

    ignore = ['junk'] 

    for i = 0L, ncalib-1L do begin
       file = strcompress(calibinfo[i].file)
       skip = 0
       for j = 0L, n_elements(ignore)-1L do if strmatch(file,'*'+ignore[j]+'*') then skip = 1
       if (skip eq 0L) then begin
          if strmatch(calibinfo[i].object,'*flat*sky*',/fold) then root = 'twilight'
          if strmatch(calibinfo[i].object,'*bias*',/fold) then root = 'bias'
          if strmatch(calibinfo[i].object,'*std*',/fold) then root = repstr(calibinfo[i].target,'-Stetson','')
          if strmatch(calibinfo[i].object,'*bias*',/fold) then $
            outfile = 'a.'+root+'_'+strtrim(calibinfo[i].date,2)+'_Q'+strtrim(calibinfo[i].quadrant,2)+'.fits' else $
              outfile = 'a.'+root+'_'+strtrim(calibinfo[i].date,2)+'_'+$
            strtrim(calibinfo[i].filter,2)+'_Q'+strtrim(calibinfo[i].quadrant,2)+'.fits'
          if strmatch(calibinfo[i].object,'*std*',/fold) then $
            caliboutpath = outpath else caliboutpath = outpath+'calib/'
          if keyword_set(unpack) then begin
             spawn, '/bin/cp -vf '+caliblist[i]+' '+caliboutpath+outfile ; NOTE difference relative to data
             im = mrdfits(caliboutpath+outfile,0,hdr,/silent)
;            hdr = headfits(caliboutpath+outfile)
             hdr = hdr[where(strcompress(hdr,/remove) ne '')]
             hdr = hdr[where((strcompress(hdr,/remove) ne '') and (strmatch(hdr,'*HIERARCH*',/fold) eq 0B) and $
               (strmatch(hdr,'*COMMENT*',/fold) eq 0B) and (strmatch(hdr,'*HISTORY*',/fold) eq 0B))]
;            sxdelpar, hdr, 'BSCALE' ; DO NOT DELETE THESE KEYWORDS
;            sxdelpar, hdr, 'BZERO'
             sxdelpar, hdr, 'FILTER1' & sxdelpar, hdr, 'FILTER2' & sxdelpar, hdr, 'FILTER3' & sxdelpar, hdr, 'FILTER4'
             sxaddpar, hdr, 'FILTER', strtrim(calibinfo[i].filter,2), ' Response Function', after='AIRMASS'
             sxaddpar, hdr, 'QUADRANT', 'Q'+strtrim(calibinfo[i].quadrant,2), ' Quadrant Number', after='FILTER'
             sxaddpar, hdr, 'PIXSCALE', float(calibinfo[i].pixscale), ' Pixel scale (arcsec/pixel)', after='QUADRANT'
             sxaddpar, hdr, 'SEEING', float(calibinfo[i].seeing), ' FWHM seeing (arcsec)', after='PIXSCALE'
             sxaddpar, hdr, 'GAIN', float(calibinfo[i].gain), ' Gain (e-/ADU)', after='SEEING'
             sxaddpar, hdr, 'RDNOISE', float(calibinfo[i].rdnoise), ' Read-noise (e-)', after='GAIN'
             mwrfits, im, caliboutpath+outfile, hdr, /create
;            modfits, caliboutpath+outfile, 0, hdr ; jm08jul18nyu - modfits is broken
          endif else splog, file+' --> '+caliboutpath+outfile
       endif
    endfor

return
end

pro unpack_vimos, unpack=unpack
; jm07jan15nyu - written

    rootpath = vimos_path(/original)
    outpath1 = vimos_path(/dec03)+'raw/'
    outpath2 = vimos_path(/feb06)+'raw/'

    splog, '###########################################################################'
    unpack_vimos_single, outpath=outpath1, datapath=rootpath+'072.A-0367A_01/2003-12-23/raw/', $
      calibpath=rootpath+'072.A-0367A_01/2003-12-23/calib/', unpack=unpack
    splog, '###########################################################################'
    unpack_vimos_single, outpath=outpath1, datapath=rootpath+'072.A-0367A_01/2003-12-26/raw/', $
      calibpath=rootpath+'072.A-0367A_01/2003-12-26/calib/', unpack=unpack
    splog, '###########################################################################'
    unpack_vimos_single, outpath=outpath1, datapath=rootpath+'072.A-0367B_01/2003-12-26/raw/', $
      calibpath=rootpath+'072.A-0367B_01/2003-12-26/calib/', unpack=unpack

    splog, '###########################################################################'
    unpack_vimos_single, outpath=outpath2, datapath=rootpath+'076.B-0362B/202848/sci_raw/', $
      calibpath=rootpath+'076.B-0362B/GEN_CALIB/raw/', /recentdata, unpack=unpack
    splog, '###########################################################################'
    unpack_vimos_single, outpath=outpath2, datapath=rootpath+'076.B-0362A/202845/sci_raw/', $
      calibpath=rootpath+'076.B-0362A/GEN_CALIB/raw/', /recentdata, unpack=unpack

return
end
    
