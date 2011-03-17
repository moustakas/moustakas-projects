;+
; NAME:
;       LDSS3_REDUX
;
; PURPOSE:
;       Wrapper script on reducing all the LDSS3/SG1120 imaging. 
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;       preproc    - make master dark, bias, and flat frames
;       ccdproc    - process the science images (overscan, trim,
;                    bias-subtract, flat-field, and make inverse
;                    variance maps) 
;
; OUTPUTS: 
;       Reduced FITS files written to LDSS3_PATH(/reduced).
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       Before this you should run UNPACK_LDSS3.
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Nov, NYU - written
;       jm07jannyu - developed more
;       jm09mar23nyu - MAKE_BADPIX keyword added 
;
; Copyright (C) 2006-2007, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro ldss3_redux, make_badpix=make_badpix, preproc=preproc, $
  ccdproc=ccdproc

    datapath = ldss3_path(/feb06)+'raw/'
    sg1120_outpath = ldss3_path(/feb06)+'sg1120/'
    standards_outpath = ldss3_path(/feb06)+'standards/'

    overscan = [2032,4110,2159,4185]                            ; for IM_MKBIAS and IM_MKFLAT
    biassec = [overscan[0],overscan[2],overscan[1],overscan[3]] ; for LDSS3_CCDPROC

; choose a generous trim to accommodate both chip 1 and 2    
    
    crop = [502,2031,760,3499]  ; for IM_MKBIAS and IM_MKFLAT
;   crop = [692,2031,799,3428]  ; for IM_MKBIAS and IM_MKFLAT
    scale = [700,1100,1100,1500] ; for IM_MKFLAT (cropped coordinates!!)
    trimsec = crop              ; for LDSS3_CCDPROC
    
    bandpass = ['g','r']
    chip = ['c1','c2']
    gain = 0.7
    rdnoise = 4.5
    satvalue = 50000.0          ; [ADU]

; make the master bad pixel masks
    if keyword_set(make_badpix) then begin
       print & splog, 'Generating the bad pixel masks'
       ldss3_make_badpix, datapath=datapath, trimsec=trimsec, /wfits
    endif
       
; do the pre-processing
    if keyword_set(preproc) then begin
       for i = 0L, n_elements(chip)-1L do begin

; combine the "fast" bias frames; do not trim, since we need the
; overscan region for later; trim on output

          print & splog, 'Generating the master bias frame for for chip '+chip[i]
          biases = file_search(datapath+'a.????_bias_fast_'+chip[i]+'.fits')
          im_mkbias, '', datapath+'bias_'+chip[i]+'.fits', biases, 0, crop=crop, overscan=overscan

; combine the darks; this needs to appear *after* the call to make the
; master bias frame, since it depends on it

          print & splog, 'Generating the master dark frame for chip '+chip[i]
          darks = file_search(datapath+'a.????_dark_fast_'+chip[i]+'.fits')
          im_mkdark, '', datapath+'dark_'+chip[i]+'.fits', darks, 0, $
            datapath+'bias_'+chip[i]+'.fits', crop=crop, overscan=overscan

; test code: generate the flat-fields from the data themselves;
; exclude images with vignetting

          print & splog, 'Generating the test r-band flat-field for chip '+chip[i]+' from the data.'
          flats = file_search(datapath+'a.????_sg1120_?_r_'+chip[i]+'.fits')
          if (i eq 0L) then bad = where(strmatch(flats,'a.2017*') or strmatch(flats,'a.2018*') or $
            strmatch(flats,'a.2020*'),comp=keep) else $
              bad = where(strmatch(flats,'a.2018*') or strmatch(flats,'a.2033*'),comp=keep)
          im_mkflat, '', datapath+'flat_fromdata_r_'+chip[i]+'.fits', flats[keep], 0, crop=crop, $
            overscan=overscan, scale=scale

          print & splog, 'Generating the test g-band flat-field for chip '+chip[i]+' from the data.'
          flats = file_search(datapath+'a.????_sg1120_?_g_'+chip[i]+'.fits')
          if (i eq 0L) then bad = where(strmatch(flats,'a.2016*') or strmatch(flats,'a.2019*') or $
            strmatch(flats,'a.2021*'),comp=keep) else $
              bad = where(strmatch(flats,'a.2019*') or strmatch(flats,'a.2032*'),comp=keep)
          im_mkflat, '', datapath+'flat_fromdata_g_'+chip[i]+'.fits', flats[keep], 0, crop=crop, $
            overscan=overscan, scale=scale

; generate the master g- and r-band flat-fields          
          
          print & splog, 'Generating the r-band flat-field for chip '+chip[i]
          flats = file_search(datapath+'a.????_twilight_r_'+chip[i]+'.fits')
          im_mkflat, '', datapath+'flat_r_'+chip[i]+'.fits', flats, 0, crop=crop, $
            overscan=overscan, scale=scale

          print & splog, 'Generating the g-band flat-field for chip '+chip[i]
          flats = file_search(datapath+'a.????_twilight_g_'+chip[i]+'.fits')
          im_mkflat, '', datapath+'flat_g_'+chip[i]+'.fits', flats, 0, crop=crop, $
            overscan=overscan, scale=scale

       endfor

; more testing of the flat-fields; write out the ratio of the flats
       flat = mrdfits(datapath+'flat_r_c1.fits',0,/silent)
       dataflat = mrdfits(datapath+'flat_fromdata_r_c1.fits',0,h,/silent)
       mwrfits, flat/dataflat, datapath+'flat_dividedby_dataflat_r_c1.fits', h, /create
       
       flat = mrdfits(datapath+'flat_r_c2.fits',0,/silent)
       dataflat = mrdfits(datapath+'flat_fromdata_r_c2.fits',0,h,/silent)
       mwrfits, flat/dataflat, datapath+'flat_dividedby_dataflat_r_c2.fits', h, /create
       
       flat = mrdfits(datapath+'flat_g_c1.fits',0,/silent)
       dataflat = mrdfits(datapath+'flat_fromdata_g_c1.fits',0,h,/silent)
       mwrfits, flat/dataflat, datapath+'flat_dividedby_dataflat_g_c1.fits', h, /create
       
       flat = mrdfits(datapath+'flat_g_c2.fits',0,/silent)
       dataflat = mrdfits(datapath+'flat_fromdata_g_c2.fits',0,h,/silent)
       mwrfits, flat/dataflat, datapath+'flat_dividedby_dataflat_g_c2.fits', h, /create

    endif

; process each bandpass and chip separately
    if keyword_set(ccdproc) then begin
       for ib = 0L, n_elements(bandpass)-1L do begin
          for ic = 0L, n_elements(chip)-1L do begin
             splog, '###########################################################################'
             splog, 'SG1120 '+bandpass[ib]+'-band, '+chip[ic]
             splog, '###########################################################################'
;            imagelist = file_search(datapath+'a.203[2-3]_sg1120_?_'+bandpass[ib]+'_'+chip[ic]+'.fits')
             imagelist = file_search(datapath+'a.????_sg1120_?_'+bandpass[ib]+'_'+chip[ic]+'.fits')
             badpixlist = repstr(imagelist,'.fits','.badpix.fits')
             im_ccdproc, imagelist, outpath=sg1120_outpath, prefix=prefix, gain=gain, $
               rdnoise=rdnoise, biassec=biassec, trimsec=trimsec, satvalue=satvalue, $
               badpixfile=badpixlist, /crrej, /wfits, biasfile=datapath+'bias_'+chip[ic]+'.fits', $
               flatfile=datapath+'flat_'+bandpass[ib]+'_'+chip[ic]+'.fits', /grow_satmask
;              darkfile=datapath+'dark_'+chip[ic]+'.fits', weight_suffix='.weight'
;            ldss3_init_wcs, sg1120_outpath+prefix+file_basename(imagelist), weight_suffix='.weight', /apodize, /wfits
          endfor
; pack everything into a MEF and delete temporary files
          chip1list = file_search(sg1120_outpath+prefix+'a.????_sg1120_?_'+bandpass[ib]+'_c1.fits',count=nimage)
          chip2list = file_search(sg1120_outpath+prefix+'a.????_sg1120_?_'+bandpass[ib]+'_c2.fits')
          for ii = 0L, nimage-1L do begin
             outfile_image = strmid(chip1list[ii],0,strpos(chip1list[ii],'_c1.fits'))+'.fits'
             outfile_weight = repstr(outfile_image,'.fits','.weight.fits')
             outfile_rms = repstr(outfile_image,'.fits','.rms.fits')
             outfile_flag = repstr(outfile_image,'.fits','.flag.fits')

             splog, 'Generating MEF: '+file_basename(outfile_image)
             niceprint, '   '+file_basename([chip1list[ii],chip2list[ii]])
             imagecube = im_fits_cube([chip1list[ii],chip2list[ii]])
             weightcube = im_fits_cube(repstr([chip1list[ii],chip2list[ii]],'.fits','.weight.fits'))
             rmscube = im_fits_cube(repstr([chip1list[ii],chip2list[ii]],'.fits','.rms.fits'))
             flagcube = im_fits_cube(repstr([chip1list[ii],chip2list[ii]],'.fits','.flag.fits'))

             mwrfits, 0, outfile_image, /create
             mwrfits, 0, outfile_weight, /create
             mwrfits, 0, outfile_rms, /create
             mwrfits, 0, outfile_flag, /create
             for icc = 0L, 1L do begin
                case icc of
                   0L: begin ; reverse chip 1
                      mwrfits, float(reverse(imagecube[icc].image)), outfile_image, $
                        imagecube[icc].header[where(strcompress(imagecube[icc].header,/remove) ne '')]
                      mwrfits, float(reverse(weightcube[icc].image)), outfile_weight, $
                        weightcube[icc].header[where(strcompress(weightcube[icc].header,/remove) ne '')]
                      mwrfits, float(reverse(rmscube[icc].image)), outfile_rms, $
                        rmscube[icc].header[where(strcompress(rmscube[icc].header,/remove) ne '')]
                      mwrfits, fix(reverse(flagcube[icc].image)), outfile_flag, $
                        flagcube[icc].header[where(strcompress(flagcube[icc].header,/remove) ne '')]
                   end
                   1L: begin
                      mwrfits, float(imagecube[icc].image), outfile_image, $
                        imagecube[icc].header[where(strcompress(imagecube[icc].header,/remove) ne '')]
                      mwrfits, float(weightcube[icc].image), outfile_weight, $
                        weightcube[icc].header[where(strcompress(weightcube[icc].header,/remove) ne '')]
                      mwrfits, float(rmscube[icc].image), outfile_rms, $
                        rmscube[icc].header[where(strcompress(rmscube[icc].header,/remove) ne '')]
                      mwrfits, fix(flagcube[icc].image), outfile_flag, $
                        flagcube[icc].header[where(strcompress(flagcube[icc].header,/remove) ne '')]
                   end
                endcase
                rmfile, imagecube[icc].fname
                rmfile, weightcube[icc].fname
                rmfile, rmscube[icc].fname
                rmfile, flagcube[icc].fname
             endfor 
          endfor 
       endfor 
; RELEGATED!
;      for ib = 0L, n_elements(bandpass)-1L do begin
;         for ic = 0L, n_elements(chip)-1L do begin
;            splog, '###########################################################################'
;            splog, 'Standards '+bandpass[ib]+'-band, '+chip[ic]
;            splog, '###########################################################################'
;            imagelist = file_search(datapath+'a.????_sdss_1_'+bandpass[ib]+'_'+chip[ic]+'.fits')
;            im_ccdproc, imagelist, outpath=standards_outpath, prefix=prefix, gain=gain, rdnoise=rdnoise, $
;              biassec=biassec, trimsec=trimsec, biasfile=datapath+'bias_'+chip[ic]+'.fits', $
;              flatfile=datapath+'flat_'+bandpass[ib]+'_'+chip[ic]+'.fits', $ ; darkfile=datapath+'dark_'+chip[ic]+'.fits', 
;              /wfits, weight_suffix='.weight'
;            ldss3_init_wcs, standards_outpath+prefix+file_basename(imagelist), weight_suffix='.weight', /apodize, /wfits
;         endfor
;      endfor
    endif 

return
end
