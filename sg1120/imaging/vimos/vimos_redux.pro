;+
; NAME:
;       VIMOS_REDUX
;
; PURPOSE:
;       Wrapper script on reducing all the VIMOS/SG1120 imaging.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;       make_badpix - generate the quadrant-specific bad pixel mask
;       preproc    - make master dark, bias, and flat frames
;       ccdproc    - process the science images (overscan, trim,
;                    bias-subtract, flat-field, and make inverse
;                    variance maps) 
;
; OUTPUTS: 
;       Reduced FITS files written to VIMOS_PATH(/reduced).
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       Before this you should run UNPACK_VIMOS.
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Nov, NYU - written
;       jm07jannyu - developed more
;       jm08julnyu - added MAKE_BADPIX keyword; additional tweaks 
;
; Copyright (C) 2006-2008, John Moustakas
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

pro vimos_redux, bandpass=bandpass, dec03=dec03, feb06=feb06, $
  make_badpix=make_badpix, preproc=preproc, ccdproc=ccdproc

    if ((not keyword_set(dec03)) and (not keyword_set(feb06))) or $
      (keyword_set(dec03) and keyword_set(feb06)) then begin
       splog, 'Must specify either DEC03 or FEB06 keyword!'
       return
    endif
    
    datapath = vimos_path(dec03=dec03,feb06=feb06)+'raw/'
    sg1120_outpath = vimos_path(dec03=dec03,feb06=feb06)+'sg1120/'
    standards_outpath = vimos_path(dec03=dec03,feb06=feb06)+'standards/'

    if keyword_set(dec03) then suffix = '2003' else suffix = '2006'
    
;   info = vimos_forage(file_search(datapath+'????_*.fits'))
;   struct_print, info

    overscan = [0,200,49,2300]
    biassec = [overscan[0],overscan[2],overscan[1],overscan[3]]
    crop = [0,2147,0,2439]
    scale = [900,1300,1000,1400]
    trimsec = crop

; http://www.eso.org/observing/dfo/quality/VIMOS/qc/scrflat_qc1.html]    
; http://www.eso.org/observing/etc/doc/ut3/vimos/vimos-etc-um.html#heading52

    if (not keyword_set(bandpass)) then bandpass = ['B','V','R']
    quadrant = ['Q1','Q2','Q3','Q4']
;   allgain = [0.58,0.54,0.51,0.56] ; [ADU/e]
    allgain = [1.70,1.86,1.95,1.80] ; [e/ADU]
    allrdnoise = [5.15,4.0,3.9,4.4]
    satvalue = 50000.0          ; [ADU]
;   satvalue = 220000.0         ; [e]

; make the master bad pixel masks
    if keyword_set(make_badpix) then begin

       print & splog, 'Generating the generic bad pixel masks for each quadrant'
       vimos_make_badpix, nx=crop[1]+1, ny=crop[3]+1, suffix=suffix, $
         trimsec=trimsec, outpath=datapath, /wfits

       for iq = 0L, n_elements(quadrant)-1L do begin
          print & splog, 'Generating customized bad pixel masks for quandrant '+quadrant[iq]
          vimos_make_custom_badpix, datapath=datapath, trimsec=trimsec, quadrant=quadrant[iq], /wfits
       endfor

    endif
    
; do the pre-processing
    
    if keyword_set(preproc) then begin

       for iq = 0L, n_elements(quadrant)-1L do begin
          print & splog, 'Generating the master bias frame for quadrant '+quadrant[iq]+'.'
          biases = file_search(datapath+'calib/a.bias*_'+quadrant[iq]+'*.fits')
          im_mkbias, '', datapath+'bias_'+quadrant[iq]+'.fits', biases, 0, crop=crop, overscan=overscan
 
          for ib = 0L, n_elements(bandpass)-1L do begin
             print & splog, 'Generating the '+bandpass[ib]+'-band flat-field for quadrant '+quadrant[iq]+'.'
             flats = file_search(datapath+'calib/a.twilight*_'+bandpass[ib]+'_'+quadrant[iq]+'.fits')
             im_mkflat, '', datapath+'flat_'+bandpass[ib]+'_'+quadrant[iq]+'.fits', flats, 0, crop=crop, $
               overscan=overscan, scale=scale
          endfor
       endfor

    endif

; process each bandpass and quadrant separately

    if keyword_set(ccdproc) then begin
       for ib = 0L, n_elements(bandpass)-1L do begin
          for iq = 0L, n_elements(quadrant)-1L do begin
             splog, '###########################################################################'
             splog, 'SG1120 '+bandpass[ib]+'-band, '+quadrant[iq]
             splog, '###########################################################################'
             imagelist = file_search(datapath+'a.sg1120*_'+bandpass[ib]+'_'+quadrant[iq]+'.fits')
             badpixlist = repstr(imagelist,'.fits','.badpix.fits')
             im_ccdproc, imagelist, outpath=sg1120_outpath, prefix=prefix, gain=allgain[iq], $
               rdnoise=allrdnoise[iq], biassec=biassec, trimsec=trimsec, satvalue=satvalue, $
               badpixfile=badpixlist, /crrej, /wfits, biasfile=datapath+'bias_'+quadrant[iq]+'.fits', $
               flatfile=datapath+'flat_'+bandpass[ib]+'_'+quadrant[iq]+'.fits', /grow_satmask
;              badpixfile=datapath+'badpix_'+quadrant[iq]+'.fits', weight_suffix='.weight'
;            vimos_apodize, sg1120_outpath+prefix+file_basename(imagelist), weight_suffix='.weight', /wfits
          endfor
; pack everything into a MEF and delete temporary files
          allimagelist = file_search(sg1120_outpath+prefix+'a.sg1120*_'+bandpass[ib]+'_Q?.fits')
          imagelist = file_search(sg1120_outpath+prefix+'a.sg1120*_'+bandpass[ib]+'_Q1.fits',count=nimage)
          for ii = 0L, nimage-1L do begin

             quadrantlist = file_search(sg1120_outpath+strmid(file_basename(imagelist[ii]),0,$
               strpos(file_basename(imagelist[ii]),':',/reverse_search))+':??.???_'+$
               bandpass[ib]+'_Q?.fits')
             if (n_elements(quadrantlist) ne 4L) then message, 'Problem here!'
             quadrantlist = quadrantlist[sort(strmid(quadrantlist,6,2,/reverse))]
             outfile_image = strmid(quadrantlist[0],0,strpos(quadrantlist[0],'_Q1.fits'))+'.fits'
             outfile_weight = repstr(outfile_image,'.fits','.weight.fits')
             outfile_rms = repstr(outfile_image,'.fits','.rms.fits')
             outfile_flag = repstr(outfile_image,'.fits','.flag.fits')

             print & splog, 'Generating MEF: '+file_basename(outfile_image)
             niceprint, '   '+file_basename(quadrantlist)
             imagecube = im_fits_cube(quadrantlist)
             weightcube = im_fits_cube(repstr(quadrantlist,'.fits','.weight.fits'))
             rmscube = im_fits_cube(repstr(quadrantlist,'.fits','.rms.fits'))
             flagcube = im_fits_cube(repstr(quadrantlist,'.fits','.flag.fits'))

             mwrfits, 0, outfile_image, /create
             mwrfits, 0, outfile_weight, /create
             mwrfits, 0, outfile_rms, /create
             mwrfits, 0, outfile_flag, /create
             for iquad = 0L, 3L do begin
                mwrfits, float(imagecube[iquad].image), outfile_image, imagecube[iquad].header
                mwrfits, float(weightcube[iquad].image), outfile_weight, weightcube[iquad].header
                mwrfits, float(rmscube[iquad].image), outfile_rms, rmscube[iquad].header
                mwrfits, long(flagcube[iquad].image), outfile_flag, flagcube[iquad].header

                rmfile, imagecube[iquad].fname
                rmfile, weightcube[iquad].fname
                rmfile, rmscube[iquad].fname
                rmfile, flagcube[iquad].fname
             endfor 
          endfor   
       endfor  
; now do the standards; don't bother packing them into MEF
; because we need to solve for the zeropoint on a quadrant-by-quadrant
; basis anyway
       for ib = 0L, n_elements(bandpass)-1L do begin
          for iq = 0L, n_elements(quadrant)-1L do begin
             splog, '###########################################################################'
             splog, 'Standards '+bandpass[ib]+'-band, '+quadrant[iq]
             splog, '###########################################################################'
             imagelist = file_search(datapath+'a.[PG,SA]*_'+bandpass[ib]+'_'+quadrant[iq]+'.fits')
             im_ccdproc, imagelist, outpath=standards_outpath, prefix=prefix, gain=allgain[iq], $
               rdnoise=allrdnoise[iq], biassec=biassec, trimsec=trimsec, satvalue=satvalue, $
               /crrej, /wfits, biasfile=datapath+'bias_'+quadrant[iq]+'.fits', $
               flatfile=datapath+'flat_'+bandpass[ib]+'_'+quadrant[iq]+'.fits', $
               badpixfile=datapath+'badpix_'+quadrant[iq]+'.fits', /grow_satmask
          endfor 
       endfor  
    endif 

return
end
