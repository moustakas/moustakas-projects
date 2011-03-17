pro vimos_fit_scattered_light, dec03=dec03, feb06=feb06, ccdproc=ccdproc, skyobjsub=skyobjsub
; jm08aug25nyu - model the scattered light in Q3

    if ((not keyword_set(dec03)) and (not keyword_set(feb06))) or $
      (keyword_set(dec03) and keyword_set(feb06)) then begin
       splog, 'Must specify either DEC03 or FEB06 keyword!'
       return
    endif
    
    datapath = vimos_path(dec03=dec03,feb06=feb06)+'raw/'
    outpath = datapath+'scattered_light/'
    sexpath = sg1120_path(/sex)
    sexconfig = sexpath+'default.sex'
    sexparam = sexpath+'sg1120.sex.param'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

    bandpass = 'B'
    
    if keyword_set(ccdproc) then begin

       overscan = [0,200,49,2300]
       biassec = [overscan[0],overscan[2],overscan[1],overscan[3]]
       trimsec = [0,2147,0,2439]

       imagelist = file_search(datapath+'a.sg1120*_'+bandpass+'_Q3.fits')
       im_ccdproc, imagelist, outpath=outpath, prefix=prefix, gain=0.51, rdnoise=3.9, $
         biassec=biassec, trimsec=trimsec, biasfile=datapath+'bias_Q3.fits', $
         flatfile=datapath+'flat_'+bandpass+'_Q3.fits', /crrej, /wfits, $
         badpixfile=datapath+'badpix_Q3.fits', weight_suffix='.weight'
       
    endif

; subtract the sky and the objects    
       
    if keyword_set(skyobjsub) then begin

       imagelist = file_search(outpath+'ra.sg1120*_'+bandpass+'_Q3.fits')
       weightlist = repstr(imagelist,'.fits','.weight.fits')
       noskyobjlist = outpath+file_basename(repstr(imagelist,'.fits','.noskyobj.fits'))

;-CATALOG_NAME '+catlist[ii]+' -CATALOG_TYPE FITS_LDAC'+$            

       t0 = systime(1)
       for ii = 0L, n_elements(imagelist)-1L do begin
          print, 'Sky-subtracting '+imagelist[ii]
          seeing = string(sxpar(headfits(imagelist[ii]),'SEEING'),format='(G0.0)')
          spawn, 'sex '+imagelist[ii]+' -c '+sexconfig+' -CATALOG_TYPE NONE'+$
            ' -PARAMETERS_NAME '+sexparam+' -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+$
            ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+weightlist[ii]+' -WEIGHT_THRESH 0 -SEEING_FWHM '+seeing+$
            ' -VERBOSE_TYPE NORMAL -NTHREADS 4 -CHECKIMAGE_TYPE FILTERED -CHECKIMAGE_NAME '+noskyobjlist[ii], /sh
       endfor
       splog, 'Total time to run SE = ', (systime(1)-t0)/60.0, ' minutes.'
       
    endif

    flist = file_search(outpath+'ra.sg1120*_'+bandpass+'_Q3*.fits')
    cube = im_fits_cube(flist)
    im_avgclip, cube.image, im, /norm, thresh=2.5
;   im = total(cube.image,3)
    writefits, outpath+'test.fits', im

stop
    
;   im_mkflat, '', datapath+'test.fits', flist[0:2], 0
    
    
stop    

return
end
    
