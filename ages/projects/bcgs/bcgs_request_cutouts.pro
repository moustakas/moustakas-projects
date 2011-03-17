pro bcgs_request_cutouts
; jm10jul23ucsd
    
    bcgspath = ages_path(/projects)+'bcgs/'
    outpath = bcgspath+'thumbs/'

    bcgs = mrdfits(bcgspath+'bcgs_photometry_v3.fits.gz',1)
    ngal = n_elements(bcgs)
    bcgsid = string(lindgen(ngal),format='(I2.2)')
    
    ra = strtrim(im_dec2hms(bcgs.i_alpha_j2000/15D,/col),2)
    dec = repstr(strtrim(im_dec2hms(bcgs.i_delta_j2000,/col),2),'+','')
    rawidth = '0.75' ; =45 arcsec

    band = 'I'
;   band = ['Bw','R','I','K']
    nband = n_elements(band)
    for iband = 0, nband-1 do begin
       fitsname = outpath+'bcg'+bcgsid+'_thumb_'+band[iband]+'.fits'

       scriptname = outpath+'cut_'+band[iband]+'.wget'
       openw, lun, scriptname, /get_lun
       for ii = 0, ngal-1 do begin
          printf, lun, 'wget "http://archive.noao.edu/ndwfs/cutout.php?ra='+$
            ra[ii]+'&dec='+dec[ii]+'&rawidth='+rawidth+'&decwidth=INDEF&fcsystem=J2000&filters='+band[iband]+'" '+$
            '-O "'+fitsname[ii]+'"'
       endfor
       free_lun, lun
       spawn, 'chmod +x '+scriptname, /sh
       spawn, scriptname, /sh
    endfor
    
return
end
    
