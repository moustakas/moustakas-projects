pro wise_sdss_readspec, getinfo=getinfo, getspec=getspec
; jm11apr13ucsd - get the SDSS spectra of WISE galaxies

    path = wise_path()
    if keyword_set(getinfo) then begin
       cat = mrdfits(path+'atlas_wise_nearby.fits.gz',1)
       these = where((cat.isdss ne -1) and (cat.inside ne -1),ngal)
       info = mrdfits(path+'sdss_atlas.fits.gz',1,row=these)

       glactc, info.ra, info.dec, 2000.0, gl, gb, 1, /degree
       ebv_mw = dust_getval(gl,gb,/interp,/noloop)
       
       out = struct_trimtags(info,select=['plate','fiberid',$
         'mjd','z','z_err','vdisp','vdisp_err','zwarning','sn_median'])
       out = struct_addtags(cat[these],out)
       out = struct_addtags(out,replicate({ebv_mw: 0.0},ngal))
       out.ebv_mw = ebv_mw
       im_mwrfits, out, 'wise_readspec_info.fits', /clobber
    endif
       
; now pass the FITS table to NYU and get the spectra using readspec;
; split into combinations of PLATE+MJD
    if keyword_set(getspec) then begin
       out = mrdfits('~/wise_readspec_info.fits.gz',1)
       ngal = n_elements(out)

       allid = string(out.plate,format='(I4.4)')+string(out.mjd,format='(I5)')
       id = allid[uniq(allid,sort(allid))]
       nid = n_elements(id)

       nmaxpix = 4000
       spec = replicate({crval1: 0D, cd1_1: 0D, flux: fltarr(nmaxpix), $
         invvar: fltarr(nmaxpix)},ngal)
       for ii = 0, nid-1 do begin
          these = where(id[ii] eq allid)
          platestr = strmid(id[ii],0,4)
          mjdstr = strmid(id[ii],4,5)
          
          platefile = file_search(getenv('SPECTRO_DATA')+'/'+platestr+'/'+$
            'spPlate-'+platestr+'-'+mjdstr+'.fits*',count=nfile)
          if (nfile ne 1) then message, 'Problem here!'

          splog, 'Reading '+platefile
          flux1 = mrdfits(platefile,0,hdr,row=out[these].fiberid-1,/silent)
          invvar1 = mrdfits(platefile,1,row=out[these].fiberid-1,/silent)

          npix = sxpar(hdr,'NAXIS1')
          spec[these].crval1 = sxpar(hdr,'CRVAL1')
          spec[these].cd1_1 = sxpar(hdr,'CD1_1')
          spec[these].flux[0:npix-1] = flux1*1D-17
          spec[these].invvar[0:npix-1] = invvar1/1D-17^2
       endfor
       
;   readspec, out.plate, out.fiberid, mjd=out.mjd, flux=flux, $
;     invvar=invvar, loglam=loglam, wave=wave, sky=sky, /silent

       spec = struct_addtags(out,spec)
       im_mwrfits, spec, 'wise_sdssspec.fits', /clobber
    endif
       
return
end
    
    
