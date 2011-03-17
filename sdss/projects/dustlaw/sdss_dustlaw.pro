pro sdss_dustlaw, gather=gather, deredshift=deredshift, $
  twomass=twomass, hammer=hammer, merge=merge
; jm10feb02ucsd - build the sample to measure the dust law above the
; plane 

    common com_dust, allphot, allspec
    
    outpath = sdss_path(/projects)+'dustlaw/'

; --------------------------------------------------
; gather the relevant data and write out    
    if keyword_set(gather) then begin
       if (n_elements(allphot) eq 0) then begin
          allphot = mrdfits(vagc_name('object_sdss_imaging'),1)
          allspec = mrdfits(vagc_name('object_sdss_spectro'),1)
;         alltwomass = mrdfits(vagc_name('object_twomass'),1)
       endif
       stars = where($
         (allspec.sdss_spectro_tag ne -1) and $
         (strtrim(allspec.class,2) eq 'STAR') and $
         (allspec.sn_median gt 3.0) and $
         (allspec.zwarning eq 0),nstar)
;      stars = stars[0:100] & nstar = n_elements(stars)
       phot = allphot[stars]
       spec = allspec[stars]
;      twomass = alltwomass[stars]
       splog, 'Nstar = ', nstar

; to minimize I/O overhead in readspec, sort by the combination of
; plate+MJD       
       srt = sort(string(spec.plate,format='(I3.3)')+'-'+$
         string(spec.mjd,format='(I5.5)'))
;      srt = sort(strtrim(spec.subclass,2))
       phot = phot[srt]
       spec = spec[srt]

       chunksize = 10000L ; spectra per chunk
       nchunk = ceil(nstar/float(chunksize))
       chunkfiles = 'stars_specdata_chunk_'+string(lindgen(nchunk)+1L,$
         format='(I3.3)')+'.fits'

; add the chunk filename and write out       
       spec = struct_addtags(temporary(spec),$
         replicate({chunkfile: '', chunkfile_rest: ''},nstar))
       for ichunk = 0L, nchunk-1L do begin
          i1 = ichunk*chunksize
          i2 = ((ichunk*chunksize+chunksize)<nstar)-1L
          nthese = i2-i1+1L
          these = lindgen(nthese)+i1
          spec[these].chunkfile = chunkfiles[ichunk]+'.gz'
          spec[these].chunkfile_rest = repstr(spec[these].chunkfile,'_specdata','_specdata_rest')
       endfor
       im_mwrfits, phot, outpath+'stars_photinfo.fits', /clobber
       im_mwrfits, spec, outpath+'stars_specinfo.fits', /clobber

; grab the spectra in chunks
;      for ichunk = 0, 9 do begin
;      for ichunk = 10, 19 do begin
;      for ichunk = 20, 29 do begin
       for ichunk = 0L, nchunk-1L do begin
          splog, 'Chunk = '+string(ichunk+1,format='(I2.2)')
          i1 = ichunk*chunksize
          i2 = ((ichunk*chunksize+chunksize)<nstar)-1L
          nthese = i2-i1+1L
          these = lindgen(nthese)+i1

          readspec, spec[these].plate, spec[these].fiberid, $
            mjd=spec[these].mjd, flux=flux, invvar=invvar, $
            loglam=loglam, /align
          npix = n_elements(lambda)
          data = {flux: flux, invvar: invvar, loglam: loglam, $
            object_position: these}
          im_mwrfits, data, outpath+chunkfiles[ichunk], /clobber
       endfor
    endif    

; --------------------------------------------------
; deredshift and interpolate onto a common wavelength scale
    if keyword_set(deredshift) then begin
; define the output wavelength scale (3810-9200 A)
       newloglam = dindgen((alog10(9200.0D)-alog10(3810.0D))/1D-4)*1D-4+alog10(3810.0D)
       npix = n_elements(newloglam)
       
; now read and parse the data
       info = mrdfits(outpath+'stars_specinfo.fits.gz',1)
       allchunkfiles = strtrim(info.chunkfile,2)
       chunkfiles = allchunkfiles[uniq(allchunkfiles,sort(allchunkfiles))]
       nchunk = n_elements(chunkfiles)

;; junk code!
;       for ichunk = 0, nchunk-1 do begin
;          splog, 'Chunk = '+string(ichunk+1,format='(I2.2)')
;          these = where(chunkfiles[ichunk] eq allchunkfiles,nthese)
;          rawdata = mrdfits(outpath+chunkfiles[ichunk],1)
;          out = struct_addtags(rawdata,{object_position: these})
;          im_mwrfits, out, repstr(outpath+chunkfiles[ichunk],'.gz',''), /clobber
;       endfor

;      for ichunk = 8, 9 do begin
;      for ichunk = 27, 29 do begin
;      for ichunk = 1, 9 do begin
;      for ichunk = 10, 19 do begin
;      for ichunk = 20, 29 do begin
       for ichunk = 0, nchunk-1 do begin
          splog, 'Chunk = '+string(ichunk+1,format='(I2.2)')
          these = where(chunkfiles[ichunk] eq allchunkfiles,nthese)
          rawdata = mrdfits(outpath+chunkfiles[ichunk],1)
;         newloglam = rawdata.loglam
;         npix = n_elements(newloglam)
          data = {loglam: newloglam, flux: fltarr(npix,nthese), $
            invvar: fltarr(npix,nthese), $
            object_position: rawdata.object_position}
          for ii = 0L, nthese-1L do begin
             print, format='("Object ",I5.5,"/",I5.5,A10,$)', $
               ii, nthese-1, string(13b)
             loglam = rawdata.loglam - alog10(info[these[ii]].z+1.0D) ; deredshift
             flux = rawdata.flux[*,ii]*(1.0+info[these[ii]].z)
             invvar = rawdata.invvar[*,ii]/(1.0+info[these[ii]].z)^2.0
             combine1fiber, loglam, flux, invvar, newloglam=newloglam, $
               newflux=newflux, newivar=newivar
             data.flux[*,ii] = newflux
             data.invvar[*,ii] = newivar
          endfor
; write out
          outfile = repstr(chunkfiles[ichunk],'_specdata','_specdata_rest')
          outfile = repstr(outfile,'.fits.gz','.fits')
          im_mwrfits, data, outpath+outfile, /clobber
       endfor
    endif

; --------------------------------------------------
; parse the 2MASS/PSC input/output
    if keyword_set(twomass) then begin
       photinfo = mrdfits(outpath+'stars_photinfo.fits.gz',1)
       nstar = n_elements(photinfo)
;      iquery_2mass, photinfo.ra, photinfo.dec, 'stars_twomass_input.tbl'    
;      tbl = im_read_tbl(outpath+'stars_twomass_output.tbl')
;      im_mwrfits, tbl, outpath+'stars_twomass_output.fits', /clobber
       tmass = mrdfits(outpath+'stars_twomass_output.fits.gz',1)

       spherematch, photinfo.ra, photinfo.dec, tmass.ra, $
         tmass.dec, 1.0/3600.0, m1, m2
       out = im_empty_structure(tmass[0],$
         empty_value=-999.0,ncopies=nstar)
       out[m1] = tmass[m2]

       im_mwrfits, out, outpath+'stars_twomass.fits', /clobber
    endif

;; --------------------------------------------------
;; classify
;    if keyword_set(hammer) then begin
;       data = mrdfits('stars_specdata_rest.fits.gz',1)
;       im_hammer, data, out
;    endif

; --------------------------------------------------
; merge the most useful info into one structure
    if keyword_set(merge) then begin
       specinfo = mrdfits(outpath+'stars_specinfo.fits.gz',1)
       photinfo = mrdfits(outpath+'stars_photinfo.fits.gz',1)
       twomass = mrdfits(outpath+'stars_twomass.fits.gz',1)
       nstar = n_elements(specinfo)

;      morephot = replicate({maggies: fltarr(8), ivarmaggies: fltarr(8)},nstar)
       morephot = replicate({ugriz: fltarr(5), ugriz_err: fltarr(5), $
         jhk: fltarr(3), jhk_err: fltarr(3)},nstar)

       aboff = [-0.036, 0.012, 0.010, 0.028, 0.040]
       for ii = 0, 4 do begin
          good = where((photinfo.psfflux[ii] gt 0.0) and $
            (photinfo.psfflux_ivar[ii] gt 0.0),ngood)
          morephot[good].ugriz[ii] = maggies2mag(photinfo[good].psfflux[ii]*1E-9,$
            ivarmaggies=photinfo[good].psfflux_ivar[ii]*1E18,magerr=magerr)
          morephot[good].ugriz[ii] = morephot[good].ugriz[ii] + aboff[ii]
          morephot[good].ugriz_err[ii] = magerr
       endfor

       tags = ['j','h','k']
       for ii = 0, 2 do begin
          ftag = tag_indx(twomass,tags[ii]+'_m')
          utag = tag_indx(twomass,tags[ii]+'_msigcom')
          good = where((twomass.(ftag) gt 0.0) and $
            (twomass.(utag) gt 0.0),ngood)
          morephot[good].jhk[ii] = twomass[good].(ftag)
          morephot[good].jhk_err[ii] = twomass[good].(utag)
       endfor
          
;      sdss_to_maggies, smaggies, ivarsmaggies, calib=photinfo, flux='psf'
;      twomass_to_maggies_psc, twomass, tmaggies, ivartmaggies
;      morephot.maggies[0:4] = smaggies
;      morephot.ivarmaggies[0:4] = ivarsmaggies
;      morephot.maggies[5:7] = tmaggies
;      morephot.ivarmaggies[5:7] = ivartmaggies

       spectags = ['ra','dec','plate','mjd','fiberid','class','subclass',$
         'z','z_err','vdisp','vdisp_err','sn_median','chunkfile_rest']
       phottags = ['extinction']
       out = struct_trimtags(specinfo,select=spectags)
       out = struct_addtags(temporary(out),$
         struct_trimtags(photinfo,select=phottags))

       out = struct_addtags(temporary(out),morephot)
       im_mwrfits, out, outpath+'stars_info.fits', /clobber
stop
    endif

return
end
