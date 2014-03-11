pro build_desi_tkrs_sample
; jm13jun15siena - build a sample of TKRS spectra for DESI

    rootpath = getenv('IM_PROJECTS_DIR')+'/desi/'
    dataroot = 'tkrs'

; read the quality file
    qualfile = rootpath+dataroot+'/'+'list.all.id.rcqual.yrange'
    readcol, qualfile, id, quality, prob, minrow, maxrow, $
      sizepix, sizearcsec, format='L,I,I,X,F,F,F,F', /silent
    info = replicate({id: 0L, quality: 0, prob: 0, minrow: 0, $
      maxrow: 0, sizepix: 0, sizearcsec: 0.0},n_elements(id))
    info.id = id
    info.quality = quality
    info.prob = prob
    info.minrow = minrow
    info.maxrow = maxrow
    info.sizepix = sizepix
    info.sizearcsec = sizearcsec

; keep just Q>=3 objects at z~1    
    good = where(info.quality ge 3,ngood)
    
    zcat = read_tkrs(/kcorr) ; includes k-corrections and ZCAT
    match, zcat.id, info[good].id, m1, m2
    zcat = zcat[m1]
    info = info[good[m2]]

    these = where(zcat.z gt 0.9 and zcat.z lt 1.1,ngal)
    zcat = zcat[these]
    info = info[these]

; add info about the rotation curves
    info = struct_addtags(info,replicate({prefix: '', $
      rotcurvefile: '', fitdatafile: ''},ngal))
    info.prefix = string(zcat.id,format='(I7.7)')
    info.rotcurvefile = dataroot+'/'+file_basename(file_search(rootpath+$
      dataroot+'/'+info.prefix+'_*.rotcurv'))
    info.fitdatafile = dataroot+'/'+file_basename(file_search(rootpath+$
      dataroot+'/'+info.prefix+'_*.fitdata'))

; read all the fitdata files
    for ii = 0, ngal-1 do begin
       readcol, rootpath+info[ii].fitdatafile, id, flux, sigma, $
         sigma_cor, vsyst, radius, radius_err, vrot, vrot_err, $
         sigma2d, sigma2d_err, ngood, chi2, /silent
       if ii eq 0 then begin
          fitdata = replicate({$
            id: 0L,$
            flux: 0.0, $
            sigma: 0.0, $
            sigma_cor: 0.0, $
            vsyst: 0.0, $
            radius: 0.0, $
            radius_err: 0.0, $
            vrot: 0.0, $
            vrot_err: 0.0, $
            sigma2d: 0.0, $
            sigma2d_err: 0.0, $
            ngood: 0.0, $
            chi2: 0.0},ngal)
       endif
       fitdata[ii].id = id 
       fitdata[ii].flux = flux
       fitdata[ii].sigma = sigma
       fitdata[ii].sigma_cor = sigma_cor
       fitdata[ii].vsyst = vsyst
       fitdata[ii].radius = radius
       fitdata[ii].radius_err = radius_err
       fitdata[ii].vrot = vrot
       fitdata[ii].vrot_err = vrot_err
       fitdata[ii].sigma2d = sigma2d
       fitdata[ii].sigma2d_err = sigma2d_err
       fitdata[ii].ngood = ngood
       fitdata[ii].chi2 = chi2
    endfor

; pack it all into one structure
    info = struct_addtags(info,struct_trimtags(fitdata,except='id'))

    outroot = rootpath+dataroot
    im_mwrfits, zcat, outroot+'_zcat.fits', /clobber
    im_mwrfits, info, outroot+'_info.fits', /clobber
    
stop    
return
end
