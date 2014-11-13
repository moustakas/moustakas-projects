pro decals_merge_tractor

    common com_sdss, specmatch

    dr10dir = getenv('IM_ARCHIVE_DIR')+'/data/sdss/dr10/'
    if n_elements(specmatch) eq 0L then specmatch = mrdfits(dr10dir+'specmatch-dr10.fits.gz',1)

    tractordir = '/global/data/desi/decam/release/edr/tractor/'
    allfile = file_search(tractordir+'tractor-??????.fits',count=nfile)

;   for ii = 0L, nfile-1 do begin
;      cat1 = mrdfits(allfile[ii],1,hdr,/silent)
;      if max(cat1.ra) gt 360.0 or abs(max(cat1.dec)) gt 90 then print, file_basename(allfile[ii])
;   endfor       

;   for ii = 0L, 10 do begin
    for ii = 0L, nfile-1 do begin
       cat1 = mrdfits(allfile[ii],1,hdr,/silent)
       spherematch, specmatch.plug_ra, specmatch.plug_dec, (cat1.ra<360)>0, (cat1.dec<90)>(-90), $
         1D/3600.0, m1, m2
       nmatch = n_elements(m1)*(m1[0] ne -1)
       splog, 'Brick '+strtrim(cat1[0].brickid,2)+', '+strtrim(nmatch,2)+' matches'

       if nmatch gt 0L then begin
          out1 = struct_addtags(cat1[m2],struct_trimtags(specmatch[m1],$
            select=['plate','mjd','fiberid','class','subclass','z','z_err',$
            'vdisp','vdisp_err','zwarning','wise_*nanomaggies*']))
          if n_elements(out) eq 0L then out = out1 else out = [temporary(out),out1]
       endif
    endfor

    im_mwrfits, out, '/global/data/desi/decam/isedfit/decals_edr.fits', /clobber
    
stop    

return
end
