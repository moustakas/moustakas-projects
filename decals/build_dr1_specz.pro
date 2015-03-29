pro merge_tractor_dr1, merge=merge, specz=specz
; jm15mar11siena - merge the DR1 tractor catalogs    

    common com_sdss, sdss_coord

    dr1dir = getenv('DECALS_DIR')+'/'
    isedfit_dir = getenv('IM_ARCHIVE_DIR')+'/decam/isedfit/dr1/'
    sdssdr12dir = '/home/work/data/sdss/dr12/'
    outfile = dr1dir+'dr1-specz.fits'

    if n_elements(sdss_coord) eq 0L then sdss_coord = mrdfits(sdssdr12dir+$
      'photoPosPlate-dr12.fits',1,columns=['ra','dec'])

    allbrick = file_basename(file_search(dr1dir+'tractor/*',/test_dir,count=nbrick))

    tall = systime(1)
;   for ii = 93, 94 do begin
    for ii = 0L, nbrick-1 do begin
       delvarx, cat
       catfile = file_search(dr1dir+'tractor/'+allbrick[ii]+'/tractor-*.fits',count=ncat)
       for ic = 0L, ncat-1 do begin
          print, format='("Brick ",I0,"/",I0," Cat ",I0,"/",I0, A10,$)', $
            ii+1, nbrick, ic+1, ncat, string(13b)
          cat1 = mrdfits(catfile[ic],1,/silent)
          cat1 = cat1[where(cat1.brick_primary eq 'T')]
          cat1 = struct_trimtags(temporary(cat1),except='sdss_*')
          if n_elements(cat) eq 0L then cat = cat1 else cat = [temporary(cat),temporary(cat1)]
       endfor          
       ngal = n_elements(cat)
          
; match against the SDSS/DR10       
       spherematch, sdss_coord.ra, sdss_coord.dec, cat.ra, $
         cat.dec, 1D/3600.0, m1, m2
       nmatch = n_elements(m1)*(m1[0] ne -1)
       splog, 'BigBrick '+strtrim(allbrick[ii],2)+', '+strtrim(nmatch,2)+' matches'
          
       if nmatch gt 0L then begin
; sort to preserve the order of the Tractor catalogs
          srt = sort(m2)
          m1 = m1[srt]
          m2 = m2[srt]
          tractor1 = cat[m2]
             
; build the SDSS indices
          if n_elements(tractor) eq 0L then begin
             tractor = temporary(tractor1)
             sdssindx = m1
          endif else begin
             tractor = [temporary(tractor),temporary(tractor1)]
             sdssindx = [sdssindx,m1]
          endelse
       endif
    endfor
    splog, 'Time to match (min) ', (systime(1)-tall)/60.0

; get the ancillary data
    spec = mrdfits(sdssdr12dir+'specObj-dr12.fits',1,rows=sdssindx)
    phot = mrdfits(sdssdr12dir+'photoPosPlate-dr12.fits',1,rows=sdssindx)

; generate the output structure with the following format:
;    HDU #0: Tractor outputs (perhaps excluding the SDSS* tags as duplicative to the below)
;    HDU #1: SDSS spectro (from the specObj file)
;    HDU #2: SDSS imaging (from the photoPosPlate file)

    splog, 'Writing '+outfile
    mwrfits, tractor, outfile, /create
    mwrfits, spec, outfile
    mwrfits, phot, outfile
    
stop    

return
end
