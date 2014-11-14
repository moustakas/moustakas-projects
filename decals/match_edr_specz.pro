pro match_edr_specz
; jm14nov13siena - match the EDR Tractor catalogs to the SDSS
; photometry (DR10), SDSS spectroscopy (DR10), and unWISE catalogs

; after this routine run RECALIBRATE_EDR_SDSS
    
; read the coordinates from the photo-dr10 file and keep them in a
; common block for testing purposes
    common com_sdss, sdss_coord

    dr10dir = '/home/data/archive/data/sdss/dr10/'
    topdir = '/global/data/desi/decam/release/edr/'
    unwisedir = dr10dir

;   dr10dir = getenv('HOME')+'/' ; hack!
;   topdir = '/global/project/projectdirs/cosmo/work/decam/release/edr/'
;   unwisedir = '/project/projectdirs/cosmo/data/unwise/unwise-phot/'
    
    if n_elements(sdss_coord) eq 0L then sdss_coord = mrdfits(dr10dir+$
      'photoPosPlate-dr10.fits',1,columns=['ra','dec'])

    allcat = file_search(topdir+'tractor/updated-tractor/tractor-??????.fits',count=ncat)

; loop on each Tractor catalog
    tall = systime(1)
    for ii = 0L, ncat-1 do begin
       cat1 = mrdfits(allcat[ii],1,hdr,/silent)
       ngal = n_elements(cat1)

; match against the SDSS/DR10       
       spherematch, sdss_coord.ra, sdss_coord.dec, cat1.ra, $
         cat1.dec, 1D/3600.0, m1, m2
       
       nmatch = n_elements(m1)*(m1[0] ne -1)
       splog, 'Brick '+strtrim(cat1[0].brickid,2)+', '+strtrim(nmatch,2)+' matches'
       
       if nmatch gt 0L then begin

; sort to preserve the order of the Tractor catalogs
          srt = sort(m2)
          m1 = m1[srt]
          m2 = m2[srt]
          tractor1 = struct_trimtags(cat1[m2],except='sdss_*')

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
    spec = mrdfits(dr10dir+'specObj-dr10.fits',1,rows=sdssindx)
    phot = mrdfits(dr10dir+'photoPosPlate-dr10.fits',1,rows=sdssindx)
    unwise = struct_trimtags(mrdfits(unwisedir+'specmatch-dr10.fits',1,$
      rows=sdssindx),select='*wise*')

; generate the output structure with the following format:
;    HDU #0: Tractor outputs (perhaps excluding the SDSS* tags as duplicative to the below)
;    HDU #1: SDSS spectro (from the specObj file)
;    HDU #2: SDSS imaging (from the photoPosPlate file)
;    HDU #3: WISE imaging (the WISE* tags in Dustin's file above)

    outfile = topdir+'edr-specz-dr10-oldcalib.fits'
    splog, 'Writing '+outfile
    mwrfits, tractor, outfile, /create
    mwrfits, spec, outfile
    mwrfits, phot, outfile
    mwrfits, unwise, outfile
    
stop    

return
end
    
