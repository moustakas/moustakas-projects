pro fix_tractor_edr
; jm14nov12siena - fix the Tractor catalogs for the EDR

    topdir = '/global/project/projectdirs/cosmo/work/decam/release/edr/'
    outdir = topdir+'tractor/updated-tractor/'
    allcat = file_search(topdir+'tractor/tractor-??????.fits',count=ncat)

; reddening factors from Schlafly & Finkbeiner (2014, private
; communication) for the DECam filter curves
    red_fac = [3.995,3.214,2.165,1.592,1.211,1.064] ; =ugrizY
;   red_fac = [5.155,3.237,2.176,1.595,1.217,1.058] ; =ugrizY
    for ii = 0L, ncat-1 do begin
       cat1 = mrdfits(allcat[ii],1,hdr,/silent)
       ngal = n_elements(cat1)
       
; fix OBJID by sorting on blob and then assigning sequential numbers 
       cat1 = cat1[sort(cat1.blob)]
       cat1.objid = lindgen(ngal)
       
; add the DECam extinction values
       cat1 = struct_addtags(cat1,replicate({decam_nexp: intarr(6), $
         extinction: fltarr(6)},ngal))
       euler, cat1.ra, cat1.dec, ll, bb, 1
       cat1.extinction = red_fac # dust_getval(ll,bb,/interp,/noloop)

; fix the out-of-range ra,dec values
       cat1.ra = (cat1.ra>0D)<360D
       cat1.dec = (cat1.dec>(-90D))<90D

; get the number of exposures per band
       gfile = topdir+'coadd/nexp-'+strtrim(cat1[0].brickid,2)+'-g.fits.gz'
       rfile = topdir+'coadd/nexp-'+strtrim(cat1[0].brickid,2)+'-r.fits.gz'
       zfile = topdir+'coadd/nexp-'+strtrim(cat1[0].brickid,2)+'-z.fits.gz'
       if file_test(gfile) eq 0 or file_test(rfile) eq 0 or file_test(zfile) eq 0 then stop
       nexp_g = mrdfits(gfile,0,/silent)
       nexp_r = mrdfits(rfile,0,/silent)
       nexp_z = mrdfits(zfile,0,/silent)

       cat1.decam_nexp[1] = nexp_g[round(cat1.x),round(cat1.y)]
       cat1.decam_nexp[2] = nexp_r[round(cat1.x),round(cat1.y)]
       cat1.decam_nexp[4] = nexp_z[round(cat1.x),round(cat1.y)]
       
; write out the updated catalog
       outfile = outdir+file_basename(allcat[ii])
       splog, 'Writing ', outfile
       mwrfits, cat1, outfile, /create
    endfor

return
end
    
