pro build_sg1120_sdss_refcat, sdss
; jm09may04nyu - build a (scamp-compatible) astrometric and
;   photometric reference catalog for SG1120

; http://cas.sdss.org/dr7/en/tools/search/sql.asp
; SELECT
;   p.run, p.rerun, p.camcol, p.field, p.obj,
;   p.type, p.ra, p.dec, p.raerr, p.decerr, p.petror90_g, p.petror90_r, 
;   p.psfmag_u,p.psfmag_g,p.psfmag_r,p.psfmag_i,p.psfmag_z,
;   p.psfmagerr_u,p.psfmagerr_g,p.psfmagerr_r,p.psfmagerr_i,p.psfmagerr_z,
;   p.extinction_u,p.extinction_g,p.extinction_r,p.extinction_i,p.extinction_z
;   FROM fGetNearbyObjEq(170.0,-12.05,15) n, PhotoPrimary p
;   WHERE n.objID=p.objID  AND p.psfmag_r BETWEEN 10 AND 30

    sexpath = sg1120_path(/sex)
    sdssfile = 'sg1120_sdss_dr7'
    
    readcol, sexpath+sdssfile+'.csv', run, obj, type, $
      ra, dec, raerr, decerr, rad_g, rad_r, $
      u, g, r, i, z, uerr, gerr, rerr, ierr, zerr, $
      udust, gdust, rdust, idust, zdust, $
      format='L,I,I,D,D,D,D,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,', $
      delimiter=',', /silent, skipline=1
;;; pack into a SDSS-style structure that can be parsed by
;;; SDSS_TO_MAGGIES
;;    sdss = {ra: 0.0D, dec: 0.0D, raerr: 0.0D, decerr: 0.0D, $
;;      radius_g: 0.0, radius_r: 0.0, type: -1, $
;;      psfmag_u: 0.0, psfmag_g: 0.0, psfmag_r: 0.0, psfmag_i: 0.0, psfmag_z: 0.0, $
;;      psfmagerr_u: 0.0, psfmagerr_g: 0.0, psfmagerr_r: 0.0, psfmagerr_i: 0.0, psfmagerr_z: 0.0, $
;;      extinction_u: 0.0, extinction_g: 0.0, extinction_r: 0.0, extinction_i: 0.0, extinction_z: 0.0}
;;    sdss = replicate(sdss,n_elements(r))
;;    sdss.ra = ra
;;    sdss.dec = dec
;;    sdss.raerr = raerr/3600.0D   ; [deg]
;;    sdss.decerr = decerr/3600.0D ; [deg]
;;    sdss.radius_g = rad_g
;;    sdss.radius_r = rad_r
;;    sdss.type = type
;;    sdss.psfmag_u = u
;;    sdss.psfmag_g = g
;;    sdss.psfmag_r = r
;;    sdss.psfmag_i = i
;;    sdss.psfmag_z = z
;;    sdss.psfmagerr_u = uerr
;;    sdss.psfmagerr_g = gerr
;;    sdss.psfmagerr_r = rerr
;;    sdss.psfmagerr_i = ierr
;;    sdss.psfmagerr_z = zerr
;;    sdss.extinction_u = udust
;;    sdss.extinction_g = gdust
;;    sdss.extinction_r = rdust
;;    sdss.extinction_i = idust
;;    sdss.extinction_z = zdust
;;    im_mwrfits, sdss, sexpath+sdssfile+'.fits'
;;
;;    sdss_to_maggies, maggies, ivarmaggies, cas=sdss, flux='psf'
;;    rmag = maggies2mag(reform(maggies[2,*]),magerr=rmagerr,$
;;      ivarmaggies=reform(ivarmaggies[2,*]))

   sdss = {ra: 0.0D, dec: 0.0D, raerr: 0.0D, decerr: 0.0D, $
     radius_g: 0.0, radius_r: 0.0, type: -1, $
     u: 0.0, g: 0.0, r: 0.0, i: 0.0, z: 0.0, $
     uerr: 0.0, gerr: 0.0, rerr: 0.0, ierr: 0.0, zerr: 0.0, $
     udust: 0.0, gdust: 0.0, rdust: 0.0, idust: 0.0, zdust: 0.0}
   sdss = replicate(sdss,n_elements(r))
   sdss.ra = ra
   sdss.dec = dec
   sdss.raerr = raerr/3600.0D   ; [deg]
   sdss.decerr = decerr/3600.0D ; [deg]
   sdss.radius_g = rad_g
   sdss.radius_r = rad_r
   sdss.type = type
; add the K-correct AB offsets (see k_abfix), but *do not* correct for
; extinction 
   sdss.u = u + (-0.036)
   sdss.g = g + (+0.012)
   sdss.r = r + (+0.010)
   sdss.i = i + (+0.028)
   sdss.z = z + (+0.040)
   sdss.uerr = uerr
   sdss.gerr = gerr
   sdss.rerr = rerr
   sdss.ierr = ierr
   sdss.zerr = zerr
   sdss.udust = udust
   sdss.gdust = gdust
   sdss.rdust = rdust
   sdss.idust = idust
   sdss.zdust = zdust
   im_mwrfits, sdss, sexpath+sdssfile+'.fits'
    
; write a DS9 region file of *all* sources
    write_ds9_regionfile, sdss.ra, sdss.dec, color='red', $
      symbol='circle', file=sexpath+sdssfile+'.reg'

; build a scamp-compatible reference catalog    
    write_scamp_catalog, sdss.ra, sdss.dec, sdss.r, ra_err=sdss.raerr, $
      dec_err=sdss.decerr, mag_err=sdss.rerr, outcat=outcat, $
      outfile=sexpath+sdssfile+'_refcat.cat'

return
end
    
