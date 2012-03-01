pro cmdm_to_maggies, cat, maggies, ivar, filterlist=filterlist

    ngal = n_elements(cat)    

    filterlist = 'subaru_suprimecam_'+['B','V','Rc','Ic','z']+'.par'
    nbands = n_elements(filterlist)

    tags = ['B','V','R','I','Z']
    errtags = 'd'+['B','V','R','I','Z']

; construct maggies and ivarmaggies in each band       
    maggies = dblarr(nbands,ngal)
    ivar = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(cat[0],tags[ib])
       utag = tag_indx(cat[0],errtags[ib])

       good = where(cat.(utag) gt 0.0,ngood)
       if (ngood ne 0L) then begin
          magerr = cat[good].(utag)
          maggies[ib,good] = 10.0^(-0.4*cat[good].(ftag))
          ivar[ib,good] = 1.0/(0.4*alog(10.0)*(maggies[ib,good]*magerr))^2
       endif
    endfor
    
    if (keyword_set(nominerror) eq 0) then begin
       minerr = replicate(0.02,nbands)
       k_minerror, maggies, ivar, minerr
    endif

return   
end


pro cmdm_kcorrect, kcorr, clobber=clobber
; jm12feb27ucsd - compute K-corrections

;   kcorr = mrdfits(cmdmpath+'cmdm_kcorrect.fits.gz',1)
;   kcorrect_qaplot, kcorr, psfile=cmdmpath+'qaplot_cmdm_kcorrect.ps', $
;     in_filterlist=cmdm_filterlist(), /clobber, vname='default.nolines'
    
    path = clash_path()+'projects/cmdm/'
    h100 = 0.7
    vname = 'default'

    kcorrfile = path+'macs1206_mass.fits'
;   kcorrfile = path+'cmdm_kcorrect.fits'
    if im_file_test(kcorrfile+'.gz',clobber=clobber) then return

    cl = 'macs1206'
    zcat1 = read_clash_catalog(cl,/redshift)
    cat1 = rsex(path+'data/'+'MACSJ1206_Subaru.cat')
    spherematch, cat1.ra, cat1.dec, zcat1.ra, zcat1.dec, 1D/3600, m1, m2

    gal = where(zcat1[m2].z gt 0.4 and zcat1[m2].z lt 0.49,ngal)
    cat = cat1[m1[gal]]
    zcat = zcat1[m2[gal]]

    cmdm_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist
    nfilt = n_elements(filterlist)

    kcorr = {$
      id:                                0L,$
      ra:                                0D,$
      dec:                               0D,$
      z:                             -999.0,$
      maggies:                fltarr(nfilt),$
      ivarmaggies:            fltarr(nfilt),$
      bestmaggies:            fltarr(nfilt),$
      mass:                          -999.0,$
      coeffs:                     fltarr(5),$
      chi2:                          -999.0,$

      ugriz_absmag:         fltarr(5)-999.0,$
      ugriz_absmag_ivar:    fltarr(5)-999.0,$
      ugriz_kcorrect:       fltarr(5)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.id = cat.id
    kcorr.ra = cat.ra
    kcorr.dec = cat.dec
    kcorr.z = zcat.z
    kcorr.maggies = maggies
    kcorr.ivarmaggies = ivarmaggies

; compute k-corrections
    out_filterlist = sdss_filterlist()
    kcorrect = im_kcorrect(zcat.z,maggies,ivarmaggies,$
      filterlist,out_filterlist,band_shift=0.0,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=absmag,$
      ivarabsmag=absmag_ivar,uvflux=uvflux,$;clineflux=cflux,uvflux=uvflux,$
      /silent,vname=vname,h100=h100) ; AB, band_shift=0.1
    
    kcorr.bestmaggies = bestmaggies
    kcorr.mass = alog10(mass) ; h=0.7, Chabrier
    kcorr.coeffs = coeffs
    kcorr.chi2 = chi2
    kcorr.ugriz_absmag      = absmag
    kcorr.ugriz_absmag_ivar = absmag_ivar
    kcorr.ugriz_kcorrect    = kcorrect

    im_mwrfits, kcorr, kcorrfile, /clobber

; write out an ascii file for Doron
    openw, lun, path+'macs1206_mass.dat', /get_lun
    printf, lun, '## Preliminary stellar mass catalog for galaxies with 0.4<zspec<0.49 '
    printf, lun, '## in the MACS1206 field.  '
    printf, lun, '## J. Moustakas '+im_today()
    printf, lun, '# 1 ID [from MACSJ1206_Subaru.cat]'
    printf, lun, '# 2 ra [degree]'
    printf, lun, '# 3 dec [degree]'
    printf, lun, '# 4 z'
    printf, lun, '# 5 mass [log_{10}]'
    struct_print, struct_trimtags(kcorr,select=['id','ra','dec','z','mass']), $
      ddigit=10, fdigit=6, lun=lun, /no_head
    free_lun, lun

stop
    
return
end
