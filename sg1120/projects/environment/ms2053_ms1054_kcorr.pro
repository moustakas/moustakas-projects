pro hst_to_maggies, cat, maggies, ivar, ebv=ebv, filterlist=filterlist

    names = 'm'+['606','814']+'_mag_aper'
    errnames = 'm'+['606','814']+'_magerr_aper'
    filterlist = ['wfpc2_f606w.par','wfpc2_f814w.par']
    minerrors = [0.02,0.02]

    vega2ab = k_vega2ab(filterlist=filterlist,/kurucz,/silent)
    
    filt = im_filterspecs(filterlist=filterlist)
    kl = k_lambda(filt.weff,/odonnell,R_V=3.1)

    euler, cat.m814_ra, cat.m814_dec, gl, gb, 1
    ebv = dust_getval(gl,gb,/interp)

    nband = n_elements(names)
    maggies = dblarr(nband,n_elements(cat))
    ivar = dblarr(nband,n_elements(cat))
    for iband = 0L, nband-1 do begin
       itag_m = tag_indx(cat[0],names[iband]+'')
       itag_msig = tag_indx(cat[0],errnames[iband])
       good = where((cat.(itag_m) gt 0.0) and $
         (cat.(itag_msig) ge 0.0) and (cat.(itag_m) lt 90.0) and $
         (cat.(itag_msig) lt 90.0) and (cat.(itag_m) gt -90.0) and $
         (cat.(itag_msig) gt -90.0),count)
       if (count gt 0) then begin
          mag = cat[good].(itag_m) - kl[iband]*ebv[good] + vega2ab[iband]
          magerr = cat[good].(itag_msig) 
          maggies[iband,good] = 10.0^(-0.4*mag)
          notzero = where((maggies[iband,good] gt 0.0),nnotzero)
          if (nnotzero ne 0L) then ivar[iband,good[notzero]] = $
            1.0/(0.4*alog(10.0)*(maggies[iband,good[notzero]]*magerr[notzero])^2)
       endif
    endfor

    k_minerror, maggies, ivar, minerrors

return
end

;function hst_kcorrect, cat
;
;    sfhpath = getenv('ISEDFIT_SFHGRID_DIR')+'/basemodels/bc03/'
;
;    ngal = n_elements(cat)
;    hst_to_maggies, cat, maggies, ivarmaggies, $
;      filterlist=in_filterlist, ebv=ebv
;    out_filterlist = 'bessell_'+['U','B','V']+'.par'
;    noutband = n_elements(out_filterlist)
;
;; do the fitting on a grid of tau
;    tau = [im_array(0.0,10.0,1.0),100]
;;   tau = [im_array(0.0,10.0,0.5),100]
;    taufile = 'chab_Z0.02_tau_'+tau2string(tau)+'.fits.gz'
;    ntau = n_elements(tau)
;
;    age = getage(cat.redshift) ; z_f=infinity
;    dist = 10.0*3.085678D18 ; 10 pc fiducial distance [cm]
;
;    tauchi2 = fltarr(ntau,ngal)
;    absmag = fltarr(noutband,ntau,ngal)
;    for ii = 0, ntau-1 do begin    
;       sed = mrdfits(sfhpath+taufile[ii],1)
;       get_element, sed.age/1E9, age, indx
;       restwave = sed.wave
;       
;       restflux = 3.826D33*sed.flux[*,indx]/(4*!dpi*dist^2.0) ; [erg/s/cm2/A/M_sun]
;
;       kk = im_simple_kcorrect(cat.redshift,maggies,ivarmaggies,$
;         in_filterlist,out_filterlist,restwave,restflux,absmag=absmag1,$
;         chi2=chi2,scale=scale,band_shift=0.0,/vega,/silent)
;         
;       absmag[*,ii,*] = absmag1
;       tauchi2[ii,*] = chi2
;    endfor
;
;    out = {$
;      ebv:        0.0, $
;      tau:        0.0, $
;      age:        0.0, $
;      chi2:       0.0, $
;      ubv_absmag: fltarr(noutband)}
;    out = struct_addtags(struct_trimtags(cat,select=['runid','run',$
;      'hst_id','keck_id','redshift','q']),replicate(out,ngal))
;    out.ebv = ebv
;    out.age = age
;    
;    for jj = 0, ngal-1 do begin
;       out[jj].chi2 = min(tauchi2[*,jj],mindx)
;       out[jj].tau = tau[mindx]
;       out[jj].ubv_absmag = absmag[*,mindx,jj]
;    endfor
;
;return, out
;end    

pro ms2053_ms1054_kcorr
; jm10jul25ucsd - compute rest-frame quantities for Vy's MS1054
; and MS2053 MIPS studies

    path = sg1120_path()+'projects/ms2053_ms1054/'
    filt = 'bessell_'+['U','B','V']+'.par'
    nfilt = n_elements(filt)
    
;; #########################
;; MS1054
;    hst1054 = rsex(path+'hst1054z.cat')
;    good = where((hst1054.redshift gt 0.0) and (hst1054.m606_mag_aper gt 0) and $
;      (hst1054.m814_mag_aper gt 0))
;;   kcorr1054 = hst_kcorrect(hst1054[good])
;
;    hst_to_maggies, hst1054[good], maggies, ivarmaggies, $
;      filterlist=filterlist, ebv=ebv
;    kcorr = im_kcorrect(hst1054[good].redshift,maggies,ivarmaggies,$
;      filterlist,filt,band_shift=0.0,absmag=absmag,ivarabsmag=absmag_ivar,$
;      chi2=chi2,/vega,/silent,vname=vname)
;    out = struct_addtags(struct_trimtags(hst1054,select=[$
;       'hst_id','keck_id','redshift','q']),replicate({ubv_absmag: fltarr(nfilt)-999.0, $
;      kcorr: fltarr(nfilt)-999.0, chi2: -999.0},n_elements(hst1054)))
;    out[good].ubv_absmag = absmag
;    out[good].kcorr = kcorr
;    out[good].chi2 = chi2
;
;    im_mwrfits, out, path+'kcorr_hst1054.fits', /clobber
    
; #########################
; MS2053
    hst2053 = rsex(path+'hst2053z.cat')
    good = where((hst2053.redshift gt 0.0) and (hst2053.m606_mag_aper gt 0) and $
      (hst2053.m814_mag_aper gt 0))

    hst_to_maggies, hst2053[good], maggies, ivarmaggies, $
      filterlist=filterlist, ebv=ebv
;   readcol, '2053_table.cat', hstid, keckid, z, m1, i1, m2, i2, format='X,X,A,A,F,F,F,F,F', comment='#'
;   match, strtrim(hst2053[good].hst_id,2)+strtrim(hst2053[good].keck_id,2), $
;     strtrim(hstid,2)+strtrim(keckid,2), j1, j2
;   niceprint, maggies[0,j1], m1[j2], maggies[1,j1], m2[j2]

    kcorr = im_kcorrect(hst2053[good].redshift,maggies,ivarmaggies,$
      filterlist,filt,band_shift=0.0,absmag=absmag,ivarabsmag=absmag_ivar,$
      chi2=chi2,/vega,/silent,vname=vname)
    out = struct_addtags(struct_trimtags(hst2053,select=[$
       'hst_id','keck_id','redshift','q']),replicate({ubv_absmag: fltarr(nfilt)-999.0, $
      kcorr: fltarr(nfilt)-999.0, chi2: -999.0},n_elements(hst2053)))
    out[good].ubv_absmag = absmag
    out[good].kcorr = kcorr
    out[good].chi2 = chi2

    im_mwrfits, out, path+'kcorr_hst2053.fits', /clobber
    
return
end
    
