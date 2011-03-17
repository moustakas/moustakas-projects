pro hst_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist

    names = 'mauto'+['606','814']
    filterlist = ['wfpc2_f606w.par','wfpc2_f814w.par']
    minerrors = [0.02,0.02]

    vega2ab = k_vega2ab(filterlist=filterlist,/kurucz,/silent)
    
    filt = im_filterspecs(filterlist=filterlist)
    kl = k_lambda(filt.weff,/odonnell,R_V=3.1)

;   euler, 15.0D*hms2dec('10:54:00'), hms2dec('-03:21:00'), gl, gb, 1
;   ebv1 = dust_getval(gl,gb,/interp)
;   ebv = replicate(ebv1,n_elements(cat))
    ebv = cat.ebv

    nband = n_elements(names)
    maggies = dblarr(nband,n_elements(cat))
    ivarmaggies = dblarr(nband,n_elements(cat))
    for iband = 0L, nband-1 do begin
       itag_m = tag_indx(cat[0],names[iband]+'')
       itag_msig = tag_indx(cat[0],names[iband]+'_err')
       indx = where((cat.(itag_m) gt 0.0) and $
         (cat.(itag_msig) ge 0.0) and (cat.(itag_m) lt 90.0) and $
         (cat.(itag_msig) lt 90.0) and (cat.(itag_m) gt -90.0) and $
         (cat.(itag_msig) gt -90.0),count)
       if (count gt 0) then begin
          maggies[iband,indx] = 10.0^(-0.4*(cat[indx].(itag_m)-$
            kl[iband]*ebv[indx]+vega2ab[iband]))
          sig = cat[indx].(itag_msig)> 0.001
          notzero = where((maggies[iband,indx] gt 0.0),nnotzero)
          if (nnotzero ne 0L) then ivarmaggies[iband,indx[notzero]] = $
            1.0/(0.4*alog(10.0)*(maggies[iband,indx[notzero]]*sig)^2)
       endif
    endfor

    k_minerror, maggies, ivarmaggies, minerrors

return
end

function best_kcorrect, cat

    sfhpath = getenv('ISEDFIT_SFHGRID_DIR')+'/basemodels/bc03/'

    ngal = n_elements(cat)
    hst_to_maggies, cat, obsmaggies, ivarobsmaggies, $
      filterlist=in_filterlist
    out_filterlist = 'bessell_'+['U','B','V']+'.par'
    noutband = n_elements(out_filterlist)
    
; do the fitting on a grid of tau
    tau = [im_array(0.0,10.0,0.5),100]
;   tau = [im_array(0.0,10.0,1.0),100] ; 0.5)
    taufile = 'salp_Z0.02_tau_'+tau2string(tau)+'.fits.gz'
    ntau = n_elements(tau)

    age = getage(cat.z) ; z_f=infinity
    dist = 10.0*3.085678D18 ; 10 pc fiducial distance [cm]

    tauchi2 = fltarr(ntau,ngal)
    absmag = fltarr(noutband,ntau,ngal)
    for ii = 0, ntau-1 do begin    
    
       sed = mrdfits(sfhpath+taufile[ii],1)
       get_element, sed.age/1E9, age, indx
       restwave = sed.wave
       
       restflux = 3.826D33*sed.flux[*,indx]/(4*!dpi*dist^2.0) ; [erg/s/cm2/A/M_sun]

       kk = im_simple_kcorrect(cat.z,obsmaggies,ivarobsmaggies,$
         in_filterlist,out_filterlist=out_filterlist,restwave,restflux,$
         absmag=absmag1,chi2=chi2,rmaggies=rmaggies,scale=scale,/silent,$
         band_shift=0.0,vega=1)
       absmag[*,ii,*] = absmag1
       tauchi2[ii,*] = chi2

    endfor

    out = {$
      ebv:        0.0, $
      tau:        0.0, $
      age:        0.0, $
      chi2:       0.0, $
      ubv_absmag: fltarr(noutband)}
    out = struct_addtags(struct_trimtags(cat,$
      select=['name','z']),replicate(out,ngal))
    out.ebv = cat.ebv
    out.age = age
    
    for jj = 0, ngal-1 do begin
       out[jj].chi2 = min(tauchi2[*,jj],mindx)
       out[jj].tau = tau[mindx]
       out[jj].ubv_absmag = absmag[*,mindx,jj]
    endfor

return, out
end    

pro sg1120_field_kcorrect, parse_cl1358=parse_cl1358
; jm09jun27nyu - compute rest-frame quantities for the field galaxies
; in cl1054 and cl2053 for Vy's environment paper

;    E       -5
;    E/S0    -4
;    S0      -2
;    S0/a     0
;    Sa       1
;    Sb       3
;    Sc       5
;    Sd       7
;    Merger  99
;    
;    anything w/fabtype<-10 was not classified.

    path = sg1120_path()+'projects/environment/'

    if keyword_set(parse_cl1358) then begin
       cat = rsex(path+'cl1358.original.cat')
       ngal = n_elements(cat)
       out = {name: 0L, ra: 0.0D, dec: 0.0D, z: 0.0, $
         absv: 0.0, bmv: 0.0, massb: 0.0, $
         mauto606: 0.0, mauto606_err: 0.0, $
         mauto814: 0.0, mauto814_err: 0.0}
       out = replicate(out,ngal)
       struct_assign, cat, out, /nozero
       out.name = cat.id
       out.mauto814 = cat.imag + 0.098 ; Imag needs a zeropoint offset!
       out.mauto606 = out.mauto814+cat.vmi ; V-I color already has the zeropoint offset!
       out.mauto814_err = 0.0
       out.mauto606_err = 0.0
       out.z = float(0.328);+randomn(seed,ngal)*0.067/2.0
       wsex, out, outfile='cl1358.cat'
       return
    endif
    
    ra = ['13:58:00','10:54:00','20:53:00']
    dec = ['+62:00:00','-03:21:00','-04:00:00']
    infile = path+['cl1358.cat','hst1054.cat','hst2053.cat']
    outfile = repstr(infile,'.cat','.kcorr.fits')
    for ii = 0, 2 do begin
       cat = rsex(infile[ii])
       euler, 15.0D*hms2dec(ra[ii]), hms2dec(dec[ii]), gl, gb, 1
       ebv = dust_getval(gl,gb,/interp)
       addebv = replicate({ebv: ebv},n_elements(cat))
       cat = struct_addtags(cat,addebv)
       out = best_kcorrect(cat)
;      wsex, out, outfile=outfile[ii]
       im_mwrfits, out, outfile[ii]
    endfor

stop
    
    
    cc = rsex('cl1358.original.cat')
    out = mrdfits('cl1358.kcorr.fits.gz',1)

    ml = 1.737*(out.ubv_absmag[1]-out.ubv_absmag[2])-0.942
    mass = ml + 0.4*(5.45-out.ubv_absmag[1])

;   niceprint, cc.massb, mass
;   niceprint, cc.bmv, out.ubv_absmag[1]-out.ubv_absmag[2]
;   niceprint, cc.absv, out.ubv_absmag[2]
    ss = im_stats(cc.absv-out.ubv_absmag[2],/ver)
    ss = im_stats(cc.bmv-(out.ubv_absmag[1]-out.ubv_absmag[2]),/ver)
    ss = im_stats(cc.massb-mass,/ver)

; test the catalog
    dm = 41.18 ; dmodulus(0.328)
    ebv = 0.014
    kl = [2.81588,1.93973]

    mv1 = (cc.imag+0.098) + 0.27*cc.vmi + 0.63 - dm
    mv2 = ((cc.imag+0.098)-ebv*kl[1]) + 0.27*(cc.vmi-ebv*(kl[0]-kl[1])) + 0.63 - dm
    print, cc.absv-mv1
    print, cc.absv-mv2

    
stop    
    
return
end
    
