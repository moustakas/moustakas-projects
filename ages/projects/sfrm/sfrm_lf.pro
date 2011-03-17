function fit_luminosity_function, mr, weight, binsize=binsize, $
  binmr=binmr, phi=phi, errphi=phierr, parinfo=parinfo, $
  nofit=nofit, plotlimit=plotlimit, mrlimit=mrlimit
; fit the stellar mass function

; round the bin centers to two decimal points
;   histmin = fix((max(mr)-binsize/2.0)*100.0)/100.0+binsize/2.0

; build the LF for the plot, setting the weight to zero for objects
; below PLOTLIMIT

    phi = im_hist1d(mr,weight*(mr lt plotlimit),$
      binsize=binsize,obin=binmr,binedge=0,$
      h_err=phierr,histmin=histmin)
    phi = phi/binsize
    phierr = phierr/binsize
    
; fit the mr function; remake the mr function, setting the weight
; to zero for objects less mrive than MRLIMIT
    if keyword_set(nofit) then fit = 0.0 else begin
       fitphi1 = im_hist1d(mr,weight*(mr lt mrlimit),$;*(mr gt -22.5+5*alog10(0.7)),$
         binsize=binsize,obin=fitmr1,binedge=0,$
         h_err=fitphierr1,histmin=histmin)
       good = where(fitphi1 ne 0.0)
       fitphi = fitphi1[good]/binsize
       fitphierr = fitphierr1[good]/binsize
       fitmr = fitmr1[good]

       lf_fit_schechter, fitmr, fitphi, fitphierr, $
         fit, parinfo=parinfo

;      ploterror, binmr, phi, phierr, ps=10, xsty=3, ysty=3, $
;        yrange=[1E-5,1E-1], /ylog
;      oploterror, fitmr, fitphi, fitphierr, ps=10, $
;        color=djs_icolor('red'), errcolor=djs_icolor('red')
;      djs_oplot, fitmr, mf_schechter(10^fitmr,fit), color='yellow'
;      djs_oplot, fitmr, mf_schechter(10^fitmr,parinfo[0].value,$
;        parinfo[1].value,parinfo[2].value), color='green'
    endelse
return, fit
end

pro sfrm_lf, jackknife=jackknife, debug=debug
; jm10feb05ucsd - for comparison, measure the 0.1r-band luminosity
; function in AGES

    sfrmpath = ages_path(/projects)+'sfrm/'
    parent = read_sfrm_sample()

    binsize = 0.25
    mraxis1 = reverse(im_array(-25,-12,0.02))
    xrange = [-18.0,-23.0]
;   yrange = [-5,-1]
    yrange = [1E-5,2E-2]

; requisite parameters for computing Vmax in each redshift bin; the
; evolution parameters are following Eisenstein+09 
    if keyword_set(iauto) then itagname = 'i_mag_auto' else itagname = 'i_tot'
    itag = tag_indx(parent,itagname)

    vname = 'default.nolines'
    area = 7.9*!dtor^2 ; [sr]
    ibright = 15.0 ; [Vega]
    ifaint = 19.95 ; [Vega]
    ifilt = 'ndwfs_I.par'
    h100 = 1.0 ; 0.7
    
;   q0 = 0.0 & q1 = 0.0 & qz0 = 0.0 ; no evolution
    q0 = 1.6 & q1 = 0.0 & qz0 = 0.1 ; moderate luminosity evolution

; parameters of the fitting; constrain M* to be >10^8 M_sun and <10^11
; M_sun, and also constrain the faint-end slope
    parinfo = {value: 0.0D, limited: [0,0], $
      limits: [0.0D,0.0D], fixed: 0}
    parinfo = replicate(parinfo,3)

    parinfo[0].value = 5D-3
    parinfo[0].limited = 0
    parinfo[0].limits = [1D-4,0.0]

    parinfo[1].value = -20.0
    parinfo[1].limited = 0
    parinfo[1].limits = [-21.0,-19.0]

    parinfo[2].value = -1.05D ; alpha
    parinfo[2].fixed = 1

;   plot, mraxis1, lf_schechter(mraxis1,parinfo[0].value,parinfo[1].value,parinfo[2].value), $
;     /ylog, xr=[-12,-24]
    
; read the predefined redshift bins and the absolute magnitude limits 
    zbins = sfrm_zbins(nzbins)
    limits = mrdfits(sfrmpath+'sfrm_limits.fits.gz',1,/silent)

; initialize the output data structures
    result = {mrlimit: 0.0, plotlimit: 0.0, $
      mstar: 0.0, phistar: 0.0, alpha: 0.0,$
      mstar_err: 0.0, phistar_err: 0.0, alpha_err: 0.0, $
      mstar_jackknife_err: 0.0, phistar_jackknife_err: 0.0, $
      alpha_jackknife_err: 0.0}
    result = replicate(result,nzbins)
    result.mrlimit = limits.mrlim
    result.plotlimit = result.mrlimit+0.5 
    result = struct_addtags(zbins,result)

    data = {ngal: 0L, nbins: 0, mr: fltarr(50)-999.0, $
      mrerr: fltarr(50)-999.0, phi: fltarr(50)-999.0, $
      phierr: fltarr(50)-999.0}
    data = replicate(data,nzbins)
    data = struct_addtags(zbins,data)
    
; now measure the mr function in each redshift bin
    for ii = 0, nzbins-1 do begin
       these = where((parent.z gt zbins[ii].zlo) and $
         (parent.z lt zbins[ii].zup),nthese)
       sample = parent[these]

; compute Vmax       
       zmin = sample.zmin>zbins[ii].zlo
       zmax = sample.zmax<zbins[ii].zup
       vmax = (area/3.0)*(lf_comvol(zmax)-lf_comvol(zmin))/h100^3.0 ; h=0.7
       vol = (area/3.0)*((lf_comvol(zbins[ii].zup)-$
         lf_comvol(zbins[ii].zlo)))[0]/h100^3.0 ; h=0.7
       weight = sample.final_weight/vmax

;; compute Vmax and the final weight; use the targetting magnitude 
;       splog, 'Computing Vmax ---'
;       vvmax = im_calc_vmax(sample.i_obs,sample.coeffs,$
;         ifilt,area,ibright,ifaint,actual_z=sample.z,$
;         zbins[ii].zlo,zbins[ii].zup,q0=q0,q1=q1,qz0=qz0,$
;         vname=vname,/silent)
       
;       vv = mrdfits(ages_path(/analysis)+'catalog.Vmax.v3.fits.gz',1,rows=sample.ages_id)
;       res = im_zminzmax(sample.z,sample.i_obs,sample.coeffs,$
;         filter=ifilt,bright=ibright,faint=ifaint,vname=vname)
;       niceprint, res.zmin, vvmax.zmin, res.zmax, vvmax.zmax
;       plot, res.zmax, vvmax.zmax, ps=4, xr=[0.04,0.16], yr=[0.04,0.16]
;       djs_oplot, res.zmax, vv.zcut2, ps=4, color='red'
;       weight = sample.final_weight/vvmax.vmax

       kk = mrdfits(ages_path(/analysis)+'catalog.kcorr.v3.fits.gz',1,rows=sample.ages_id)
       mr = kk.m_r01
;      mr = sample.ugriz_absmag[2] ; absolute magnitude
stop
;; test
;       kk = mrdfits(ages_path(/analysis)+'catalog.kcorr.v3.fits.gz',1,rows=sample.ages_id)
;       vv = mrdfits(ages_path(/analysis)+'catalog.Vmax.v3.fits.gz',1,rows=sample.ages_id)
;
;       plot, vvmax.zmax, vv.zcut2<zbins[ii].zup[0], ps=4
       
; fit the luminosity function
       splog, 'Fitting the LF ---'
       fit = fit_luminosity_function(mr,weight,binsize=binsize,$
         parinfo=parinfo,plotlimit=result[ii].plotlimit,$
         mrlimit=result[ii].mrlimit,binmr=binmr,$
         phi=phi,errphi=phierr)
       nbins = n_elements(phi)
       data[ii].nbins = nbins
       data[ii].mr[0:nbins-1] = binmr
       data[ii].mrerr[0:nbins-1] = replicate(binsize/2.0,nbins)
       data[ii].phi[0:nbins-1] = phi
       data[ii].phierr[0:nbins-1] = phierr
       data[ii].ngal = total(mr gt result[ii].mrlimit)

       result[ii] = im_struct_assign(fit,result[ii],/nozero)
       
; use jack-knife to get the errors on the parameters 
       if keyword_set(jackknife) then begin
          allfield = sample.field
          field = allfield[uniq(allfield,sort(allfield))]
          nfield = n_elements(field)
          for ff = 0, nfield-1 do begin
             keep = where(allfield ne field[ff],ngal)
             jfit1 = fit_luminosity_function(mr[keep],weight[keep],$
               binsize=binsize,parinfo=parinfo,plotlimit=result[ii].plotlimit,$
               mrlimit=result[ii].mrlimit,binmr=binmr,$
               phi=phi,errphi=phierr)
             if (ff eq 0) then jfit = jfit1 else jfit = [jfit,jfit1]
          endfor
; for the jackknife error equation see
; http://www.physics.utah.edu/~detar/phycs6730/handouts/jackknife/jackknife
          result[ii].mstar_jackknife_err = sqrt((nfield-1)*total((jfit.mstar-fit.mstar)^2)/float(nfield))
          result[ii].phistar_jackknife_err = sqrt((nfield-1)*total((jfit.phistar-fit.phistar)^2)/float(nfield))
          result[ii].alpha_jackknife_err = sqrt((nfield-1)*total((jfit.alpha-fit.alpha)^2)/float(nfield))
       endif

; plot it!       
       if keyword_set(debug) then begin
          ploterror, binmr, phi, replicate(binsize/2.0,nbins), phierr, $
            yrange=yrange, xrange=xrange, xsty=1, ysty=1, /ylog, $
            xtitle=textoidl('log (M/M_{'+sunsymbol()+'})'), $
            ytitle=textoidl('\Phi'), psym=symcat(15), symsize=1.2
          legend, 'z='+strtrim(string(zbins[ii].zlo,format='(F12.2)'),2)+$
            '-'+strtrim(string(zbins[ii].zup,format='(F12.2)'),2), $
            /right, /top, box=0
          djs_oplot, mraxis1, lf_schechter(mraxis1,fit), color='red'
          djs_oplot, result[ii].mrlimit*[1,1], 10^!y.crange, color='cyan'
          cc = get_kbrd(1)
       endif
    endfor

; write out
    im_mwrfits, result, sfrmpath+'lf_fit_all.fits', /clobber
    im_mwrfits, data, sfrmpath+'lf_data_all.fits', /clobber
    
return
end
    
