function do_vmax_mffit, mfdata, parinfo=parinfo, double_parinfo=double_parinfo, $
  triple_parinfo=triple_parinfo, modified_parinfo=modified_parinfo, $
  plus_parinfo=plus_parinfo, mffit_double=mffit_double, mffit_triple=mffit_triple, $
  mffit_modified=mffit_modified, mffit_plus=mffit_plus, debug=debug, quiet=quiet, $
  title=title, uniform_weights=uniform_weights, minfitmass=minfitmass
; fit a Schechter function to the 1/Vmax weighted histogram; get the
; errors by Monte Carlo (not sure why, but the formal MPFIT errors are
; way too small, and in the SDSS sample the weighted fit is not
; particularly good) 

    maxfitmass = 12.5
    maxis1 = range(8.0,12.5,50)
    maxis2 = range(8.0,11.5,50)
    if (n_elements(minfitmass) eq 0) then minfitmass = 7.0

    gdlim = where(mfdata.limit eq 1)
    gdfit = where(mfdata.limit eq 1 and mfdata.meanmass lt maxfitmass and $
      mfdata.meanmass gt minfitmass,ngdfit)
    if keyword_set(uniform_weights) then begin
;      splog, 'Hack! - Uniform weights and no errors!'
       mffit = mf_fit_schechter(mfdata.meanmass[gdfit],mfdata.phi[gdfit],$
         mfdata.phi[gdfit]*0.05,parinfo=parinfo,quiet=1)
    endif else begin
       mffit = mf_fit_schechter(mfdata.meanmass[gdfit],mfdata.phi[gdfit],$
         mfdata.phierr[gdfit],parinfo=parinfo,quiet=1)
    endelse

; scale the errors
    mffit.phistar_err *= sqrt(mffit.chi2_dof)
    mffit.logmstar_err *= sqrt(mffit.chi2_dof)
    mffit.alpha_err *= sqrt(mffit.chi2_dof)

    mffit = struct_addtags(mffit,{unweighted: keyword_set(uniform_weights), $
      phistar_cv_err: 0.0, logmstar_cv_err: 0.0, alpha_cv_err: 0.0})

; render a QAplot
    if keyword_set(debug) then begin
;      niceprint, mfdata.meanmass[gdlim], mfdata.phi[gdlim], mfdata.number[gdlim]
       nofit = where(mfdata.number gt 0 and (mfdata.limit eq 0 or $
         mfdata.meanmass gt maxfitmass or mfdata.meanmass lt minfitmass))
       djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=[8.0,13.0], $
         yrange=[-7.0,-1], title=title
       djs_oplot, minfitmass*[1,1], !y.crange, line=0
       djs_oploterr, mfdata.meanmass[gdfit], alog10(mfdata.phi[gdfit]), $
         yerr=mfdata.phierr[gdfit]/mfdata.phi[gdfit]/alog(10), $
         psym=symcat(16), symsize=1.5
       djs_oploterr, mfdata.meanmass[nofit], alog10(mfdata.phi[nofit]), $
         yerr=mfdata.phierr[nofit]/mfdata.phi[nofit]/alog(10), psym=symcat(9), $
         symsize=1.5, color=im_color('yellow')
       djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit)), $
         color=im_color('orange'), thick=3
       cc = get_kbrd(1)
    endif

return, mffit
end

function redbaryons_binsize, bin=bin, minmass=minmass, maxmass=maxmass
    if n_elements(bin) eq 0 then binsize = 0.1 else binsize = bin
    minmass = 10.0-binsize/2 ; shift by a half-bin
    maxmass = 13.0
return, binsize
end

function mge_1d_phimodel, sol, linearmass=linearmass
    sigma = sol[1,*]
    surf = sol[0,*]/(sqrt(2*!pi)*sigma)

    if n_elements(linearmass) eq 0 then linearmass = $
      range(min(sigma)/2.0,max(sigma)*2.0,100,/log)
    phimodel = linearmass*0.0
    
    for jj = 0, n_elements(linearmass)-1 do phimodel[jj] = $
      total(surf*exp(-0.5*(linearmass[jj]/sigma)^2))
return, phimodel
end

pro build_redbaryons_csmf, debug=debug, clobber=clobber
; jm13aug27siena - build the CSMF for satellites and centrals 

    common com_redmapper, centrals, satellites
    
    datapath = redmapper_path(/redbaryons)

; binning parameters
    binsize = redbaryons_binsize(bin=binsize,$
      minmass=minmass,maxmass=maxmass)
    sat_masslimit = 10.0 ; 10.0 ; 10.5
    cen_masslimit = 10.9 ; 10.0 ; 11.0
    
; read the catalogs
    if n_elements(centrals) eq 0L then read_redmapper, $
      centrals=centrals, satellites=satellites

    zbins = redbaryons_zbins(nzbins)
    for iz = 0, nzbins-1 do begin
       title = 'z='+string(zbins[iz].zlo,format='(F4.2)')+'-'+$
         string(zbins[iz].zup,format='(F4.2)')

       isat = where(satellites.z ge zbins[iz].zlo and satellites.z lt zbins[iz].zup,nsat)
       icen = where(centrals.z ge zbins[iz].zlo and centrals.z lt zbins[iz].zup,ncen)

       satweight = replicate(1.0/float(nsat),nsat)
       cenweight = replicate(1.0/float(ncen),ncen)
;      satweight = replicate(1.0/zbins[iz].vol,nsat)*0.0+1.0
;      cenweight = replicate(1.0/zbins[iz].vol,ncen)*0.0+1.0
       sat_mfdata1 = im_mf_vmax(satellites[isat].mstar_avg,satweight,$
         binsize=binsize,minmass=minmass,maxmass=maxmass,masslimit=sat_masslimit)
       cen_mfdata1 = im_mf_vmax(centrals[icen].mstar_avg,cenweight,$
         binsize=binsize,minmass=minmass,maxmass=maxmass,masslimit=cen_masslimit)

; fit!
       sat_good = where(sat_mfdata1.number ge 1 and $
         sat_mfdata1.limit eq 1,comp=sat_ignore)
       cen_good = where(cen_mfdata1.number ge 1 and $
         cen_mfdata1.limit eq 1,comp=cen_ignore)
;      mge_fit_1d, 10D^sat_mfdata1.mass[sat_good], sat_mfdata1.phi[sat_good], $
;        ngauss=10, sol=sat_sol, /quiet, outer_slope=-5.0, /negative
;      mge_fit_1d, 10D^cen_mfdata1.mass[cen_good], cen_mfdata1.phi[cen_good], $
;        ngauss=10, sol=cen_sol, /quiet, outer_slope=-5.0, /negative
       
;      mffit1 = do_vmax_mffit(sat_mfdata1,parinfo=parinfo,debug=debug,$
;        uniform_weights=1,title=title)

; integrate numerically       
       sat_rho1 = integrate_mf(sat_mfdata1)
       cen_rho1 = integrate_mf(cen_mfdata1)

; debugging plot
       if keyword_set(debug) then begin
          djs_plot, [0], [0], /nodata, /ylog, yr=[1E-10,1E-3], title=title, $
            xsty=1, ysty=1, xrange=[9,12.5]

          oploterror, sat_mfdata1.mass[sat_good], sat_mfdata1.phi[sat_good], $
            sat_mfdata1.phierr[sat_good], psym=symcat(16), yr=[1E-10,1E-3], $
            color=im_color('forest green'), errcolor=im_color('forest green')
          oploterror, sat_mfdata1.mass[sat_ignore], sat_mfdata1.phi[sat_ignore], $
            sat_mfdata1.phierr[sat_ignore], psym=symcat(9), yr=[1E-10,1E-3], $
            color=im_color('yellow'), errcolor=im_color('yellow')
;         djs_oplot, sat_mfdata1.mass[sat_good], mge_1d_phimodel(sat_sol,$
;           linearmass=10D^sat_mfdata1.mass[sat_good])
            
          oploterror, cen_mfdata1.mass[cen_good], cen_mfdata1.phi[cen_good], $
            cen_mfdata1.phierr[cen_good], psym=symcat(15), yr=[1E-10,1E-3], $
            color=im_color('firebrick'), errcolor=im_color('firebrick')
          oploterror, cen_mfdata1.mass[cen_ignore], cen_mfdata1.phi[cen_ignore], $
            cen_mfdata1.phierr[cen_ignore], psym=symcat(6), yr=[1E-10,1E-3], $
            color=im_color('tomato'), errcolor=im_color('tomato')
;         djs_oplot, cen_mfdata1.mass[cen_good], mge_1d_phimodel(cen_sol,$
;           linearmass=10D^cen_mfdata1.mass[cen_good])
          cc = get_kbrd(1)
       endif
       
; pack it in!       
       if iz eq 0 then begin
          sat_mfdata = sat_mfdata1
          cen_mfdata = cen_mfdata1
          sat_rho = sat_rho1
          cen_rho = cen_rho1
       endif else begin
          sat_mfdata = [sat_mfdata,sat_mfdata1]
          cen_mfdata = [cen_mfdata,cen_mfdata1]
          sat_rho = [sat_rho,sat_rho1]
          cen_rho = [cen_rho,cen_rho1]
       endelse
    endfor

; write out
    outfile = datapath+'csmf_sat.fits'
    im_mwrfits, sat_mfdata, outfile, clobber=clobber, /nogzip
    im_mwrfits, sat_rho, outfile, /append, /gzip
    
    outfile = datapath+'csmf_cen.fits'
    im_mwrfits, cen_mfdata, outfile, clobber=clobber, /nogzip
    im_mwrfits, cen_rho, outfile, /append, /gzip
    
return
end
    
