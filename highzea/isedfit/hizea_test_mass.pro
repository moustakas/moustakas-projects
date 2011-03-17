pro hizea_test_mass, dogrid=dogrid
; jm09feb02nyu 

    path = '~/home/research/projects/hizea/'
    
    lsun = 3.826D33
    light = 2.99792458D18 ; speed of light [A/s]

    filterlist = ['sdss_'+['u0','g0','r0','i0','z0'],$
      'spitzer_irac_'+['ch1','ch2','ch3','ch4']]+'.par'
    filtinfo = im_filterspecs(filterlist=filterlist,/verbose)
    nfilt = n_elements(filterlist)

    readcol, path+'SDSSJ231122.58-083923.7.flux', weff, flux, $
      ferr, filter, format='F,F,F,A', /silent, comment='#'
    zobj = 0.73

    flam2maggies = filtinfo.weff^2.0/light*10.0^(0.4*48.6)
    maggies = flux[2:10]*flam2maggies
    errmaggies = ferr[2:10]*flam2maggies
    ivarmaggies = 1.0/errmaggies^2.0

    mab = -2.5*alog10(maggies)
    errmab = 2.5*errmaggies/maggies/alog(10.0)
    
    taufiles = file_search('$BC03_SFHGRID_DIR/fits/salp_m62_tau_??.?Gyr.fits.gz')
;   taufiles = file_search('$BC03_SFHGRID_DIR/measure/salp_m62_tau_??.?Gyr.info.fits')
    ntau = n_elements(taufiles)

    if keyword_set(dogrid) then begin
       for itau = 0L, ntau-1L do begin

          splog, 'Reading '+taufiles[itau]
          fits = mrdfits(taufiles[itau],1,/silent)
          lambda = k_lambda_to_edges(fits.wave*(1+zobj))
          these = where(fits.age/1D9 lt getage(zobj),nage)
          vmatrix = lsun*fits.flux[*,these]/(1+zobj)/(4.0*!dpi*dluminosity(zobj,/cm)^2)

          if (itau eq 0L) then begin
             out = replicate({tau: 0.0, age: fits.age[these], $
               chi2: dblarr(nage), mass: dblarr(nage), $
               modelmaggies: fltarr(nfilt,nage)},ntau)
          endif
          out[itau].tau = fits.tau

          for iage = 0L, nage-1L do begin
             out[itau].modelmaggies[*,iage] = k_project_filters(lambda,$
               vmatrix[*,iage],filterlist=filterlist,/silent)
             out[itau].mass[iage] = total(ivarmaggies*maggies*out[itau].modelmaggies[*,iage])/$
               total(ivarmaggies*out[itau].modelmaggies[*,iage]^2.0)
             out[itau].chi2[iage] = total(ivarmaggies*(maggies-out[itau].mass[iage]*out[itau].modelmaggies[*,iage])^2.0)
;         print, mass[itau,iage], chi2[itau,iage]
          endfor
       endfor
       mwrfits, out, path+'hizea_chi2grid.fits', /create
    endif else out = mrdfits(path+'hizea_chi2grid.fits',1)

; minimize

    tau = out.tau
    age = out[0].age
    ntau = n_elements(tau)
    nage = n_elements(age)
    
    finalchi2 = min(out.chi2,chi2indx)
    indx = array_indices([nage,ntau],chi2indx,/dim)
    ageindx = indx[0] & tauindx = indx[1]
    mstar = out[tauindx].mass[ageindx]
    finalage = age[ageindx]
    finaltau = tau[tauindx]

;   niceprint, mstar*out[tauindx].modelmaggies[*,ageindx], maggies

; test

    nmodel = ntau*nage
    vmaggies = rebin(reform(maggies,nfilt,1),nfilt,nmodel)
    vivarmaggies = rebin(reform(ivarmaggies,nfilt,1),nfilt,nmodel)
    vmodelmaggies = reform(out.modelmaggies,nfilt,nmodel)

    vmass = total(reform((ivarmaggies*maggies),1,nfilt)#vmodelmaggies,1)/$
      total(reform(ivarmaggies,1,nfilt)#vmodelmaggies^2.0,1)
    vchi2 = total(vivarmaggies*(vmaggies-rebin(reform(vmass,1,nmodel),nfilt,nmodel)*$
      vmodelmaggies)^2.0,1)
    print, min(vchi2,jj) & print, vmass[jj]
    
stop    
    
; make a plot    

    allsed = mrdfits(taufiles[tauindx],1,/silent)
    sedwave = allsed.wave*(1+zobj)
    get_element, allsed.age, out[0].age, these

;   zfactor = 1.0
    zfactor = 1.0/((1+zobj)*(4.0*!dpi*dluminosity(zobj,/cm)^2))
    sedflam = mstar*lsun*allsed.flux[*,these[ageindx]]*zfactor ; [erg/s/cm2/A]
    
    sedfnu = sedflam*sedwave^2.0/light
    sedmab = -2.5*alog10(sedfnu)-48.6

    dfpsplot, path+'SDSSJ2311_mass.ps', /color, /landscape
    im_plotfaves, /post
    djs_plot, [0], [0], /nodata, xrange=[0.1,20.0], yrange=[21.5,18], $
      /xsty, /ysty, /xlog, ytitle='m_{AB}', xtitle='Observed Wavelength (\mu'+'m)', $
      charsize=1.8
    label = ['\chi^{2} = '+strtrim(string(chi2min,format='(F12.1)'),2),$
      'log (M_{*}/M_{\odot}) = '+string(alog10(mstar),format='(F5.2)'),$
      '\tau = '+string(out[tauindx].tau,format='(F4.1)')+' Gyr',$
      't = '+string(out[tauindx].age[ageindx]/1E6,format='(F5.1)')+' Myr',$
      'Z/Z_{\odot} = 1']
    im_legend, label, /left, /top, box=0, charsize=1.8
    djs_oplot, sedwave/1E4, sedmab, line=0, color='green'
    plotsym, 8, 2.0, fill=1
    oploterror, filtinfo.weff/1E4, mab, filtinfo.fwhm/1E4, errmab, $
      psym=8, errthick=4.0
    im_plotfaves
    dfpsclose

    
    
stop    

return
end
    
