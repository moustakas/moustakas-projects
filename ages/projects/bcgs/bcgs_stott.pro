function get_jk, sfhfile, zform=zform
; get observed frame J and Ks-band photometry

    filt = 'twomass_'+['J','Ks']+'.par'
    v2ab = k_vega2ab(filterlist=filt,/kurucz,/silent)
    ksolar = k_solar_magnitudes(filterlist='twomass_Ks.par',/silent)
    nfilt = n_elements(filt)
    
    dz = 0.05 & minz = 0.01 & maxz = 5.0
    zz = findgen((maxz-minz)/dz+1)*dz+minz
    nz = n_elements(zz)

    sfhpath = getenv('ISEDFIT_SFHGRID_DIR')+'/basemodels/bc03/'
    sed = mrdfits(sfhpath+sfhfile,1)
    restwave = sed.wave

    age = getage(zz)-getage(zform) ; galaxy age [Gyr]
    goodages = where(age gt 0.0,nage)
    age = age[goodages]
    redshift = zz[goodages]

    vmatrix = interpolate(sed.flux,findex(sed.age/1D9,age)) ; [erg/s/cm^2/M_sun] at 10 pc
       
; project the desired filters onto all the SEDs vs redshift
    out = replicate({jk_ab: 0.0, jk_vega: 0.0, age: 0.0, $
      z: 0.0, kcorr: 0.0, ml_k: 0.0},nage)
    out.age = age
    out.z = redshift
    for iz = 0, nage-1 do begin
       wave = restwave*(1+redshift[iz])
       jk = reform(k_project_filters(k_lambda_to_edges(wave),$
         reform(vmatrix[*,iz]/(1+redshift[iz])),filterlist=filt,/silent))
       out[iz].jk_ab = reform(-2.5*alog10(jk[0]/jk[1])) ; J-Ks, observed
       krest = reform(k_project_filters(k_lambda_to_edges(restwave),$ ; rest-frame
         reform(vmatrix[*,iz]),filterlist='twomass_Ks.par',/silent))
       out[iz].kcorr = +2.5*alog10(krest/jk[1])
       out[iz].ml_k = 10^(-0.4*ksolar)/krest ; K-band M/L 
;      djs_plot, restwave, vmatrix[*,iz], xsty=3, ysty=3, xrange=[3E3,3E4], /xlog
;      splog, age[iz], redshift[iz], out[iz].ml_k
;      cc = get_kbrd(1)
    endfor
    out.jk_vega = out.jk_ab - (v2ab[0]-v2ab[1]) ; AB --> Vega

return, out
end

function get_mass, kmag_vega, z, info
; m_R = M_Q + DM(z) + K_QR(z) (all in AB)     

    kv2ab = k_vega2ab(filterlist='twomass_Ks.par',/kurucz,/silent)
    ksolar = k_solar_magnitudes(filterlist='twomass_Ks.par',/silent)

    nobj = n_elements(kmag_vega)
    out = {z: 0.0, kmag: 0.0, kabs: 0.0, mass: 0.0}
    out = replicate(out,nobj)
    out.z = z
    out.kmag = kmag_vega
    
    out.kabs = (kmag_vega + kv2ab) - dmodulus(z) - $ ; AB
      interpol(info.kcorr,info.z,z)
    out.mass = interpol(info.ml_k,info.z,z)*10^(-0.4*(out.kabs-ksolar))
    
return, out
end

pro bcgs_stott
; jm10jun28ucsd - reproduce the Stott+ results

    use the right filter curves!


    
    kv2ab = k_vega2ab(filterlist='twomass_Ks.par',/kurucz,/silent)

    path = ages_path(/projects)+'bcgs/'
    bootes = rsex(path+'bcgs_sample.sex')
    good = where(bootes.good)
    bootes = bootes[good]
    ised = mrdfits(path+'BwRIJHKsirac_bc03_chab_calzetti_sfhgrid02.fits.gz',1,rows=good)
;   phot = mrdfits(path+'bcgs_photometry.fits.gz',1,rows=good)

    stott08 = mrdfits(path+'08stott.fits.gz',1)
    stott10 = rsex('10stott.sex')
;   stott = stott[where(stott.z lt 0.4)]

; figure out the right formation redshift for the CSP model
;   zf = 8.2
;   age = range(0.0,getage(0.0)-getage(zf),100)
;   sfr = exp(-age/1.0)
;   print, im_integral(age,sfr,getage(zf),getage(5.0))/im_integral(age,sfr,0.0,getage(zf))
    
    zf2 = get_jk('chab_Z0.02_tau_00.0Gyr.fits.gz',zform=2.0)
    zf5 = get_jk('chab_Z0.02_tau_00.0Gyr.fits.gz',zform=5.0)
    csp = get_jk('chab_Z0.02_tau_01.0Gyr.fits.gz',zform=8.0)

; get absolute magnitudes and stellar masses, for each set of models
    mstott08_zf2 = get_mass(stott08.kmag,stott08.z,zf2)
    mstott08_zf5 = get_mass(stott08.kmag,stott08.z,zf5)
    mstott08_csp = get_mass(stott08.kmag,stott08.z,csp)

    mstott10_zf2 = get_mass(stott10.kmag,stott10.z,zf2)
    mstott10_zf5 = get_mass(stott10.kmag,stott10.z,zf5)
    mstott10_csp = get_mass(stott10.kmag,stott10.z,csp)

    mbootes_zf2 = get_mass(bootes.kmag,bootes.z,zf2)
    mbootes_zf5 = get_mass(bootes.kmag,bootes.z,zf5)
    mbootes_csp = get_mass(bootes.kmag,bootes.z,csp)

; make some plots    
    im_plotconfig, 0, pos, psfile=path+'stott.ps', $
      xmargin=[1.3,0.4], width=6.8, height=6.8
; J-K vs redshift
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xrange=[0.03,2.5], /xlog, yrange=[0.5,2.8], $
      xtitle='Redshift', ytitle='J-K_{s} (Vega)', $
      position=pos
    djs_oplot, zf2.z, zf2.jk_vega, line=5, thick=5
    djs_oplot, zf5.z, zf5.jk_vega, line=0, thick=5
    djs_oplot, csp.z, csp.jk_vega, line=3, thick=5, color='blue'
    djs_oplot, bootes.z, bootes.jmag-bootes.kmag, psym=6, color='orange', symsize=1.5
    djs_oplot, stott08.z, stott08.jmag-stott08.kmag, psym=7, color='red', symsize=1.5
    im_legend, ['SSP (z_{f}=2)','SSP (z_{f}=5)','\tau=1 Gyr (z_{f}=8)'], /left, $
      /top, box=0, line=[5,0,3], thick=5, pspacing=1.2, $
      color=['','','blue']
    im_legend, ['Stott+08','Bootes'], /left, /bottom, $
      box=0, psym=[7,6], color=['red','orange'], symsize=1.5

; mass vs mass
; # results in Stott+10 vs my analysis
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xrange=[10.6,12.4], yrange=[10.6,12.4], $
      xtitle='log_{10} (M/M_{'+sunsymbol()+'}) [Stott+10]', $
      ytitle='log_{10} (M/M_{'+sunsymbol()+'}) [Simple BC03 method]', $
      position=pos
    djs_oplot, !x.crange, !y.crange, line=0, thick=5
    djs_oplot, alog10(stott10.mass*1D12), alog10(mstott10_zf2.mass), $
      psym=symcat(16), color=fsc_color('forest green',101), symsize=1.5
    djs_oplot, alog10(stott10.mass*1D12), alog10(mstott10_zf5.mass), $
      psym=symcat(15), color=fsc_color('dodger blue',101), symsize=1.5
    djs_oplot, alog10(stott10.mass*1D12), alog10(mstott10_csp.mass), $
      psym=symcat(14), color=fsc_color('firebrick',101), symsize=1.5
    im_legend, ['SSP (z_{f}=2)','SSP (z_{f}=5)','\tau=1 Gyr (z_{f}=8)'], /left, $
      /top, box=0, psym=[16,15,14], color=['forest green','dodger blue','firebrick'], $
      symsize=1.5

; # Bootes sample: my analysis vs fancy SED-fitting
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xrange=[10.6,12.4], yrange=[10.6,12.4], $
      xtitle='log_{10} (M/M_{'+sunsymbol()+'}) [Bootes, SED-fitting]', $
      ytitle='log_{10} (M/M_{'+sunsymbol()+'}) [Bootes, Simple BC03 method]', $
      position=pos
    djs_oplot, !x.crange, !y.crange, line=0, thick=5
    djs_oplot, ised.mass50, alog10(mstott10_zf2.mass), $
      psym=symcat(16), color=fsc_color('forest green',101), symsize=1.5
    djs_oplot, ised.mass50, alog10(mstott10_zf5.mass), $
      psym=symcat(15), color=fsc_color('dodger blue',101), symsize=1.5
    djs_oplot, ised.mass50, alog10(mstott10_csp.mass), $
      psym=symcat(14), color=fsc_color('firebrick',101), symsize=1.5
    im_legend, ['SSP (z_{f}=2)','SSP (z_{f}=5)','\tau=1 Gyr (z_{f}=8)'], /left, $
      /top, box=0, psym=[16,15,14], color=['forest green','dodger blue','firebrick'], $
      symsize=1.5

; mass vs redshift for each sample and set of models
; ## SSP, zf=2
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xrange=[0.03,2.5], /xlog, yrange=[10.3,12.7], $
      xtitle='Redshift', ytitle='log_{10} (M/M_{'+sunsymbol()+'})', $
      position=pos
    im_legend, ['Stott+08','Stott+10','Bootes'], /left, /bottom, $
      box=0, psym=[7,16,6], color=['red','blue','orange'], symsize=1.5
    im_legend, 'SSP (z_{f}=2)', /left, /top, box=0
    djs_oplot, mstott08_zf2.z, alog10(mstott08_zf2.mass), psym=7, color='red', symsize=1.5
    djs_oplot, mstott10_zf2.z, alog10(mstott10_zf2.mass), psym=symcat(16), color='blue', symsize=1.5
    djs_oplot, mbootes_zf2.z, alog10(mbootes_zf2.mass), psym=6, color='orange', symsize=1.5

; ## SSP, zf=5
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xrange=[0.03,2.5], /xlog, yrange=[10.3,12.7], $
      xtitle='Redshift', ytitle='log_{10} (M/M_{'+sunsymbol()+'})', $
      position=pos
    im_legend, ['Stott+08','Stott+10','Bootes'], /left, /bottom, $
      box=0, psym=[7,16,6], color=['red','blue','orange'], symsize=1.5
    im_legend, 'SSP (z_{f}=5)', /left, /top, box=0
    djs_oplot, mstott08_zf5.z, alog10(mstott08_zf5.mass), psym=7, color='red', symsize=1.5
    djs_oplot, mstott10_zf5.z, alog10(mstott10_zf5.mass), psym=symcat(16), color='blue', symsize=1.5
    djs_oplot, mbootes_zf5.z, alog10(mbootes_zf5.mass), psym=6, color='orange', symsize=1.5

; ## CSP, zf=8
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xrange=[0.03,2.5], /xlog, yrange=[10.3,12.7], $
      xtitle='Redshift', ytitle='log_{10} (M/M_{'+sunsymbol()+'})', $
      position=pos
    im_legend, ['Stott+08','Stott+10','Bootes'], /left, /bottom, $
      box=0, psym=[7,16,6], color=['red','blue','orange'], symsize=1.5
    im_legend, '\tau=1 Gyr (z_{f}=8)', /left, /top, box=0
    djs_oplot, mstott08_csp.z, alog10(mstott08_csp.mass), psym=7, color='red', symsize=1.5
    djs_oplot, mstott10_csp.z, alog10(mstott10_csp.mass), psym=symcat(16), color='blue', symsize=1.5
    djs_oplot, mbootes_csp.z, alog10(mbootes_csp.mass), psym=6, color='orange', symsize=1.5

    im_plotconfig, /psclose

stop    
    
return
end
;
;    mstott08 = [$
;      get_mass(stott08.kmag,stott08.z,zf2),$
;      get_mass(stott08.kmag,stott08.z,zf5),$
;      get_mass(stott08.kmag,stott08.z,csp)]
;    mstott10 = [$
;      get_mass(stott10.kmag,stott10.z,zf2),$
;      get_mass(stott10.kmag,stott10.z,zf5),$
;      get_mass(stott10.kmag,stott10.z,csp)]
;    mbootes = [$
;      get_mass(bootes.kmag,bootes.z,zf2),$
;      get_mass(bootes.kmag,bootes.z,zf5),$
;      get_mass(bootes.kmag,bootes.z,csp)]
