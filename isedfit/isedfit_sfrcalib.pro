pro isedfit_sfrcalib
; jm13aug11siena - look at the SFR calibrations based on different
; models 
    
    ssppath = getenv('ISEDFIT_SSP_DIR')+'/'

    imf = 'salp'
    sps = ['fsps_v2.4_miles','bc03_stelib','pegase','maraston05']
    label = ['FSPS (v2.4)','BC03','PEGASE-HR',$
      'Maraston+05']
    color = ['black','firebrick','dodger blue','forest green']
    psym = [16,15,14,17]
    line = [0,5,3,2]
    nsps = n_elements(sps)

    zsun = 0.019
    light = im_light(/ang)
    dist = 10.0*3.085678D18 ; fiducial distance [10 pc in cm]
    
    time = range(0.01D,0.1D,10,/log) ; output time vector [Gyr]
;   time = range(0.01D,0.1D,30) ; output time vector [Gyr]
    ntime = n_elements(time)

    psfile = getenv('IM_PROJECTS_DIR')+'/isedfit/qa_sfrcalib.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, $
      xmargin=[1.5,0.2], width=6.8

; Q(H^0) calibration    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=[-0.05,2.8], yrange=[52.5,53.4], $
      xtitle='Stellar Metallicity Z/Z_{'+sunsymbol()+'}', $
      ytitle='log [Q(H^{0}) / SFR] (s^{-1} / M_{'+sunsymbol()+'} yr^{-1})'
    for ss = 0, nsps-1 do begin
       info = mrdfits(ssppath+'info_'+sps[ss]+'_'+imf+'.fits.gz',1,/silent)
       nZ = n_elements(info.Zmetal)
       sfrcalib = fltarr(ntime,nZ)
       for iZ = 0, nZ-1 do begin
          ssp = mrdfits(ssppath+sps[ss]+'/'+strtrim(info.sspfile[iZ],2),1,/silent)
          flux = isedfit_convolve_sfh(ssp.flux,age=ssp.age,tau=100D,$
            nlyc=ssp.nlyc,time=time,sfh=sfh,cspnlyc=nlyc)
          sfrcalib[*,iZ] = nlyc-alog10(sfh)
       endfor
       djs_oplot, info.Zmetal/zsun, sfrcalib[ntime-1,*], $
         psym=-symcat(psym[ss]), color=im_color(color[ss]), $
         line=line[ss], symsize=2.0, thick=6
    endfor
    k98 = -alog10(1.08D-53)
    plots, 1.0, k98, psym=symcat(6,thick=6), /data, symsize=1.5
    plots, 1.0, k98, psym=symcat(7,thick=6), /data, symsize=1.5
    xyouts, 1.2, k98+0.05, 'Kennicutt (1998)', charsize=1.4, $
      align=0.5, /data
    im_legend, label, /right, /top, box=0, color=color, $
      line=line, pspacing=1.9, margin=0, psym=-psym, thick=6, $
      charsize=1.6
    
; L_nu calibration    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=[-0.05,2.8], yrange=[27.5999,28.0], $
      xtitle='Stellar Metallicity Z/Z_{'+sunsymbol()+'}', $
      ytitle='log [L_{\nu}(0.15-0.28\mu'+'m) / SFR] (erg s^{-1} Hz^{-1} / M_{'+$
      sunsymbol()+'} yr^{-1})'
    for ss = 0, nsps-1 do begin
       info = mrdfits(ssppath+'info_'+sps[ss]+'_'+imf+'.fits.gz',1,/silent)
       nZ = n_elements(info.Zmetal)
       sfrcalib = fltarr(ntime,nZ)
       for iZ = 0, nZ-1 do begin
          ssp = mrdfits(ssppath+sps[ss]+'/'+strtrim(info.sspfile[iZ],2),1,/silent)
          flux = isedfit_convolve_sfh(ssp.flux,age=ssp.age,tau=100D,time=time,sfh=sfh)
          fnu = flux*0.0
          for tt = 0, ntime-1 do fnu[*,tt] = sfh[tt]/(flux[*,tt]*$
            ssp.wave^2/light*4.0*!dpi*dist^2)
          ww = where(ssp.wave gt 1500.0 and ssp.wave lt 2800.0)
          for tt = 0, ntime-1 do sfrcalib[tt,iZ] = -alog10(djs_mean(fnu[ww,tt]))
       endfor
       djs_oplot, info.Zmetal/zsun, sfrcalib[ntime-1,*], $
         psym=-symcat(psym[ss]), color=im_color(color[ss]), $
         line=line[ss], symsize=2.0, thick=6
    endfor
    k98 = -alog10(1.4D-28)
    plots, 1.0, k98, psym=symcat(6,thick=6), /data, symsize=1.5
    plots, 1.0, k98, psym=symcat(7,thick=6), /data, symsize=1.5
    xyouts, 1.2, k98+0.01, 'Kennicutt (1998)', charsize=1.4, $
      align=0.5, /data
    im_legend, label, /right, /top, box=0, color=color, $
      line=line, pspacing=1.9, margin=0, psym=-psym, thick=6, $
      charsize=1.6
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    

return
end
