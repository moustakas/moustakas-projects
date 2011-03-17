pro qaplot_mz_selection
; jm10may08ucsd - QAplots to help select our selection parameters

    mzpath = ages_path(/projects)+'mz/'

INCOMPLETE CODE!
    
    ppxf = read_ages(/ppxf)
    phot = read_ages(/photo)
    phot = phot[ppxf.ages_id]

    psfile = mzpath+'qaplots/mz_selection.ps'
    im_plotconfig, 0, pos, psfile=psfile

       test = where(phot.main_weight and (phot.z gt sample_zmin) and $
         (phot.z lt sample_zmax) and (phot.i_mag_auto gt 0.0) and $
         (phot.i_mag_auto lt 20.0))
       losnr = where(ppxf[test].continuum_snr lt 3.0,comp=hisnr)
       im_plothist, phot[test[hisnr]].i_mag_auto, bin=0.02, yr=[0,25], xr=[19,20.5]
       im_plothist, phot[test[losnr]].i_mag_auto, bin=0.02, /over, /fill
       djs_oplot, 19.8*[1,1], !y.crange, color='red'

       djs_plot, ppxf[test].continuum_snr, phot[test].i_mag_auto, ps=6, sym=0.2, $
         /xlog, xsty=3, ysty=3, yr=[18,20], xr=[0.1,50]
       djs_oplot, ppxf[these].continuum_snr, phot[these].i_mag_auto, $
         ps=6, sym=0.2, color='red'
    
    
    for ii = 0, nzbins-1 do begin
       these = where((parent.z gt zbins[ii].zlo) and $
         (parent.z lt zbins[ii].zup),nthese)
       sample = parent[these]

; compute Vmax       
       zmin = sample.zmin_noevol>zbins[ii].zlo
       zmax = sample.zmax_noevol<zbins[ii].zup
       vmax = (area/3.0)*(lf_comvol(zmax)-lf_comvol(zmin))/h100^3.0 ; h=0.7
       vol = (area/3.0)*((lf_comvol(zbins[ii].zup)-$
         lf_comvol(zbins[ii].zlo)))[0]/h100^3.0 ; h=0.7

       neg = where(vmax le 0,nneg)
       if (nneg ne 0) then stop
       
       hi = where(sample.k_chi2 gt 100,comp=lo)
       
       djs_plot, sample[lo].ugriz_absmag[2], vmax[lo]/vol[lo], psym=6, $
         xtitle='M_{0.1r}', ytitle='V/V_{max}', sym=0.4, xsty=1, ysty=3, $
         xrange=[-15,-25.5], yrange=[0,1.05]
       djs_oplot, sample[hi].ugriz_absmag[2], vmax[hi]/vol[hi], $
         psym=6, symsize=0.4, color='red'
       legend, 'z='+strtrim(string(zbins[ii].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[ii].zup,format='(F12.2)'),2), $
         /left, /top, box=0, margin=0
    endfor    
    im_plotconfig, /psclose, psfile=psfile, /gzip

; emission-line selection    
    
; diagnostic plots       
;      plot, ppxf.h_beta_ew[0], ppxf.h_beta[0]/ppxf.h_beta[1], $
;        ps=4, /xlog, /ylog, xr=[0.1,500], yr=[0.1,1E3]
;      djs_oplot, ewhbcut1*[1,1], 10^!y.crange, line=0, color='red'

       ww = where((ppxf.h_beta_ew[0] gt ewhbcut1) and $
         (ppxf.oii_3727_ew[1] le 0.0))
       data = ppxf[ww] & specfit = read_ages_gandalf_specfit(data)
       qaplot_ages_gandalf_specfit, data, specfit, psfile='junk.ps'

    
    
return
end
    
