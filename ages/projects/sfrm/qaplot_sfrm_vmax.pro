pro plot_sfrm_vmax
; jm10feb15ucsd - build some QAplots demonstrating that my Vmax values
;   are reasonable 

    sfrmpath = ages_path(/projects)+'sfrm/'
    paperpath = ages_path(/papers)+'sfrm/'

    parent = read_sfrm_sample()
    zbins = sfrm_zbins(nzbins)

    area = ages_survey_area()*!dtor^2 ; [sr]
    h100 = 0.7

    psfile = sfrmpath+'qaplots/sfrm_vmax.ps'
    im_plotconfig, 0, pos, psfile=psfile

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

return
end
    
