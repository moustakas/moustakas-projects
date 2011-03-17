pro test_lf_code
; jm10feb06ucsd - try to reproduce Eisenstein's LFs using *his*
;   calculations of Vmax and rest-frame quantities

    spectweight1 = mrdfits(ages_path(/analysis)+'catalog.spectweight.fits.gz',1)
    keep = where(spectweight1.main_weight,ngal)
    spectweight = spectweight1[keep]
    kcorr = mrdfits(ages_path(/analysis)+'catalog.kcorr.v3.fits.gz',1,rows=keep)
    vall = mrdfits(ages_path(/analysis)+'catalog.Vmax.v3.fits.gz',1,rows=keep)

    zbins = sfrm_zbins(nzbins)
    for ii = 0, nzbins-1 do begin
       these = where((kcorr.z gt zbins[ii].zlo) and $
         (kcorr.z lt zbins[ii].zup),nthese)
       sample = kcorr[these]
       mr = sample.m_r01

       weight = spectweight[these].spec_weight*spectweight[these].target_weight*$
         spectweight[these].fiber_weight/vall[these].vmax2
       im_plothist, mr, weight=weight, xx, yy, h_err=yyerr, /noplot
       ploterror, xx, yy, yyerr, psym=8, xr=[-17,-23], /ylog, yr=[0.1,1E6], sym=3

stop

    endfor
stop    

return
end
    
