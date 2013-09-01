pro clash_spitzer_plots
; jm13aug28siena - build some QAplots

    isedfit_dir = getenv('IM_RESEARCH_DIR')+'/projects/clash/spitzer/'

    prefix = 'spitzer'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    ff = clash_filterlist(weff=weff,/useirac,short=short)
    nfilt = n_elements(ff)

    rr = read_isedfit(isedfit_paramfile,isedfit_dir=isedfit_dir)
    ngal = n_elements(rr)

    psfile = isedfit_dir+'qa_clash_spitzer.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=[0.1,5.0], yrange=[-4,4], /xlog, $
      xtitle='Rest Wavelength \lambda (\mu'+'m)', $
      ytitle='Residuals (Data minus Model, mag)'
    for ii = 0, ngal-1 do begin
       ww = where(rr[ii].maggies gt 0.0,nww)
       irac = where(strmatch(short[ww],'*ch[1-2]*'),comp=hst)
       lambda = weff[ww]/(1.0+rr[ii].z)/1D4
       resid = -2.5*alog10(rr[ii].maggies[ww]/rr[ii].bestmaggies[ww])
       djs_oplot, [lambda[hst]], [resid[hst]], psym=symcat(16), $
         color=im_color('dodger blue'), symsize=0.7
       djs_oplot, [lambda[irac]], [resid[irac]], psym=symcat(15), $
         color=im_color('tan'), symsize=0.7
       if ii eq 0L then begin
          bigresid = resid
          biglambda = lambda
       endif else begin
          bigresid = [bigresid,resid]
          biglambda = [biglambda,lambda]
       endelse
    endfor
    mm = im_medxbin(alog10(biglambda),bigresid,0.05,minpts=15,/ver)
    djs_oplot, 10D^mm.medx, mm.medy, line=0, thick=8
    djs_oplot, 10D^mm.medx, mm.quant75, line=5, thick=8
    djs_oplot, 10D^mm.medx, mm.quant25, line=5, thick=8

    im_legend, ['HST','IRAC'], /right, /top, box=0, $
      color=['dodger blue','tan'], psym=[16,15]
    im_plotconfig, psfile=psfile, /psclose, /pdf
    

return
end
    
    
