pro ages_ppxf_specphot_tweak, clobber=clobber
; jm09nov13ucsd - take the output from AGES_PPXF_PRETWEAK and derive a
;   low-order correction to the spectrophotometry for each PASS

    version = ages_version(/ppxf_specfit)
    base_specfitpath = ages_path(/ppxf)
    specfitpath = ages_path(/spec1d)+'fluxed/pretweak/'+version+'/' 

    allfiles = file_search(specfitpath+'ppxf_???.fits.gz')
    nfile = n_elements(allfiles)
    
    nbin = 40
    out = replicate({pass: 0, nobj: 0, wave: fltarr(nbin)-999.0, $
      tweak: fltarr(nbin)-999.0},nfile)
    
    medbin = 150.0 ; [A]
    xrange = [3650,9250]
    yrange = [-1,1]
    flag = 1E-10

    outfile = base_specfitpath+'ppxf_specphot_tweak_'+version+'.fits'
    if file_test(outfile+'.gz',/reg) and (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+outfile+' exists; use /CLOBBER'
       return
    endif
    
    psfile = repstr(outfile,'.fits','.ps')
    im_plotconfig, 8, pos, psfile=psfile, charsize=1.5, $
      ymargin=[0.8,1.1]
;   for ii = 29, 30 do begin
    for ii = 0, nfile-1 do begin
       splog, format='("PASS ",I0,"/",I0)', ii+1, nfile

;      splog, 'Reading '+allfiles[ii]
       result = mrdfits(allfiles[ii],1)
       nobj = n_elements(result)

       index = where((result.chi2 lt 1.0) and (result.ebv lt 0.3) and $
         (result.zabs gt 0.05) and (result.zabs lt 0.5),nobj)

; all in the observed frame       
       good = where((result[index].wave gt 0.0))
       wave = (result[index].wave)[good]
       data = (result[index].flux)[good]
       model = (result[index].bestfit)[good]
       resid = 2.5*alog10((data/(model+(model eq 0.0))*(model ne 0.0))>flag)
       weight = (resid*0.0+1.0)*(resid gt 2.5*alog10(flag)) ; 0=bad
;      zero = where(weight eq 0.0) & if zero[0] ne -1 then stop

       med = im_medxbin(exp(wave),resid,medbin,weight=weight)
       nmedbins = n_elements(med)

; pack into a structure
       out[ii].pass = result[0].pass
       out[ii].nobj = nobj
       out[ii].wave[0:nmedbins-1] = med.xbin
       out[ii].tweak[0:nmedbins-1] = med.medy

; QAplot       
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xrange=xrange, yrange=yrange, xtitle='Observed Wavelength (\AA)', $
         ytitle='+2.5 log_{10} (Observed Spectrum / Best-Fit Model) (mag)', $
         title='Pass '+string(result[0].pass,format='(I0)')
       djs_oplot, !x.crange, [1,1], line=0, color='grey'
       djs_oplot, 8500*[1,1], !y.crange, line=5, color='grey'
       djs_oplot, med.xbin, med.medy, line=0, color='red'
       djs_oplot, med.xbin, med.quant75, line=5, color='red'
       djs_oplot, med.xbin, med.quant25, line=5, color='red'
       im_legend, 'N='+string(nobj,format='(I0)'), /right, /top, $
         box=0, charsize=1.6
       
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

; write out
    im_mwrfits, out, outfile, /clobber
    
return
end
    
