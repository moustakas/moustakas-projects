pro ages_kcorrect_zptoffset, compute=compute, qaplot=qaplot

    iopath = ages_path(/analysis)
    niter = 10
    
; --------------------------------------------------
; compute the zeropoint offsets
    if keyword_set(compute) then begin
    
; read the AGES photometry
       maggies = read_ages_photometry(ivarmaggies=ivarmaggies,$
         index=index,filterlist=filterlist,redshift=redshift,$
         nozpoffset=0)
    
;      index = index[900:1000]
       maggies = maggies[*,index]
       ivarmaggies = ivarmaggies[*,index]
       redshift = redshift[index]
       ngal = n_elements(redshift)
       nfilt = n_elements(filterlist)

; toss out filters that are not being used: galex and irac ch2-4 
       keep = where((total(ivarmaggies gt 0,2) ne 0.0),nfilt)
       filterlist = filterlist[keep]
       maggies = maggies[keep,*]
       ivarmaggies = ivarmaggies[keep,*]
       
; demand photometry in all NFILT bandpasses
       good = where(total(maggies gt 0.0,1) eq nfilt,ngal)
       redshift = redshift[good]
       maggies = maggies[*,good]
       ivarmaggies = ivarmaggies[*,good]
       
; compute the zeropoint shifts relative to the reference band
       refindx = where(strmatch(filterlist,'*_I*',/fold))

; initialize the input/output arrays       
       out = {iter: 0, filter: filterlist, ngal: ngal, $
         zptoffset: fltarr(nfilt)}
       out = replicate(out,niter)
       out.iter = intarr(niter)
       in_maggies = maggies
       in_ivarmaggies = ivarmaggies
       
;; test!
;       test = mrdfits(iopath+'ages_zptoffsets.fits.gz',1)
;       zptoffset = total(test.zptoffset,2)
;       for ii = 0, nfilt-1 do begin
;          factor = 10.0^(-0.4*zptoffset[ii])
;          in_maggies[ii,*] = in_maggies[ii,*]*factor
;          in_ivarmaggies[ii,*] = in_ivarmaggies[ii,*]/factor^2.0
;       endfor
       
       t0 = systime(1)
       for iter = 0, niter-1 do begin
          delvarx, kcorr_qafile
          if (iter eq 0) then kcorr_qafile = iopath+'qa_iterfirst.ps'
          if (iter eq niter-1) then kcorr_qafile = iopath+'qa_iterlast.ps'
          kcorr = im_kcorrect(redshift,in_maggies,in_ivarmaggies,$
            filterlist,sdss_filterlist(),band_shift=0.1,/silent,$
            coeffs=coeffs,bestmaggies=bestmaggies,chi2=chi2,$
            vname='default.nolines',psfile=kcorr_qafile)
          
          bestmag = -2.5*alog10(bestmaggies)
          mag = maggies2mag(in_maggies,ivarmaggies=in_ivarmaggies,magerr=errmag)
          ivarmag = 1.0/(errmag+(errmag eq 0.0))^2.0*(errmag ne 0.0)
          
; compute the zeropoint tweak using MPFIT and convert to magnitude 
          for ii = 0, nfilt-1 do begin
             good = where((mag[ii,*] gt 0.0) and (chi2 lt 100.0))
;            im_plothist, bestmag[ii,good]-mag[ii,good], bin=0.01
             out[iter].zptoffset[ii] = median(bestmag[ii,good]-mag[ii,good])
             if (abs(out[iter].zptoffset[ii]) lt 0.005) then out[iter].zptoffset[ii] = 0.0 
          endfor
          out[iter].zptoffset = out[iter].zptoffset - out[iter].zptoffset[refindx[0]]
             
; print some info       
          splog, '##################################################'
          splog, 'Iteration '+string(iter+1,format='(I2.2)')
          for ii = 0, nfilt-1 do splog, filterlist[ii], $
            out[iter].zptoffset[ii], format='(A30,F12.4)'
          
          if (iter eq 0) then begin
             splog, format='("Time for the first iteration = '+$
               '",G0," minutes")', (systime(1)-t0)/60.0
             print
          endif

; apply the derived offsets
          zptoffset = total(out.zptoffset,2) ; cumulative zeropoint shift
          for ii = 0, nfilt-1 do begin
             factor = 10.0^(-0.4*zptoffset[ii])
             in_maggies[ii,*] = maggies[ii,*]*factor
             in_ivarmaggies[ii,*] = ivarmaggies[ii,*]/factor^2.0
          endfor
       endfor
       
; write out    
       outfile = iopath+'ages_zptoffsets.fits'
       im_mwrfits, out, outfile, /clobber
    endif 
    
; ---------------------------------------------------------------------------
; build a QAplot
    if keyword_set(qaplot) then begin
       psfile = iopath+'ages_zptoffsets.ps'
       im_plotconfig, 0, pos, psfile=psfile, charsize=1.6, $
         xmargin=[1.3,0.2], ymargin=[0.6,0.9]
       
       allcolor = ['blue','red','dark green','orange','purple',$
         'cyan','magenta','dark red','brown']

       splog, 'Reading '+iopath+'ages_zptoffsets.fits.gz'
       out = mrdfits(iopath+'ages_zptoffsets.fits.gz',1)
       allfilter = strtrim(out[0].filter,2)
       nfilt = n_elements(allfilter)
          
       color = allcolor[0:nfilt-1]
       
       djs_plot, [0], [0], /nodata, xsty=1, ysty=3, $
         xrange=[0.1,niter+0.9], yrange=[-0.31,0.31], $
         ytitle='Cumulative Zeropoint Shift (model minus data, mag)', $
         xtitle='Iteration Number', position=pos
       djs_oplot, !x.crange, [0,0], line=1, thick=4
       for ii = 0, nfilt-1 do djs_oplot, findgen(niter)+1, $
         total(out.zptoffset[ii],/cumulative), psym=-symcat(16), symsize=1.5, $
         color=color[ii]
       legend, repstr(allfilter,'.par',''), /left, /bottom, $
         box=0, charsize=1.2, textcolor=djs_icolor(color)
       legend, 'N='+string(out[0].ngal,format='(I0)'), /left, /top, $
         box=0, charsize=1.2
       im_plotconfig, psfile=psfile, /psclose, /gzip
       spawn, 'rsync -auv '+psfile+'.gz ~', /sh
    endif
       
return
end
    
