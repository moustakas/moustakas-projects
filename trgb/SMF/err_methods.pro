; determine the one-sigma error of the bootstrapped TRGB detections
; and generate a plot of the errors
pro err_methods, infobase, lf, iter, trgblin, trgblinboot, trgblog, trgblogboot, linerr, logerr

; ----------------------------------------------------------------------
; calculate the error using the cumulative distribution function
; ----------------------------------------------------------------------

	x1 = long(0.1827*iter) ; 18.27%
        x2 = long(0.50*iter)   ; 50%
	linsrt = sort(trgblinboot)
	logsrt = sort(trgblogboot)

        linerr = (trgblinboot[linsrt])[x2] - (trgblinboot[linsrt])[x1]
        logerr = (trgblogboot[logsrt])[x2] - (trgblogboot[logsrt])[x1]
        
; ----------------------------------------------------------------------
; calculate the errors by iterative sigma clipping
; ----------------------------------------------------------------------
        
;        djs_iterstat, trgblinboot, sigrej=3.0, maxiter=15, mean=trgbmean, median=trgbmed, $
;          sigma=trgblinerr, mask=tmask
;        goodtrgb = where(tmask eq 0B,gtcount)
;        print & print, 'A total of '+strn(iter-gtcount)+' linear TRGB detections were rejected at 3-sigma.'
;        print

;        djs_iterstat, trgblogboot, sigrej=3.0, maxiter=15, mean=trgbmean, median=trgbmed, $
;          sigma=trgblogerr, mask=tmask
;        goodtrgb = where(tmask eq 0B, gtcount)
;        print & print, 'A total of '+strn(iter-gtcount)+' logarithmic TRGB detections were rejected at 3-sigma.'
;        print

; ----------------------------------------------------------------------
; calculate the errors by the root-median-square version of standard deviation 
; ----------------------------------------------------------------------

;        trgblinerr = stdev2(trgblinboot,mean(trgblinboot))
;        trgblogerr = stdev2(trgblogboot,mean(trgblogboot))

;        trgblinerr = trgblinerr/0.67
;        trgblogerr = trgblogerr/0.67

; ----------------------------------------------------------------------
        
        window, 0, xs=450, ys=450
        plothist, trgblinboot, bin=0.01, color=7, $ ; xrange=[trgblin-0.5,trgblin+0.5]
          ytit='Number', xtit='I', ysty=3, title='Distribution of Results', line=0, xsty=3
        plothist, trgblogboot, color=3, /overplot, bin=0.01, line=2

        path = trgb_datapath()
        ps_open, path[3]+'errors/'+infobase.object+'_smf_error', /ps_fonts
        device, /times, /inches
        plothist, trgblinboot, bin=0.01, color=7, $ ; xrange=[trgblin-0.5,trgblin+0.5]
          ytit='Number', xtit='I', ysty=3, title=infobase.truename, line=0, xsty=3
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=0, color=2, thick=3
        plothist, trgblogboot, color=3, /overplot, bin=0.01, line=2
        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, color=5, thick=3
        legend, ['Linear','Log'], lines=[0,2], box=0, charsize=1.5, color=1
        ps_close
        
return
end

