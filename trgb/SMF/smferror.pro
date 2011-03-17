; determine the one-sigma error of the bootstrapped TRGB detections
; and generate a plot of the errors
pro smferror, infobase, lf, iter, trgblin, trgblinboot, trgblog, trgblogboot, linerr, logerr

;	colortable1

; ----------------------------------------------------------------------
; calculate the error using the cumulative distribution function
; ----------------------------------------------------------------------

	x1 = long(0.1827*iter) ; 18.27%
        x2 = long(0.50*iter)   ; 50%
	linsrt = sort(trgblinboot)
	logsrt = sort(trgblogboot)

        linerr = (trgblinboot[linsrt])[x2] - (trgblinboot[linsrt])[x1]
        logerr = (trgblogboot[logsrt])[x2] - (trgblogboot[logsrt])[x1]
        
;       window, 0, xs=450, ys=450
;       plothist, trgblinboot, bin=0.01, color=7, $ ; xrange=[trgblin-0.5,trgblin+0.5]
;         ytit='Number', xtit='I', ysty=3, title='Distribution of Results', line=0, xsty=3
;       plothist, trgblogboot, color=3, /overplot, bin=0.01, line=2

return
end

