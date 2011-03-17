; jm00apr17ucb

; determine the error budget in detecting the trgb using bootstrap resampling

pro trgb_err, objname

	on_error, 2 	; return to user
        npar = N_params()

        if npar ne 1 then begin
            print, 'Syntax: trgb_err, object_name'
            return
        endif
        
	spawn, ['cls']
        spawn, ['pwd'], datapath

        read_data, objname+'_IVmags.dat', data, ncol=8	; read the data

	nstars = n_elements(data[*,0])        
        print & print, 'There are '+strn(nstars)+' in this starlist.' & print

        imags = data[*,3]	; I magnitudes
        ierr = data[*,4]	; I magnitude errors

; determine the nominal trgb detection (with the original starlist)

        trgb_lfunction, objname, data, trgb_nom

; now create new starlists to find the rms scatter in the edge detection

        iter = 100L	; number of iterations
        trgb = fltarr(iter)

        for k = 0L, iter-1L do begin
            bootindx = floor(randomu(seed,nstars)*nstars) ; resample the starlist uniformly
            
; perturb the magnitude values to take into account photometric measurement errors

            i_perturbed = imags[bootindx] + ierr[bootindx]*randomn(nseed,nstars)
            bootdata = [ [data[bootindx,0:2]],[i_perturbed],[data[bootindx,4:7]] ]

            trgb_lfunction, objname, bootdata, edge
            trgb[k] = edge[0]
        endfor

; plot the residuals

        trgb_resid = trgb-trgb_nom[0]

        plotsym, 0, 1, /fill
        window, 0, xs=500, ys=500
        plot, trgb_resid, ps=8, xsty=3, ysty=3, $
          xtit='Number of Trials', ytit='True_TRGB-Random_TRGB', $
          tit='Residual Plot of TRGB Detection', yminor=3
        oplot, [!x.crange[0],!x.crange[1]], [0,0], line=2

        print
        print, 'Mean   : ', mean(trgb_resid)
        print, 'Median : ', median(trgb_resid)
        print, 'StDev  : ', stdev(trgb_resid)
        print

return
end
