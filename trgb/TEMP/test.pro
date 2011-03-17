pro test

	objname = 'SextansB'
        halo = 1
        trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core
        imags = data[3,*] - infobase.a_i   ; apply extinction correction
        color = data[7,*] - infobase.e_v_i ; color excess correction
        nstars = n_elements(data[0,*])

        binsize = 0.05
	trgb_dacosta, gcrgb     ; globular cluster rgb
        xy2image, color, imags, cmdim, binsize=binsize
        xy2image, gcrgb[1].color, gcrgb[1].imag, gcim, binsize=binsize

        imcompare, reverse(gcim,2), reverse(cmdim,2)
        
;        window, 2, xs=350, ys=350
;        display, cmdim
        

stop
return
end
