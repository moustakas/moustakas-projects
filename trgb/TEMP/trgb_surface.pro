; jm00june25ucb

pro trgb_surface, objname, hst=hst

	colortable1

;       trgb_readata, objname, datapath, data, infobase, hst=hst

;       readfast, '/deepscr1/ioannis/trgb/ugc07577/ugc07577_starlist.dat', data, skip=2, ncols=5
        readfast, '/deep3/marc/trgb/data/22dec97/sextansb/SextansB_IVmags.dat', data, ncol=8, skip=2

        srt = sort(data[3,*])	; sort by magnitude
        imags = (data[3,srt] > 18.) < 30.
        xcenter = data[1,srt]
        ycenter = data[2,srt]

        array = fltarr([[xcenter,ycenter],imags]

        t1 = systime(/seconds)
        window, 0, xs=450, ys=450
        surf = min_curve_surf(imags,xcenter,ycenter)
        surface, surf
        print, systime(/seconds)-t1

        result = tri_surf(imags,xcenter,ycenter,xvalues=xcenter,yvalues=ycenter,missing=median(imags))
        
        stop
        
return
end
