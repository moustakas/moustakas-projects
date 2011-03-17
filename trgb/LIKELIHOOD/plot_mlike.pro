; jm00oct19uofa
; generate a plot of the maximum likelihood results
pro plot_mlike, lnlike, lf, error_data, infobase, trgb, trgbwidth, beta, nbeta, alpha=alpha

	!x.ticklen=0.02
        p_thick = !p.thick & chars = !p.charsize
        charth = !p.charthick & xth = !x.thick & yth = !y.thick
        !p.thick = 2.0 & !p.charsize=1.8 & !p.charthick=2.0 & !x.thick=2 & !y.thick=2

        colortable1
        device, decompose=0

	temp = max(lnlike,maxindx)

        x1 = maxindx mod lf.nhist
        z1 = maxindx/(lf.nhist*lf.nhist)
        y1 = maxindx/(lf.nhist*nbeta)
        
; regenerate the best fit
        
        trgb_lfunction, lf.mags, lf.merr, lfplot, binsize=0.02, $
          minmag=lf.minmag, maxmag=lf.maxmag
        
        gbest = trgb_lfmodel(error_data,lfplot.mag-trgb,trgbwidth,beta,$
                             alpha=alpha,/nosmooth)
        gnorm = total(lfplot.hist)/total(gbest)

        window, 0, xs=500, ys=500
        plot, lfplot.mag, alog10(float(lfplot.hist)>1.), color=7, xtit='I', ytit='log N', $
          xsty=3, ysty=3, title=infobase.truename, ps=10
        oplot, lfplot.mag, alog10(gnorm*gbest), line=0, color=5, thick=3.0
        oplot, [trgb,trgb], [!y.crange[0],!y.crange[1]], line=2, color=2, thick=2.5
;       legend, ['TRGB '+strn(trgb,format='(F6.2)')+' mag', $
;                'Strength '+strn(10^(-trgbwidth),format='(F6.2)'), $
;                'Beta  '+strn(beta,format='(F6.2)')], $
;         charsize = 1.2, textcolor=[1,1,1], box=0, /clear
        
;        path = trgb_datapath()
;        ps_open, path[3]+infobase.object+'_lhood', /ps_fonts
;        device, /times, /inches
;        plot, [min(lfplot.mag),max(lfplot.mag)], [0,max(alog10(float(lfplot.hist)>1.))*1.1], /nodata, $
;          line=0, yminor=3, xminor=3, xsty=3, ysty=3, $
;          yrange=[0,max(alog10(float(lfplot.hist)>1.))*1.1], xtickname=replicate(' ',10), $
;          color=1, title=infobase.truename, xtit='I', ytit='log N(m)'
;        oplot, lfplot.mag, alog10(float(lfplot.hist)>1.), color=7
;        oplot, lfplot.mag, alog10(gnorm*gbest), line=0, color=5, thick=1.5
;        oplot, [trgb,trgb], [!y.crange[0],!y.crange[1]], line=2, color=2, thick=2.5
;        ps_close

        !x.ticklen=0
        !p.thick = p_thick & !p.charsize=chars & !p.charthick=charth
        !x.thick = xth & !y.thick=yth

return
end

