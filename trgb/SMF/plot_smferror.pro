pro plot_smferror, objname

	!x.ticklen=0.02
        p_thick = !p.thick & chars = !p.charsize
        charth = !p.charthick & xth = !x.thick & yth = !y.thick
        !p.thick = 2.0 & !p.charsize=1.9 & !p.charthick=2.0 & !x.thick=2 & !y.thick=2

        colortable1
        device, decomposed=0

        trgb_readata, objname, datapath, data, infobase

        path = trgb_datapath()
        smfile = path[4]+'SMF/'+objname+'_smf.txt'

        readcol, smfile, lin, log, format='F,F', /silent

        plothist, lin, bin=0.01, color=7, $ ; xrange=[trgblin-0.5,trgblin+0.5]
          ytit='Number', xtit='I', ysty=3, title=infobase.truename, line=0, xsty=3
;       oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=0, color=2, thick=3
        plothist, log, color=3, /overplot, bin=0.01, line=2
;       oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, color=5, thick=3
;       legend, ['Linear','Log'], lines=[0,2], box=0, charsize=1.5, color=1

stop

;        ps_open, path[3]+'SMF/'+objname+'_smf_error', /ps_fonts
;        device, /times, /inches
;        plothist, trgblinboot, bin=0.01, color=7, $ ; xrange=[trgblin-0.5,trgblin+0.5]
;          ytit='Number', xtit='I', ysty=3, title=infobase.truename, line=0, xsty=3
;        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=0, color=2, thick=3
;        plothist, trgblogboot, color=3, /overplot, bin=0.01, line=2
;        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, color=5, thick=3
;        legend, ['Linear','Log'], lines=[0,2], box=0, charsize=1.5, color=1
;        ps_close

        !x.ticklen=0
        !p.thick = p_thick & !p.charsize=chars & !p.charthick=charth
        !x.thick = xth & !y.thick=yth

return
end
