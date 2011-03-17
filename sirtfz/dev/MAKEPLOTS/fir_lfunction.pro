pro fir_lfunction
;+
; NAME:
;	FIR_LF
;
; PURPOSE:
;	Generate a plot of the 60 mu (far-IR) IRAS Bright Galaxy
;	Sample luminosity function from Soifer & Neugebauer 1991
;	(Guiderdoni et al 1998).
;
; PROCEDURES USED:
;	COLORTABLE1, SIRTF_DATAPATH(), READFAST, JMWINDOW, PLOTSYM,
;	TEXSTRINGS(), PLOTFAVES
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 September 15, U of A
;-

    plotfaves, charthick=1.8
    colortable1
    s = texstrings()

    path = sirtf_datapath() & path=path[0]
    readfast, path+'data/fir_lf.dat', lf, skip=4, ncols=3

    jmwindow, 0.4, yratio=0.5, windx=0, title='FIR Luminosity Function'
    plotsym, 0, 1, /fill
    plot, lf[0,*], lf[1,*], ps=8, color=7, xsty=3, ysty=3, $
      xrange=[6,15], yrange=[-9,1], title='IRAS 60 '+s.mu+' Luminosity Function', $
      xtitle='log '+s.nu+'L'+textoidl('_{\nu}')+' (L'+textoidl('_{'+s.sun+'}')+')', $
      ytitle='log '+s.cphi+' ('+s.nu+'L'+textoidl('_{\nu}')+') Mpc'+textoidl('^{-3}')+' mag'+textoidl('^{-1}')

    plotfaves, /restore
    
return
end
