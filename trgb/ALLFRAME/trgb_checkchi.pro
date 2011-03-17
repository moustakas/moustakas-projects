pro trgb_checkchi, objname, ncol, mini=mini
;+
; NAME:
;	TRGB_CHECKCHI
;
; PURPOSE:
;	Check the distribution of CHI values from the ALLFRAME
;	photometry. 
;
; INPUTS:
;	objname : galaxy name
;	ncol	: number of columns in the .raw photometry file
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; COMMON BLOCKS:
;
; PROCEDURES USED:
;	READFAST
;
; MODIFICATION HISTORY:
;	Bryan Mendez, 2000 April, UCB
;	generalized & documented, jm00may26ucb	
;-

	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]

        if keyword_set(mini) then filename = datapath+'/'+objname+'_mini.raw' else $
          filename = datapath+'/'+objname+'.raw'

	readfast, filename, data, skip=3, ncols=ncol

        chi = data[ncol-2,*]
        imag = data[3,*]

        window, 0, xs=450, ys=450
        plothist, chi, bin=0.1, xr=[0,2], xtit = 'Chi', thick=2, $
          xthick=2, ythick=2, charthick=2, charsize=1.5, ytit='Number', $
          tit=strupcase(objname), xsty=3, ysty=3

        window, 2, xs=450, ys=450
        plot, imag, chi, ps=3, xr=[5,20], yr=[0,10], $
          xtitle = 'I magnitude', ytitle = 'Chi', thick=2, $
          xthick=2, ythick=2, charthick=2, charsize=1.5, $
          tit=strupcase(objname), xsty=3, ysty=3

        print, 'Mean Chi   = ', mean(chi)
        print, 'Median Chi = ', median(chi)
        print, 'StDev Chi  = ', stdev(chi)

return
end
