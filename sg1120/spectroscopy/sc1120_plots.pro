pro sc1120_plots, sc1120dust, sc1120nodust, postscript=postscript, hiiregions=hiiregions, $
  literature=literature, local=local, sdssdust=sdssdust, sdssnodust=sdssnodust, $
  fullsample=fullsample, cleanpng=cleanpng, bptmodels=bptmodels, paper=paper, talk=talk, $
  _extra=extra
; jm05jan31uofa
; jm04dec02uofa

; sc1120_plots, sc1120dust, sc1120nodust, /bptmodels, /postscript, /literature, /local, /paper, /hiiregions, /talk
    
    if keyword_set(postscript) then postthick = 8.0 else postthick = 2.0

    htmlbase = 'sc1120'

    html_path = sc1120_path(/web)
    pspath = html_path+htmlbase+'/'
;   paperpath = sc1120_path(/papers)+'FIG_PAPER1/'

    snrcut = 1.0

    if (n_elements(sc1120dust) eq 0L) then sc1120dust = read_sc1120(linefitnodust=sc1120nodust)
    
    grids = read_kewley_grids(Z=Z,U=U)

    if keyword_set(bptmodels) then models = kewley_bpt_lines(/kauffmann,_extra=extra)
    if keyword_set(hiiregions) then hii = read_hii_regions()

    if keyword_set(local) then begin
       atlasdust = read_integrated(/hiionly,/snrcuts,/silent)
       nfgsdust = read_nfgs(/hiionly,/snrcuts,/silent)
    endif

; clean up old PS and PNG files

    if keyword_set(cleanpng) then begin
       splog, 'Deleting all PNG and PS files in '+pspath
       spawn, ['rm -f '+pspath+'*.png*'], /sh
       spawn, ['rm -f '+pspath+'*.ps'], /sh
    endif

    if (not keyword_set(postscript)) then im_window, 0, /square

; initialize plotting variables

    @'xyrange_sc1120'

; ------------------------------------------------------------
; EW(Ha) vs EW(OII)
; ------------------------------------------------------------

    psname = 'sc1120_ewha_vs_ewoii.ps'
    im_openclose, pspath+psname, postscript=postscript

    lineratio, sc1120dust, 'H_ALPHA_EW', '', 'OII_3727_EW', '', $
      x, xerr, y, yerr, index=indx, nindex=nindx, snrcut=snrcut
    
    xtitle = 'log EW(H\alpha)  [\AA]' 
    ytitle = 'log EW([O II] \lambda3727)  [\AA]'
    
    xrange = ewharange
    yrange = ewharange

    if keyword_set(local) then begin
 
       lineratio, atlasdust, 'H_ALPHA_EW', '', 'OII_3727_EW', '', xatlas, xerratlas, $
         yatlas, yerratlas, index=indxatlas, nindex=nindxatlas, snrcut=snrcut
 
       lineratio, nfgsdust, 'H_ALPHA_EW', '', 'OII_3727_EW', '', xnfgs, xerrnfgs, $
         ynfgs, yerrnfgs, index=indxnfgs, nindex=nindxnfgs, snrcut=snrcut
 
       xlocal = [xatlas,xnfgs]
       xerrlocal = [xerratlas,xerrnfgs]
 
       ylocal = [yatlas,ynfgs]
       yerrlocal = [yerratlas,yerrnfgs]
       
    endif
    
    sc1120_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, /top, $
      xlocal=xlocal, ylocal=ylocal, xerrlocal=xerrlocal, yerrlocal=yerrlocal
    
    im_openclose, postscript=postscript, /close        

; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) then $
      im_ps2html, htmlbase, html_path=html_path, cleanpng=0, _extra=extra

;; --------------------------------------------------    
;; SELECT PLOTS FOR THE PAPER HERE
;; --------------------------------------------------    
;
;    if keyword_set(paper) then begin
;
;       splog, 'Writing paper plots to '+paperpath+'.'
;       paperplots = []
;
;       for k = 0L, n_elements(paperplots)-1L do $
;         spawn, ['/bin/cp -f '+html_path+htmlbase+'/'+paperplots[k]+' '+paperpath], /sh
;
;    endif

    
stop    
    
return
end    
