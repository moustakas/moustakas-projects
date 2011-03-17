pro sfroii, atlas1, postscript=postscript

; read the data       

    if (n_elements(atlas1) eq 0L) then atlas1 = read_integrated()

    lsun = 3.826D33             ; [erg/s]
    haconst = 7.9D-42           ; K98 conversion L(IR) --> SFR(IR)
    haoii = 1.0                 ; assumed intrinsic ratio!
    Mpc2cm = 3.086D24  ; [cm/Mpc]
    irconst = 4.5D-44            ; K98 conversion L(IR) --> SFR(IR)

; plotting variables    
    
    postthick1 = 2.0
    postthick2 = 2.0
    textcolor1 = 'white'

;   pspath = getenv('PAPERSPATH')+'/literature/'
    pspath = './'
    if keyword_set(postscript) then begin
       postthick1 = 4.0
       postthick2 = 3.0
       textcolor1 = 'black'
    endif

    charsize_0 = 1.0
    charsize_1 = 1.1
    charsize_2 = 1.2
    charsize_3 = 1.3
    charsize_4 = 1.4
    charsize_5 = 1.5
    charsize_6 = 1.6
    charsize_7 = 1.7
    charsize_8 = 1.8
    charsize_9 = 1.9
    singlecharsize_0 = 2.0
    singlecharsize_1 = 2.1
    singlecharsize_2 = 2.2
    singlecharsize_3 = 2.3
    singlecharsize_4 = 2.4
    singlecharsize_5 = 2.5
    charsize_30 = 3.0

    armexten = '.ps'
    
; ---------------------------------------------------------------------------    
; SFR([OII]) vs SFR(IR)
; ---------------------------------------------------------------------------    

    psname = 'sfroii'
    if keyword_set(postscript) or keyword_set(pdf) then begin
       psfile = pspath+psname+armexten
       if keyword_set(postscript) then splog, 'Writing '+psfile
    endif else delvarx, psfile

    xmargin = [1.1,0.3] & ymargin = [0.3,1.0]
    xspace = 0.0 & yspace = 0.0
    width = 4.5*[1,1] & height = 4.5
    xpage = total(height)+total(ymargin)+total(yspace)
    ypage = total(width)+total(xmargin)+total(xspace)

    xspace1 = xspace & yspace1 = yspace
    xmargin1 = xmargin & ymargin1 = ymargin
    width1 = width & height1 = height
    xpage1 = xpage & ypage1 = ypage
    
    arm_plotconfig, landscape=1, nx=2, ny=1, xmargin=xmargin1, $
      ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
      height=height1, coord=pos, xpage=xpage1, ypage=ypage1, $
      psfile=psfile, /writeover, bw=0l;, /show
    cleanplot, /silent

    xtitle1 = 'log \psi([O II])_{uncor} (M_{\odot} yr^{-1})'
    ytitle1 = 'log \psi(IR) (M_{\odot} yr^{-1})'
    xrange1 = [-3.0,2.9]
    yrange1 = [-3.0,2.9]

    xtitle2 = 'log \psi([O II])_{MK06} (M_{\odot} yr^{-1})'
    ytitle2 = 'log \psi(IR) (M_{\odot} yr^{-1})'
    xrange2 = [-3.0,2.9]
    yrange2 = [-3.0,2.9]

    cc = iclassification(atlas1,snrcut_class=5.0,/kauffmann)
    sf = where(cc.bpt_nii_kauffmann_class eq 4)
    agn = where(cc.bpt_nii_kauffmann_class eq 2)
    unk = where(cc.bpt_nii_kauffmann_class eq 0)
    atlas = atlas1[sf]

    indx = where((atlas.rc3_b_lum gt -900.0) and (atlas.oii_3727[0]/atlas.oii_3727[1] gt 3.0) and $
      (atlas.ir_lum gt -900.0))

    oii_lum = alog10(atlas[indx].oii_3727[0]*(4.0*!dpi*(atlas[indx].distance*Mpc2cm)^2.0))
    oii_uncor = oii_lum + alog10(haconst) + alog10(haoii)
    oii_best = alog10(oii_sfr(atlas[indx].rc3_b_lum,oii_lum))
    ir = atlas[indx].ir_lum + alog10(lsun) + alog10(irconst)

    djs_plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange1, $
      xtitle=textoidl(xtitle1), ytitle=textoidl(ytitle1), xsty=1, ysty=1, $
      charsize=charsize_8, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,0]
    djs_oplot, !x.crange, !y.crange, line=5, thick=2.0
    im_symbols, 108, psize=1.0, fill=1, color=fsc_color('dodger blue',102)
    djs_oplot, oii_uncor[sf], ir[sf], ps=8
    im_symbols, 106, psize=1.0, fill=1, color=fsc_color('forest green',101)
    djs_oplot, oii_uncor[agn], ir[agn], ps=8

    djs_plot, [0], [0], /nodata, /noerase, xrange=xrange1, yrange=yrange1, $
      xtitle=textoidl(xtitle2), ytitle='', ytickname=replicate(' ',10), xsty=1, ysty=1, $
      charsize=charsize_8, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,1]
    djs_oplot, !x.crange, !y.crange, line=5, thick=2.0
    im_symbols, 108, psize=1.0, fill=1, color=fsc_color('dodger blue',102)
    djs_oplot, oii_best[sf], ir[sf], ps=8
    im_symbols, 106, psize=1.0, fill=1, color=fsc_color('forest green',101)
    djs_oplot, oii_best[agn], ir[agn], ps=8
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    

return
end
    
