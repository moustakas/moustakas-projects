pro mzplot_cmdevol, ps=ps
; jm09apr23nyu - CMD evolution

; read the data    
    ageskcorr = read_mz_emline_sample(/mzhiiplus_ancillary)
    agesohdust = read_mz_log12oh_sample()

    mzpath = ages_path(/projects)+'mz/'
    pspath = ages_path(/papers)+'mz/FIG_MZ/'

    encap = 1 ; default is to make EPS files
    if keyword_set(encap) then begin
       ps = 1
       suffix = '.eps' 
    endif else suffix = '.ps'

    mgrange1 = [-16.1,-24.5]
    grrange1 = [0.0,1.3]

    mbrange1 = [-16.1,-24.5]
    ubrange1 = [0.21,1.59]

    mgtitle1 = textoidl('M_{0.1g} - 5 log (h_{70})')
    grtitle1 = textoidl('^{0.1}(g - r)')
    mbtitle1 = textoidl('M_{B} - 5 log (h_{70})')
    ubtitle1 = textoidl('(U - B)')

; M_g vs g-r (AB) division between red and blue galaxies
    mgaxis1 = im_array(-24.0,-17.0,0.01)
    grseq1 = poly(mgaxis1,[0.46,-0.022])

; M_B vs U-B (AB) division between red and blue galaxies
    uv2ab = k_vega2ab(filterlist='bessell_U.par',/silent,/kurucz)
    bv2ab = k_vega2ab(filterlist='bessell_B.par',/silent,/kurucz)
    
    mbaxis1 = im_array(-24.0,-20.0,0.01)
    ubseq1 = poly(mbaxis1+21.52,[0.454-0.25+(uv2ab-bv2ab),-0.032])
    ubseq2 = ubseq1*0.0+0.65
    ubseq3 = ubseq1*0.0+0.79
    
; ---------------------------------------------------------------------------    
; M_B vs U-B CMD evolution

; 2x3 panel plot
    if keyword_set(ps) then psfile = pspath+'cmd_mb_vs_ub'+suffix
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.7, $
      xmargin=[1.1,0.2], width=[3.6,3.6]

; zbin0    
    zmin = 0.05 & zmax = 0.15
    aa0 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    mzages_hogg_scatterplot, aa0.mb_ab, aa0.ub_ab, position=pos[*,0], $
      xstyle=1, ystyle=1, xtitle='', ytitle=ubtitle1, $
      xrange=mbrange1, yrange=ubrange1, xtickname=replicate(' ',10)
    djs_oplot, mbaxis1, ubseq1, line=0, color='blue', thick=6.0
    djs_oplot, mbaxis1, ubseq2, line=5, color='blue', thick=6.0
    djs_oplot, mbaxis1, ubseq3, line=5, color='blue', thick=6.0
    im_legend, '0.05<z<0.15', /left, /top, box=0, charsize=1.5, margin=0

; zbin1
    zmin = 0.15 & zmax = 0.25
    aa1 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    mzages_hogg_scatterplot, aa1.mb_ab, aa1.ub_ab, /noerase, position=pos[*,1], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=mbrange1, yrange=ubrange1, xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    djs_oplot, mbaxis1, ubseq1, line=0, color='blue', thick=6.0
    djs_oplot, mbaxis1, ubseq2, line=5, color='blue', thick=6.0
    djs_oplot, mbaxis1, ubseq3, line=5, color='blue', thick=6.0
    im_legend, '0.15<z<0.25', /left, /top, box=0, charsize=1.5, margin=0

; zbin2
    zmin = 0.25 & zmax = 0.35
    aa2 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    mzages_hogg_scatterplot, aa2.mb_ab, aa2.ub_ab, /noerase, position=pos[*,2], $
      xstyle=1, ystyle=1, xtitle='', ytitle=ubtitle1, $
      xrange=mbrange1, yrange=ubrange1, xtickname=replicate(' ',10)
    djs_oplot, mbaxis1, ubseq1, line=0, color='blue', thick=6.0
    djs_oplot, mbaxis1, ubseq2, line=5, color='blue', thick=6.0
    djs_oplot, mbaxis1, ubseq3, line=5, color='blue', thick=6.0
    im_legend, '0.25<z<0.35', /left, /top, box=0, charsize=1.5, margin=0

; zbin3
    zmin = 0.35 & zmax = 0.45
    aa3 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    mzages_hogg_scatterplot, aa3.mb_ab, aa3.ub_ab, /noerase, position=pos[*,3], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=mbrange1, yrange=ubrange1, xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    djs_oplot, mbaxis1, ubseq1, line=0, color='blue', thick=6.0
    djs_oplot, mbaxis1, ubseq2, line=5, color='blue', thick=6.0
    djs_oplot, mbaxis1, ubseq3, line=5, color='blue', thick=6.0
    im_legend, '0.35<z<0.45', /left, /top, box=0, charsize=1.5, margin=0

; zbin4
    zmin = 0.45 & zmax = 0.55
    aa4 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    mzages_hogg_scatterplot, aa4.mb_ab, aa4.ub_ab, /noerase, position=pos[*,4], $
      xstyle=1, ystyle=1, xtitle=mbtitle1, ytitle=ubtitle1, $
      xrange=mbrange1, yrange=ubrange1
    djs_oplot, mbaxis1, ubseq1, line=0, color='blue', thick=6.0
    djs_oplot, mbaxis1, ubseq2, line=5, color='blue', thick=6.0
    djs_oplot, mbaxis1, ubseq3, line=5, color='blue', thick=6.0
    djs_oplot, -20*[1,1], !y.crange, line=1, thick=10
    im_legend, '0.45<z<0.55', /left, /top, box=0, charsize=1.5, margin=0

; zbin5
    zmin = 0.55 & zmax = 0.75
    aa5 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    mzages_hogg_scatterplot, aa5.mb_ab, aa5.ub_ab, /noerase, position=pos[*,5], $
      xstyle=1, ystyle=1, xtitle=mbtitle1, ytitle='', $
      xrange=mbrange1, yrange=ubrange1, ytickname=replicate(' ',10)
    djs_oplot, mbaxis1, ubseq1, line=0, color='blue', thick=6.0
    djs_oplot, mbaxis1, ubseq2, line=5, color='blue', thick=6.0
    djs_oplot, mbaxis1, ubseq3, line=5, color='blue', thick=6.0
    djs_oplot, -20*[1,1], !y.crange, line=1, thick=10
    im_legend, '0.55<z<0.75', /left, /top, box=0, charsize=1.5, margin=0

    im_plotconfig, /psclose

; ---------------------------------------------------------------------------    
; M_g vs g-r CMD evolution

; 2x3 panel plot
    if keyword_set(ps) then psfile = pspath+'cmd_mg_vs_gr'+suffix
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.7, $
      xmargin=[1.1,0.2], width=[3.6,3.6]

; zbin0    
    zmin = 0.05 & zmax = 0.15
    aa0 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    mzages_hogg_scatterplot, aa0.mr_ab, aa0.gr_ab, position=pos[*,0], $
      xstyle=1, ystyle=1, xtitle='', ytitle=grtitle1, $
      xrange=mgrange1, yrange=grrange1, xtickname=replicate(' ',10)
    djs_oplot, mgaxis1, grseq1, line=0, color='blue', thick=6.0
    im_legend, '0.05<z<0.15', /left, /top, box=0, charsize=1.5, margin=0

; zbin1
    zmin = 0.15 & zmax = 0.25
    aa1 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    mzages_hogg_scatterplot, aa1.mr_ab, aa1.gr_ab, /noerase, position=pos[*,1], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=mgrange1, yrange=grrange1, xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    djs_oplot, mgaxis1, grseq1, line=0, color='blue', thick=6.0
    im_legend, '0.15<z<0.25', /left, /top, box=0, charsize=1.5, margin=0

; zbin2
    zmin = 0.25 & zmax = 0.35
    aa2 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    mzages_hogg_scatterplot, aa2.mr_ab, aa2.gr_ab, /noerase, position=pos[*,2], $
      xstyle=1, ystyle=1, xtitle='', ytitle=grtitle1, $
      xrange=mgrange1, yrange=grrange1, xtickname=replicate(' ',10)
    djs_oplot, mgaxis1, grseq1, line=0, color='blue', thick=6.0
    im_legend, '0.25<z<0.35', /left, /top, box=0, charsize=1.5, margin=0

; zbin3
    zmin = 0.35 & zmax = 0.45
    aa3 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    mzages_hogg_scatterplot, aa3.mr_ab, aa3.gr_ab, /noerase, position=pos[*,3], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=mgrange1, yrange=grrange1, xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    djs_oplot, mgaxis1, grseq1, line=0, color='blue', thick=6.0
    im_legend, '0.35<z<0.45', /left, /top, box=0, charsize=1.5, margin=0

; zbin4
    zmin = 0.45 & zmax = 0.55
    aa4 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    mzages_hogg_scatterplot, aa4.mr_ab, aa4.gr_ab, /noerase, position=pos[*,4], $
      xstyle=1, ystyle=1, xtitle=mgtitle1, ytitle=grtitle1, $
      xrange=mgrange1, yrange=grrange1
    djs_oplot, mgaxis1, grseq1, line=0, color='blue', thick=6.0
    im_legend, '0.45<z<0.55', /left, /top, box=0, charsize=1.5, margin=0

; zbin5
    zmin = 0.55 & zmax = 0.75
    aa5 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    mzages_hogg_scatterplot, aa5.mr_ab, aa5.gr_ab, /noerase, position=pos[*,5], $
      xstyle=1, ystyle=1, xtitle=mgtitle1, ytitle='', $
      xrange=mgrange1, yrange=grrange1, ytickname=replicate(' ',10)
    djs_oplot, mgaxis1, grseq1, line=0, color='blue', thick=6.0
    im_legend, '0.55<z<0.75', /left, /top, box=0, charsize=1.5, margin=0

    im_plotconfig, /psclose

return
end
    
