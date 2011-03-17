function sfrd_frac, sfr, inzaxis, outzaxis, plotzaxis
; jm09apr17nyu - compute the cumulative fraction of stellar mass (or
; metals) formed (synthesized) versus time/redshift

    lookback = getage(0.0)-getage(outzaxis) ; lookback time
    ntime = n_elements(lookback)
    
    sfr = interpol(sfr,inzaxis,outzaxis)
    sfrd_total = im_integral(lookback,sfr*1D9,$ ; [M_sun]
      min(lookback),max(lookback)) 
    sfrd_time = im_integral(lookback,sfr*1D9,$  ; [M_sun]
      lookback,replicate(max(lookback),ntime))
    sfrd_frac = sfrd_time/sfrd_total

; extrapolate to z<0.1    
;   plotzaxis = im_array(0.0,0.6,0.01)
    plotzaxis = [0.0,outzaxis]
    notzero = where(sfrd_frac gt 0.0)
    sfrd_frac = interpol(1.0*(sfrd_frac[notzero]),$
      outzaxis[notzero],plotzaxis)
    
return, sfrd_frac
end

function zformpeg, peg, zform=pegzform, peg_zaxis=peg_zaxis
; jm09apr17nyu - attach a formation redshift to a Pegase model

    ppeg = peg
    rev = reverse(sort(ppeg.age)) ; reverse the time array!
    ppeg = ppeg[rev]

    pegtform = getage(pegzform)
    peg_zaxis = getredshift(ppeg.age/1D3+pegtform)
    peg_lookback = getage(0.0)-getage(peg_zaxis)

    peg_good = where(peg_zaxis gt 0.0)
    peg_zaxis = peg_zaxis[peg_good]
    peg_lookback = peg_lookback[peg_good]
    ppeg = ppeg[peg_good]
    
return, ppeg
end
    
pro osu_mzsdss_hogg_scatterplot, x, y, _extra=extra, keynote=keynote, $
  hackit=hackit
; mass-metallicity plot wrapper on HOGG_SCATTERPLOT
    if keyword_set(keynote) then begin
       loadct, 1, /silent
       acolor = 'white'
       ccolor = 'black'
       outcolor = 'black'
       outsymsize = 0.15
    endif else begin
       loadct, 0, /silent
       acolor = 'black'
       ccolor = 'black'
       outcolor = 'grey'
       outsymsize = 0.3
    endelse
;   plot, [0], [0], /nodata, xsty=5, ysty=5, _extra=extra, $
;     ytitle='', xtitle='', xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    eextra = struct_trimtags(extra,except=['POSITION'])
    pos = extra.position
    if keyword_set(hackit) then $
      pos[2] = pos[2]-0.015
    hogg_scatterplot, x, y, /outliers, outpsym=symcat(16), $
      outsymsize=outsymsize, outcolor=outcolor, _extra=eextra, $
      color=djs_icolor(acolor), ccolor=djs_icolor(ccolor), $
      position=pos
; overplot black axes
    eextra = struct_trimtags(extra,except=['XTITLE','YTITLE'])
    plot, [0], [0], /noerase, xsty=1, ysty=1, _extra=eextra, $
      color=djs_icolor(ccolor), ytitle='', xtitle='', $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    if keyword_set(keynote) then loadct, 0, /silent
return
end

pro osu_mzages_hogg_scatterplot, x, y, _extra=extra, keynote=keynote
; mass-metallicity plot wrapper on HOGG_SCATTERPLOT
    if keyword_set(keynote) then begin
       loadct, 3, /silent
       acolor = 'white'
       ccolor = 'black'
    endif else begin
       loadct, 0, /silent
       acolor = 'black'
       ccolor = 'black'
    endelse
    hogg_scatterplot, x, y, /outliers, outpsym=symcat(16), $
      outsymsize=0.4, outcolor='grey', _extra=extra, $
      color=djs_icolor(acolor), ccolor=djs_icolor(ccolor)
; overplot black axes
    eextra = struct_trimtags(extra,except=['XTITLE','YTITLE'])
    plot, [0], [0], /noerase, xsty=1, ysty=1, _extra=eextra, $
      color=djs_icolor(ccolor), ytitle='', xtitle='', $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    if keyword_set(keynote) then loadct, 0, /silent
return
end

pro osu_ages_oplot, x, y, _extra=extra, keynote=keynote
    if keyword_set(keynote) then cc = 'cyan' else cc = 'grey'
    djs_oplot, x, y, psym=symcat(16), symsize=0.4, $
      color=cc, _extra=extra
return
end

pro talk_09apr_osu, keynote=keynote
; jm09apr09nyu - plots for my colloquium

    pspath = '~/home/research/talks/2009/09apr_osu/'
    if keyword_set(keynote) then begin
       pspath = pspath+'keynote/'
       textcolor = 'white'
    endif else begin
       textcolor = 'black'
    endelse

    litpath = getenv('PAPERSPATH')+'/literature/data/'

    red, h100=0.7
    mzpath = ages_path(/projects)+'mz/'

    ageskcorr = read_mz_emline_sample(/mzhiiplus_ancillary)
    agesohdust = read_mz_log12oh_sample()
    sdsskcorr = read_mz_emline_sample(/mzhii_ancillary,/sdss)
    sdssohnodust = read_mz_log12oh_sample(/sdss,/nodust)

    sdssages_noevol0 = read_mz_sdss2ages_sample(/zbin1)
    sdssages_noevol1 = read_mz_sdss2ages_sample(/zbin2)
    sdssages_noevol2 = read_mz_sdss2ages_sample(/zbin3)
    sdssages_noevol3 = read_mz_sdss2ages_sample(/zbin4)
    sdssages_noevol4 = read_mz_sdss2ages_sample(/zbin5)
    sdssages_noevol5 = read_mz_sdss2ages_sample(/zbin6)

    sdssages_levol0 = read_mz_sdss2ages_sample(/evolve,/zbin1)
    sdssages_levol1 = read_mz_sdss2ages_sample(/evolve,/zbin2)
    sdssages_levol2 = read_mz_sdss2ages_sample(/evolve,/zbin3)
    sdssages_levol3 = read_mz_sdss2ages_sample(/evolve,/zbin4)
    sdssages_levol4 = read_mz_sdss2ages_sample(/evolve,/zbin5)
    sdssages_levol5 = read_mz_sdss2ages_sample(/evolve,/zbin6)

;   sdssparent = read_sdss_vagc_mpa(/ispec,sample='dr7',$
;     letter='bsafe',poststr='32')

; read the output from FIT_MZLZEVOL
    ohevol = mrdfits(mzpath+'ohevol.fits.gz',1)

; representative Pegase models
    splog, 'Reading the Pegase models'
    pegpath = getenv('PEGASE_HR_SFHGRID_DIR')+'/measure/'
    peg1 = mrdfits(pegpath+'salp_tau_001.0Gyr.info.fits',1,silent=0)
    peg2 = mrdfits(pegpath+'salp_tau_001.0Gyr.info.fits',1,silent=0)
    peg3 = mrdfits(pegpath+'salp_tau_003.0Gyr.info.fits',1,silent=0)
    ipeg3 = mrdfits(pegpath+'salp_kennlaw_0.10_infall_003.0Gyr.info.fits',1,silent=0)
    pegconst = mrdfits(pegpath+'salp_tau_999.0Gyr.info.fits',1,silent=0)

    zf = 5.0
    zf2_peg1 = zformpeg(peg1,zform=zf,peg_zaxis=zf2_peg_zaxis)
    zf2_peg2 = zformpeg(peg2,zform=zf)
    zf2_peg3 = zformpeg(peg3,zform=zf)
    zf2_ipeg3 = zformpeg(ipeg3,zform=zf)
    zf2_pegconst = zformpeg(pegconst,zform=zf)

    zf = 1.5
    zf1_peg1 = zformpeg(peg1,zform=zf,peg_zaxis=zf1_peg_zaxis)
    zf1_peg2 = zformpeg(peg2,zform=zf)
    zf1_peg3 = zformpeg(peg3,zform=zf)
    zf1_ipeg3 = zformpeg(ipeg3,zform=zf)
    zf1_pegconst = zformpeg(pegconst,zform=zf)
    
    massrange1 = [8.1,12.0]
    massrange2 = [8.35,11.75]
    mgrange1 = [-16.1,-23.9]
    ohrange1 = [8.35,9.37]
    ohrange2 = [8.5,9.3]

    masstitle1 = textoidl('log (M_{*}/M'+sunsymbol()+')')
    mgtitle1 = textoidl('M_{0.1g} - 5 log (h_{70})')
    ohtitle1 = textoidl('12 + log (O/H)')

    linestyle1 = 0 & linecolor1 = 'black'
    linestyle2 = 5  & linecolor2 = 'dodger blue' ; 'firebrick'
    linestyle3 = 0 & linecolor3 = 'black'
    linestyle4 = 3 & linecolor4 = 'black'
    linestyle5 = 5 & linecolor5 = 'dodger blue' ; 'firebrick'

    timelabel1 = [1.0,3.0,5.0,7.0] ; [Gyr]
    zaxis1 = im_array(0.0,0.8,0.01)

    ps = 1
    suffix = '.ps'

; ------------------------------------------------------------
; luminosity and mass evolution for blue galaxies (literature
; compilation) 

    if keyword_set(ps) then psfile = pspath+'z_vs_mr_mass_lit'+suffix
    im_plotconfig, 6, pos, psfile=psfile, charsize=2.0, $
      xmargin=[1.5,0.1], ymargin=[1.0,1.1], keynote=keynote

    h100 = 0.7
    B2g01 = +0.0759 ; ^{0.1}g = B+0.0759+0.0620*[(B-V)-0.5870] [AB, Blanton & Roweis]
    B2r01 = -0.6429 ; ^{0.1}r = B-0.6429-1.0845*[(B-V)-0.5870] [AB, Blanton & Roweis]
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par',/kurucz,/silent))[0]

    xtitle = 'Redshift'
    ytitle1 = textoidl('M_{0.1r}^{*}')
    ytitle2 = textoidl('log (M^{*}/M'+sunsymbol()+')')

    xrange = [-0.02,0.82]
    yrange1 = [-20.2,-22.4]
    yrange2 = [10.6,11.29]

; luminosity function colors from the literature
    if keyword_set(keynote) then begin
       w06color   = 'orange'      & w06psym = 16 & w06symsize = 2.4
       f06color   = 'dodger blue' & f06psym = 14 & f06symsize = 2.8
       e07color   = 'firebrick'   & e07psym = 15 & e07symsize = 2.0
       i09color   = 'firebrick'   & i09psym = 15 & i09symsize = 2.0
       b06color   = 'dodger blue' & b06psym = 14 & b06symsize = 2.8
    endif else begin
       w06color   = 'orange'    & w06psym = 16 & w06symsize = 2.4
       f06color   = 'navy'      & f06psym = 14 & f06symsize = 2.8
       e07color   = 'firebrick' & e07psym = 15 & e07symsize = 2.0
       i09color   = 'firebrick' & i09psym = 15 & i09symsize = 2.0
       b06color   = 'navy'      & b06psym = 14 & b06symsize = 2.8
    endelse

; luminosity evolution plot    
    djs_plot, [0], [0], /nodata, ysty=1, xsty=9, xrange=xrange, $
      yrange=yrange1, xtitle='', ytitle=ytitle1, position=pos[*,0], $
      xtickname=replicate(' ',10)
    axis, /xaxis, xsty=1, xrange=xrange, xtickv=getredshift(getage(0.0)-timelabel1), $
      xtitle='Lookback Time (Gyr)', xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)')

; overplot the line for 1.5 mag/z of luminosity evolution    
    djs_oplot, zaxis1, poly(zaxis1-0.1,[-20.95,-1.5]), line=0, thick=6.0

; -------------------------
; Eisenstein et al. 2009 [AGES]; Omega_0=0.3, Omega_lamba=0.7, h=1.0;
; AB; alpha fixed at -1.10 for all blue galaxies; evolving color cut
; based on the (u0.1-r0.1) color: A=u0.1-r0.1+0.08(M_r0.1+20);  

    z_e07    = [0.1,0.2,0.3,0.4,0.5,0.65]
    zerr_e07 = [0.05,0.05,0.05,0.05,0.05,0.1]

    mstar_e07 = [-19.95,-20.39,-20.43,-20.57,-20.86,-20.83] + 5.0*alog10(h100) ; h=1-->h=0.7
    mstarerr_e07 = [0.12,0.07,0.04,0.05,0.11,0.08]

    oploterror, z_e07, mstar_e07, zerr_e07, mstarerr_e07, $
      errthick=6.0, symsize=e07symsize, psym=symcat(e07psym), $
      color=fsc_color(e07color,1E8), errcolor=fsc_color(e07color,1E8)

; -------------------------
; Willmer et al. 2006 [DEEP2]; "minimal" model adopted, as
; recommended; Omega_0=0.3, Omega_lamba=0.7, h=0.7; Vega; alpha=-1.30
; fixed for all blue galaxies, based on a fit to all galaxies at
; 0.2<z<0.6 in COMBO-17 (Faber et al. 2005); color cut: (U-B) =
; -0.032(M_B+21.52)+0.454-0.25; drop the last bin

    z_w06 = [0.3,0.5,0.7,0.9,1.1];,1.3]
    zerr_w06 = [0.1,0.1,0.1,0.1,0.1];,0.1]

    mstar_w06 = [-20.36,-20.72,-21.15,-21.21,-21.38] + Bvega2ab + B2r01 ; Vega-->AB; B-->r0.1 ; ,-21.86]
    mstarerr_up_w06 = [0.13,0.05,0.07,0.01,0.04];,0.07]
    mstarerr_lo_w06 = [0.11,0.07,0.07,0.03,0.05];,0.08]
    mstarerr_w06 = mstar_w06*0.0
    for i = 0, n_elements(mstar_w06)-1L do mstarerr_w06[i] = mean([mstarerr_up_w06[i],mstarerr_lo_w06[i]])
    
    oploterror, z_w06, mstar_w06, zerr_w06, mstarerr_up_w06, /hi, symsize=w06symsize, $
      errthick=6.0, psym=symcat(w06psym), color=fsc_color(w06color,1E8), errcolor=fsc_color(w06color,1E8)
    oploterror, z_w06, mstar_w06, zerr_w06, mstarerr_lo_w06, /lo, symsize=w06symsize, $
      errthick=6.0, psym=symcat(w06psym), color=fsc_color(w06color,1E8), errcolor=fsc_color(w06color,1E8)

; -------------------------
; Faber et al. 2007 [COMBO17]; Omega_0=0.3, Omega_lamba=0.7, h=0.7;
; Vega; alpha=-1.3 fixed for all blue galaxies;

    z_f06 = [0.3,0.5,0.7,0.9,1.1]
    zerr_f06 = [0.1,0.1,0.1,0.1,0.1]

    mstar_f06 = [-20.74,-21.10,-21.30,-21.10,-21.25] + Bvega2ab + B2r01 ; Vega-->AB; B-->r0.1
    mstarerr_f06 = [0.20,0.15,0.16,0.17,0.18]

    oploterror, z_f06, mstar_f06, zerr_f06, mstarerr_f06, symsize=f06symsize, $
      errthick=6.0, psym=symcat(f06psym), color=fsc_color(f06color,1E8), errcolor=fsc_color(f06color,1E8)
    
; legend
    label = ['Willmer+06','Faber+07','Eisenstein+09']
    im_legend, label, /left, /top, box=0, psym=[w06psym,f06psym,e07psym], $
      symsize=[w06symsize,f06symsize,e07symsize], color=[w06color,f06color,e07color], $
      spacing=2.0, charsize=1.7, textcolor=textcolor
    im_legend, 'Q = 1.5 mag/z', /right, /bottom, box=0, $
      line=0, thick=6.0, charsize=1.7, textcolor=textcolor, $
      pspacing=1.4, margin=0, color=textcolor

; mass evolution plot    
    djs_plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, xrange=xrange, $
      yrange=yrange2, xtitle=xtitle, ytitle=ytitle2, position=pos[*,1], $
      yminor=3

; overplot the line for 1.5 mag/z of luminosity evolution    
    djs_oplot, zaxis1, poly(zaxis1-0.1,[10.87,0.25]), line=5, thick=6.0

; -------------------------
; Ilbert+09 [COSMOS]; Omega_0=0.3, Omega_lamba=0.7, h=0.7; Chabrier
; IMF; see ILBERT_MASS_FUNCTION for details

    z_i09 = [0.3,0.5,0.7,0.9,1.1,1.35,1.75]
    zerr_i09 = [0.1,0.1,0.1,0.1,0.1,0.15,0.25]
    
; intermediate-activity star-forming galaxies    
;   mstar_i09 = [11.01,10.95,10.89,10.78,10.73,10.76,10.86]
;   mstarerr_i09 = [0.03,0.03,0.02,0.02,0.02,0.02,0.02]
; high-activity star-forming galaxies    
;   mstar_i09 = [10.26,10.30,10.39,10.42,10.42,10.49,10.86]
;   mstarerr_i09 = [0.05,0.03,0.02,0.02,0.02,0.01,0.5]
; all star-forming galaxies    
;   mstar_i09 = [11.01,10.96,10.92,10.83,10.80,10.81,10.94]
;   mstarerr_i09 = mstar_i09*0.0+0.1
; all star-forming galaxies - from Ilbert
    mstar_i09 = [10.98,10.97,10.94,10.89,10.86,10.79,10.89]
    mstarerr_i09 = [0.029,0.026,0.022,0.017,0.017,0.015,0.015]

    oploterror, z_i09, mstar_i09, zerr_i09, mstarerr_i09, symsize=i09symsize, $
      errthick=6.0, psym=symcat(i09psym), color=fsc_color(i09color,1E8), errcolor=fsc_color(i09color,1E8)

; -------------------------
; Borch+06 [COMBO-17]; Omega_0=0.3, Omega_lamba=0.7, h=0.7; Kroupa+93
; IMF (same as Chabrier) 
    
    z_b06 = [0.1,0.3,0.5,0.7,0.9]
    zerr_b06 = [0.1,0.1,0.1,0.1,0.1]

    mstar_b06 = [10.80,10.86,10.93,11.07,11.00]
    mstarerr_b06 = [0.10,0.10,0.19,0.12,0.07,0.25]

    oploterror, z_b06, mstar_b06, zerr_b06, mstarerr_b06, symsize=b06symsize, $
      errthick=6.0, psym=symcat(b06psym), color=fsc_color(b06color,1E8), errcolor=fsc_color(b06color,1E8)

; fit to both samples
    xx = [z_b06,z_i09] & yy = [mstar_b06,mstar_i09] & yyerr = [mstarerr_b06,mstarerr_i09]
    ww = where(xx lt 0.8)
;   print, linfit(xx[ww],yy[ww],measure_errors=yyerr[ww])
    
; legend
    label = ['Borch+06','Ilbert+09']
    im_legend, label, /left, /top, box=0, psym=[b06psym,i09psym], $
      symsize=[b06symsize,i09symsize], color=[b06color,i09color], $
      spacing=2.0, charsize=1.7, textcolor=textcolor
    im_legend, 'R = 0.25 dex/z', /right, /bottom, box=0, $
      line=5, thick=6.0, charsize=1.7, textcolor=textcolor, $
      pspacing=1.4, margin=0, color=textcolor
    
    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif

stop    
    
; ------------------------------------------------------------
; redshift vs fraction of stars/metals formed
    
    if keyword_set(ps) then psfile = pspath+'z_vs_ohfraction'+suffix
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, width=6.5, $
      height=6.5, xmargin=[1.6,0.4], ymargin=[1.0,1.0], thick=6.0, keynote=keynote

    xrange = [0.0,0.77]
    yrange = [0.74,1.0499]
;   yrange = [73.0,104.0] ; per cent
;   yrange = [60.0,109.0] ; per cent
    xtitle = 'Redshift'
    ytitle = 'Fraction of Metals Synthesized by z=0.1'
    
; coefficients from plot_redshift_vs_sfr_mass_density    
    sfrd_zaxis = findgen((6.0-0.1)/0.01+1)*0.01+0.1 ; for integral
    sfrd_lookback = getage(0.0)-getage(sfrd_zaxis)  ; lookback time
    
    coeff = [0.010464392D,0.11303816D,3.4534319D,3.7109016D]
    sfrd = (coeff[0]+coeff[1]*sfrd_zaxis)/(1.0+(sfrd_zaxis/coeff[2])^coeff[3])
;   plot, alog10(1+sfrd_zaxis), alog10(sfrd), xsty=3, ysty=3, yr=[-2.3,-0.3]
    
; what fraction of the stars (and therefore metals) were formed by
; various redshifts; integrate the model; convert the SFRD to
; M_sun/Gyr; the factor of 1/64 (1.6%) is the IMF-weighted yield based
; on the Woosley & Weaver (1995) models, assuming a Salpeter IMF
; 0.1-125 M_sun (Conti et al. 2003; Pettini et al. 2006); note that
; Madau et al. (1996) used 1/42 (2.4%) for a Salpeter 0.1-100 IMF;
; also note that the exact value is ~independent of the IMF; and in
; detail, of course, the factor drops out because we are doing
; fractions here
    sfrd_frac = sfrd_frac(sfrd,sfrd_zaxis,sfrd_zaxis,plot_zaxis)

; interpolate the Pegase models and then integrate; note that I'm not
; dividing by the co-moving volume, but I could
    zf2_peg1_sfrd_frac = sfrd_frac(zf2_peg1.sfr,zf2_peg_zaxis,sfrd_zaxis)
    zf1_peg1_sfrd_frac = sfrd_frac(zf1_peg1.sfr,zf1_peg_zaxis,sfrd_zaxis)

    zf2_peg3_sfrd_frac = sfrd_frac(zf2_peg3.sfr,zf2_peg_zaxis,sfrd_zaxis)
    zf1_peg3_sfrd_frac = sfrd_frac(zf1_peg3.sfr,zf1_peg_zaxis,sfrd_zaxis)

    zf2_pegconst_sfrd_frac = sfrd_frac(zf2_pegconst.sfr,zf2_peg_zaxis,sfrd_zaxis)
    zf1_pegconst_sfrd_frac = sfrd_frac(zf1_pegconst.sfr,zf1_peg_zaxis,sfrd_zaxis)

; now make the plot    
    djs_plot, [0], [0], /nodata, ysty=1, xsty=9, xrange=xrange, $
      yrange=yrange, xtitle=xtitle, ytitle=ytitle, position=pos
    axis, /xaxis, xsty=1, xrange=xrange, xtickv=getredshift(getage(0.0)-timelabel1), $
      xtitle='Lookback Time (Gyr)', xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)')
    djs_oplot, plot_zaxis, sfrd_frac, line=0, thick=12, color=djs_icolor(textcolor)

; overplot the Pegase models
; tau=1
;   djs_oplot, plot_zaxis, zf2_peg1_sfrd_frac, line=5, $
;     color=fsc_color('dodger blue',148), thick=8
    djs_oplot, plot_zaxis, zf1_peg1_sfrd_frac, line=5, $
      color=fsc_color('dodger blue',148), thick=8

; tau=3    
;   djs_oplot, plot_zaxis, zf2_peg3_sfrd_frac, line=3, $
;     color=fsc_color('orange',149), thick=8
    djs_oplot, plot_zaxis, zf1_peg3_sfrd_frac, line=3, $
      color=fsc_color('orange',149), thick=8

; tau=inf    
;   djs_oplot, plot_zaxis, zf2_pegconst_sfrd_frac, line=1, $
;     color=fsc_color('purple',150), thick=10
    djs_oplot, plot_zaxis, zf1_pegconst_sfrd_frac, line=1, $
      color=fsc_color('yellow',150), thick=10

; and finally the data    
    color1 = 'firebrick' & psym1 = 6 & symsize1 = 3.0 & line1 = 3
    color2 = 'royal blue' & psym2 = 9 & symsize2 = 3.0 & line2 = 5
    color3 = 'forest green' & psym3 = 4 & symsize3 = 3.5 & line3 = 0
    color4 = 'orange'       & psym3 = 7 & symsize3 = 3.5 & line3 = 4

    zerr = [0.05,0.05,0.05,0.05,0.05,0.1] ; no error bar on the first point
    z1 = ohevol.z
    z2 = ohevol.z-0.015
    z3 = ohevol.z+0.015
    z4 = ohevol.z+0.015

    oh1 = 1.0*10^ohevol.ldlogoh_noevol_cor
    oh2 = 1.0*10^ohevol.ldlogoh_levol_cor
    oh3 = 1.0*10^ohevol.mdlogoh_noevol_cor
    oh4 = 1.0*10^ohevol.mdlogoh_levol_cor

    oh1_err = ohevol.ldlogoh_noevol_cor_err*oh1*alog(10.0)
    oh2_err = ohevol.ldlogoh_levol_cor_err*oh2*alog(10.0)
    oh3_err = ohevol.mdlogoh_noevol_cor_err*oh3*alog(10.0)
    oh4_err = ohevol.mdlogoh_levol_cor_err*oh4*alog(10.0)

;   oh1_err[0] = 0.0
;   oh2_err[0] = 0.0
;   oh3_err[0] = 0.0
    
; plot the mean values
    zz = z1
    moh = zz*0.0
    moh_err = zz*0.0
    weighted = 0
    for ii = 0, n_elements(zz)-1L do begin
       if weighted then begin
          moh[ii] = im_weighted_mean([oh2[ii],oh3[ii]],$
            [oh2_err[ii],oh3_err[ii]],wsigma=jj)
          moh_err[ii] = jj
       endif else begin
          oo = [oh1[ii],oh2[ii],oh3[ii],oh4[ii]]
;         oo = [oh2[ii],oh3[ii]]
          moh[ii] = mean(oo)
          moh_err[ii] = stddev(oo)
       endelse
    endfor

;   for ii = 0, n_elements(zz)-1L do im_oplot_box, 1.9*zerr[ii], $
;     2*moh_err[ii], 0.0, xoffset=zz[ii], yoffset=moh[ii]
    oploterror, zz, moh, zerr, moh_err, psym=-symcat(6,thick=10), $
      symsize=5.0, color=djs_icolor('cyan'), errcolor=djs_icolor('cyan'), $
      thick=8.0, errthick=10

; now plot the individual points    
    plotpoints = 0
    if plotpoints then begin
       oploterror, z1, oh1, zerr, oh1_err, psym=-symcat(psym1,thick=8), $
         symsize=symsize1, line=line1, color=fsc_color(color1,91), $
         errcolor=fsc_color(color1,91)
       oploterror, z2, oh2, zerr, oh2_err, psym=-symcat(psym2,thick=8), $
         symsize=symsize2, color=fsc_color(color2,92), errcolor=fsc_color(color2,92), line=line2
       oploterror, z3, oh3, zerr, oh3_err, psym=-symcat(psym3,thick=8), $
         symsize=symsize3, color=fsc_color(color3,93), errcolor=fsc_color(color3,93), line=line3
       oploterror, z4, oh4, zerr, oh4_err, psym=-symcat(psym4,thick=8), $
         symsize=symsize4, color=fsc_color(color4,94), errcolor=fsc_color(color4,94), line=line4
    endif

; label the Pegase curves
    xyouts, 0.52, 1.01, textoidl('\tau=1 Gyr'), align=0.5, $
      charsize=2.0, color=fsc_color('dodger blue',10)
    xyouts, 0.52, 0.96, textoidl('\tau=3 Gyr'), align=0.5, $
      charsize=2.0, color=fsc_color('orange',11), $
      orientation=-20
    xyouts, 0.3, 0.83, textoidl('\tau=!7y!6'), align=0.5, $
      charsize=2.2, color=fsc_color('yellow',12), $
      orientation=-60
;   xyouts, 0.6, 0.79, textoidl('!MI!N'+'d\rho_{*}(t)'), $
;     align=0.5, charsize=2.0, color=fsc_color('black',13)
;   xyouts, 0.6, 0.77, textoidl('!MI!N'+'dt\prime'+'(d\rho_{*}/dt\prime)'), $
    xyouts, 0.605, 0.78, textoidl('!MI!N'+'dt'+'(d\rho_{*}/dt)'), $
      align=0.5, charsize=1.9, color=fsc_color(textcolor,13)
    
; legend    
    label = ['\tau=1 Gyr','\tau=3 Gyr','\tau=const']
    cc = ['dodger blue','orange','yellow']
    line = [5,3,1]
    spacing = 2.2

;   label = ['!MI!N'+'dt'+'\rho_{*}(t)','\tau=1 Gyr','\tau=3 Gyr','\tau=const']
;   cc = ['','dodger blue','orange','purple']
;   line = [0,5,3,1]
;   spacing = 2.0
;;  pp = [0.01,!y.crange[0]+1.5]
    
;   im_legend, label, /left, /bottom, box=0, pspacing=1.3, color=cc, $
;     textcolor=cc, line=line, spacing=spacing, margin=0, charsize=1.7, $
;     position=pp, thick=10

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif

; ------------------------------------------------------------
; redshift vs absolute magnitude and stellar mass, with the
; completeness curves overlaid

    agesparent = read_mz_parent_sample()
    ageskcorr = read_mz_emline_sample(/mzhiiplus_ancillary)
    sdsskcorr = read_mz_emline_sample(/mzhii_ancillary,/sdss)

    psfile = pspath+'z_vs_mg_mass.ps'
    im_plotconfig, 6, pos, psfile=psfile, height=3.0*[1,1], xmargin=[1.3,0.2], $
      ymargin=[0.8,1.0], charsize=1.8, width=5.5, keynote=keynote

; compute the limiting absolute magnitude and stellar mass

; desired redshift grid
    sspzform = 5.0
    tauzform = 2.0
    csfzform = 2.0
    dz = 0.02 & minz = 0.0 & maxz = 1.0
    zz = (findgen((maxz-minz)/dz+1)*dz+minz)>0.01
    nz = n_elements(zz)

; specify the models    
    modelspath = getenv('ISEDFIT_SFHGRID_DIR')+'/basemodels/bc03/'
    sspfits = 'chab_Z0.02_tau_00.0Gyr.fits.gz'
    taufits = 'chab_Z0.02_tau_03.0Gyr.fits.gz'
    csffits = 'chab_Z0.02_tau_100Gyr.fits.gz'

    ssp = mrdfits(modelspath+sspfits,1,/silent)
    tau = mrdfits(modelspath+taufits,1,/silent)
    csf = mrdfits(modelspath+csffits,1,/silent)
    wave = csf.wave

; interpolate the models at the appropriate ages, given ZFORM and the
; desired redshift grid
    sspflux = interpolate(ssp.flux,findex(ssp.age/1D9,getage(zz)-getage(sspzform)))
    tauflux = interpolate(tau.flux,findex(tau.age/1D9,getage(zz)-getage(tauzform)))
    csfflux = interpolate(csf.flux,findex(csf.age/1D9,getage(zz)-getage(csfzform)))

; now given the apparent magnitude (I=20), compute the rest-frame
; luminosity and stellar mass at each redshift
    in_filterlist = 'ndwfs_I.par'
    out_filterlist = 'sdss_g0.par'
    sdss_band_shift = 0.1
    rsolarmag = k_solar_magnitudes(filterlist=out_filterlist,$
      band_shift=sdss_band_shift,/silent)

    I20 = 20.0 ; 19.95
    Ivega2ab = k_vega2ab(filterlist=in_filterlist,/kurucz,/silent)
    maggies = reform(10^(-0.4*(I20+Ivega2ab)),1,1)

    mgvz = replicate({z: 0.0, ssp_Mg: 0.0, tau_Mg: 0.0, csf_Mg: 0.0, $
      ssp_mass: 0.0, tau_mass: 0.0, csf_mass: 0.0},nz)
    mgvz.z = zz
    for ii = 0L, nz-1L do begin
; SSP
       kk = im_simple_kcorrect(zz[ii],maggies,maggies*0.0+1.0,$
         in_filterlist,wave,sspflux[*,ii],band_shift=sdss_band_shift,$
         out_filterlist=out_filterlist,absmag=abs)
       mgvz[ii].ssp_Mg = abs
       mgvz[ii].ssp_mass = 10^(-0.4*(abs-rsolarmag))
; CSF
       kk = im_simple_kcorrect(zz[ii],maggies,maggies*0.0+1.0,$
         in_filterlist,wave,csfflux[*,ii],band_shift=sdss_band_shift,$
         out_filterlist=out_filterlist,absmag=abs)
       mgvz[ii].csf_Mg = abs
       mgvz[ii].csf_mass = 10^(-0.4*(abs-rsolarmag))
; TAU
       kk = im_simple_kcorrect(zz[ii],maggies,maggies*0.0+1.0,$
         in_filterlist,wave,tauflux[*,ii],band_shift=sdss_band_shift,$
         out_filterlist=out_filterlist,absmag=abs)
       mgvz[ii].tau_Mg = abs
       mgvz[ii].tau_mass = 10^(-0.4*(abs-rsolarmag))
    endfor

; some other previously calculated limits    
    limits = mrdfits(mzpath+'mz_limits.fits.gz',1)

; now make the plot    
    z = ageskcorr.z
    mg = ageskcorr.ugriz_absmag[1]
    gr = ageskcorr.ugriz_absmag[1]-ageskcorr.ugriz_absmag[2]
    mass = ageskcorr.isedfit_mass
    weight = ageskcorr.spec_weight
    
    xtitle = 'Redshift'
    ytitle1 = textoidl('M_{0.1g} - 5 log (h_{70})')
    ytitle2 = textoidl('log (M_{*}/M'+sunsymbol()+')')

    xrange = [0.03,0.77]
    yrange1 = [-15.0,-23.2]
    yrange2 = [7.9,11.7]

; redshift vs ^{0.1} M_r
    plot, [0], [0], /nodata, ysty=1, xsty=9, xrange=xrange, yrange=yrange1, $
      xtitle='', ytitle=ytitle1, xtickname=replicate(' ',10), $
      position=pos[*,0], yminor=4
    axis, /xaxis, xsty=1, xrange=xrange, xtitle='Lookback Time (Gyr)', $
      xtickv=getredshift(getage(0.0)-timelabel1), $
      xticks=n_elements(timelabel1)-1L, xtickname=string(timelabel1,format='(I0)')
    osu_ages_oplot, z, mg, keynote=keynote
    djs_oplot, mgvz.z, mgvz.csf_Mg, line=0, color=fsc_color('firebrick',101), thick=8
; z vs stellar mass
    plot, [0], [0], /nodata, /noerase, ysty=1, xsty=9, $
      xrange=xrange, yrange=yrange2, xtitle=xtitle, ytitle=ytitle2, $
      position=pos[*,1], yminor=5
    osu_ages_oplot, z, mass, keynote=keynote
    djs_oplot, mgvz.z, alog10(mgvz.csf_mass), line=0, $
      color=fsc_color('firebrick',101), thick=8
    legend, textoidl(['\psi(t)=const']), /right, /bottom, box=0, $
      pspacing=1.4, color=fsc_color('firebrick',100), line=0, $
      charsize=1.5, thick=8, textcolor=!p.color

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif

; ------------------------------------------------------------
; delta-log(O/H) with redshift

    if keyword_set(ps) then psfile = pspath+'dlogohevol'+suffix
    im_plotconfig, 6, pos, psfile=psfile, charsize=2.0, $
      height=4.0*[1,1], ymargin=[0.8,1.0], xmargin=[1.4,0.1], thick=10, $
      keynote=keynote

    ztitle1 = 'Redshift'
    dohtitle1 = textoidl('\Delta'+'log(O/H) from !8L-Z!6')
    dohtitle2 = textoidl('\Delta'+'log(O/H) from !8M-Z!6')
    zrange = [0.0,0.75]
    dohrange = [-0.58,0.13]
    zz = im_array(0.0,1.0,0.02)

    color1 = 'black' & color2 = 'red' & color3 = 'forest green' & color4 = 'royal blue'
    line1 = 5       & line2 = 0      & line3 = 5         & line4 = 0
    psym1 = 6       & psym2 = 15     & psym3 = 4         & psym4 = 14
    symsize1 = 3.0  & symsize2 = 3.0 & symsize3 = 4.0    & symsize4 = 4.0
    if keyword_set(keynote) then color1 = 'white' else color1 = 'black' 

; LZ    
    djs_plot, [0], [0], /nodata, ysty=1, xsty=9, position=pos[*,0], $
      xrange=zrange, yrange=dohrange, xtitle='', ytitle=dohtitle1, $
      xtickname=replicate(' ',10)
    axis, /xaxis, xsty=1, xrange=xrange, xtickv=getredshift(getage(0.0)-timelabel1), $
      xtitle='Lookback Time (Gyr)', xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)')
    djs_oplot, zz, zz*0.0, line=1, thick=3
;   cc = poly_fit([0.1,0.65],[0.0,-alog10(2.0)],1) ; factor of 2 to z=0.65
;   djs_oplot, zz, poly(zz,cc), line=1, thick=3
;   cc = poly_fit([0.1,0.65],[0.0,-alog10(4.0)],1) ; factor of 4
;   djs_oplot, zz, poly(zz,cc), line=1, thick=3
; Q=0
    oploterror, ohevol.z, ohevol.ldlogoh_noevol, ohevol.ldlogoh_noevol_err, $
      psym=-symcat(psym1,thick=8), color=fsc_color(color1,91), errcolor=fsc_color(color1,91), $
      symsize=symsize1, line=line1
    oploterror, ohevol.z, ohevol.ldlogoh_noevol_cor, ohevol.ldlogoh_noevol_cor_err, $
      psym=-symcat(psym2), color=fsc_color(color2,92), errcolor=fsc_color(color2,92), $
      symsize=symsize2, line=line2
; Q=1.5
    oploterror, ohevol.z, ohevol.ldlogoh_levol, ohevol.ldlogoh_levol_err, $
      psym=-symcat(psym3,thick=8), color=fsc_color(color3,93), errcolor=fsc_color(color3,93), $
      symsize=symsize3, line=line3
    oploterror, ohevol.z, ohevol.ldlogoh_levol_cor, ohevol.ldlogoh_levol_cor_err, $
      psym=-symcat(psym4), color=fsc_color(color4,94), errcolor=fsc_color(color4,94), $
      symsize=symsize4, line=line4

    im_legend, ['Q=0, Observed','Q=0, Corrected','Q=1.5, Observed','Q=1.5, Corrected'], $
      /left, /bottom, box=0, charsize=1.8, color=[color1,color2,color3,color4], $
      linestyle=[line1,line2,line3,line4], psym=[psym1,psym2,psym3,psym4], $
      symthick=7, symsize=0.6*[symsize1,symsize2,symsize3,symsize4], textcolor=textcolor

; MZ
    djs_plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, position=pos[*,1], $
      xrange=zrange, yrange=dohrange, xtitle=ztitle1, ytitle=dohtitle2
    djs_oplot, zz, zz*0.0, line=1, thick=3

    oploterror, ohevol.z, ohevol.mdlogoh, ohevol.mdlogoh_err, $
      psym=-symcat(psym1,thick=8), color=fsc_color(color1,91), errcolor=fsc_color(color1,91), $
      symsize=symsize1, line=line1
    oploterror, ohevol.z, ohevol.mdlogoh_levol_cor, ohevol.mdlogoh_levol_cor_err, $
      psym=-symcat(psym2), color=fsc_color(color2,92), errcolor=fsc_color(color2,92), $
      symsize=symsize2, line=line2

    im_legend, ['Observed','Corrected'], $
      /left, /bottom, box=0, charsize=1.8, color=[color1,color2], $
      linestyle=[line1,line2], psym=[psym1,psym2], $
      symthick=7, symsize=0.7*[symsize1,symsize2], textcolor=textcolor

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif

; ------------------------------------------------------------
; luminosity evolution for blue galaxies (literature compilation)
    psfile = pspath+'redshift_vs_mg_lit.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, $
      xmargin=[1.6,0.4], width=6.5, height=6.5, ymargin=[1.0,1.0], $
      keynote=keynote

    h100 = 0.7
    B2g01 = +0.0759 ; ^{0.1}g = B+0.0759+0.0620*[(B-V)-0.5870] [AB, Blanton & Roweis]
    B2r01 = -0.6429 ; ^{0.1}r = B-0.6429-1.0845*[(B-V)-0.5870] [AB, Blanton & Roweis]
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par',/kurucz,/silent))[0]

; luminosity function colors from the literature
    if keyword_set(keynote) then begin
       w06color   = 'orange'       & w06sym   = 16 & w06psize = 2.0
       f06color   = 'dodger blue'  & f06sym   = 14 & f06psize = 3.0
       e07color   = 'firebrick'    & e07sym   = 15 & e07psize = 2.0
       b06color   = 'dodger blue'  & b06sym   = 34 & b06psize = 3.0
    endif else begin
       w06color   = 'orange'       & w06sym   = 16 & w06psize = 2.0
       f06color   = 'navy'         & f06sym   = 14 & f06psize = 3.0
       e07color   = 'firebrick'    & e07sym   = 15 & e07psize = 2.0
       b06color   = 'navy'         & b06sym   = 34 & b06psize = 3.0
    endelse

    xtitle = 'Redshift'
    ytitle = 'M_{0.1r}^{*} for Blue Galaxies'
    xrange = [-0.02,0.85]
    yrange = [-20.2,-22.5]

    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, $
      xsty=9, ysty=1, xrange=xrange, yrange=yrange, position=pos
    axis, /xaxis, xthick=postthick1, xsty=1, xrange=xrange, $
      xtitle='Lookback Time (Gyr)', xtickv=getredshift(getage(0.0)-timelabel1), $
      xticks=n_elements(timelabel1)-1L, xtickname=string(timelabel1,format='(I0)')

; overplot the line for 1.5 mag/z of luminosity evolution    
    djs_oplot, zaxis1, poly(zaxis1-0.1,[-20.95,-1.5]), line=0

; -------------------------
; Eisenstein et al. 2009 [AGES]; Omega_0=0.3, Omega_lamba=0.7, h=1.0;
; AB; alpha fixed at -1.10 for all blue galaxies; evolving color cut
; based on the (u0.1-r0.1) color: A=u0.1-r0.1+0.08(M_r0.1+20);  
    z_e07    = [0.1,0.2,0.3,0.4,0.5,0.65]
    zerr_e07 = [0.05,0.05,0.05,0.05,0.05,0.1]

    mstar_e07 = [-19.95,-20.39,-20.43,-20.57,-20.86,-20.83] + 5.0*alog10(h100) ; h=1-->h=0.7
    mstarerr_e07 = [0.12,0.07,0.04,0.05,0.11,0.08]

; -------------------------
; Willmer et al. 2006 [DEEP2]; "minimal" model adopted, as
; recommended; Omega_0=0.3, Omega_lamba=0.7, h=0.7; Vega; alpha=-1.30
; fixed for all blue galaxies, based on a fit to all galaxies at
; 0.2<z<0.6 in COMBO-17 (Faber et al. 2005); color cut: (U-B) =
; -0.032(M_B+21.52)+0.454-0.25; drop the last bin
    z_w06 = [0.3,0.5,0.7,0.9,1.1];,1.3]
    zerr_w06 = [0.1,0.1,0.1,0.1,0.1];,0.1]

    mstar_w06 = [-20.36,-20.72,-21.15,-21.21,-21.38] + Bvega2ab + B2r01 ; Vega-->AB; B-->r0.1 ; ,-21.86]
    mstarerr_up_w06 = [0.13,0.05,0.07,0.01,0.04];,0.07]
    mstarerr_lo_w06 = [0.11,0.07,0.07,0.03,0.05];,0.08]
    mstarerr_w06 = mstar_w06*0.0
    for i = 0, n_elements(mstar_w06)-1L do mstarerr_w06[i] = mean([mstarerr_up_w06[i],mstarerr_lo_w06[i]])
    
; -------------------------
; Faber et al. 2007 [COMBO17]; Omega_0=0.3, Omega_lamba=0.7, h=0.7;
; Vega; alpha=-1.3 fixed for all blue galaxies;
    z_f06 = [0.3,0.5,0.7,0.9,1.1]
    zerr_f06 = [0.1,0.1,0.1,0.1,0.1]

    mstar_f06 = [-20.74,-21.10,-21.30,-21.10,-21.25] + Bvega2ab + B2r01 ; Vega-->AB; B-->r0.1
    mstarerr_f06 = [0.20,0.15,0.16,0.17,0.18]

    oploterror, z_f06, mstar_f06, zerr_f06, mstarerr_f06, symsize=3.8, $
      psym=symcat(f06sym), color=fsc_color(f06color,1E8), errcolor=fsc_color(f06color,1E8)
    
    oploterror, z_w06, mstar_w06, zerr_w06, mstarerr_up_w06, /hi, symsize=3.5, $
      psym=symcat(w06sym), color=fsc_color(w06color,1E8), errcolor=fsc_color(w06color,1E8)
    oploterror, z_w06, mstar_w06, zerr_w06, mstarerr_lo_w06, /lo, symsize=3.5, $
      psym=symcat(w06sym), color=fsc_color(w06color,1E8), errcolor=fsc_color(w06color,1E8)

    oploterror, z_e07, mstar_e07, zerr_e07, mstarerr_e07, symsize=3.5, $
      psym=symcat(e07sym), color=fsc_color(e07color,1E8), errcolor=fsc_color(e07color,1E8)

; legend
;   label = ['AGES+SDSS','COMBO-17','DEEP2']
    label = ['Willmer et al. 2006','Faber et al. 2007','Eisenstein et al. 2009']
    im_legend, label, /left, /top, box=0, charsize=charsize_6, $
      psym=[w06sym,f06sym,e07sym], symsize=[2.4,2.8,2.0], $
      spacing=2.0, color=[w06color,f06color,e07color], textcolor=textcolor

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif

; ---------------------------------------------------------------------------    
; MZ evolution - AGES

; read the output from FIT_MZLZEVOL
    mzevol_ages = mrdfits(mzpath+'mzevol.fits.gz',1)
    pivotmass = mzevol_ages[0].pivotmass ; the same for all

; 2x3 panel plot
    if keyword_set(ps) then psfile = pspath+'mzevol'+suffix
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.8, $
      xmargin=[1.1,0.2], width=[3.6,3.6], thick=8, keynote=keynote

; zbin0    
    zmin = 0.05 & zmax = 0.15
    aa0 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzages_hogg_scatterplot, aa0.mass, aa0.oh, position=pos[*,0], $
      xstyle=1, ystyle=1, xtitle='', ytitle=ohtitle1, $
      xrange=massrange1, yrange=ohrange1, xtickname=replicate(' ',10), $
      keynote=keynote
    oplot_mzfit, mzevol_ages[0].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle1, linecolor=linecolor1
    im_legend, '0.05<z<0.15', /left, /top, box=0, charsize=1.5, margin=0

; zbin1
    zmin = 0.15 & zmax = 0.25
    aa1 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzages_hogg_scatterplot, aa1.mass, aa1.oh, /noerase, position=pos[*,1], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=massrange1, yrange=ohrange1, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      keynote=keynote
    oplot_mzfit, mzevol_ages[0].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle1, linecolor=linecolor1
    oplot_mzfit, mzevol_ages[1].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle2, linecolor=linecolor2
    im_legend, '0.15<z<0.25', /left, /top, box=0, charsize=1.5, margin=0

; zbin2
    zmin = 0.25 & zmax = 0.35
    aa2 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzages_hogg_scatterplot, aa2.mass, aa2.oh, /noerase, position=pos[*,2], $
      xstyle=1, ystyle=1, xtitle='', ytitle=ohtitle1, $
      xrange=massrange1, yrange=ohrange1, xtickname=replicate(' ',10), $
      keynote=keynote
    oplot_mzfit, mzevol_ages[0].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle1, linecolor=linecolor1
    oplot_mzfit, mzevol_ages[2].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle2, linecolor=linecolor2
    im_legend, '0.25<z<0.35', /left, /top, box=0, charsize=1.5, margin=0

; zbin3
    zmin = 0.35 & zmax = 0.45
    aa3 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzages_hogg_scatterplot, aa3.mass, aa3.oh, /noerase, position=pos[*,3], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=massrange1, yrange=ohrange1, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      keynote=keynote
    oplot_mzfit, mzevol_ages[0].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle1, linecolor=linecolor1
    oplot_mzfit, mzevol_ages[3].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle2, linecolor=linecolor2
    im_legend, '0.35<z<0.45', /left, /top, box=0, charsize=1.5, margin=0

; zbin4
    zmin = 0.45 & zmax = 0.55
    aa4 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzages_hogg_scatterplot, aa4.mass, aa4.oh, /noerase, position=pos[*,4], $
      xstyle=1, ystyle=1, xtitle=masstitle1, ytitle=ohtitle1, $
      xrange=massrange1, yrange=ohrange1, $
      keynote=keynote
    oplot_mzfit, mzevol_ages[0].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle1, linecolor=linecolor1
    oplot_mzfit, mzevol_ages[4].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle2, linecolor=linecolor2
    im_legend, '0.45<z<0.55', /left, /top, box=0, charsize=1.5, margin=0

; zbin5
    zmin = 0.55 & zmax = 0.75
    aa5 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzages_hogg_scatterplot, aa5.mass, aa5.oh, /noerase, position=pos[*,5], $
      xstyle=1, ystyle=1, xtitle=masstitle1, ytitle='', $
      xrange=massrange1, yrange=ohrange1, ytickname=replicate(' ',10), $
      keynote=keynote
    oplot_mzfit, mzevol_ages[0].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle1, linecolor=linecolor1
    oplot_mzfit, mzevol_ages[5].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle2, linecolor=linecolor2
    im_legend, '0.55<z<0.75', /left, /top, box=0, charsize=1.5, margin=0

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif
    
; ---------------------------------------------------------------------------    
; gLZ evolution - AGES

; read the output from FIT_MZLZEVOL
    glzevol_ages = mrdfits(mzpath+'lzevol_g.fits.gz',1)
    pivotmag = glzevol_ages[0].pivotmag ; the same for all

; 2x3 panel plot
    if keyword_set(ps) then psfile = pspath+'glzevol'+suffix
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.8, $
      xmargin=[1.1,0.2], width=[3.6,3.6], thick=8, keynote=keynote

; zbin0    
    zmin = 0.05 & zmax = 0.15
    aa0 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzages_hogg_scatterplot, aa0.mg_ab, aa0.oh, position=pos[*,0], $
      xstyle=1, ystyle=1, xtitle='', ytitle=ohtitle1, $
      xrange=mgrange1, yrange=ohrange1, xtickname=replicate(' ',10), $
      keynote=keynote
    oplot_lzfit, glzevol_ages[0].lzcoeff_local_noevol, pivotmag, band='g', $
      linestyle=linestyle3, linecolor=linecolor3
    im_legend, '0.05<z<0.15', /left, /top, box=0, charsize=1.5, margin=0
    
; zbin1
    zmin = 0.15 & zmax = 0.25
    aa1 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzages_hogg_scatterplot, aa1.mg_ab, aa1.oh, /noerase, position=pos[*,1], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=mgrange1, yrange=ohrange1, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      keynote=keynote
    oplot_lzfit, glzevol_ages[1].lzcoeff_local_noevol, pivotmag, band='g', $
      linestyle=linestyle3, linecolor=linecolor3
;   oplot_lzfit, glzevol_ages[1].lzcoeff_local_levol, pivotmag, band='g', $
;     linestyle=linestyle4, linecolor=linecolor4
    oplot_lzfit, glzevol_ages[1].lzcoeff_levol, pivotmag, band='g', $
      linestyle=linestyle5, linecolor=linecolor5
    im_legend, '0.15<z<0.25', /left, /top, box=0, charsize=1.5, margin=0

; zbin2
    zmin = 0.25 & zmax = 0.35
    aa2 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzages_hogg_scatterplot, aa2.mg_ab, aa2.oh, /noerase, position=pos[*,2], $
      xstyle=1, ystyle=1, xtitle='', ytitle=ohtitle1, $
      xrange=mgrange1, yrange=ohrange1, xtickname=replicate(' ',10), $
      keynote=keynote
    oplot_lzfit, glzevol_ages[2].lzcoeff_local_noevol, pivotmag, band='g', $
      linestyle=linestyle3, linecolor=linecolor3
;   oplot_lzfit, glzevol_ages[2].lzcoeff_local_levol, pivotmag, band='g', $
;     linestyle=linestyle4, linecolor=linecolor4
    oplot_lzfit, glzevol_ages[2].lzcoeff_levol, pivotmag, band='g', $
      linestyle=linestyle5, linecolor=linecolor5
    im_legend, '0.25<z<0.35', /left, /top, box=0, charsize=1.5, margin=0

; zbin3
    zmin = 0.35 & zmax = 0.45
    aa3 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzages_hogg_scatterplot, aa3.mg_ab, aa3.oh, /noerase, position=pos[*,3], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=mgrange1, yrange=ohrange1, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      keynote=keynote
    oplot_lzfit, glzevol_ages[3].lzcoeff_local_noevol, pivotmag, band='g', $
      linestyle=linestyle3, linecolor=linecolor3
;   oplot_lzfit, glzevol_ages[3].lzcoeff_local_levol, pivotmag, band='g', $
;     linestyle=linestyle4, linecolor=linecolor4
    oplot_lzfit, glzevol_ages[3].lzcoeff_levol, pivotmag, band='g', $
      linestyle=linestyle5, linecolor=linecolor5
    im_legend, '0.35<z<0.45', /left, /top, box=0, charsize=1.5, margin=0

; zbin4
    zmin = 0.45 & zmax = 0.55
    aa4 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzages_hogg_scatterplot, aa4.mg_ab, aa4.oh, /noerase, position=pos[*,4], $
      xstyle=1, ystyle=1, xtitle=mgtitle1, ytitle=ohtitle1, $
      xrange=mgrange1, yrange=ohrange1, $
      keynote=keynote
    oplot_lzfit, glzevol_ages[4].lzcoeff_local_noevol, pivotmag, band='g', $
      linestyle=linestyle3, linecolor=linecolor3
;   oplot_lzfit, glzevol_ages[4].lzcoeff_local_levol, pivotmag, band='g', $
;     linestyle=linestyle4, linecolor=linecolor4
    oplot_lzfit, glzevol_ages[4].lzcoeff_levol, pivotmag, band='g', $
      linestyle=linestyle5, linecolor=linecolor5
    im_legend, '0.45<z<0.55', /left, /top, box=0, charsize=1.5, margin=0

; zbin5
    zmin = 0.55 & zmax = 0.75
    aa5 = mzlz_grab_info(agesohdust,ageskcorr,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzages_hogg_scatterplot, aa5.mg_ab, aa5.oh, /noerase, position=pos[*,5], $
      xstyle=1, ystyle=1, xtitle=mgtitle1, ytitle='', $
      xrange=mgrange1, yrange=ohrange1, ytickname=replicate(' ',10), $
      keynote=keynote
    oplot_lzfit, glzevol_ages[5].lzcoeff_local_noevol, pivotmag, band='g', $
      linestyle=linestyle3, linecolor=linecolor3
;   oplot_lzfit, glzevol_ages[5].lzcoeff_local_levol, pivotmag, band='g', $
;     linestyle=linestyle4, linecolor=linecolor4
    oplot_lzfit, glzevol_ages[5].lzcoeff_levol, pivotmag, band='g', $
      linestyle=linestyle5, linecolor=linecolor5
    im_legend, '0.55<z<0.75', /left, /top, box=0, charsize=1.5, margin=0

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif

; ---------------------------------------------------------------------------    
; MZ evolution - SDSS2AGES/levol

    linecolor2 = 'firebrick'

; read the output from FIT_MZLZEVOL
    mzevol_levol = mrdfits(mzpath+'mzevol.fits.gz',3)
    pivotmass = mzevol_levol[0].pivotmass ; the same for all

; 2x3 panel plot
    if keyword_set(ps) then psfile = pspath+'mzevol_levol'+suffix
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.8, $
      xmargin=[1.1,0.2], width=[3.6,3.6], thick=8, keynote=keynote

; zbin0    
    zmin = 0.05 & zmax = 0.15
    ss0 = mzlz_grab_info(sdssages_levol0,sdssages_levol0,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzsdss_hogg_scatterplot, ss0.mass, ss0.oh, position=pos[*,0], $
      xstyle=1, ystyle=1, xtitle='', ytitle=ohtitle1, $
      xrange=massrange1, yrange=ohrange1, xtickname=replicate(' ',10), $
      keynote=keynote
    oplot_mzfit, mzevol_levol[0].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle1, linecolor=linecolor1
    im_legend, '0.05<z<0.15', /left, /top, box=0, charsize=1.5, margin=0

; zbin1
    zmin = 0.15 & zmax = 0.25
    ss1 = mzlz_grab_info(sdssages_levol1,sdssages_levol1,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzsdss_hogg_scatterplot, ss1.mass, ss1.oh, /noerase, position=pos[*,1], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=massrange1, yrange=ohrange1, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      keynote=keynote
    oplot_mzfit, mzevol_levol[0].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle1, linecolor=linecolor1
    oplot_mzfit, mzevol_levol[1].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle2, linecolor=linecolor2
    im_legend, '0.15<z<0.25', /left, /top, box=0, charsize=1.5, margin=0

; zbin2
    zmin = 0.25 & zmax = 0.35
    ss2 = mzlz_grab_info(sdssages_levol2,sdssages_levol2,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzsdss_hogg_scatterplot, ss2.mass, ss2.oh, /noerase, position=pos[*,2], $
      xstyle=1, ystyle=1, xtitle='', ytitle=ohtitle1, $
      xrange=massrange1, yrange=ohrange1, xtickname=replicate(' ',10), $
      keynote=keynote
    oplot_mzfit, mzevol_levol[0].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle1, linecolor=linecolor1
    oplot_mzfit, mzevol_levol[2].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle2, linecolor=linecolor2
    im_legend, '0.25<z<0.35', /left, /top, box=0, charsize=1.5, margin=0

; zbin3
    zmin = 0.35 & zmax = 0.45
    ss3 = mzlz_grab_info(sdssages_levol3,sdssages_levol3,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzsdss_hogg_scatterplot, ss3.mass, ss3.oh, /noerase, position=pos[*,3], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=massrange1, yrange=ohrange1, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      keynote=keynote
    oplot_mzfit, mzevol_levol[0].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle1, linecolor=linecolor1
    oplot_mzfit, mzevol_levol[3].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle2, linecolor=linecolor2
    im_legend, '0.35<z<0.45', /left, /top, box=0, charsize=1.5, margin=0

; zbin4
    zmin = 0.45 & zmax = 0.55
    ss4 = mzlz_grab_info(sdssages_levol4,sdssages_levol4,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzsdss_hogg_scatterplot, ss4.mass, ss4.oh, /noerase, position=pos[*,4], $
      xstyle=1, ystyle=1, xtitle=masstitle1, ytitle=ohtitle1, $
      xrange=massrange1, yrange=ohrange1, $
      keynote=keynote
    oplot_mzfit, mzevol_levol[0].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle1, linecolor=linecolor1
    oplot_mzfit, mzevol_levol[4].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle2, linecolor=linecolor2
    im_legend, '0.45<z<0.55', /left, /top, box=0, charsize=1.5, margin=0

; zbin5
    zmin = 0.55 & zmax = 0.75
    ss5 = mzlz_grab_info(sdssages_levol5,sdssages_levol5,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzsdss_hogg_scatterplot, ss5.mass, ss5.oh, /noerase, position=pos[*,5], $
      xstyle=1, ystyle=1, xtitle=masstitle1, ytitle='', $
      xrange=massrange1, yrange=ohrange1, ytickname=replicate(' ',10), $
      keynote=keynote
    oplot_mzfit, mzevol_levol[0].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle1, linecolor=linecolor1
    oplot_mzfit, mzevol_levol[5].mzcoeff, pivotmass, /nopsym, $
      linestyle=linestyle2, linecolor=linecolor2
    im_legend, '0.55<z<0.75', /left, /top, box=0, charsize=1.5, margin=0

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif

; ---------------------------------------------------------------------------    
; gLZ evolution - SDSS2AGES/levol

    linecolor5 = 'firebrick'

; read the output from FIT_MZLZEVOL
    glzevol_levol = mrdfits(mzpath+'lzevol_g.fits.gz',3)
    pivotmag = glzevol_levol[0].pivotmag ; the same for all

; 2x3 panel plot
    if keyword_set(ps) then psfile = pspath+'glzevol_levol'+suffix
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.8, $
      xmargin=[1.1,0.2], width=[3.6,3.6], thick=8, keynote=keynote

; zbin0
    zmin = 0.05 & zmax = 0.15
    ss0 = mzlz_grab_info(sdssages_levol0,sdssages_levol0,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzsdss_hogg_scatterplot, ss0.mg_ab, ss0.oh, position=pos[*,0], $
      xstyle=1, ystyle=1, xtitle='', ytitle=ohtitle1, $
      xrange=mgrange1, yrange=ohrange1, xtickname=replicate(' ',10), $
      keynote=keynote
    oplot_lzfit, glzevol_levol[0].lzcoeff_local_noevol, pivotmag, band='g', $
      linestyle=linestyle3, linecolor=linecolor3
    im_legend, '0.05<z<0.15', /left, /top, box=0, charsize=1.5, margin=0
    
; zbin1
    zmin = 0.15 & zmax = 0.25
    ss1 = mzlz_grab_info(sdssages_levol1,sdssages_levol1,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzsdss_hogg_scatterplot, ss1.mg_ab, ss1.oh, /noerase, position=pos[*,1], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=mgrange1, yrange=ohrange1, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      keynote=keynote
    oplot_lzfit, glzevol_levol[1].lzcoeff_local_noevol, pivotmag, band='g', $
      linestyle=linestyle3, linecolor=linecolor3
;   oplot_lzfit, glzevol_levol[1].lzcoeff_local_levol, pivotmag, band='g', $
;     linestyle=linestyle4, linecolor=linecolor4
    oplot_lzfit, glzevol_levol[1].lzcoeff_levol, pivotmag, band='g', $
      linestyle=linestyle5, linecolor=linecolor5
    im_legend, '0.15<z<0.25', /left, /top, box=0, charsize=1.5, margin=0

; zbin2
    zmin = 0.25 & zmax = 0.35
    ss2 = mzlz_grab_info(sdssages_levol2,sdssages_levol2,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzsdss_hogg_scatterplot, ss2.mg_ab, ss2.oh, /noerase, position=pos[*,2], $
      xstyle=1, ystyle=1, xtitle='', ytitle=ohtitle1, $
      xrange=mgrange1, yrange=ohrange1, xtickname=replicate(' ',10), $
      keynote=keynote
    oplot_lzfit, glzevol_levol[2].lzcoeff_local_noevol, pivotmag, band='g', $
      linestyle=linestyle3, linecolor=linecolor3
;   oplot_lzfit, glzevol_levol[2].lzcoeff_local_levol, pivotmag, band='g', $
;     linestyle=linestyle4, linecolor=linecolor4
    oplot_lzfit, glzevol_levol[2].lzcoeff_levol, pivotmag, band='g', $
      linestyle=linestyle5, linecolor=linecolor5
    im_legend, '0.25<z<0.35', /left, /top, box=0, charsize=1.5, margin=0

; zbin3
    zmin = 0.35 & zmax = 0.45
    ss3 = mzlz_grab_info(sdssages_levol3,sdssages_levol3,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzsdss_hogg_scatterplot, ss3.mg_ab, ss3.oh, /noerase, position=pos[*,3], $
      xstyle=1, ystyle=1, xtitle='', ytitle='', $
      xrange=mgrange1, yrange=ohrange1, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      keynote=keynote, /hackit
    oplot_lzfit, glzevol_levol[3].lzcoeff_local_noevol, pivotmag, band='g', $
      linestyle=linestyle3, linecolor=linecolor3
;   oplot_lzfit, glzevol_levol[3].lzcoeff_local_levol, pivotmag, band='g', $
;     linestyle=linestyle4, linecolor=linecolor4
    oplot_lzfit, glzevol_levol[3].lzcoeff_levol, pivotmag, band='g', $
      linestyle=linestyle5, linecolor=linecolor5
    im_legend, '0.35<z<0.45', /left, /top, box=0, charsize=1.5, margin=0

; zbin4
    zmin = 0.45 & zmax = 0.55
    ss4 = mzlz_grab_info(sdssages_levol4,sdssages_levol4,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzsdss_hogg_scatterplot, ss4.mg_ab, ss4.oh, /noerase, position=pos[*,4], $
      xstyle=1, ystyle=1, xtitle=mgtitle1, ytitle=ohtitle1, $
      xrange=mgrange1, yrange=ohrange1, $
      keynote=keynote
    oplot_lzfit, glzevol_levol[4].lzcoeff_local_noevol, pivotmag, band='g', $
      linestyle=linestyle3, linecolor=linecolor3
;   oplot_lzfit, glzevol_levol[4].lzcoeff_local_levol, pivotmag, band='g', $
;     linestyle=linestyle4, linecolor=linecolor4
    oplot_lzfit, glzevol_levol[4].lzcoeff_levol, pivotmag, band='g', $
      linestyle=linestyle5, linecolor=linecolor5
    im_legend, '0.45<z<0.55', /left, /top, box=0, charsize=1.5, margin=0

; zbin5
    zmin = 0.55 & zmax = 0.75
    ss5 = mzlz_grab_info(sdssages_levol5,sdssages_levol5,/ewunity,zmin=zmin,zmax=zmax)
    osu_mzsdss_hogg_scatterplot, ss5.mg_ab, ss5.oh, /noerase, position=pos[*,5], $
      xstyle=1, ystyle=1, xtitle=mgtitle1, ytitle='', $
      xrange=mgrange1, yrange=ohrange1, ytickname=replicate(' ',10), $
      keynote=keynote
    oplot_lzfit, glzevol_levol[5].lzcoeff_local_noevol, pivotmag, band='g', $
      linestyle=linestyle3, linecolor=linecolor3
;   oplot_lzfit, glzevol_levol[5].lzcoeff_local_levol, pivotmag, band='g', $
;     linestyle=linestyle4, linecolor=linecolor4
    oplot_lzfit, glzevol_levol[5].lzcoeff_levol, pivotmag, band='g', $
      linestyle=linestyle5, linecolor=linecolor5
    im_legend, '0.55<z<0.75', /left, /top, box=0, charsize=1.5, margin=0

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif
    
return
end

