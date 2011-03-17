pro oplot_lzevol, lzevol, absmagaxis=absmagaxis, ohrange=ohrange, $
  plot_evolved=plot_evolved, plot_local_levol=plot_local_levol, $
  thick=thick
; overplot the evolving LZ relation semi-intelligently
    
    if (n_elements(thick) eq 0L) then thick = 5.0

    ohoffset = 0.09
; local relation
    inrange = where((poly(absmagaxis-lzevol.pivot,lzevol.lzcoeff_local) gt ohrange[0]+ohoffset) and $
      (poly(absmagaxis-lzevol.pivot,lzevol.lzcoeff_local) lt ohrange[1]-ohoffset))
    djs_oplot, absmagaxis[inrange], (poly(absmagaxis-lzevol.pivot,lzevol.lzcoeff_local))[inrange], $
      thick=thick, line=0;, color=fsc_color('',1E8)
; shifted local relation
    if keyword_set(plot_evolved) then begin
       inrange = where((poly(absmagaxis-lzevol.pivot,lzevol.lzcoeff) gt ohrange[0]+ohoffset) and $
         (poly(absmagaxis-lzevol.pivot,lzevol.lzcoeff) lt ohrange[1]-ohoffset))
       djs_oplot, absmagaxis[inrange], (poly(absmagaxis-lzevol.pivot,lzevol.lzcoeff))[inrange], $
         thick=thick, line=5, color=fsc_color('blue',1E8)
    endif
; evolved local LZ relation
    if keyword_set(plot_local_levol) then begin
       inrange = where((poly(absmagaxis-lzevol.pivot,lzevol.lzcoeff_levol) gt ohrange[0]+ohoffset) and $
         (poly(absmagaxis-lzevol.pivot,lzevol.lzcoeff_levol) lt ohrange[1]-ohoffset))
       djs_oplot, absmagaxis[inrange], (poly(absmagaxis-lzevol.pivot,lzevol.lzcoeff_levol))[inrange], $
         thick=thick, line=3, color=fsc_color('royal blue',1E8)
    endif

return
end
    
pro oplot_mzevol, mzevol, massaxis=massaxis, plot_evolved=plot_evolved, thick=thick
; overplot the evolving MZ relation semi-intelligently

    if (n_elements(thick) eq 0L) then thick = 5.0
    
    djs_oplot, massaxis, poly(massaxis-mzevol.pivot,mzevol.mzcoeff_local), $
      thick=thick, line=0;, color=fsc_color(color,1E8)
; shifted local relation
    if keyword_set(plot_evolved) then begin
       djs_oplot, massaxis, poly(massaxis-mzevol.pivot,mzevol.mzcoeff), $
         thick=thick, line=5, color=fsc_color('blue',1E8)
    endif

return
end
    
pro talk_08nov_upenn, postscript=postscript, pdf=pdf
; jm08nov03nyu - plots for my colloquium talk

    common ages_common, ages, ancillary

; default plotting variables - note that POSTSCRIPT over-rides PDF

    im_plotfaves, pdf=pdf, postscript=postscript

    postthick1 = 2.0
    postthick2 = 2.0
    postthick3 = 2.0
    postthick4 = 2.0
    charsize1 = 2.0
    charsize2 = 2.0
    charsize3 = 1.5
    charsize4 = 1.8
    textcolor1 = 'white'
    axis_color = 'white'
    xcolortable = 0

    mzpath = ages_path(/projects)+'mz/'
    pspath = getenv('PAPERSPATH')+'/literature/'
    datapath2 = pspath+'data/'
    datapath1 = '~/home/research/talks/2008/08nov_upenn/'

    if keyword_set(pdf) then begin
       postscript = 0L
       pspath = datapath1+'keynote/' ; for keynote presentations
       textcolor1 = 'white'
       axis_color = 'black'
       xcolortable = 3

       charsize1 = 2.2
       charsize2 = 1.6
       charsize3 = 1.2
       charsize4 = 1.4
       postthick1 = 8.0
       postthick2 = 6.0
       postthick3 = 4.0
       postthick4 = 10.0
       erbcolor = 'forest green'
       liucolor_loz = 'orange'
       liucolor_hiz = 'red'
       shapleycolor_loz = 'orange'
       shapleycolor_hiz = 'red'
       maiercolor = liucolor_hiz ; same color
       linethick1 = 10.0
       linethick2 = 8.0
       symthick1 = 8.0
    endif
    if keyword_set(postscript) then begin
       pspath = datapath1
       textcolor1 = 'black'
       axis_color = 'black'
       xcolortable = 3

       postthick1 = 4.0
       postthick2 = 3.0
       postthick3 = 3.0
       postthick4 = 8.0
       linethick1 = 8.0
       linethick2 = 6.0
       symthick1 = 6.0
    endif

    loadct, xcolortable, /silent

    if keyword_set(pdf) then begin
       agespsize1 = 0.4 & agessym1 = 16 & agescolor1 = 'grey'
       agespsize2 = 0.5 & agessym2 = 16 & agescolor2 = 'cyan'
       sdsspsize1 = 0.2 & sdsssym1 = 16 & sdsscolor1 = 'grey'
    endif else begin
       agespsize1 = 0.4 & agessym1 = 16 & agescolor1 = 'grey'
       agespsize2 = 0.4 & agessym2 = 16 & agescolor2 = 'grey'
       sdsspsize1 = 0.2 & sdsssym1 = 16 & sdsscolor1 = 'grey'
    endelse

    if keyword_set(pdf) then begin
       dlogohsym = 9     & dlogohpsize = 4.2    & dlogohcolor = 'cyan'
       dlogohcorsym = 15 & dlogohcorpsize = 3.5 & dlogohcorcolor = 'red'
    endif else begin
       dlogohsym = 4     & dlogohpsize = 4.2    & dlogohcolor = 'dark green'
       dlogohcorsym = 15 & dlogohcorpsize = 3.2 & dlogohcorcolor = 'red'
    endelse
       
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

    timelabel1 = [1.0,3.0,5.0,7.0] ; [Gyr]
    zaxis1 = findgen(((0.8)-(0.0))/0.01+1)*0.01+(0.0)
    ohaxis = findgen((10.0-5.0)/0.01+1)*0.01+5.0

    minmass1 = 9.4
    maxmass1 = 11.9 ; 11.2

    bigmgaxis1 = findgen(((-10.0)-(-30))/0.01+1)*0.01+(-30.0)
    bigmraxis1 = findgen(((-10.0)-(-30))/0.01+1)*0.01+(-30.0)
;   bigmassaxis1 = findgen((maxmass1-minmass1)/0.01+1)*0.01+minmass1
    bigmassaxis1 = findgen((maxmass1-8.9)/0.01+1)*0.01+8.9
;   bigmassaxis1 = findgen((11.2-9.0)/0.01+1)*0.01+9.0
    
    mrtitle1 = textoidl('M_{0.1r}')
    mgtitle1 = textoidl('M_{0.1g}')
    masstitle1 = textoidl('log (M_{*}/M'+sunsymbol()+')')
    ohtitle1 = textoidl('12 + log (O/H)')
;   ohtitle1 = textoidl('12 + log (O/H)_{EW}')
    ohtitle2 = textoidl('12 + log (O/H)')

    light = 2.99792458D10       ; speed of light [cm/s]

; read some of the data we're going to need

    mz_ancillary = read_ages_mz_sample(/mzhiiplus_ancillary)
    mz_ispec = read_ages_mz_sample(/mzhiiplus_ispec)
    
    ageskcorr = read_ages_mz_sample(/mzhiiplus_ancillary)
    agesohdust = read_ages_mz_sample(/mzhiiplus_log12oh)

; luminosity function colors from the literature

    if keyword_set(pdf) then begin
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

; some constants

    h100 = 0.7
    mrstar = -20.44 + 5.0*alog10(h100) ; from Blanton+03
    mstar = 10.648 - alog10(0.7)       ; from Baldry+08 (diet Salpeter --> Salpeter)

    mgsun = (k_solar_magnitudes(filterlist='sdss_g0.par',/silent))[0]
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par',/kurucz,/silent))[0]
    B2r01 = -0.6429 ; ^{0.1}r = B-0.6429-1.0845*[(B-V)-0.5870] [AB, Blanton & Roweis]
    B2g01 = +0.0759 ; ^{0.1}g = B+0.0759+0.0620*[(B-V)-0.5870] [AB, Blanton & Roweis]

; read the LZ/MZ files (see MZPLOTS)

    lzevol_mr_file = mzpath+'ages_lzevol_mr.fits'
    lzevol_mg_file = mzpath+'ages_lzevol_mg.fits'
    mzevol_file = mzpath+'ages_mzevol.fits'

    lzevol_mr_noevol_file = mzpath+'ages_lzevol_mr_noevol.fits'
    lzevol_mg_noevol_file = mzpath+'ages_lzevol_mg_noevol.fits'
    mzevol_noevol_file = mzpath+'ages_mzevol_noevol.fits'

    lzevol_mr_levol_file = mzpath+'ages_lzevol_mr_levol.fits'
    lzevol_mg_levol_file = mzpath+'ages_lzevol_mg_levol.fits'
    mzevol_levol_file = mzpath+'ages_mzevol_levol.fits'
    
    splog, 'Reading '+lzevol_mr_file
    lzevol_mr = mrdfits(lzevol_mr_file,1,/silent)
    splog, 'Reading '+lzevol_mr_noevol_file
    lzevol_mr_noevol = mrdfits(lzevol_mr_noevol_file,1,/silent)
    splog, 'Reading '+lzevol_mr_levol_file
    lzevol_mr_levol = mrdfits(lzevol_mr_levol_file,1,/silent)

    splog, 'Reading '+lzevol_mg_file
    lzevol_mg = mrdfits(lzevol_mg_file,1,/silent)
    splog, 'Reading '+lzevol_mg_noevol_file
    lzevol_mg_noevol = mrdfits(lzevol_mg_noevol_file,1,/silent)
    splog, 'Reading '+lzevol_mg_levol_file
    lzevol_mg_levol = mrdfits(lzevol_mg_levol_file,1,/silent)

    splog, 'Reading '+mzevol_file
    mzevol = mrdfits(mzevol_file,1,/silent)
    splog, 'Reading '+mzevol_noevol_file
    mzevol_noevol = mrdfits(mzevol_noevol_file,1,/silent)
    splog, 'Reading '+mzevol_levol_file
    mzevol_levol = mrdfits(mzevol_levol_file,1,/silent)

; read some representative Pegase models

    splog, 'Reading the Pegase models'
    pegpath = getenv('PEGASE_HR_SFHGRID_DIR')+'/MEASURE/'
    peg1 = mrdfits(pegpath+'salp_tau_001.0Gyr.info.fits',1,silent=0)
    peg2 = mrdfits(pegpath+'salp_tau_001.0Gyr.info.fits',1,silent=0)
    peg3 = mrdfits(pegpath+'salp_tau_003.0Gyr.info.fits',1,silent=0)
    pegconst = mrdfits(pegpath+'salp_tau_999.0Gyr.info.fits',1,silent=0)

    rev = reverse(sort(peg1.age)) ; reverse the time array!
    peg1 = peg1[rev]
    peg2 = peg2[rev]
    peg3 = peg3[rev]
    pegconst = pegconst[rev]
    
    peg_mgalaxy = 10D10
    pegzform = 1.5
    pegtform = getage(pegzform)
    peg_zaxis = getredshift(peg1.age/1D3+pegtform)
    peg_lookback = getage(0.0)-getage(peg_zaxis)
    
    peg_good = where(peg_zaxis gt 0.0)
    peg_zaxis = peg_zaxis[peg_good]
    peg_lookback = peg_lookback[peg_good]
    peg1 = peg1[peg_good]
    peg2 = peg2[peg_good]
    peg3 = peg3[peg_good]
    pegconst = pegconst[peg_good]

; ---------------------------------------------------------------------------
; MZ: AGES, SDSS-AGES (noevol), SDSS-AGES (levol), zbins
; ---------------------------------------------------------------------------

    xpage = 8.5 & ypage = 9.5
    pagemaker, nx=3, ny=4, xspace=0, yspace=0.0, width=2.4*[1,1,1], $
      height=2.0*[1,1,1,1], xmargin=[1.0,0.3], ymargin=[0.4,1.1], $
      xpage=xpage, ypage=ypage, position=pos, /normal

; -------------------------
; define the data    
; -------------------------

; AGES    
    
; 0.01<z<0.20

    indx_zbin1 = where((agesohdust.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(agesohdust.r23branch_ew_alpha_gr_kk04,2) eq 'U') and $
      (ageskcorr.z gt 0.01) and (ageskcorr.z lt 0.20),nindx_zbin1)

    mr_zbin1 = ageskcorr[indx_zbin1].ugriz_absmag[2]
    mg_zbin1 = ageskcorr[indx_zbin1].ugriz_absmag[1]
    mass_zbin1 = ageskcorr[indx_zbin1].mass + im_convert_imf(/from_chabrier)
    oh_zbin1 = agesohdust[indx_zbin1].zstrong_ew_alpha_gr_12oh_kk04
    oherr_zbin1 = agesohdust[indx_zbin1].zstrong_ew_alpha_gr_12oh_kk04_err
    weight_zbin1 = ageskcorr[indx_zbin1].spec_weight

; 0.20<z<0.40

    indx_zbin2 = where((agesohdust.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(agesohdust.r23branch_ew_alpha_gr_kk04,2) eq 'U') and $
      (ageskcorr.z gt 0.20) and (ageskcorr.z lt 0.40),nindx_zbin2)

    mr_zbin2 = ageskcorr[indx_zbin2].ugriz_absmag[2]
    mg_zbin2 = ageskcorr[indx_zbin2].ugriz_absmag[1]
    mass_zbin2 = ageskcorr[indx_zbin2].mass + im_convert_imf(/from_chabrier)
    oh_zbin2 = agesohdust[indx_zbin2].zstrong_ew_alpha_gr_12oh_kk04
    oherr_zbin2 = agesohdust[indx_zbin2].zstrong_ew_alpha_gr_12oh_kk04_err
    weight_zbin2 = ageskcorr[indx_zbin2].spec_weight

; 0.40<z<0.60

    indx_zbin3 = where((agesohdust.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(agesohdust.r23branch_ew_alpha_gr_kk04,2) eq 'U') and $
      (ageskcorr.z gt 0.40) and (ageskcorr.z lt 0.60),nindx_zbin3)

    mr_zbin3 = ageskcorr[indx_zbin3].ugriz_absmag[2]
    mg_zbin3 = ageskcorr[indx_zbin3].ugriz_absmag[1]
    mass_zbin3 = ageskcorr[indx_zbin3].mass + im_convert_imf(/from_chabrier)
    oh_zbin3 = agesohdust[indx_zbin3].zstrong_ew_alpha_gr_12oh_kk04
    oherr_zbin3 = agesohdust[indx_zbin3].zstrong_ew_alpha_gr_12oh_kk04_err
    weight_zbin3 = ageskcorr[indx_zbin3].spec_weight

; 0.60<z<0.80

    indx_zbin4 = where((agesohdust.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(agesohdust.r23branch_ew_alpha_gr_kk04,2) eq 'U') and $
      (ageskcorr.z gt 0.60) and (ageskcorr.z lt 0.80),nindx_zbin4)

    mr_zbin4 = ageskcorr[indx_zbin4].ugriz_absmag[2]
    mg_zbin4 = ageskcorr[indx_zbin4].ugriz_absmag[1]
    mass_zbin4 = ageskcorr[indx_zbin4].mass + im_convert_imf(/from_chabrier)
    oh_zbin4 = agesohdust[indx_zbin4].zstrong_ew_alpha_gr_12oh_kk04
    oherr_zbin4 = agesohdust[indx_zbin4].zstrong_ew_alpha_gr_12oh_kk04_err
    weight_zbin4 = ageskcorr[indx_zbin4].spec_weight

; SDSS/AGES - no luminosity evolution

    sdssagesohdust_noevol_zbin1 = read_sdssages_mz_sample(/ohdust_zbin1,evolve=0)
    sdssagesohdust_noevol_zbin2 = read_sdssages_mz_sample(/ohdust_zbin2,evolve=0)
    sdssagesohdust_noevol_zbin3 = read_sdssages_mz_sample(/ohdust_zbin3,evolve=0)
    sdssagesohdust_noevol_zbin4 = read_sdssages_mz_sample(/ohdust_zbin4,evolve=0)

; 0.01<z<0.20

    indx_noevol_zbin1 = where((sdssagesohdust_noevol_zbin1.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(sdssagesohdust_noevol_zbin1.r23branch_ew_alpha_gr_kk04,2) eq 'U'),nindx_noevol_zbin1)

    mr_noevol_zbin1 = sdssagesohdust_noevol_zbin1[indx_noevol_zbin1].ugriz_absmag[2]
    mg_noevol_zbin1 = sdssagesohdust_noevol_zbin1[indx_noevol_zbin1].ugriz_absmag[1]
    mass_noevol_zbin1 = sdssagesohdust_noevol_zbin1[indx_noevol_zbin1].mass + im_convert_imf(/from_chabrier)
    oh_noevol_zbin1 = sdssagesohdust_noevol_zbin1[indx_noevol_zbin1].zstrong_ew_alpha_gr_12oh_kk04
    oherr_noevol_zbin1 = sdssagesohdust_noevol_zbin1[indx_noevol_zbin1].zstrong_ew_alpha_gr_12oh_kk04_err

; 0.20<z<0.40

    indx_noevol_zbin2 = where((sdssagesohdust_noevol_zbin2.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(sdssagesohdust_noevol_zbin2.r23branch_ew_alpha_gr_kk04,2) eq 'U'),nindx_noevol_zbin2)

    mr_noevol_zbin2 = sdssagesohdust_noevol_zbin2[indx_noevol_zbin2].ugriz_absmag[2]
    mg_noevol_zbin2 = sdssagesohdust_noevol_zbin2[indx_noevol_zbin2].ugriz_absmag[1]
    mass_noevol_zbin2 = sdssagesohdust_noevol_zbin2[indx_noevol_zbin2].mass + im_convert_imf(/from_chabrier)
    oh_noevol_zbin2 = sdssagesohdust_noevol_zbin2[indx_noevol_zbin2].zstrong_ew_alpha_gr_12oh_kk04
    oherr_noevol_zbin2 = sdssagesohdust_noevol_zbin2[indx_noevol_zbin2].zstrong_ew_alpha_gr_12oh_kk04_err

; 0.40<z<0.60

    indx_noevol_zbin3 = where((sdssagesohdust_noevol_zbin3.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(sdssagesohdust_noevol_zbin3.r23branch_ew_alpha_gr_kk04,2) eq 'U'),nindx_noevol_zbin3)

    mr_noevol_zbin3 = sdssagesohdust_noevol_zbin3[indx_noevol_zbin3].ugriz_absmag[2]
    mg_noevol_zbin3 = sdssagesohdust_noevol_zbin3[indx_noevol_zbin3].ugriz_absmag[1]
    mass_noevol_zbin3 = sdssagesohdust_noevol_zbin3[indx_noevol_zbin3].mass + im_convert_imf(/from_chabrier)
    oh_noevol_zbin3 = sdssagesohdust_noevol_zbin3[indx_noevol_zbin3].zstrong_ew_alpha_gr_12oh_kk04
    oherr_noevol_zbin3 = sdssagesohdust_noevol_zbin3[indx_noevol_zbin3].zstrong_ew_alpha_gr_12oh_kk04_err
    
; 0.60<z<0.80

    indx_noevol_zbin4 = where((sdssagesohdust_noevol_zbin4.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(sdssagesohdust_noevol_zbin4.r23branch_ew_alpha_gr_kk04,2) eq 'U'),nindx_noevol_zbin4)

    mr_noevol_zbin4 = sdssagesohdust_noevol_zbin4[indx_noevol_zbin4].ugriz_absmag[2]
    mg_noevol_zbin4 = sdssagesohdust_noevol_zbin4[indx_noevol_zbin4].ugriz_absmag[1]
    mass_noevol_zbin4 = sdssagesohdust_noevol_zbin4[indx_noevol_zbin4].mass + im_convert_imf(/from_chabrier)
    oh_noevol_zbin4 = sdssagesohdust_noevol_zbin4[indx_noevol_zbin4].zstrong_ew_alpha_gr_12oh_kk04
    oherr_noevol_zbin4 = sdssagesohdust_noevol_zbin4[indx_noevol_zbin4].zstrong_ew_alpha_gr_12oh_kk04_err

; SDSS/AGES - with luminosity evolution
    
    sdssagesohdust_levol_zbin1 = read_sdssages_mz_sample(/ohdust_zbin1,evolve=1)
    sdssagesohdust_levol_zbin2 = read_sdssages_mz_sample(/ohdust_zbin2,evolve=1)
    sdssagesohdust_levol_zbin3 = read_sdssages_mz_sample(/ohdust_zbin3,evolve=1)
    sdssagesohdust_levol_zbin4 = read_sdssages_mz_sample(/ohdust_zbin4,evolve=1)

; 0.01<z<0.20

    indx_levol_zbin1 = where((sdssagesohdust_levol_zbin1.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(sdssagesohdust_levol_zbin1.r23branch_ew_alpha_gr_kk04,2) eq 'U'),nindx_levol_zbin1)

    mr_levol_zbin1 = sdssagesohdust_levol_zbin1[indx_levol_zbin1].ugriz_absmag[2]
    mg_levol_zbin1 = sdssagesohdust_levol_zbin1[indx_levol_zbin1].ugriz_absmag[1]
    mass_levol_zbin1 = sdssagesohdust_levol_zbin1[indx_levol_zbin1].mass + im_convert_imf(/from_chabrier)
    oh_levol_zbin1 = sdssagesohdust_levol_zbin1[indx_levol_zbin1].zstrong_ew_alpha_gr_12oh_kk04
    oherr_levol_zbin1 = sdssagesohdust_levol_zbin1[indx_levol_zbin1].zstrong_ew_alpha_gr_12oh_kk04_err

; 0.20<z<0.40

    indx_levol_zbin2 = where((sdssagesohdust_levol_zbin2.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(sdssagesohdust_levol_zbin2.r23branch_ew_alpha_gr_kk04,2) eq 'U'),nindx_levol_zbin2)

    mr_levol_zbin2 = sdssagesohdust_levol_zbin2[indx_levol_zbin2].ugriz_absmag[2]
    mg_levol_zbin2 = sdssagesohdust_levol_zbin2[indx_levol_zbin2].ugriz_absmag[1]
    mass_levol_zbin2 = sdssagesohdust_levol_zbin2[indx_levol_zbin2].mass + im_convert_imf(/from_chabrier)
    oh_levol_zbin2 = sdssagesohdust_levol_zbin2[indx_levol_zbin2].zstrong_ew_alpha_gr_12oh_kk04
    oherr_levol_zbin2 = sdssagesohdust_levol_zbin2[indx_levol_zbin2].zstrong_ew_alpha_gr_12oh_kk04_err

; 0.40<z<0.60

    indx_levol_zbin3 = where((sdssagesohdust_levol_zbin3.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(sdssagesohdust_levol_zbin3.r23branch_ew_alpha_gr_kk04,2) eq 'U'),nindx_levol_zbin3)

    mr_levol_zbin3 = sdssagesohdust_levol_zbin3[indx_levol_zbin3].ugriz_absmag[2]
    mg_levol_zbin3 = sdssagesohdust_levol_zbin3[indx_levol_zbin3].ugriz_absmag[1]
    mass_levol_zbin3 = sdssagesohdust_levol_zbin3[indx_levol_zbin3].mass + im_convert_imf(/from_chabrier)
    oh_levol_zbin3 = sdssagesohdust_levol_zbin3[indx_levol_zbin3].zstrong_ew_alpha_gr_12oh_kk04
    oherr_levol_zbin3 = sdssagesohdust_levol_zbin3[indx_levol_zbin3].zstrong_ew_alpha_gr_12oh_kk04_err
    
; 0.60<z<0.80

    indx_levol_zbin4 = where((sdssagesohdust_levol_zbin4.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(sdssagesohdust_levol_zbin4.r23branch_ew_alpha_gr_kk04,2) eq 'U'),nindx_levol_zbin4)

    mr_levol_zbin4 = sdssagesohdust_levol_zbin4[indx_levol_zbin4].ugriz_absmag[2]
    mg_levol_zbin4 = sdssagesohdust_levol_zbin4[indx_levol_zbin4].ugriz_absmag[1]
    mass_levol_zbin4 = sdssagesohdust_levol_zbin4[indx_levol_zbin4].mass + im_convert_imf(/from_chabrier)
    oh_levol_zbin4 = sdssagesohdust_levol_zbin4[indx_levol_zbin4].zstrong_ew_alpha_gr_12oh_kk04
    oherr_levol_zbin4 = sdssagesohdust_levol_zbin4[indx_levol_zbin4].zstrong_ew_alpha_gr_12oh_kk04_err
    
; -------------------------
; now make the plots!    
; -------------------------

    mrrange1 = [-16.2,-25.5]
    mgrange1 = [-15.0,-23.9]
    massrange1 = [8.3,12.1]
    ohrange1 = [8.3,9.39]

; LZ relation    
    
    if keyword_set(pdf) then loadct, xcolortable, /silent ; change the color table
    psname = 'ages_sdssages_mr_vs_12oh_zbins'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)
    
; 1a) 0.01<z<0.20

    im_hogg_scatterplot, mr_zbin1, oh_zbin1, position=pos[*,0], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle='', ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.01<z<0.2', /left, /top, box=0, charsize=charsize_2, margin=0
    oplot_lzevol, lzevol_mr[0], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      plot_local_levol=0, plot_evolved=0, thick=linethick1

    im_hogg_scatterplot, mr_noevol_zbin1, oh_noevol_zbin1, /noerase, position=pos[*,1], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle='', ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.01<z<0.2', /left, /top, box=0, charsize=charsize_2, margin=0
    oplot_lzevol, lzevol_mr_noevol[0], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      plot_local_levol=0, plot_evolved=0, thick=linethick1

    im_hogg_scatterplot, mr_levol_zbin1, oh_levol_zbin1, /noerase, position=pos[*,2], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle='', ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.01<z<0.2', /left, /top, box=0, charsize=charsize_2, margin=0
    oplot_lzevol, lzevol_mr_levol[0], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      plot_local_levol=0, plot_evolved=0, thick=linethick1

; 1b) 0.2<z<0.4

    im_hogg_scatterplot, mr_zbin2, oh_zbin2, /noerase, position=pos[*,3], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle='', ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.2<z<0.4', /left, /top, box=0, charsize=charsize_2, margin=0
    oplot_lzevol, lzevol_mr[1], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      /plot_evolved, thick=linethick1

    im_hogg_scatterplot, mr_noevol_zbin2, oh_noevol_zbin2, /noerase, position=pos[*,4], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle='', ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.2<z<0.4', /left, /top, box=0, charsize=charsize_2, margin=0
    oplot_lzevol, lzevol_mr_noevol[1], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      /plot_evolved, thick=linethick1

    im_hogg_scatterplot, mr_levol_zbin2, oh_levol_zbin2, /noerase, position=pos[*,5], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle='', ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.2<z<0.4', /left, /top, box=0, charsize=charsize_2, margin=0
    oplot_lzevol, lzevol_mr_levol[1], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      /plot_evolved, thick=linethick1

; 1c) 0.4<z<0.6

    im_hogg_scatterplot, mr_zbin3, oh_zbin3, /noerase, position=pos[*,6], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle='', ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.4<z<0.6', /left, /top, box=0, charsize=charsize_2, margin=0
    oplot_lzevol, lzevol_mr[2], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      /plot_evolved, thick=linethick1

    im_hogg_scatterplot, mr_noevol_zbin3, oh_noevol_zbin3, /noerase, position=pos[*,7], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle='', ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.4<z<0.6', /left, /top, box=0, charsize=charsize_2, margin=0
    oplot_lzevol, lzevol_mr_noevol[2], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      /plot_evolved, thick=linethick1

    im_hogg_scatterplot, mr_levol_zbin3, oh_levol_zbin3, /noerase, position=pos[*,8], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle='', ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.4<z<0.6', /left, /top, box=0, charsize=charsize_2, margin=0
    oplot_lzevol, lzevol_mr_levol[2], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      /plot_evolved, thick=linethick1

; 1d) 0.6<z<0.8

    im_hogg_scatterplot, mr_zbin4, oh_zbin4, /noerase, position=pos[*,9], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle=mrtitle1, ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.6<z<0.8', /left, /top, box=0, charsize=charsize_2, margin=0
    oplot_lzevol, lzevol_mr[3], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      /plot_evolved, thick=linethick1

    im_hogg_scatterplot, mr_noevol_zbin4, oh_noevol_zbin4, /noerase, position=pos[*,10], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle=mrtitle1, ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.6<z<0.8', /left, /top, box=0, charsize=charsize_2, margin=0
    oplot_lzevol, lzevol_mr_noevol[3], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      /plot_evolved, thick=linethick1

    im_hogg_scatterplot, mr_levol_zbin4, oh_levol_zbin4, /noerase, position=pos[*,11], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle=mrtitle1, ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.6<z<0.8', /left, /top, box=0, charsize=charsize_2, margin=0
    oplot_lzevol, lzevol_mr_levol[3], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      /plot_evolved, thick=linethick1

; y-title

    xyouts, pos[0,0]-0.08, pos[1,0], textoidl(ohtitle1), charsize=charsize_7, $
      align=0.5, orientation=90, /normal, color=fsc_color(textcolor1,100)
    xyouts, pos[0,6]-0.08, pos[1,6], textoidl(ohtitle1), charsize=charsize_7, $
      align=0.5, orientation=90, /normal, color=fsc_color(textcolor1,100)
    
; big-title

    xyouts, (pos[2,0]-pos[0,0])/2.0+pos[0,0], pos[3,0]+0.01, textoidl('AGES'), $
      charsize=charsize_4, align=0.5, /normal, color=fsc_color(textcolor1,100)
    xyouts, (pos[2,1]-pos[0,1])/2.0+pos[0,1], pos[3,0]+0.01, textoidl('SDSS-AGES; Q=0'), $
      charsize=charsize_4, align=0.5, /normal, color=fsc_color(textcolor1,100)
    xyouts, (pos[2,2]-pos[0,2])/2.0+pos[0,2], pos[3,0]+0.01, textoidl('SDSS-AGES; Q=1.5'), $
      charsize=charsize_4, align=0.5, /normal, color=fsc_color(textcolor1,100)
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    
    
; MZ relation

    if keyword_set(pdf) then loadct, xcolortable, /silent ; change the color table
    psname = 'ages_sdssages_mass_vs_12oh_zbins'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

; 1a) 0.01<z<0.20

    im_hogg_scatterplot, mass_zbin1, oh_zbin1, position=pos[*,0], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle='', ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.01<z<0.2', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol[0], massaxis=bigmassaxis1, plot_evolved=0

    im_hogg_scatterplot, mass_noevol_zbin1, oh_noevol_zbin1, /noerase, position=pos[*,1], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle='', ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.01<z<0.2', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol_noevol[0], massaxis=bigmassaxis1, plot_evolved=0

    im_hogg_scatterplot, mass_levol_zbin1, oh_levol_zbin1, /noerase, position=pos[*,2], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle='', ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.01<z<0.2', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol_levol[0], massaxis=bigmassaxis1, plot_evolved=0

; 1b) 0.2<z<0.4

    im_hogg_scatterplot, mass_zbin2, oh_zbin2, /noerase, position=pos[*,3], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle='', ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.2<z<0.4', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol[1], massaxis=bigmassaxis1, plot_evolved=1

    im_hogg_scatterplot, mass_noevol_zbin2, oh_noevol_zbin2, /noerase, position=pos[*,4], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle='', ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.2<z<0.4', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol_noevol[1], massaxis=bigmassaxis1, plot_evolved=1

    im_hogg_scatterplot, mass_levol_zbin2, oh_levol_zbin2, /noerase, position=pos[*,5], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle='', ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.2<z<0.4', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol_levol[1], massaxis=bigmassaxis1, plot_evolved=1

; 1c) 0.4<z<0.6

    im_hogg_scatterplot, mass_zbin3, oh_zbin3, /noerase, position=pos[*,6], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle='', ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.4<z<0.6', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol[2], massaxis=bigmassaxis1, plot_evolved=1

    im_hogg_scatterplot, mass_noevol_zbin3, oh_noevol_zbin3, /noerase, position=pos[*,7], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle='', ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.4<z<0.6', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol_noevol[2], massaxis=bigmassaxis1, plot_evolved=1

    im_hogg_scatterplot, mass_levol_zbin3, oh_levol_zbin3, /noerase, position=pos[*,8], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle='', ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.4<z<0.6', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol_levol[2], massaxis=bigmassaxis1, plot_evolved=1

; 1d) 0.6<z<0.8

    im_hogg_scatterplot, mass_zbin4, oh_zbin4, /noerase, position=pos[*,9], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle=masstitle1, ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.6<z<0.8', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol[3], massaxis=bigmassaxis1, plot_evolved=1

    im_hogg_scatterplot, mass_noevol_zbin4, oh_noevol_zbin4, /noerase, position=pos[*,10], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle=masstitle1, ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.6<z<0.8', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol_noevol[3], massaxis=bigmassaxis1, plot_evolved=1

    im_hogg_scatterplot, mass_levol_zbin4, oh_levol_zbin4, /noerase, position=pos[*,11], /outliers, $
      outpsym=symcat(sdsssym1), outsymsize=sdsspsize1, outcolor=sdsscolor1, $
      xtitle=masstitle1, ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.6<z<0.8', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol_levol[3], massaxis=bigmassaxis1, plot_evolved=1

; y-title

    xyouts, pos[0,0]-0.08, pos[1,0], textoidl(ohtitle1), charsize=charsize_7, $
      align=0.5, orientation=90, /normal, color=fsc_color(textcolor1,100)
    xyouts, pos[0,6]-0.08, pos[1,6], textoidl(ohtitle1), charsize=charsize_7, $
      align=0.5, orientation=90, /normal, color=fsc_color(textcolor1,100)
    
; big-title

    xyouts, (pos[2,0]-pos[0,0])/2.0+pos[0,0], pos[3,0]+0.01, textoidl('AGES'), $
      charsize=charsize_4, align=0.5, /normal, color=fsc_color(textcolor1,100)
    xyouts, (pos[2,1]-pos[0,1])/2.0+pos[0,1], pos[3,0]+0.01, textoidl('SDSS-AGES; Q=0'), $
      charsize=charsize_4, align=0.5, /normal, color=fsc_color(textcolor1,100)
    xyouts, (pos[2,2]-pos[0,2])/2.0+pos[0,2], pos[3,0]+0.01, textoidl('SDSS-AGES; Q=1.5'), $
      charsize=charsize_4, align=0.5, /normal, color=fsc_color(textcolor1,100)
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ---------------------------------------------------------------------------
; AGES - LZ and MZ in 4 redshift bins
; ---------------------------------------------------------------------------
    
    mrrange1 = [-16.2,-25.5]
    mgrange1 = [-15.0,-23.9]
    massrange1 = [8.3,12.1]
    ohrange1 = [8.3,9.39]
    
    xpage = 8.5 & ypage = 11.0
    pagemaker, nx=2, ny=4, xspace=0, yspace=0.0, width=[3.5,3.5], height=2.375*[1,1,1,1], $
      xmargin=[1.1,0.4], ymargin=[0.2,1.3], xpage=xpage, ypage=ypage, $
      position=pos, /normal

    if keyword_set(pdf) then loadct, xcolortable, /silent ; change the color table
    psname = 'ages_mr_mass_vs_12oh_zbins'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)
    
; don't change the "1" if you're still using loadct!!    

; -------------------------
; define the data    
; -------------------------
    
; 0.01<z<0.20

    indx_zbin1 = where((agesohdust.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(agesohdust.r23branch_ew_alpha_gr_kk04,2) eq 'U') and $
      (ageskcorr.z gt 0.01) and (ageskcorr.z lt 0.20),nindx_zbin1)

    mr_zbin1 = ageskcorr[indx_zbin1].ugriz_absmag[2]
    mg_zbin1 = ageskcorr[indx_zbin1].ugriz_absmag[1]
    mass_zbin1 = ageskcorr[indx_zbin1].mass + im_convert_imf(/from_chabrier)
    oh_zbin1 = agesohdust[indx_zbin1].zstrong_ew_alpha_gr_12oh_kk04
    oherr_zbin1 = agesohdust[indx_zbin1].zstrong_ew_alpha_gr_12oh_kk04_err
    weight_zbin1 = ageskcorr[indx_zbin1].spec_weight

; 0.20<z<0.40

    indx_zbin2 = where((agesohdust.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(agesohdust.r23branch_ew_alpha_gr_kk04,2) eq 'U') and $
      (ageskcorr.z gt 0.20) and (ageskcorr.z lt 0.40),nindx_zbin2)

    mr_zbin2 = ageskcorr[indx_zbin2].ugriz_absmag[2]
    mg_zbin2 = ageskcorr[indx_zbin2].ugriz_absmag[1]
    mass_zbin2 = ageskcorr[indx_zbin2].mass + im_convert_imf(/from_chabrier)
    oh_zbin2 = agesohdust[indx_zbin2].zstrong_ew_alpha_gr_12oh_kk04
    oherr_zbin2 = agesohdust[indx_zbin2].zstrong_ew_alpha_gr_12oh_kk04_err
    weight_zbin2 = ageskcorr[indx_zbin2].spec_weight

; 0.40<z<0.60

    indx_zbin3 = where((agesohdust.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(agesohdust.r23branch_ew_alpha_gr_kk04,2) eq 'U') and $
      (ageskcorr.z gt 0.40) and (ageskcorr.z lt 0.60),nindx_zbin3)

    mr_zbin3 = ageskcorr[indx_zbin3].ugriz_absmag[2]
    mg_zbin3 = ageskcorr[indx_zbin3].ugriz_absmag[1]
    mass_zbin3 = ageskcorr[indx_zbin3].mass + im_convert_imf(/from_chabrier)
    oh_zbin3 = agesohdust[indx_zbin3].zstrong_ew_alpha_gr_12oh_kk04
    oherr_zbin3 = agesohdust[indx_zbin3].zstrong_ew_alpha_gr_12oh_kk04_err
    weight_zbin3 = ageskcorr[indx_zbin3].spec_weight

; 0.60<z<0.80

    indx_zbin4 = where((agesohdust.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(agesohdust.r23branch_ew_alpha_gr_kk04,2) eq 'U') and $
      (ageskcorr.z gt 0.60) and (ageskcorr.z lt 0.80),nindx_zbin4)

    mr_zbin4 = ageskcorr[indx_zbin4].ugriz_absmag[2]
    mg_zbin4 = ageskcorr[indx_zbin4].ugriz_absmag[1]
    mass_zbin4 = ageskcorr[indx_zbin4].mass + im_convert_imf(/from_chabrier)
    oh_zbin4 = agesohdust[indx_zbin4].zstrong_ew_alpha_gr_12oh_kk04
    oherr_zbin4 = agesohdust[indx_zbin4].zstrong_ew_alpha_gr_12oh_kk04_err
    weight_zbin4 = ageskcorr[indx_zbin4].spec_weight

; -------------------------
; now make the plot!    
; -------------------------
    
; 0.01<z<0.20

    im_hogg_scatterplot, mr_zbin1, oh_zbin1, weight=weight_zbin1, position=pos[*,0], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle='', ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.01<z<0.2', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_lzevol, lzevol_mr[0], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      plot_local_levol=0, plot_evolved=0, thick=linethick1
    
    im_hogg_scatterplot, mass_zbin1, oh_zbin1, /noerase, position=pos[*,1], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle='', ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.01<z<0.2', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol[0], massaxis=bigmassaxis1, plot_evolved=0, thick=linethick1

; 0.2<z<0.4

    im_hogg_scatterplot, mr_zbin2, oh_zbin2, /noerase, position=pos[*,2], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle='', ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.2<z<0.4', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_lzevol, lzevol_mr[1], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      /plot_evolved, thick=linethick1

    im_hogg_scatterplot, mass_zbin2, oh_zbin2, /noerase, position=pos[*,3], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle='', ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.2<z<0.4', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol[1], massaxis=bigmassaxis1, plot_evolved=1, thick=linethick1

; 0.4<z<0.6

    im_hogg_scatterplot, mr_zbin3, oh_zbin3, /noerase, position=pos[*,4], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle='', ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.4<z<0.6', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_lzevol, lzevol_mr[2], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      /plot_evolved, thick=linethick1

    im_hogg_scatterplot, mass_zbin3, oh_zbin3, /noerase, position=pos[*,5], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle='', ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.4<z<0.6', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol[2], massaxis=bigmassaxis1, plot_evolved=1, thick=linethick1

; 0.6<z<0.8

    im_hogg_scatterplot, mr_zbin4, oh_zbin4, /noerase, position=pos[*,6], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle=mrtitle1, ytitle='', xrange=mrrange1, yrange=ohrange1, $
      charsize=charsize_7, $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.6<z<0.8', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_lzevol, lzevol_mr[3], absmagaxis=bigmraxis1, ohrange=ohrange1, $
      /plot_evolved, thick=linethick1

    im_hogg_scatterplot, mass_zbin4, oh_zbin4, /noerase, position=pos[*,7], /outliers, $
      outpsym=symcat(agessym1), outsymsize=agespsize1, outcolor=agescolor1, $
      xtitle=masstitle1, ytitle='', xrange=massrange1, yrange=ohrange1, $
      charsize=charsize_7, ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(axis_color,101)
    legend, '0.6<z<0.8', /left, /top, box=0, charsize=charsize_4, margin=0
    oplot_mzevol, mzevol[3], massaxis=bigmassaxis1, plot_evolved=1, thick=linethick1

; y-title

    xyouts, pos[0,0]-0.08, pos[1,0], textoidl(ohtitle1), charsize=charsize_7, $
      align=0.5, orientation=90, /normal, color=fsc_color(textcolor1,100)
    xyouts, pos[0,4]-0.08, pos[1,4], textoidl(ohtitle1), charsize=charsize_7, $
      align=0.5, orientation=90, /normal, color=fsc_color(textcolor1,100)

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
       loadct, xcolortable, /silent ; restore the color table
    endif

; ------------------------------------------------------------
; Redshift vs Delta [log(O/H)]
; ------------------------------------------------------------

    qz0 = 0.1 ; all evolutionary measurements are relative to z=0.1
    
    psname = 'redshift_vs_dlogoh'
    xpage = 8.5 & ypage = 8.8
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    pagemaker, nx=2, ny=2, xspace=0, yspace=0.0, width=3.3*[1,1], height=3.3*[1,1], $
      xmargin=[1.5,0.4], ymargin=[1.1,1.1], xpage=xpage, ypage=ypage, $
      position=pos, /normal

; make the plot    
    
    xtitle = 'Redshift'
    ytitle = '<\Delta'+'log(O/H)>'

    xrange = [-0.04,0.89] ; [-0.02,0.85]
    yrange = [-0.57,0.15]

; fixed luminosity, Q = 0
    
    plot, [0], [0], /nodata, xtitle='', ytitle=textoidl(ytitle), charsize=charsize_7, $
      xsty=9, ysty=1, xrange=xrange, position=pos[*,0], yrange=yrange, xtickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100)
    axis, /xaxis, xsty=1, xrange=xrange, charsize=charsize_7, xtitle='Lookback Time (Gyr)', $
      xtickv=getredshift(getage(0.0)-timelabel1), xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)'), color=fsc_color(textcolor1,100)
    djs_oplot, !x.crange, [0,0], line=5, color=fsc_color(textcolor1,100)
    im_legend, textoidl('L-Z; Q = 0'), /left, /bottom, margin=0, box=0, charsize=charsize_7, $
      textcolor=textcolor1

    djs_oplot, lzevol_mr.z, lzevol_mr.ldlogoh_cor_mean, psym=-symcat(dlogohcorsym,thick=symthick1), $
      color=fsc_color(dlogohcorcolor,1E8), symsize=dlogohcorpsize, thick=linethick2
    djs_oplot, lzevol_mr.z, lzevol_mr.ldlogoh_mean, psym=-symcat(dlogohsym,thick=symthick1), $
      color=fsc_color(dlogohcolor,1E8), symsize=dlogohpsize, thick=linethick2

; fixed luminosity, Q = 1.5
    
    plot, [0], [0], /nodata, /noerase, xtitle='', ytitle='', charsize=charsize_7, $
      xsty=9, ysty=1, xrange=xrange, position=pos[*,1], yrange=yrange, ytickname=replicate(' ',10), $
      xtickname=replicate(' ',10), color=fsc_color(textcolor1,100)
    axis, /xaxis, xsty=1, xrange=xrange, charsize=charsize_7, xtitle='Lookback Time (Gyr)', $
      xtickv=getredshift(getage(0.0)-timelabel1), xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)'), color=fsc_color(textcolor1,100)
    djs_oplot, !x.crange, [0,0], line=5, color=fsc_color(textcolor1,100)
    im_legend, textoidl('L-Z; Q = 1.5'), /left, /bottom, margin=0, box=0, charsize=charsize_7, $
      textcolor=textcolor1

    djs_oplot, lzevol_mr.z, lzevol_mr.ldlogoh_levol_cor_mean, psym=-symcat(dlogohcorsym,thick=symthick1), $
      color=fsc_color(dlogohcorcolor,1E8), symsize=dlogohcorpsize, thick=linethick2
    djs_oplot, lzevol_mr.z, lzevol_mr.ldlogoh_levol_mean, psym=-symcat(dlogohsym,thick=symthick1), $
      color=fsc_color(dlogohcolor,1E8), symsize=dlogohpsize, thick=linethick2

; fixed stellar mass, Q=0
    
    plot, [0], [0], /nodata, /noerase, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), charsize=charsize_7, $
      xsty=9, ysty=1, xrange=xrange, position=pos[*,2], yrange=yrange, color=fsc_color(textcolor1,100)
    djs_oplot, !x.crange, [0,0], line=5, color=fsc_color(textcolor1,100)
    im_legend, textoidl('M-Z; Q=0'), /left, /bottom, margin=0, box=0, charsize=charsize_7, $
      textcolor=textcolor1

    djs_oplot, mzevol.z, mzevol.mdlogoh_cor_mean, psym=-symcat(dlogohcorsym,thick=symthick1), $
      color=fsc_color(dlogohcorcolor,1E8), symsize=dlogohcorpsize, thick=linethick2
    djs_oplot, mzevol.z, mzevol.mdlogoh_mean, psym=-symcat(dlogohsym,thick=symthick1), $
      color=fsc_color(dlogohcolor,1E8), symsize=dlogohpsize, thick=linethick2

; fixed stellar mass, Q=1.5
    
    plot, [0], [0], /nodata, /noerase, xtitle=textoidl(xtitle), ytitle='', charsize=charsize_7, $
      xsty=9, ysty=1, xrange=xrange, position=pos[*,3], yrange=yrange, ytickname=replicate(' ',10), $
      color=fsc_color(textcolor1,100)
    djs_oplot, !x.crange, [0,0], line=5, color=fsc_color(textcolor1,100)
    im_legend, textoidl('M-Z; Q=1.5'), /left, /bottom, margin=0, box=0, charsize=charsize_7, $
      textcolor=textcolor1

    djs_oplot, mzevol.z, mzevol.mdlogoh_levol_cor_mean, psym=-symcat(dlogohcorsym,thick=symthick1), $
      color=fsc_color(dlogohcorcolor,1E8), symsize=dlogohcorpsize, thick=linethick2
    djs_oplot, mzevol.z, mzevol.mdlogoh_mean, psym=-symcat(dlogohsym,thick=symthick1), $
      color=fsc_color(dlogohcolor,1E8), symsize=dlogohpsize, thick=linethick2

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ------------------------------------------------------------
; redshift vs fraction of stars/metals formed
; ------------------------------------------------------------
    
    psname = 'redshift_vs_metal_fraction'
    xpage = 8.5 & ypage = 8.5
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.6,0.4], ymargin=[1.0,1.0], xpage=xpage, ypage=ypage, $
      position=pos, /normal

    xrange = [0.0,0.85]
    yrange = [83.0,104.0] ; percent

    xtitle = 'Redshift'
    ytitle = 'Percentage of Metals Synthesized by z=0.1'
    
; coefficients from plot_redshift_vs_sfr_mass_density    

    sfrd_zaxis = findgen((6.0-0.1)/0.01+1)*0.01+0.1 ; for integral
    sfrd_lookback = getage(0.0)-getage(sfrd_zaxis)  ; lookback time
    plot_sfrd_zaxis = [0.0,sfrd_zaxis]
    
    coeff = [0.010464392D,0.11303816D,3.4534319D,3.7109016D]
    sfrd = (coeff[0]+coeff[1]*sfrd_zaxis)/(1.0+(sfrd_zaxis/coeff[2])^coeff[3])
;   plot, alog10(1+sfrd_zaxis), alog10(sfrd), xsty=3, ysty=3, yr=[-2.3,-0.3]
    
; what fraction of the stars (and therefore metals) were formed by
; various redshifts; integrate the model; convert the SFRD to
; M_sun/Gyr; the factor of 1/64 (1.6%) is the IMF-weighted yield based
; on the Woosley & Weaver (1995) models, assuming a Salpeter IMF
; 0.1-125 M_sun (Conti et al. 2003; Pettini et al. 2006); note that
; Madau et al. (1996) used 1/42 (2.4%) for a Salpeter 0.1-100 IMF;
; also note that the exact value is ~independent of the IMF

    sfrd_total = im_integral(sfrd_lookback,sfrd*1D9,$              ; M_sun/Mpc^3
      min(sfrd_lookback),max(sfrd_lookback))                       
    sfrd_time = im_integral(sfrd_lookback,sfrd*1D9,sfrd_lookback,$ ; M_sun/Mpc^3
      replicate(max(sfrd_lookback),n_elements(sfrd_lookback)))
    sfrd_frac = sfrd_time/sfrd_total

    metal_total = sfrd_total/64.0 ; M_sun/Mpc^3
    metal_time = sfrd_time/64.0   ; M_sun/Mpc^3
    metal_frac = sfrd_frac

; interpolate the Pegase models and then integrate; note that I'm not
; dividing by the co-moving volume, but I could

    peg1_sfr = interpol(peg1.sfr,peg_zaxis,sfrd_zaxis)
    peg1_sfr_total = im_integral(sfrd_lookback,peg1_sfr*1D9,$ ; M_sun
      min(sfrd_lookback),max(sfrd_lookback)) 
    peg1_sfr_time = im_integral(sfrd_lookback,peg1_sfr*1D9,$ ; M_sun
      sfrd_lookback,replicate(max(sfrd_lookback),n_elements(sfrd_lookback)))
    peg1_sfr_frac = peg1_sfr_time/peg1_sfr_total
    
    peg3_sfr = interpol(peg3.sfr,peg_zaxis,sfrd_zaxis)
    peg3_sfr_total = im_integral(sfrd_lookback,peg3_sfr*1D9,$ ; M_sun
      min(sfrd_lookback),max(sfrd_lookback)) 
    peg3_sfr_time = im_integral(sfrd_lookback,peg3_sfr*1D9,$ ; M_sun
      sfrd_lookback,replicate(max(sfrd_lookback),n_elements(sfrd_lookback)))
    peg3_sfr_frac = peg3_sfr_time/peg3_sfr_total
    
    pegconst_sfr = interpol(pegconst.sfr,peg_zaxis,sfrd_zaxis)
    pegconst_sfr_total = im_integral(sfrd_lookback,pegconst_sfr*1D9,$ ; M_sun
      min(sfrd_lookback),max(sfrd_lookback)) 
    pegconst_sfr_time = im_integral(sfrd_lookback,pegconst_sfr*1D9,$ ; M_sun
      sfrd_lookback,replicate(max(sfrd_lookback),n_elements(sfrd_lookback)))
    pegconst_sfr_frac = pegconst_sfr_time/pegconst_sfr_total
    
; now make the plot    
    
    plot, [0], [0], /nodata, ysty=1, xsty=9, xrange=xrange, yrange=yrange, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=charsize_9, position=pos[*,0], color=fsc_color(textcolor1,100)
    axis, /xaxis, xthick=postthick1, xsty=1, xrange=xrange, charsize=charsize_9, $
      xtitle='Lookback Time (Gyr)', xtickv=getredshift(getage(0.0)-timelabel1), $
      xticks=n_elements(timelabel1)-1L, xtickname=string(timelabel1,format='(I0)'), $
      color=fsc_color(textcolor1,100)
; extrapolate to z<0.1    
    notzero = where(sfrd_frac gt 0.0)
    plot_sfrd_frac1 = interpol(100.0*(sfrd_frac[notzero]),$
      sfrd_zaxis[notzero],plot_sfrd_zaxis)
;   plot_sfrd_frac1 = interpol(alog10(sfrd_frac[notzero]),$
;     sfrd_zaxis[notzero],plot_sfrd_zaxis)
    djs_oplot, plot_sfrd_zaxis, plot_sfrd_frac1, line=0, thick=linethick1, color=fsc_color(textcolor1,100)
; here is the integral to z=0.1
;   djs_oplot, sfrd_zaxis, alog10(sfrd_frac[notzero]), line=0, $
;     thick=postthick3, color=fsc_color(textcolor1,150)

; overplot the Pegase models

    notzero = where(peg1_sfr_frac gt 0.0)
    plot_peg1_sfr_frac = interpol(100.0*(peg1_sfr_frac[notzero]),$
      sfrd_zaxis[notzero],plot_sfrd_zaxis)
;   plot_peg1_sfr_frac = interpol(alog10(peg1_sfr_frac[notzero]),$
;     sfrd_zaxis[notzero],plot_sfrd_zaxis)
    djs_oplot, plot_sfrd_zaxis, plot_peg1_sfr_frac, line=5, $
      thick=linethick1, color=fsc_color('dodger blue',148)

    notzero = where(peg3_sfr_frac gt 0.0)
    plot_peg3_sfr_frac = interpol(100.0*(peg3_sfr_frac[notzero]),$
      sfrd_zaxis[notzero],plot_sfrd_zaxis)
;   plot_peg3_sfr_frac = interpol(alog10(peg3_sfr_frac[notzero]),$
;     sfrd_zaxis[notzero],plot_sfrd_zaxis)
    djs_oplot, plot_sfrd_zaxis, plot_peg3_sfr_frac, line=3, $
      thick=linethick1, color=fsc_color('orange',149)

;;    notzero = where(pegconst_sfr_frac gt 0.0)
;;    plot_pegconst_sfr_frac = interpol(100.0*(pegconst_sfr_frac[notzero]),$
;;      sfrd_zaxis[notzero],plot_sfrd_zaxis)
;;;   plot_pegconst_sfr_frac = interpol(alog10(pegconst_sfr_frac[notzero]),$
;;;     sfrd_zaxis[notzero],plot_sfrd_zaxis)
;;    djs_oplot, plot_sfrd_zaxis, plot_pegconst_sfr_frac, line=3, $
;;      thick=postthick3, color=fsc_color('purple',147)

; and finally the data    
    
    dlogoh = [ [lzevol_mr.ldlogoh_cor_median], [lzevol_mr.ldlogoh_levol_cor_median], $
      [mzevol.mdlogoh_cor_median], [mzevol.mdlogoh_levol_cor_median] ]
;   dlogoh = [ [lzevol_mr.ldlogoh_cor_mean], [lzevol_mr.ldlogoh_levol_cor_mean], [mzevol.mdlogoh_cor_mean] ]
    dlogoh_z = lzevol_mr.z
    dlogoh_z_err = [0.0,0.0,0.0,0.0] ; no error bar for the first point
;   dlogoh_z_err = [0.1,0.1,0.1,0.1] ; no error bar for the first point
    dlogoh_mean = fltarr(4) & dlogoh_mean_err = dlogoh_mean*0.0 & dlogoh_median = dlogoh_mean*0.0
    for ii = 0L, n_elements(dlogoh_z)-1L do begin
       dlogoh_median[ii] = djs_median(dlogoh[ii,*])
       dlogoh_mean[ii] = djs_mean(dlogoh[ii,*])
       dlogoh_mean_err[ii] = djsig(dlogoh[ii,*])
    endfor

; I'm cheating here by plotting the median!!    
;   oploterror, dlogoh_z, 100.0*10^dlogoh_median, dlogoh_z_err, $
;     100.0*dlogoh_mean_err*10^dlogoh_mean*alog(10.0), $
    oploterror, dlogoh_z, 100.0*10^dlogoh_mean, dlogoh_z_err, $
      100.0*dlogoh_mean_err*10^dlogoh_mean*alog(10.0), $
      psym=symcat(15), color=fsc_color('red',1E8), $
      errcolor=fsc_color('red',1E8), symsize=4.5, $
      errthick=linethick1

    im_legend, ['!9'+string(105B)+'!Xdt \rho_{*}(t)','\tau=1 Gyr','\tau=3 Gyr'], /left, /bottom, $
      box=0, charsize=charsize_9, spacing=2.5, color=[textcolor1,'dodger blue','orange'], $
      textcolor=[textcolor1,'dodger blue','orange'], line=[0,5,3], number=0.0, thick=linethick1, $
      pspacing=1.1*charsize_9

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ------------------------------------------------------------
; Luminosity evolution for blue galaxies (literature compilation)
; ------------------------------------------------------------

    qz0 = 0.1 ; all evolutionary measurements are relative to z=0.1

    psname = 'redshift_vs_mr_lit'
    xpage = 8.5 & ypage = 8.5
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.6,0.4], ymargin=[1.0,1.0], xpage=xpage, ypage=ypage, $
      position=pos, /normal

    xtitle = 'Redshift'
    ytitle = 'M_{0.1r}^{*} for Blue Galaxies'

    xrange = [-0.02,0.85]
    yrange = [-20.2,-22.5]

;   xrange = [-0.02,1.3]
;   yrange = [-20.2,-22.7]
    
    plot, [0], [0], /nodata, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, xsty=9, ysty=1, xrange=xrange, yrange=yrange, $
      position=pos[*,0], color=fsc_color(textcolor1,100)
    axis, /xaxis, xthick=postthick1, xsty=1, xrange=xrange, charsize=singlecharsize_0, $
      xtitle='Lookback Time (Gyr)', xtickv=getredshift(getage(0.0)-timelabel1), $
      xticks=n_elements(timelabel1)-1L, xtickname=string(timelabel1,format='(I0)'), $
      color=fsc_color(textcolor1,100)

; overplot the line for 1.5 mag/z of luminosity evolution    

    djs_oplot, zaxis1, poly(zaxis1-0.1,[-20.95,-1.5]), line=0, $
      thick=postthick4, color=fsc_color(textcolor1,100)
;   im_legend, textoidl('1.5 mag z^{-1}'), charsize=singlecharsize_0, line=0, $
;     thick=postthick4, color=textcolor1, /right, /bottom, box=0, $
;     textcolor=textcolor1

; -------------------------
; Eisenstein et al. 2009 [AGES]; Omega_0=0.3, Omega_lamba=0.7, h=1.0;
; AB; alpha fixed at -1.10 for all blue galaxies; evolving color cut
; based on the (u0.1-r0.1) color: A=u0.1-r0.1+0.08(M_r0.1+20);  
; -------------------------

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
; -------------------------

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
; -------------------------

    z_f06 = [0.3,0.5,0.7,0.9,1.1]
    zerr_f06 = [0.1,0.1,0.1,0.1,0.1]

    mstar_f06 = [-20.74,-21.10,-21.30,-21.10,-21.25] + Bvega2ab + B2r01 ; Vega-->AB; B-->r0.1
    mstarerr_f06 = [0.20,0.15,0.16,0.17,0.18]

; now make the plot!    
    
    oploterror, z_f06, mstar_f06, zerr_f06, mstarerr_f06, symsize=3.8, errthick=linethick1, $
      psym=symcat(f06sym), color=fsc_color(f06color,1E8), errcolor=fsc_color(f06color,1E8)
    oploterror, z_w06, mstar_w06, zerr_w06, mstarerr_up_w06, /hi, symsize=3.5, errthick=linethick1, $
      psym=symcat(w06sym), color=fsc_color(w06color,1E8), errcolor=fsc_color(w06color,1E8)
    oploterror, z_w06, mstar_w06, zerr_w06, mstarerr_lo_w06, /lo, symsize=3.5, errthick=linethick1, $
      psym=symcat(w06sym), color=fsc_color(w06color,1E8), errcolor=fsc_color(w06color,1E8)
    oploterror, z_e07, mstar_e07, zerr_e07, mstarerr_e07, symsize=3.5, errthick=linethick1, $
      psym=symcat(e07sym), color=fsc_color(e07color,1E8), errcolor=fsc_color(e07color,1E8)

; legend

;   label = ['AGES+SDSS','COMBO-17','DEEP2']
    label = ['Willmer et al. 2006','Faber et al. 2007','Eisenstein et al. 2009']
    im_legend, label, /left, /top, box=0, charsize=charsize_6, $
      psym=[w06sym,f06sym,e07sym], symsize=[2.4,2.8,2.0], $
      spacing=2.0, color=[w06color,f06color,e07color], $
      textcolor=textcolor1

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ------------------------------------------------------------
; Redshift versus (a) ^{0.1}M_r and (b) stellar mass
; ------------------------------------------------------------

    psname = 'ages_redshift_vs_mr_mass'
    xpage = 8.5 & ypage = 7.8
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    pagemaker, nx=1, ny=2, width=7.0, height=3.0*[1,1], xmargin=[1.2,0.3], $
      ymargin=[0.8,1.0], xspace=0.0, yspace=0.0, xpage=xpage, ypage=ypage, $
      position=pos, /normal

; now make the plot    
    
    z = mz_ancillary.z
    mr = mz_ancillary.ugriz_absmag[2]
    gr = mz_ancillary.ugriz_absmag[1]-mz_ancillary.ugriz_absmag[2]
    mass = mz_ancillary.mass + im_convert_imf(/from_chabrier)
    weight = mz_ancillary.spec_weight
    
    xtitle = 'Redshift'
    ytitle1 = mrtitle1
    ytitle2 = masstitle1

    xrange = [-0.02,0.85]
    yrange1 = [-14.5,-24.7]
    yrange2 = [7.5,12.5]

; redshift vs ^{0.1} M_r
    
    plot, [0], [0], /nodata, ysty=1, xsty=9, xrange=xrange, yrange=yrange1, $
      xtitle='', ytitle=ytitle1, xtickname=replicate(' ',10), $
      charsize=charsize_8, position=pos[*,0], yminor=4, $
      color=fsc_color(textcolor1,100)
    djs_oplot, z, mr, psym=symcat(agessym2), symsize=agespsize2, color=agescolor2
    djs_oplot, !x.crange, mrstar*[1,1], line=2, thick=linethick2, color='red'
    axis, /xaxis, xsty=1, xrange=xrange, charsize=charsize_8, $
      xtitle='Lookback Time (Gyr)', xtickv=getredshift(getage(0.0)-timelabel1), $
      xticks=n_elements(timelabel1)-1L, xtickname=string(timelabel1,format='(I0)'), $
      color=fsc_color(textcolor1,100)

; z vs stellar mass
    
    plot, [0], [0], /nodata, /noerase, ysty=1, xsty=9, $
      xrange=xrange, yrange=yrange2, xtitle=xtitle, ytitle=ytitle2, $
      charsize=charsize_8, position=pos[*,1], yminor=5, $
      color=fsc_color(textcolor1,100)
    djs_oplot, z, mass, psym=symcat(agessym2), symsize=agespsize2, color=agescolor2
    djs_oplot, !x.crange, mstar*[1,1], line=2, thick=linethick2, color='red'

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    

return
end

