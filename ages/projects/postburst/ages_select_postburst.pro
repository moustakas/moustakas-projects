pro ages_select_postburst, postscript=postscript
; jm07nov16nyu - write out an AGES data structure for the AGES
;                postburst/AGN project with M. Brown
; jm08augnyu - major rewrite/update

    path = ages_path(/projects)+'postburst/brown/'
    template_dir = ages_path(/specfit)
    
    ancillary = read_ages(/ancillary)
    ispec = read_ages(/ispec)
;   if (n_elements(ages) eq 0L) then ages = mrdfits(path+'ages_brown.fits.gz',1,/silent)
;   if (n_elements(sdss) eq 0L) then sdss = mrdfits(sdss_path()+'sdss_specdata_main_dr4.fits.gz',1,/silent)

    im_plotfaves, postscript=postscript

    plotsym, 0, 0.8, /fill
    snrcut = 5.0

    badplates = [106,110,209,310,311]

; AGES: EW(Ha)

; plates 612, 608, 504    
    
    indx = where($
      (ancillary.pass ne 104) and (ancillary.pass ne 113) and $
      (ancillary.pass ne 308) and $
      (ancillary.pass ne 106) and (ancillary.pass ne 110) and $
      (ancillary.pass ne 209) and (ancillary.pass ne 310) and $
      (ancillary.pass ne 311) and (ancillary.infiber_r gt 0.1) and $
      (ispec.continuum_snr gt 4.0) and $
      (ispec.minwave/(1.0+ancillary.z) le 3500.0) and $
      (ispec.nii_6584[1] ne -2.0) and $
      (ancillary.z gt 0.01) and (ancillary.z lt 0.4),nindx)
;     (ispec.h_alpha[0]/ispec.h_alpha[1] gt snrcut) and $
;     (ispec.lick_hd_a_model[0]/ispec.lick_hd_a_model[1] gt snrcut) and $
;     (ispec.lick_hg_a_model[0]/ispec.lick_hg_a_model[1] gt snrcut) and $
    
    hbalmer = (ispec[indx].lick_hd_a_model[0]+ispec[indx].lick_hg_a_model[0])/2.0
    ewha = alog10(ispec[indx].h_alpha_ew[0]>1E-10)
    ewoii = alog10(ispec[indx].oii_3727_ew[0]>1E-10)

    hda = ispec[indx].lick_hd_a[0]
    d4000 = ispec[indx].d4000_narrow[0]
;   d4000 = ispec[indx].d4000_narrow_model[0]

; do the selection

    balmermin = 3.0 ; [A]
    xselect = findgen((15.0-balmermin)/0.01)*0.01+balmermin
    yselect = 0.2*xselect ; from Mike's paper

    pb = where((ewha lt interpol(yselect,xselect,hbalmer)) and $
      (hbalmer gt balmermin),npb)
    ispec_pb = ispec[indx[pb]]
    ancillary_pb = ancillary[indx[pb]]
    splog, npb, nindx, npb/float(nindx)

; make some plots    

    if keyword_set(postscript) then begin
       psfile = 'ages_select_postburst.ps'
       dfpsplot, psfile, /square, /color
    endif
    
    djs_plot, [0], [0], /nodata, xrange=[-2,13], yrange=[-0.5,2.5], xsty=1, ysty=1, $
      charsize=2.0, xtitle='(H\gamma_{A}+H\delta_{A})/2 (\AA)', ytitle='log EW(H\alpha) (\AA)'
    djs_oplot, balmermin*[1,1], [!y.crange[0],interpol(yselect,xselect,balmermin)]
    djs_oplot, xselect, yselect
    djs_oplot, hbalmer, ewha, ps=6, sym=0.2
    djs_oplot, hbalmer[pb], ewha[pb], ps=8, color='red'
    if (not keyword_set(postscript)) then cc = get_kbrd(1)
    
    djs_plot, [0], [0], /nodata, xrange=[-2,13], yrange=[-0.5,2.5], xsty=1, ysty=1, $
      charsize=2.0, xtitle='(H\gamma_{A}+H\delta_{A})/2 (\AA)', $
      ytitle='log EW([O II] \lambda3727) (\AA)'
    djs_oplot, hbalmer, ewoii, ps=6, sym=0.2
    djs_oplot, hbalmer[pb], ewoii[pb], ps=8, color='red'
    if (not keyword_set(postscript)) then cc = get_kbrd(1)

    djs_plot, [0], [0], /nodata, xrange=[0.8,2.4], yrange=[-7,12], xsty=1, ysty=1, $
      charsize=2.0, ytitle='(H\gamma_{A}+H\delta_{A})/2 (\AA)', xtitle='D_{n}(4000)'
    djs_oplot, d4000, hbalmer, ps=6, sym=0.2
    djs_oplot, d4000[pb], hbalmer[pb], ps=8, color='red'
;   djs_oplot, d4000, hda, ps=6, sym=0.2
;   djs_oplot, d4000[pb], hda[pb], ps=8, color='red'

    if keyword_set(postscript) then begin
       dfpsclose
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf')
    endif

; write the spectra out    
stop    
    srt = reverse(sort(hbalmer[pb])) & ; srt = srt[0:1]
    ages_display_spectrum, ispec_pb[srt], /plotobswave, $
      /postscript, labeltype=6, psname='ages_postburst_all.ps', /pdf

; these were identified through visual inspection    
    
    crap = 'ages_'+[$
      '404/005',$
      '312/081',$
      '406/086',$
      '411/103',$
      '511/074',$
      '411/108',$
      '608/053',$
      '410/043',$
      '309/228',$
      '309/230',$
      '421/134',$
      '612/152',$
      '608/115',$
      '420/134',$
      '304/060',$
;     '312/038',$ ; possibly reject
      '302/154',$
      '304/008',$
      '604/282']
    rej = where_array(crap,ancillary_pb[srt].galaxy)
    final = lindgen(npb) & remove, rej, final

    ages_display_spectrum, ancillary_pb[srt[final]], /plotobswave, $
      /postscript, labeltype=6, psname='ages_postburst_final.ps', /pdf
    
    ages_display_spectrum, ancillary_pb[srt[rej]], /plotobswave, $
      /postscript, labeltype=6, psname='ages_postburst_rej.ps', /pdf

; write out the binary FITS tables of the good spectra

    if keyword_set(postscript) then begin
       mwrfits, ancillary_pb, 'ages_postburst_ancillary.fits', /create
       mwrfits, ispec_pb, 'ages_postburst_ispec.fits', /create
    endif
    
stop    
    
; ###########################################################################
; ###########################################################################

    alltemplates = 'BC03_'+['Z004','Z02','Z05']+'_salpeter_ages_v1.2_templates.fits'
    alltempinfo = [mrdfits(template_dir+alltemplates[0],1,/silent),$
      mrdfits(template_dir+alltemplates[1],1,/silent),$
      mrdfits(template_dir+alltemplates[2],1,/silent)]

    age_mass = fltarr(nindx)
    age_V = fltarr(nindx)
    for ii = 0L, nindx-1L do begin
       tempinfo = alltempinfo[ispec[indx[ii]].template_bestindx]
       ages = tempinfo.template_age/1D9                            ; template ages [yr]
       ml_V = tempinfo.template_ml_V                               ; V-band M/L
       mass_template = ispec[indx[ii]].continuum_coeff*4.0D*!dpi*$ ; stellar mass in each template [M_sun]
         (lumdist(ispec[indx[ii]].z_abs,H0=70.0,/silent)*3.085678D24)^2.0 
       age_mass[ii] = tsum(ages,ages*mass_template)/tsum(ages,mass_template)
       age_V[ii] = tsum(ages,ages*mass_template/ml_V)/tsum(ages,mass_template/ml_V)
    endfor
    
    window, 1
    djs_plot, [0], [0], /nodata, xrange=[-0.1,3.0], yrange=[-0.5,2.5], xsty=1, ysty=1, $
      charsize=2.0, xtitle='Mass-Weighted Age (Gyr)', ytitle='log EW(H\alpha) (\AA)'
;   djs_oplot, age_V, ewha, ps=3
    djs_oplot, age_mass, ewha, ps=4

    pb2 = where((age_mass lt 1.0) and (ewha lt 1.0),npb2)
;   pb2 = where((age_V lt 0.8),npb2)
;   pb2 = where((age_V lt 1.2) and (ewha lt 1.0),npb2)
    splog, npb2, nindx, npb2/float(nindx)
    
    ages_display_spectrum, ispec[indx[pb2]]
    
stop    
    
    
; SDSS: EW(Ha)
    
;   indx = where((ispec.lick_hd_a_model[0]/ispec.lick_hd_a_model[1] gt snrcut) and $
;     (ispec.lick_hg_a_model[0]/ispec.lick_hg_a_model[1] gt snrcut) and $
;     (ispec.oii_3727_ew[0]/ispec.oii_3727_ew[1] gt snrcut),nindx)
;   
;   hda = ispec[indx].lick_hd_a_model[0]
;   hga = ispec[indx].lick_hg_a_model[0]
;   ewoii = ispec[indx].oii_3727_ew[0]
;
;   djs_plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
;     charthick=postthick2, xrange=[-2,13], yrange=[-0.5,2.5], xsty=1, ysty=1, $
;     charsize=2.0, xtitle='(H\gamma_{A}+H\delta_{A})/2 (\AA)', ytitle='log EW([O II]) (\AA)'
;   djs_oplot, (hda+hga)/2.0, alog10(ewoii>1E-10), ps=8

; AGES: EW([O II])

    indx = where((ispec.lick_hd_a_model[0]/ispec.lick_hd_a_model[1] gt snrcut) and $
      (ispec.lick_hg_a_model[0]/ispec.lick_hg_a_model[1] gt snrcut) and $
      (ispec.oii_3727_ew[0]/ispec.oii_3727_ew[1] gt snrcut),nindx)
    
    hda = ispec[indx].lick_hd_a_model[0]
    hga = ispec[indx].lick_hg_a_model[0]
    ewoii = ispec[indx].oii_3727_ew[0]

    djs_plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      charthick=postthick2, xrange=[-2,13], yrange=[-0.5,2.5], xsty=1, ysty=1, $
      charsize=2.0, xtitle='(H\gamma_{A}+H\delta_{A})/2 (\AA)', ytitle='log EW([O II]) (\AA)'
    djs_oplot, (hda+hga)/2.0, alog10(ewoii>1E-10), ps=3

    pb = where(((hda+hga)/2.0 gt 5.0) and (ewoii lt 10.0),npb)
    splog, npb, nindx, npb/float(nindx)

stop
    
; write out a postscript plot of low- and high-z postbursts

    loz = where(ispec[indx[pb]].z lt 0.4,nloz,comp=hiz,ncomp=nhiz)
    splog, 'z<0.4: ', nloz, nloz/float(npb)
    splog, 'z>0.4: ', nhiz, nhiz/float(npb)
    
    psname = path+'ages_postburst_hiz.ps'
    if keyword_set(postscript) then splog, 'Writing '+psname
    ages_display_spectrum, ispec[indx[pb[hiz]]], labeltype=2L, $
      /plotobswave, postscript=postscript, psname=psname

    psname = path+'ages_postburst_loz.ps'
    if keyword_set(postscript) then splog, 'Writing '+psname
    ages_display_spectrum, ispec[indx[pb[loz]]], labeltype=2L, $
      /plotobswave, postscript=postscript, psname=psname
    
stop
    
stop    
    
return
end
    
