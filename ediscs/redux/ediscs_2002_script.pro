pro ediscs_2002_script, doplot=doplot
; jm04jun19uofa

    paramfile = ['ibatch_lss.txt','ibatch_mos_left.txt','ibatch_mos_center.txt','ibatch_mos_right.txt']
    procfile = ['objlist_lss.txt','objlist_mos_left.txt','objlist_mos_center.txt','objlist_mos_right.txt']
    skyfile = ['skylist_lss.txt','skylist_mos_left.txt','skylist_mos_center.txt','skylist_mos_right.txt']
    skyapfile = ['skyaplist_lss.txt','skyaplist_mos_left.txt','skyaplist_mos_center.txt','skyaplist_mos_right.txt']
    calibfile = ['caliblist_lss.txt','caliblist_mos_left.txt','caliblist_mos_center.txt','caliblist_mos_right.txt']
    dividefile = ['','dividelist_mos_left.txt','dividelist_mos_center.txt','dividelist_mos_right.txt']
    stdfile = 'stdlist_mos.txt'
    stdallfile = 'stdlist_mos_allstars.txt'

; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root='a.'
;   ispec_makelists, ['a.lss','a.mos*left','a.mos*center','a.mos*right'], $
;     ['lss','mos_left','mos_center','mos_right'], /all, /overwrite, /gzip
    
; ---------------------------------------------------------------------------    
; initial reductions - do not overscan subtract the MOS frames
; ---------------------------------------------------------------------------    

    iall, paramfile, procfile=procfile, doplot=doplot, /noarcfit, rednight=0
    iall, paramfile, procfile=procfile, doplot=doplot, norder_overscan=-1L, /noarcfit, rednight=[1,2,3]

; ---------------------------------------------------------------------------    
; generate the wavelength solutions    
; ---------------------------------------------------------------------------    
    
    lamppath = cwd()

; LSS    

    lampname = 'lamplines_lss.dat'
;   imake_lamplines, 'ra.lss_arclamp.fits.gz', 'lamplines_input.dat', /nolines, $
;     output_lampfile=lampname, lambda0=5615.0, dispersion=1.24

    readcol, lamppath+lampname, arc_pixel, format='X,X,X,X,F', /silent
    ibatch, paramfile[0], lampname=lampname, lamppath=lamppath, $
      intensity_cut=0.0, arc_pixel=arc_pixel, mintol=0.5, norder_arc=5, /arcfit

; MOS LEFT

    lampname = 'lamplines_mos_left.dat'
;   imake_lamplines, 'ra.mos_arclamp_left.fits.gz', 'lamplines_input.dat', /nolines, $
;     output_lampfile=lampname, lambda0=6760.0, dispersion=1.32

    readcol, lamppath+lampname, arc_pixel, format='X,X,X,X,F', /silent
    ibatch, paramfile[1], lampname=lampname, lamppath=lamppath, $
      intensity_cut=0.0, arc_pixel=arc_pixel, mintol=0.5, norder_arc=5, /arcfit

; MOS CENTER

    lampname = 'lamplines_mos_center.dat'
;   imake_lamplines, 'ra.mos_arclamp_center.fits.gz', 'lamplines_input.dat', /nolines, $
;     output_lampfile=lampname, lambda0=5465.0, dispersion=1.3

    readcol, lamppath+lampname, arc_pixel, format='X,X,X,X,F', /silent
    ibatch, paramfile[2], lampname=lampname, lamppath=lamppath, $
      intensity_cut=0.0, arc_pixel=arc_pixel, mintol=0.5, norder_arc=5, /arcfit

; MOS RIGHT

    lampname = 'lamplines_mos_right.dat'
;   imake_lamplines, 'ra.mos_arclamp_right.fits.gz', 'lamplines_input.dat', /nolines, $
;     output_lampfile=lampname, lambda0=4255.0, dispersion=1.26

    readcol, lamppath+lampname, arc_pixel, format='X,X,X,X,F', /silent
    ibatch, paramfile[3], lampname=lampname, lamppath=lamppath, $
      intensity_cut=0.0, arc_pixel=arc_pixel, mintol=0.5, norder_arc=5, /arcfit
    
; ---------------------------------------------------------------------------    
; select sky apertures    
; ---------------------------------------------------------------------------    
   
;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

; ---------------------------------------------------------------------------    
; sky subtract
; ---------------------------------------------------------------------------    

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; wavelength-calibrate and extract the LSS hot stars; generate the
; mean telluric spectrum
; ---------------------------------------------------------------------------    

    iall, paramfile, calibfile=calibfile, extfile='', sensname='', doplot=doplot, rednight=0

    readcol, calibfile[0], lsslist, format='A', comment='#', /silent
    ispec, 'w'+lsslist, specnames=hotlist, aperture=10, /tracespec, $
      traceorder=2, noplot=1, /wfits, /gzip

    hotlist = i1dnames(lsslist,aperture=10)
    
    tellfile = 'telluric_2002.fits'
    psname = 'qaplot_telluric_2002.ps'

    iconstruct_telluric, hotlist, tellfile=tellfile, psname=psname, $
      tellmethod=2L, contmethod=1L, /doplot;, /write

; ---------------------------------------------------------------------------    
; wavelength calibrate and extract the three MOS dome flats, and make
; a QA plot
; ---------------------------------------------------------------------------    

    ibatch, paramfile[1], caliblist='ra.mos_domeflat_left.fits.gz', $
      extfile='', sensname='', /calibrate
    ibatch, paramfile[2], caliblist='ra.mos_domeflat_center.fits.gz', $
      extfile='', sensname='', /calibrate
    ibatch, paramfile[3], caliblist='ra.mos_domeflat_right.fits.gz', $
      extfile='', sensname='', /calibrate

    speclist = 'wra.mos_domeflat_'+['left','center','right']+'.fits.gz'
    ispec, speclist, aperture=10.0, refrow=140, specnames=domeflatlist, /noskysub, $
      noplot=1, /wfits, /gzip

    dleft = rd1dspec(domeflatlist[[0]])
    dcenter = rd1dspec(domeflatlist[[1]])
    dright = rd1dspec(domeflatlist[[2]])

    dfpsplot, 'qaplot_mos_domeflats_2002.ps', /color, /landscape & thick = 8.0
;   window, xsize=650, ysize=350 & thick = 2.0
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=[3900,10000], $
      yrange=[1E4,2E6], charsize=2.0, charthick=thick, title='MOS Dome Flats - 2002', $
      xtitle='Wavelength [\AA]', ytitle='Counts', $
      xthick=thick, ythick=thick
    djs_oplot, dcenter.wave, dcenter.spec, ps=10, color='green', thick=2
    djs_oplot, dleft.wave, dleft.spec, ps=10, color='red', thick=2
    djs_oplot, dright.wave, dright.spec, ps=10, color='blue', thick=2
    dfpsclose & thick = 2.0
      
;   plot1dspec, domeflatlist, prefix='qaplot_', /postscript;, /normalize, normwave=6780.0
    
; ---------------------------------------------------------------------------    
; wavelength calibrate and extract the MOS standard stars, and 
; divide by the appropriate domeflats
; ---------------------------------------------------------------------------    

    domeflatlist = 'wra.mos_domeflat_'+['left','center','right']+'_010.ms.fits.gz'
    
    for i = 0L, 2L do begin

       readcol, dividefile[i+1L], dividelist, format='A', /silent, comment='#'
       ibatch, paramfile[i+1L], caliblist=dividelist, extfile='', sensname='', /calibrate

       ispec, 'w'+dividelist, /tracespec, traceorder=2, aperture=10.0, $
         specnames=specnames, /noskyshift, /noplot, /wfits, /gzip

       specnames = i1dnames(dividelist,aperture=10.0)
       ediscs_domeflat_normalize, specnames, domeflatlist[i], /wfits, /gzip

    endfor

; ---------------------------------------------------------------------------    
; next, stitch together the three (left, center, right)
; domeflat-normalized standard stars  
; ---------------------------------------------------------------------------    

    ediscs_stitch_standards, leftstdlist=dividefile[1], centerstdlist=dividefile[2], $
      rightstdlist=dividefile[3], suffix='2002', aperture=10.0, /doplot, /wfits, /gzip

; ---------------------------------------------------------------------------    
; finally, generate the sensitivity function    
; ---------------------------------------------------------------------------    

; sensitivity function including all stars

    readcol, stdallfile, stdalllist, format='A', /silent, comment='#'

    sensname = 'sensitivity_2002_allstars.fits' & senstitle = 'Sensitivity Function - All Stars - 2002'
    info = isensfunc(stdalllist,grey=1,senstitle=senstitle,stdpath='ctionewcal',$
      extfile='lasillaextinct.dat',bsorder_sens=8,sensname=sensname,$
      searchrad=300.0,/write,/doplot,slit_width=5.0)
        
; sensitivity function with outlier stars rejected
    
    readcol, stdfile, stdlist, format='A', /silent, comment='#'

    sensname = 'sensitivity_2002.fits' & senstitle = 'Sensitivity Function - 2002'
    info = isensfunc(stdlist,grey=1,senstitle=senstitle,stdpath='ctionewcal',$
      extfile='lasillaextinct.dat',bsorder_sens=8,sensname=sensname,$
      searchrad=300.0,/write,/doplot,slit_width=5.0)

; ---------------------------------------------------------------------------    
; flux-calibrate the standard stars and compare them with the
; published fluxes
; ---------------------------------------------------------------------------    

    sensname = 'sensitivity_2002.fits'
    readcol, stdallfile, stdalllist, format='A', /silent, comment='#'

    ediscs_calibrate, stdalllist, sensname, prefix='f', /wfits, /gzip
    icompare_fluxed_standards, 'f'+stdalllist, searchrad=300.0, $
      stdpath='ctionewcal', suffix='2002', /postscript
        
; ---------------------------------------------------------------------------    
; compare my sensitivity function and telluric spectrum with Pascale's
; ---------------------------------------------------------------------------    

    mysens = readfits('sensitivity_2002.fits',myhead,/silent)
    mywave = make_wave(myhead)

    psens = readfits('/d0/ioannis/ediscs/MPE_2002/sensitivity.fits',phead,/silent)
    pwave = make_wave(phead)

    dfpsplot, 'qaplot_compare_sensitivity_2002.ps', /color, /landscape & thick = 8.0
;   window, xsize=750, ysize=550 & thick = 2.0
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=minmax(mywave), $
      yrange=[29,32], charsize=2.0, charthick=thick, title='Sensitivity Function Comparison', $
      xtitle='Wavelength [\AA]', ytitle='2.5 log [Counts s^{-1} '+$
      '\AA^{-1}] / ['+flam_units()+']', xthick=thick, ythick=thick
    djs_oplot, mywave, mysens+4.2, thick=4, color='blue'
    djs_oplot, pwave, psens, line=2, thick=4, color='red'
    legend, ['Moustakas','Jablonka'], /right, /top, charsize=2, charthick=thick, $
      box=0, line=[0,2], thick=thick, color=djs_icolor(['blue','red'])
    dfpsclose & thick = 2.0

;   tell = 'telluri_2002.fits'

; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    html_path = '/home/ioannis/public_html/research/ediscs/redux/'
    ispec_webpage, '2002_stars', rootname='sra.*', html_path=html_path

; ---------------------------------------------------------------------------    
; generate tarballs of interest    
; ---------------------------------------------------------------------------    

    datapath = ediscs_path(/d2002)

    tarfile = '2002_lss2d.tar.gz'
    spawn, ['tar czvf '+datapath+tarfile+' wsra*lss_hip?????.fits.gz'], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh

    tarfile = '2002_lss1d.tar.gz'
    spawn, ['tar czvf '+datapath+tarfile+' wsra.lss*.ms.*'], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh

    tarfile = '2002_mos_domeflats2d.tar.gz'
    spawn, ['tar czvf '+datapath+tarfile+' '+strjoin('wra.mos_domeflat_'+$
      ['left','center','right']+'.fits.gz',' ')], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh

    tarfile = '2002_mos_domeflats1d.tar.gz'
    spawn, ['tar czvf '+datapath+tarfile+' wra.mos_domeflat_*.ms.fits.gz'], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh

    stars1d = strjoin(file_search('wsra.mos*.ms.fits.gz'),' ')
    stars2d = strjoin(repstr(stars1d,'_010.ms',''),' ')

    tarfile = '2002_mos_stars1d.tar.gz'
    spawn, ['tar czvf '+datapath+tarfile+' '+stars1d], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh

    tarfile = '2002_mos_stars2d.tar.gz'
    spawn, ['tar czvf '+datapath+tarfile+' '+stars2d], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh

; tar up everything

    pushd, datapath
    tarfile = 'redux02_1d.tar.gz'
    spawn, ['tar czvf '+tarfile+' 2002_lss1d.tar.gz 2002_mos_domeflats1d.tar.gz 2002_mos_stars1d.tar.gz'], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh
    popd

    pushd, datapath
    tarfile = 'redux02_2d.tar.gz'
    spawn, ['tar czvf '+tarfile+' 2002_lss2d.tar.gz 2002_mos_domeflats2d.tar.gz 2002_mos_stars2d.tar.gz'], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh
    popd

stop    

return
end
