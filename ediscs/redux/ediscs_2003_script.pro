pro ediscs_2003_script, doplot=doplot
; jm04jun19uofa

    paramfile = ['ibatch_lss.txt','ibatch_mos_left.txt','ibatch_mos_center.txt','ibatch_mos_right.txt']
    procfile = ['objlist_lss.txt','objlist_mos_left.txt','objlist_mos_center.txt','objlist_mos_right.txt']
    skyfile = ['skylist_lss.txt','skylist_mos_left.txt','skylist_mos_center.txt','skylist_mos_right.txt']
    skyapfile = ['skyaplist_lss.txt','skyaplist_mos_left.txt','skyaplist_mos_center.txt','skyaplist_mos_right.txt']
    calibfile = ['caliblist_lss.txt','caliblist_mos_left.txt','caliblist_mos_center.txt','caliblist_mos_right.txt']

;   dividefile = ['','dividelist_mos_left.txt','dividelist_mos_center.txt','dividelist_mos_right.txt']
;   stdfile = 'stdlist_mos.txt'
;   stdallfile = 'stdlist_mos_allstars.txt'

; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root='a.'
;   ispec_makelists, ['a.lss','a.mos*left','a.mos*center','a.mos*right'], $
;     ['lss','mos_left','mos_center','mos_right'], /all, /overwrite, /gzip
    
; ---------------------------------------------------------------------------    
; initial reductions
; ---------------------------------------------------------------------------    

;;;   iall, paramfile, procfile=procfile, doplot=doplot, /noarcfit

; ---------------------------------------------------------------------------    
; generate the wavelength solutions    
; ---------------------------------------------------------------------------    
    
;   imake_lamplines, 'ra.lss_arclamp.fits.gz', 'lamplines_input.dat', /nolines, $
;     output_lampfile='lamplines_lss.dat', lambda0=5289.0, dispersion=1.634
;;;   ibatch, paramfile[0], lampname='lamplines_lss.dat', lamppath=cwd(), mintol=0.25, norder_arc=8, /arcfit

;   imake_lamplines, 'ra.mos_arclamp_left.fits.gz', 'lamplines_input.dat', /nolines, $
;     output_lampfile='lamplines_mos_left.dat', lambda0=6329.0, dispersion=1.654
;;;   ibatch, paramfile[1], lampname='lamplines_mos_left.dat', lamppath=cwd(), mintol=0.25, norder_arc=10, /arcfit

;   imake_lamplines, 'ra.mos_arclamp_center.fits.gz', 'lamplines_input.dat', /nolines, $
;     output_lampfile='lamplines_mos_center.dat', lambda0=5080.0, dispersion=1.656
;;;   ibatch, paramfile[2], lampname='lamplines_mos_center.dat', lamppath=cwd(), mintol=0.25, norder_arc=5, /arcfit
    
;   imake_lamplines, 'ra.mos_arclamp_right.fits.gz', 'lamplines_input.dat', /nolines, $
;     output_lampfile='lamplines_mos_right.dat', lambda0=4100.0, dispersion=1.617
;;;   ibatch, paramfile[3], lampname='lamplines_mos_right.dat', lamppath=cwd(), mintol=0.25, norder_arc=9, /arcfit
    
; ---------------------------------------------------------------------------    
; select sky apertures    
; ---------------------------------------------------------------------------    
   
;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

; ---------------------------------------------------------------------------    
; sky subtract
; ---------------------------------------------------------------------------    

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; wavelength-calibrate and extract the LSS/MOS stars
; ---------------------------------------------------------------------------    

    iall, paramfile, calibfile=calibfile, extfile='', sensname='', doplot=doplot

    for i = 0L, n_elements(paramfile)-1L do begin
       readcol, calibfile[i], starlist, format='A', comment='#', /silent
       ispec, 'w'+starlist, aperture=10, /tracespec, traceorder=2, noplot=1, /wfits, /gzip
    endfor
    
; ---------------------------------------------------------------------------    
; for testing, wavelength calibrate the arc lamps
; ---------------------------------------------------------------------------    

    ibatch, paramfile[1], caliblist='ra.mos_arclamp_left.fits.gz', extfile='', sensname='', /calibrate
    ibatch, paramfile[2], caliblist='ra.mos_arclamp_center.fits.gz', extfile='', sensname='', /calibrate
    ibatch, paramfile[3], caliblist='ra.mos_arclamp_right.fits.gz', extfile='', sensname='', /calibrate

; ---------------------------------------------------------------------------    
; wavelength calibrate and extract the three MOS dome flats 
; ---------------------------------------------------------------------------    

    ibatch, paramfile[1], caliblist='ra.mos_domeflat_left.fits.gz', extfile='', sensname='', /calibrate
    ibatch, paramfile[2], caliblist='ra.mos_domeflat_center.fits.gz', extfile='', sensname='', /calibrate
    ibatch, paramfile[3], caliblist='ra.mos_domeflat_right.fits.gz', extfile='', sensname='', /calibrate

    speclist = 'wra.mos_domeflat_'+['left','center','right']+'.fits.gz'
    ispec, speclist, aperture=10.0, refrow=140, specnames=domeflatlist, $
      /noskysub, /noskyshift, noplot=1, /wfits, /gzip

; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    html_path = '/home/ioannis/public_html/research/ediscs/redux/'
    ispec_webpage, '2003_stars', rootname='sra.*', html_path=html_path

; ---------------------------------------------------------------------------    
; generate tarballs of interest    
; ---------------------------------------------------------------------------    

    datapath = ediscs_path(/d2003)

    tarfile = '2003_lss2d.tar.gz'
    spawn, ['tar czvf '+datapath+tarfile+' wsra*lss_hip?????.fits.gz wsra*lss*_?.fits.gz'], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh

    tarfile = '2003_lss1d.tar.gz'
    spawn, ['tar czvf '+datapath+tarfile+' wsra.lss*.ms.*'], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh

    tarfile = '2003_mos_domeflats2d.tar.gz'
    spawn, ['tar czvf '+datapath+tarfile+' '+strjoin('wra.mos_domeflat_'+$
      ['left','center','right']+'.fits.gz',' ')], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh

    tarfile = '2003_mos_domeflats1d.tar.gz'
    spawn, ['tar czvf '+datapath+tarfile+' wra.mos_domeflat_*.ms.fits.gz'], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh

    stars1d = strjoin(file_search('wsra.mos*.ms.fits.gz'),' ')
    stars2d = strjoin(repstr(stars1d,'_010.ms',''),' ')

    tarfile = '2003_mos_stars1d.tar.gz'
    spawn, ['tar czvf '+datapath+tarfile+' '+stars1d], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh

    tarfile = '2003_mos_stars2d.tar.gz'
    spawn, ['tar czvf '+datapath+tarfile+' '+stars2d], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh

; tar up everything

    pushd, datapath
    tarfile = 'redux03_1d.tar.gz'
    spawn, ['tar czvf '+tarfile+' 2003_lss1d.tar.gz 2003_mos_domeflats1d.tar.gz 2003_mos_stars1d.tar.gz'], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh
    popd

    pushd, datapath
    tarfile = 'redux03_2d.tar.gz'
    spawn, ['tar czvf '+tarfile+' 2003_lss2d.tar.gz 2003_mos_domeflats2d.tar.gz 2003_mos_stars2d.tar.gz'], /sh
    spawn, ['ln -s '+datapath+tarfile+' '+html_path+tarfile], /sh
    popd

stop
    
; ###########################################################################    
; OLD!!
; ###########################################################################    
    
; ---------------------------------------------------------------------------    
; wavelength-calibrate and extract the LSS hot stars; generate the
; mean telluric spectrum
; ---------------------------------------------------------------------------    

    iall, paramfile, calibfile=calibfile, extfile='', sensname='', doplot=doplot, rednight=0

    readcol, calibfile[0], lsslist, format='A', comment='#', /silent
    ispec, 'w'+lsslist, specnames=hotlist, aperture=10, /tracespec, $
      traceorder=2, noplot=1, /wfits, /gzip

    hotlist = i1dnames(lsslist,aperture=10)
    
    tellfile = 'telluric_2003.fits'
    psname = 'qaplot_telluric_2003.ps'

    iconstruct_telluric, hotlist, tellfile=tellfile, psname=psname, $
      tellmethod=2L, contmethod=1L, /doplot, /write

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

    dfpsplot, 'qaplot_mos_domeflats_2003.ps', /color, /landscape & thick = 8.0
;   window, xsize=650, ysize=350 & thick = 2.0
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=[3900,10000], $
      yrange=[9E3,9.5E5], charsize=2.0, charthick=thick, title='MOS Dome Flats - 2003', $
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
         specnames=specnames, /noskyshift, noplot=0, /wfits, /gzip

       specnames = i1dnames(dividelist,aperture=10.0)
       ediscs_domeflat_normalize, specnames, domeflatlist[i], /wfits, /gzip

    endfor

; ---------------------------------------------------------------------------    
; next, stitch together the three (left, center, right)
; domeflat-normalized standard stars  
; ---------------------------------------------------------------------------    

    ediscs_stitch_standards, leftstdlist=dividefile[1], centerstdlist=dividefile[2], $
      rightstdlist=dividefile[3], suffix='2003', aperture=10.0, /doplot;, /wfits, /gzip

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
; wavelength calibrate and extract the three MOS dome flats 
; ---------------------------------------------------------------------------    

    ibatch, paramfile[1], caliblist='ra.mos_domeflat_left.fits.gz', $
      extfile='', sensname='', /calibrate
    ibatch, paramfile[2], caliblist='ra.mos_domeflat_center.fits.gz', $
      extfile='', sensname='', /calibrate
    ibatch, paramfile[3], caliblist='ra.mos_domeflat_right.fits.gz', $
      extfile='', sensname='', /calibrate

    ispec, 'wra.mos_domeflat_'+['left','center','right']+'.fits.gz', $
      aperture=10.0, refrow=123, specnames=domeflatlist, /noskysub, $
      noplot=1, /wfits, /gzip
    
; ---------------------------------------------------------------------------    
; wavelength calibrate and extract the MOS standard stars.  next,
; divide by the appropriate domeflats and create the sensitivity
; function; for comparison purposes, also generate the sensitivity
; function based on *all* the stars
; ---------------------------------------------------------------------------    

    domeflatlist = 'wra.mos_domeflat_'+['left','center','right']+'_010.ms.fits.gz'
    senstitle = 'MOS '+['Left','Center','Right']+' - 2003'
    sensalllist = ['sens_mos_left_allstars.fits','sens_mos_center_allstars.fits','sens_mos_right_allstars.fits']
    senslist = ['sens_mos_left.fits','sens_mos_center.fits','sens_mos_right.fits']
    bsorder_sens = [7L,18L,15L]
    
    for i = 0L, 2L do begin

       readcol, stdfile[i+1L], stdlist, format='A', /silent, comment='#'
       readcol, stdallstarsfile[i+1L], stdalllist, format='A', /silent, comment='#'

       ibatch, paramfile[i+1L], caliblist=stdalllist, extfile='', sensname='', /calibrate

       ispec, 'w'+stdalllist, /tracespec, traceorder=2, aperture=10.0, $
         specnames=specnames, /noskyshift, noplot=1, /wfits, /gzip

       ediscs_domeflat_normalize, specnames, domeflatlist[i], /wfits, /gzip

; full sensitivity function
       
       specnames = i1dnames(stdalllist,aperture=10.0)
       info = isensfunc('d'+specnames,grey=1,senstitle=senstitle[i],stdpath='ctionewcal',$
         extfile='lasillaextinct.dat',bsorder_sens=bsorder_sens[i],sensname=sensalllist[i],$
         searchrad=400.0,/write,/doplot,slit_width=5.0)

; selected sensitivity function
       
       specnames = i1dnames(stdlist,aperture=10.0)
       info = isensfunc('d'+specnames,grey=1,senstitle=senstitle[i],stdpath='ctionewcal',$
         extfile='lasillaextinct.dat',bsorder_sens=bsorder_sens[i],sensname=senslist[i],$
         searchrad=400.0,/write,/doplot,slit_width=5.0)

    endfor

; ---------------------------------------------------------------------------    
; stitch together the three sensitivity functions to create
; sensitivity_2003.fits 
; ---------------------------------------------------------------------------    
    
    ediscs_stitch_sensitivity, suffix='2003', leftsens='sens_mos_left.fits', $
      centersens='sens_mos_center.fits', rightsens='sens_mos_right.fits', /wfits

; ---------------------------------------------------------------------------    
; wavelength-calibrate (already done above), extract one-dimensional
; spectra, normalize by the dome-flat spectrum, and flux-calibrate
; ---------------------------------------------------------------------------    

;   iall, paramfile, calibfile=calibfile, sensname='', extfile='', rednight=[1,2,3]

; sensitivity function and telluric spectrum

    sensname = 'sensitivity_2003.fits'
    domeflatlist = 'wra.mos_domeflat_'+['left','center','right']+'_010.ms.fits.gz'

    for i = 0L, 2L do begin
    
       readcol, calibfile[i+1L], caliblist, format='A', /silent, comment='#'

       ispec, 'w'+caliblist, /tracespec, traceorder=2, aperture=10.0, $
         specnames=specnames, /noskyshift, /noplot, /wfits, /gzip

       ediscs_domeflat_normalize, specnames, domeflatlist[i], /wfits, /gzip

       ediscs_calibrate, 'd'+specnames, sensname, prefix='f', /wfits, /gzip

    endfor
    
; ---------------------------------------------------------------------------    
; compare the fluxed standard stars with the published
; spectrophotometric standards
; ---------------------------------------------------------------------------    

    starlist = file_search('fdwsra*center*.ms.fits.gz')
    icompare_fluxed_standards, starlist, searchrad=400.0, stdpath='ctionewcal'
        
; ###########################################################################    
; the code below here is work in progress
; ###########################################################################    
    
; telluric-correct the LSS hot stars    
    
    hotlist = i1dnames(lsslist,aperture=10)
    idivide_telluric, hotlist, tellfile, method=2, suffix='lss', /postscript

    f = i1dnames(stdlist,aper=20.0)

    info = isensfunc(i1dnames(stdlist,aperture=20),grey=2,sensname=senslist[i],$
      stdpath='ctionewcal',senstitle=senstitle[i],/doplot,/write,$
      slit_width=1.0,searchrad=240.0)

; ---------------------------------------------------------------------------    
; flux-calibrate calibrate the MOS standard stars
; ---------------------------------------------------------------------------    
    
    sensname = 'sens_4.5_04mar.fits'
    iall, paramfile, calibfile=calibfile, tellfile=tellfile, sensname=sensname

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    


; old wavelength code



    lamppath = cwd()

    lampname = 'lamplines_lss.dat'
;   imake_lamplines, 'ra.lss_arclamp.fits.gz', 'lamplines_input.dat', result, /nolines, $
;     output_lampfile=lampname, lambda0=5290.0, dispersion=1.63

    readcol, lamppath+lampname, arc_pixel, format='X,X,X,X,F', /silent
    ibatch, paramfile[0], lampname=lampname, lamppath=lamppath, $
      intensity_cut=0.0, arc_pixel=arc_pixel, mintol=0.5, norder_arc=5, /arcfit

; MOS LEFT

    lampname = 'lamplines_mos_left.dat'
;   imake_lamplines, 'ra.mos_arclamp_left.fits.gz', 'lamplines_input.dat', result, /nolines, $
;     output_lampfile=lampname, lambda0=6329.0, dispersion=1.65

    readcol, lamppath+lampname, arc_pixel, format='X,X,X,X,F', /silent
    ibatch, paramfile[1], lampname=lampname, lamppath=lamppath, $
      intensity_cut=0.0, arc_pixel=arc_pixel, mintol=0.5, norder_arc=8, /arcfit

; MOS CENTER

    lampname = 'lamplines_mos_center.dat'
;   imake_lamplines, 'ra.mos_arclamp_center.fits.gz', 'lamplines_input.dat', result, /nolines, $
;     output_lampfile=lampname, lambda0=5080.0, dispersion=1.65

    readcol, lamppath+lampname, arc_pixel, format='X,X,X,X,F', /silent
    ibatch, paramfile[2], lampname=lampname, lamppath=lamppath, $
      intensity_cut=0.0, arc_pixel=arc_pixel, mintol=0.5, norder_arc=5, /arcfit

; MOS RIGHT

    imake_lamplines, 'ra.mos_arclamp_right.fits.gz', 'lamplines_input.dat', /nolines, $
      output_lampfile='lamplines_mos_right.dat', lambda0=4100.0, dispersion=1.617
    ibatch, paramfile[3], lampname='lamplines_mos_right.dat', lamppath=cwd(), mintol=0.25, norder_arc=9, /arcfit

    lampname = 'lamplines_mos_right.dat'
;   imake_lamplines, 'ra.mos_arclamp_right.fits.gz', 'lamplines_input.dat', result, $
;     /nolines, output_lampfile=lampname, lambda0=4309.0, dispersion=1.62

    readcol, lamppath+lampname, arc_pixel, format='X,X,X,X,F', /silent
    ibatch, paramfile[3], lampname=lampname, lamppath=lamppath, $
      intensity_cut=0.0, arc_pixel=arc_pixel, mintol=0.5, norder_arc=5, /arcfit

    ibatch, paramfile[3], lampname=lampname, lamppath=lamppath, dxcoeff=[100,0.2], $
      intensity_cut=0.0, arc_pixel=arc_pixel, mintol=2, norder_arc=5, /arcfit


return
end
