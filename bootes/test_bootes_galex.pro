pro test_bootes_galex
; jm10feb13ucsd - build a bunch of QAplots
    
    common test_galex, bootes, galex, cosmos, em_cosmos
    
    qaplotspath = ages_path(/qaplots)
    bootespath = getenv('RESEARCHPATH')+'/data/bootes/'
    cosmospath = getenv('RESEARCHPATH')+'/data/cosmos/'

; BOOTES/GALEX     
    if (n_elements(bootes) eq 0) then begin
       nuv = mrdfits(bootespath+'bootes_NUV.fits.gz',1)
       good = where((nuv.mag_psf gt 17.0) and (nuv.mag_psf lt 26.0))
       nuv = nuv[good]
       iband = mrdfits(bootespath+'bootes_I.fits.gz',1,rows=good)
       fuv = mrdfits(bootespath+'bootes_FUV.fits.gz',1,rows=good)

       nnuv = im_struct_trimtags(nuv,$
         select=tag_names(nuv),newtags='nuv_'+tag_names(nuv))
       ffuv = im_struct_trimtags(fuv,$
         select=tag_names(fuv),newtags='fuv_'+tag_names(fuv))
       iiband = im_struct_trimtags(iband,$
         select=tag_names(iband),newtags='i_'+tag_names(iband))
       bootes = struct_addtags(struct_addtags(iiband,nnuv),ffuv)

;      bootes = mrdfits(agespath+'ages_bootes.fits.gz',1,silent=0)
;      galex = mrdfits(catpath+'GALEX/ages_galex_gr4.fits.gz',1,silent=0)
    endif

;; (EM)COSMOS    
;    if (n_elements(cosmos) eq 0) then begin
;       cosmos1 = [$
;         mrdfits(cosmospath+'cosmos_b-galexmatch.fits.gz',1),$
;         mrdfits(cosmospath+'cosmos_v1a-galexmatch.fits.gz',1),$
;         mrdfits(cosmospath+'cosmos_v1b-galexmatch.fits.gz',1)]
;       em_cosmos1 = mrdfits(cosmospath+'COSMOS_GALEX_emphot_v3.fits.gz',1)
;       spherematch, cosmos1.primus_ra, cosmos1.primus_dec, $
;         em_cosmos1.ra, em_cosmos1.dec, 1.0/3600.0, m1, m2
;       cosmos = cosmos1[m1]
;       em_cosmos = em_cosmos1[m2]
;    endif

; curve-of-growth    
    iso = where((bootes.i_segflags_aper_10 eq 0) and $ ; isolated
      (bootes.nuv_mag_psf gt 0.0) and (bootes.nuv_mag_psf lt 30.0),niso)
    tagrad = ['01','02','03','04','05','06','07','08','09','10','15','20']
    radius = float(tagrad)
    nradius = n_elements(radius)
    diff = fltarr(nradius,niso)
    for ii = 0, nradius-1 do begin
       indx = tag_indx(bootes,'nuv_mag_aper_'+tagrad[ii])
       diff[ii,*] = bootes[iso].nuv_mag_psf-bootes[iso].(indx)
    endfor       
    bigradius = rebin(reform(radius,nradius,1),nradius,niso)

    psfile = qaplotspath+'bootes_galex_cog.ps'
    im_plotconfig, 0, pos, psfile=psfile
    djs_plot, bigradius, diff, psym=3, xsty=3, ysty=1, $
      xrange=[1,20], yrange=[-10,2], position=pos, $
      xtitle='Aperture Diameter (arcsec)', $
      ytitle='\Delta'+'m (PSF minus Aperture, mag)'
    im_plotconfig, psfile=psfile, /psclose, /gzip

;   hogg_scatterplot, bigradius, diff, xsty=3, ysty=1, $
;     xrange=[1,20], yrange=[-10,2], position=pos, $
;     /internal, levels=[0.5,0.75,0.9], /outliers, $
;     outcolor=djs_icolor('grey'), xtitle='Aperture Diameter (arcsec)', $
;     ytitle='\Delta'+'m (PSF minus Aperture, mag)'
;   med = im_medxbin(bigradius,diff,1.0,minx=1.0,maxx=20.0,minpts=1)
    
; build the plot    
    psfile = qaplotspath+'galex_bootes_tests.ps'
    im_plotconfig, 0, pos, psfile=psfile

; ---------------------------------------------------------------------------
; bootes PSF vs pipeline    
    im_galex_to_maggies, bootes, bootes_maggies, bootes_ivarmaggies
    im_galex_to_maggies, galex, galex_maggies, galex_ivarmaggies
    bootes_mag = maggies2mag(bootes_maggies,ivar=bootes_ivarmaggies,magerr=bootes_magerr)
    galex_mag = maggies2mag(galex_maggies,ivar=galex_ivarmaggies,magerr=galex_magerr)

; NUV
    these = where((galex_mag[1,*] gt 0.0) and (bootes_mag[1,*] gt 0.0))
    hogg_scatterplot, galex_mag[1,these], bootes_mag[1,these]-galex_mag[1,these], $
      position=pos, xrange=[18,26], yrange=[-2.9,2.9], xsty=1, ysty=1, $
      /internal, levels=[0.5,0.75,0.9], /outliers, $
      outcolor=djs_icolor('grey'), xtitle='NUV (Pipeline, AB mag)', $
      ytitle='NUV Residuals (PSF minus Pipeline, AB mag)'
    djs_oplot, !x.crange, [0,0], line=0, thick=4
    legend, 'Bootes', /right, /bottom, box=0
;   iso = where(bootes.i_segflags_aper_10 eq 0,niso) ; isolated
;   djs_oplot, galex_mag[1,iso], bootes_mag[1,iso]-galex_mag[1,iso], $
;     psym=symcat(6,thick=2), symsize=0.3, color='red'
;   im_legend, 'Isolated (10" segflags=0)', /left, /bottom, $
;     box=0, margin=0, psym=symcat(6,thick=6), color='red', symsize=1.5
; FUV
    these = where((galex_mag[0,*] gt 0.0) and (bootes_mag[0,*] gt 0.0))
    hogg_scatterplot, galex_mag[0,these], bootes_mag[0,these]-galex_mag[0,these], $
      position=pos, xrange=[19,27], yrange=[-2.9,2.9], xsty=1, ysty=1, $
      /internal, levels=[0.5,0.75,0.9], /outliers, $
      outcolor=djs_icolor('grey'), xtitle='FUV (Pipeline, AB mag)', $
      ytitle='FUV Residuals (PSF minus Pipeline, AB mag)'
    djs_oplot, !x.crange, [0,0], line=0, thick=4
    legend, 'Bootes', /right, /bottom, box=0
    
; ---------------------------------------------------------------------------
; cosmos EM vs pipeline    
    im_galex_to_maggies, cosmos, cosmos_maggies, cosmos_ivarmaggies
    im_galex_to_maggies_em, em_cosmos, em_cosmos_maggies, em_cosmos_ivarmaggies
    cosmos_mag = maggies2mag(cosmos_maggies,ivar=cosmos_ivarmaggies,magerr=cosmos_magerr)
    em_cosmos_mag = maggies2mag(em_cosmos_maggies,ivar=em_cosmos_ivarmaggies,magerr=em_cosmos_magerr)

; NUV
    these = where((cosmos_mag[1,*] gt 0.0) and (em_cosmos_mag[1,*] gt 0.0))
    hogg_scatterplot, cosmos_mag[1,these], em_cosmos_mag[1,these]-cosmos_mag[1,these], $
      position=pos, xrange=[19,26], yrange=[-2.9,2.9], xsty=1, ysty=1, $
      /internal, levels=[0.5,0.75,0.9], /outliers, $
      outcolor=djs_icolor('grey'), xtitle='NUV (Pipeline, AB mag)', $
      ytitle='NUV Residuals (EM minus Pipeline, AB mag)'
    djs_oplot, !x.crange, [0,0], line=0, thick=4
    legend, 'COSMOS', /right, /bottom, box=0
; FUV
    these = where((cosmos_mag[0,*] gt 0.0) and (em_cosmos_mag[0,*] gt 0.0))
    hogg_scatterplot, cosmos_mag[0,these], em_cosmos_mag[0,these]-cosmos_mag[0,these], $
      position=pos, xrange=[19,27], yrange=[-2.9,2.9], xsty=1, ysty=1, $
      /internal, levels=[0.5,0.75,0.9], /outliers, $
      outcolor=djs_icolor('grey'), xtitle='FUV (Pipeline, AB mag)', $
      ytitle='FUV Residuals (EM minus Pipeline, AB mag)'
    djs_oplot, !x.crange, [0,0], line=0, thick=4
    legend, 'COSMOS', /right, /bottom, box=0

    im_plotconfig, psfile=psfile, /psclose, /gzip
    
stop    
    
;; magnitude vs magnitude error
;    djs_plot, nmag, nmagerr, position=pos, psym=3, $
;      ysty=3, yr=[0.0,1.0]
;    djs_oplot, pmag, pmagerr, ps=3, color='green'

return
end
    
