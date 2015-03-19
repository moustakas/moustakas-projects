pro decals_dr1_diagnostics, debug=debug, build_info=build_info, qaplots=qaplots
; jm15mar17siena - build DR1 diagnostics

    dr1path = getenv('DECALS_DIR')+'/'
    infofile = dr1path+'qa-dr1-info.fits'
    
; ---------------------------------------------------------------------------
; build the info structure    
    if keyword_set(build_info) then begin
    
       info_template = {$
         brickid:   0L,$
         brickname: '',$
         ra:       0D,$
         dec:      0D,$
         dm_g:  -99.0,$
         dm_r:  -99.0,$
         dm_z:  -99.0,$
         dra_med:  0D,$
         ddec_med: 0D,$
         dra_sig:  0D,$
         ddec_sig: 0D}
       
;      dindx = [1,2,4]
       dindx = [0,1,2]
       sindx = [6,7,9]
       magcut = [20.0,19.5,19.0]
       band = ['g','r','z']
       nband = n_elements(band)
       
       radius = [0.5,0.75,1.0,1.5,2.0,3.5,5.0,7.0]
    
       allbrick = file_basename(file_search(dr1path+'tractor/*',/test_dir,count=nbrick))
;      for ii = 0L, 30 do begin
       for ii = 0L, nbrick-1 do begin
          catfile = file_search(dr1path+'tractor/'+allbrick[ii]+'/tractor-*.fits',count=ncat)
          info1 = replicate(info_template,ncat)
          for ic = 0L, ncat-1 do begin
             print, format='("Brick ",I0,"/",I0," Cat ",I0,"/",I0, A10,$)', $
               ii+1, nbrick, ic+1, ncat, string(13b)
             cat = mrdfits(catfile[ic],1,/silent)
             cat = cat[where(cat.brick_primary eq 'T')]
             nobj = n_elements(cat)
             
             info1[ic].brickid = cat[0].brickid
             info1[ic].brickname = cat[0].brickname
             info1[ic].ra = djs_median(cat.ra)
             info1[ic].dec = djs_median(cat.dec)
             
; compare against SDSS coordinates          
             isdss = where(cat.sdss_objid ne 0)
             info1[ic].dra_med = djs_median(cat[isdss].ra-cat[isdss].sdss_ra)*3600
             info1[ic].ddec_med = djs_median(cat[isdss].dec-cat[isdss].sdss_dec)*3600
             info1[ic].dra_sig = djsig(cat[isdss].ra-cat[isdss].sdss_ra)*3600
             info1[ic].ddec_sig = djsig(cat[isdss].dec-cat[isdss].sdss_dec)*3600
             
; compare against SDSS PSF photometry          
             decals_to_maggies, cat, maggies, ivarmaggies, /sdss, $
               /aperture, apmaggies=apmaggies, apivarmaggies=apivarmaggies, $
               /decam_grz, /shortwise
             decals_to_maggies, cat, psfmaggies, psfivarmaggies, /sdss, /psf, $
               /decam_grz, /shortwise
             
             for ib = 0, nband-1 do begin
                all = where(strtrim(cat.type,2) eq 'PSF' and cat.decam_saturated eq 'F' and $
                  maggies[dindx[ib],*] gt 0 and apmaggies[5,dindx[ib],*] gt 0 and $
                  psfmaggies[sindx[ib],*] gt 0,nall)
                these = where(strtrim(cat.type,2) eq 'PSF' and cat.decam_fracflux[dindx[ib]] lt 0.2 and $
                  cat.decam_saturated eq 'F' and maggies[dindx[ib],*] gt 10^(-0.4*magcut[ib]) and $
                  apmaggies[5,dindx[ib],*] gt 0 and psfmaggies[sindx[ib],*] gt 0,nmatch)

                if nmatch gt 40 then begin
;               if nmatch gt 100 and ib lt 2 then begin
                   dmtag = tag_indx(info1[0],'dm_'+band[ib])
                   info1[ic].(dmtag) = djs_median(-2.5*alog10($
                     apmaggies[5,dindx[ib],these]/psfmaggies[sindx[ib],these]))
                   if keyword_set(debug) then begin
                      djs_plot, -2.5*alog10(maggies[dindx[ib],all]), $
                        -2.5*alog10(apmaggies[5,dindx[ib],all]/psfmaggies[sindx[ib],all]), $
                        psym=8, symsize=0.5, xsty=3, ysty=3, yr=[-3,3], xrange=[15,26]
                      djs_oplot, -2.5*alog10(maggies[dindx[ib],these]), $
                        -2.5*alog10(apmaggies[5,dindx[ib],these]/psfmaggies[sindx[ib],these]), $
                        psym=6, color='orange'
                      djs_oplot, !x.crange, [0,0], line=0
                      djs_oplot, !x.crange, info1[ic].(dmtag)*[1,1], line=5
                      im_legend, [info1[ic].brickname,band[ib]], /left, /top, box=0, charsize=2
                      cc = get_kbrd(1)
                   endif 
                endif 
             endfor 
          endfor 
          if n_elements(info) eq 0L then info = info1 else $
            info = [temporary(info),info1]
       endfor 
stop    
       im_mwrfits, info, infofile, /clobber
    endif
       
; ---------------------------------------------------------------------------
; build some QAplots
    if keyword_set(qaplots) then begin

       info = mrdfits(infofile+'.gz',1)

; #########################
; delta-mag 6-panel plots
       dmrange = [-0.25,0.15]
       rarange = [-4,364]
       decrange = [3,35]
       ytitle = '\Delta'+'m (DECaLS 7" diameter aperture minus SDSS PSF, mag)'

       psfile = dr1path+'qa-dr1-dmstars.ps'
       im_plotconfig, 7, pos, psfile=psfile, xmargin=[1.3,0.2], $
         yspace=0.05, height=[2.8,2.8,2.8]

; g
       djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
         xrange=rarange, yrange=dmrange, xtickname=replicate(' ',10)
       good = where(info.dm_g gt -90)
       djs_oplot, info[good].ra, info[good].dm_g, psym=symcat(16), $
         color=cgcolor('dodger blue'), symsize=0.8
       djs_oplot, !x.crange, [0,0], line=0, color=cgcolor('grey')
       im_legend, 'g', /left, /top, box=0, margin=0
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
         xrange=decrange, yrange=dmrange, xtickname=replicate(' ',10), $
         ytickname=replicate(' ',10)
       djs_oplot, info[good].dec, info[good].dm_g, psym=symcat(16), $
         color=cgcolor('dodger blue'), symsize=0.8
       djs_oplot, !x.crange, [0,0], line=0, color=cgcolor('grey')

; r
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,2], xsty=1, ysty=1, $
         xrange=rarange, yrange=dmrange, xtickname=replicate(' ',10), $
         ytitle=ytitle
       good = where(info.dm_r gt -90)
       djs_oplot, info[good].ra, info[good].dm_r, psym=symcat(16), $
         color=cgcolor('forest green'), symsize=0.8
       djs_oplot, !x.crange, [0,0], line=0, color=cgcolor('grey')
       im_legend, 'r', /left, /top, box=0, margin=0
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,3], xsty=1, ysty=1, $
         xrange=decrange, yrange=dmrange, xtickname=replicate(' ',10), $
         ytickname=replicate(' ',10)
       djs_oplot, info[good].dec, info[good].dm_r, psym=symcat(16), $
         color=cgcolor('forest green'), symsize=0.8
       djs_oplot, !x.crange, [0,0], line=0, color=cgcolor('grey')

; z
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,4], xsty=1, ysty=1, $
         xrange=rarange, yrange=dmrange, xtitle='RA (degrees)'
       good = where(info.dm_z gt -90)
       djs_oplot, info[good].ra, info[good].dm_z, psym=symcat(16), $
         color=cgcolor('tomato'), symsize=0.6
       djs_oplot, !x.crange, [0,0], line=0, color=cgcolor('grey')
       im_legend, 'z', /left, /top, box=0, margin=0
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,5], xsty=1, ysty=1, $
         xrange=decrange, yrange=dmrange, xtitle='Dec (degrees)', $
         ytickname=replicate(' ',10)
       djs_oplot, info[good].dec, info[good].dm_z, psym=symcat(16), $
         color=cgcolor('tomato'), symsize=0.6
       djs_oplot, !x.crange, [0,0], line=0, color=cgcolor('grey')
       
       im_plotconfig, psfile=psfile, /psclose, /pdf

       
       

stop    
       ploterror, info.ra, info.dra_med, info.dra_sig, psym=6, xsty=3, ysty=3, /trad    

    endif
    

stop    

return
end
