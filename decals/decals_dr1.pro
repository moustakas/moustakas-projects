pro decals_dr1
; jm15mar13siena - build QAplots for the DR1 test bricks

    common com_dr1, cat

; get the merged tractor catalog (the only cut on this is
; BRICK_PRIMARY='T') 
    qa_dir = getenv('IM_PROJECTS_DIR')+'/decals/dr1/qaplots/'
    isedfit_dir = getenv('IM_PROJECTS_DIR')+'/decals/isedfit/dr1/'
    if n_elements(cat) eq 0L then cat = $
      mrdfits(isedfit_dir+'tractor-dr1.fits.gz',1)

    decals_to_maggies, cat, maggies, ivarmaggies, /sdss
    decals_to_maggies, cat, psfmaggies, psfivarmaggies, /sdss, /psf
    fix = where(finite(maggies) eq 0)
    if fix[0] ne -1L then begin
       maggies[fix] = 0
       ivarmaggies[fix] = 0
    endif
    fix = where(finite(psfmaggies) eq 0)
    if fix[0] ne -1L then begin
       psfmaggies[fix] = 0
       psfivarmaggies[fix] = 0
    endif

    mag = maggies2mag(maggies,ivarmaggies=ivarmaggies)
    psfmag = maggies2mag(psfmaggies,ivarmaggies=psfivarmaggies)

; compare SDSS/DECaLS grz fluxes
    band = ['g','r','z']
    dindx = [0,1,2]
    sindx = [6,7,9]
;   maxx = [23.5,23.5,22]
    maxx = [24,24,24]
    depth = [24.0,23.6,23.0]

    type = ['stars','galaxies']
    for tt = 0, n_elements(type)-1 do begin
       psfile = qa_dir+'qa_dr1_sdss_'+type[tt]+'.ps'
       im_plotconfig, 4, pos, psfile=psfile, xmargin=[1.2,0.3]
       for ii = 0, 2 do begin
          case tt of
             0: begin
                these = where(psfmag[dindx[ii],*] gt 0 and psfmag[sindx[ii],*] gt 0 and $
                  cat.decam_saturated ne 'T' and strtrim(cat.type,2) eq 'PSF',nobj) ; stars
                mm = psfmag[dindx[ii],these]
                dm = mm-psfmag[sindx[ii],these]
             end
             1: begin
                these = where(mag[dindx[ii],*] gt 0 and mag[sindx[ii],*] gt 0 and $
                  cat.decam_saturated ne 'T' and strtrim(cat.type,2) ne 'PSF',nobj) ; galaxies
                mm = mag[dindx[ii],these]
                dm = mm-mag[sindx[ii],these]
             end
          endcase
          med = im_medxbin(mm,dm,0.25,minx=15.0,maxx=maxx[ii],minpts=100,verbose=tt eq 0)
          if ii eq 2 then delvarx, xtickname else xtickname = replicate(' ',10)
          
          im_hogg_scatterplot, mm, dm, /outliers, /internal, xsty=1, ysty=1, $
            position=pos[*,ii], noerase=ii gt 0, xrange=[14,25], yrange=1.3*[-1,1], $
            outcolor=cgcolor('grey'), levels=[0.5,0.75,0.9,0.95,0.975,0.99], /nogrey, $
            xnpix=50, ynpix=50, contour_color=cgcolor('grey'), xtickname=xtickname
          djs_oplot, !x.crange, [0,0], line=0 ;, color=cgcolor('grey')
          djs_oplot, depth[ii]*[1,1], !y.crange, line=0 ;, color=cgcolor('grey')
          oploterror, med.medx, med.medy, med.quant75-med.medy, /hibar, $
            color=cgcolor('firebrick'), errcolor=cgcolor('firebrick'), psym=symcat(16)
          oploterror, med.medx, med.medy, med.medy-med.quant25, /lobar, $
            color=cgcolor('firebrick'), errcolor=cgcolor('firebrick'), psym=symcat(16)
          im_legend, band[ii], /left, /top, box=0, charsize=2, margin=0
       endfor
       xyouts, pos[0,1]-0.09, mean([pos[1,1],pos[3,1]]), $
         textoidl('\Delta'+'m (DECaLS - SDSS, mag)'), $
         /normal, orientation=90, align=0.5
       xyouts, mean([pos[0,2],pos[2,2]]), pos[1,2]-0.07, 'DECaLS (mag)', $
         /normal, align=0.5
       im_plotconfig, psfile=psfile, /psclose, /png
    endfor

stop    

return
end
