pro zprimus_qaplot_extract, makefits=makefits, makeplot=makeplot
; jm08jun24nyu - compare the A&B object and sky spectra of PRIMUS
;                objects with and without successful redshifts 

    rerun = '1500'
    
    caliblist = getenv('PRIMUS_DIR')+'/pro/reduction_scripts/calib_masks_rerun11.list'
    readcol, caliblist, night, allmask, format='A,A', comment='#', /silent
    keep = where(strmatch(allmask,'*vvds*',/fold) eq 0B) ; exclude VVDS masks
    allmask = allmask[keep] 
    nmask = n_elements(allmask)

    fitsfile = im_primus_path()+'qa_zprimus_calib_masks_extract.fits'
    psfile = im_primus_path()+'qa_zprimus_calib_masks_extract.ps'

    magcut = 23.5
    catacut = 0.05

    if keyword_set(makefits) then begin
       
       splog, 'Generating FITS file'

       for imask = 0L, nmask-1L do begin
;      for imask = 0L, 1L do begin

          primus_read_1dinputs, allmask[imask], rerun, $
            extract=ext1, slits=slits1, oned=oned1

; grab good object spectra and distinguish between those with
; good/crappy redshifts
          
          these = where((slits1.ztype eq 0) or (slits1.ztype eq 1) or $
            (slits1.ztype eq 2) or (slits1.ztype eq 6) and (slits1.currz gt 0.0) and $
            (slits1.mag lt magcut) and $
            (total((ext1.wave1 gt 6000) and (ext1.wave1 lt 10000),1) gt 10 and $
            total((ext1.wave2 gt 6000) and (ext1.wave2 lt 10000),1) gt 10) and $
            (ext1.badextract eq 0),nthese)
;           ((ext1.badextract and primus_flagval('SPEXTRACT', 'A_FULLMASK')) eq 0),nthese)
          if (nthese eq 0L) then message, 'This cannot be good!'
          ext = ext1[these] & slits = slits1[these] & oned = oned1[these]
          moretags = replicate({maskname: allmask[imask], mag: 0.0},nthese)
          moretags.mag = slits.mag

          good = where((abs(oned.zmin_gal[0]-slits.currz)/$
            (1.0+slits.currz) lt catacut),ngood,comp=crap,$
            ncomp=ncrap)

          splog, 'Comparing A-B object and sky extractions'
          if (ngood ne 0L) then begin
             objgood1 = struct_addtags(qaplot_extract_work($
               ext[good],sky=0),moretags[good])
             skygood1 = struct_addtags(qaplot_extract_work($
               ext[good],/sky),moretags[good])
          endif
          if (ncrap ne 0L) then begin
             objcrap1 = struct_addtags(qaplot_extract_work($
               ext[crap],sky=0),moretags[crap])
             skycrap1 = struct_addtags(qaplot_extract_work($
               ext[crap],/sky),moretags[crap])
          endif

          if (n_elements(objgood) eq 0L) and $
            (n_elements(skygood) eq 0L) then begin
             objgood = objgood1
             skygood = skygood1
          endif else begin
             if (ngood ne 0L) then begin
                objgood = [objgood,objgood1]
                skygood = [skygood,skygood1]
             endif
          endelse
          if (n_elements(objcrap) eq 0L) and $
            (n_elements(skycrap) eq 0L) then begin
             objcrap = objcrap1
             skycrap = skycrap1
          endif else begin
             if (ncrap ne 0L) then begin
                objcrap = [objcrap,objcrap1]
                skycrap = [skycrap,skycrap1]
             endif
          endelse
       endfor

; write out

       splog, 'Writing '+fitsfile
       mwrfits, objgood, fitsfile, /create
       mwrfits, skygood, fitsfile
       mwrfits, objcrap, fitsfile
       mwrfits, skycrap, fitsfile
       
    endif
       
; now make the plot

    if keyword_set(makeplot) then begin

       splog, 'Generating QA plot'
       arm_plotconfig, nx=4L, ny=4L, xmargin=[0.9,0.2], ymargin=[0.8,0.8], $
         xspace=0.1, yspace=[0.0,0.8,0.0], /landscape, coords=pos, $
         /normal, /write, xpage=xpage, ypage=ypage, psfile=psfile
       cleanplot, /silent

       splog, 'Reading '+fitsfile
       objgood = mrdfits(fitsfile,1)
       skygood = mrdfits(fitsfile,2)
       objcrap = mrdfits(fitsfile,3)
       skycrap = mrdfits(fitsfile,4)

       nccd = 8L
       ccdarray = [1,2,3,4,6,5,8,7] ; [see http://howdy.physics.nyu.edu/index.php/IMACS]
       minwave = 4200.0
       maxwave = 1E4
       binwave = 200.0
          
       for kk = 0L, 1L do begin ; object then sky
          if (kk eq 0L) then begin
             resgood = objgood
             rescrap = objcrap
             prefix = 'Object'
             yrange = 249.9*[-1,1]
          endif else begin
             resgood = skygood
             rescrap = skycrap
             prefix = 'Sky'
             yrange = 39.9*[-1,1]
          endelse
          nmask = n_elements(uniq(resgood.maskname))
          for iccd = 0L, nccd-1L do begin
             for ii = 0L, 1L do begin ; 0=z-good, 1=z-bad
                jj = iccd+ii*nccd/2+(iccd gt (nccd/2-1))*nccd/2
                if (ii eq 0L) then begin
                   result = resgood ; z-good
                   thislabel = '\delta'+'z<'+string(catacut,format='(F4.2)')
;                  thislabel = '\delta'+'z/(1+z)<'+string(catacut,format='(F4.2)')
                   thistitle = 'CCD '+string(ccdarray[iccd],format='(I0)')
                   xtickname = replicate(' ',10)
                   delvarx, xtitle
                endif else begin
                   result = rescrap ; z-crap
                   thislabel = '\delta'+'z>'+string(catacut,format='(F4.2)')
;                  thislabel = '\delta'+'z/(1+z)>'+string(catacut,format='(F4.2)')
                   thistitle = ''
                   delvarx, xtickname
                   xtitle = textoidl('Wavelength (\AA)')
                endelse
                if ((iccd+1L) eq 1L) or ((iccd+1L) eq 5L) then begin
                   delvarx, ytickname 
                   ytitle = ''
                endif else begin
                   ytickname = replicate(' ',10)
                   delvarx, ytitle
                endelse
                these = where(result.ccdnum eq ccdarray[iccd],nobj)
                stats = im_medxbin(reform(result[these].wave),$
                  weight=reform(result[these].mask)*1.0,$
                  reform(result[these].resid),binwave)
                plotsym, 8, 0.2, /fill
                hogg_scatterplot, result[these].wave, result[these].resid, $
                  weight=result[these].mask*1.0, $ darkest=30.0, $
                  xthick=4.0, ythick=4.0, charthick=3.0, $ ; xvec=xvec, yvec=yvec, $
                  charsize=1.1, xsty=1, ysty=1, xrange=[3999.0,9999.0], $
                  yrange=yrange, position=pos[*,jj], noerase=(jj gt 0L), $
                  title=thistitle, xtickinterval=2E3, xtickname=xtickname, $
                  ytickname=ytickname, xtitle=xtitle, ytitle=ytitle, $
                  /nocontour, /outliers, /internal_weight, outpsym=8, outsymsize=1.0, $
                  levels=errorf(0.5*[1.0,2.0,3.0]), outcolor=djs_icolor('default')
                djs_oplot, !x.crange, [0,0], line=0, thick=2.0
                djs_oplot, stats.binctr, stats.medy, color='red', thick=4.0, line=0
                djs_oplot, stats.binctr, stats.sigy75, color='red', thick=4.0, line=2
                djs_oplot, stats.binctr, stats.sigy25, color='red', thick=4.0, line=2

                legend, textoidl(thislabel), charsize=1.1, charthick=3.0, $
                  /right, /top, box=0, /clear, margin=0
                legend, 'N='+string(nobj,format='(I0)'), charsize=1.1, charthick=3.0, $
                  /right, /bottom, box=0, /clear, margin=0
             endfor
          endfor 

; y-title       
          xyouts, pos[0,0]-0.05, (pos[3,0]-pos[1,4])/2.0+pos[1,4], $
            prefix+' (A-B)/[0.5*(A+B)] (%)', align=0.5, charsize=1.1, $
            charthick=3.0, /normal, orientation=90.0
          xyouts, pos[0,0]-0.05, (pos[3,8]-pos[1,12])/2.0+pos[1,12], $
            prefix+' (A-B)/[0.5*(A+B)] (%)', align=0.5, charsize=1.1, $
            charthick=3.0, /normal, orientation=90.0

; big title       
          xyouts, (pos[2,3]-pos[0,0])/2.0+pos[0,0], 0.95, 'rerun '+rerun+' - '+$
            string(nmask,format='(I0)')+' calibration masks', $
            align=0.5, charsize=1.7, charthick=3.0, /normal
;         xyouts, (pos[2,3]-pos[0,0])/2.0+pos[0,0], 0.95, rerun+'/'+allmask[imask], $
;           align=0.5, charsize=1.7, charthick=3.0, /normal
          
       endfor 

       dfpsclose
       spawn, 'gzip -f '+psfile

    endif
       
stop       
       
return
end
    
