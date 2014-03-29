pro noao_decam_14apr
; jm14mar21siena - make some figures for the NOAO Large Survey
; observing proposal

    outpath = getenv('IM_PROJECTS_DIR')+'/desi/proposals/'
    photo = rsex(outpath+'photo_surveys.txt')
    nphoto = n_elements(photo)

    specz = rsex(outpath+'specz_surveys.txt')
    nspecz = n_elements(specz)
    specz = struct_addtags(specz,replicate({vol: 0.0},nspecz))
    for ii = 0, nspecz-1 do begin
       if specz[ii].zmin gt 0 then begin
          specz[ii].vol = jhnvol(specz[ii].zmin,specz[ii].zmax)*$
            specz[ii].area*3600.0/1D9 ; [Gpc^3]
       endif
    endfor
    
; ---------------------------------------------------------------------------
; photometric/spectroscopic depth vs areal coverage (deg^2)
    these = [1,2,3,4,5]
    xoff = [0,0,0,0,0]
    yoff = [-0.3,-0.35,0.2,0.2,-0.45]
    nthese = n_elements(these)

    frac = photo[these].specz_area/photo[these].area
    symsize = (frac-min(frac))/max(frac)*3.0+1

    psfile = outpath+'area_vs_depth.eps'
    im_plotconfig, 0, pos, xmargin=[1.3,0.2], charsize=2.2, $
      charthick=4.0, psfile=psfile, height=5.0
    djs_plot, [0], [0], /nodata, xrange=[19.5,25.5], yrange=[1.5,5.5], $
      /xsty, /ysty, xtitle='z-band AB Magnitude Limit', $
      ytitle='Log(Survey Area / deg^{2})', position=pos

; photometric surveys    
    for ii = 0, nthese-1 do begin
       if strtrim(photo[these[ii]].survey,2) eq 'Legacy' then begin
          xyouts, photo[these[ii]].zlim+xoff[ii], alog10(photo[these[ii]].area)+yoff[ii], $
            'SDSS/DECam!c Legacy', charsize=2.0, align=0.5
          plots, photo[these[ii]].zlim, alog10(photo[these[ii]].area), $
            psym=symcat(15), symsize=symsize[ii], color=cgcolor('dodger blue')
          plots, photo[these[ii]].zlim, alog10(photo[these[ii]].area), $
            psym=symcat(6,thick=8.0), symsize=symsize[ii], color=cgcolor('navy')
       endif else begin
          xyouts, photo[these[ii]].zlim+xoff[ii], alog10(photo[these[ii]].area)+yoff[ii], $
            photo[these[ii]].survey, charsize=1.5, align=0.5
          plots, photo[these[ii]].zlim, alog10(photo[these[ii]].area), $
            psym=symcat(16), symsize=symsize[ii], color=cgcolor('grey')
          plots, photo[these[ii]].zlim, alog10(photo[these[ii]].area), $
            psym=symcat(9,thick=3.0), symsize=symsize[ii], color=cgcolor('black')
       endelse
    endfor

    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

; ---------------------------------------------------------------------------
; number of redshifts vs comoving volume
    these = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
;               0     1    2    3    4    5    6    7     8    9    10    11    12  13      15
    xoff = [ -0.6, -0.4,+0.5,-0.5, 0.0, 0.0, 0.0,-0.50,-0.40, 0.0, -1.0, -0.9,+0.40,-0.1,+0.80,+0.1]
    yoff = [-0.05,-0.05,-0.2,+0.2,-0.3,+0.2,-0.3,-0.05,+0.15,-0.3,+0.15,+0.15,+0.15,+0.1,-0.1,-0.25]
    symsize = (specz[these].zmedian-min(specz[these].zmedian))/$
      max(specz[these].zmedian)*2.5+1
    symsize = symsize*0+1.5

    psfile = outpath+'nz_vs_volume.eps'
    im_plotconfig, 0, pos, xmargin=[1.3,0.2], charsize=2.2, psfile=psfile, height=5.0
    djs_plot, [0], [0], /nodata, xrange=[-3.5,3.5], yrange=[3.5,8], $
      /xsty, /ysty, xtitle='Log(Comoving Volume / h_{70}^{-3} Gpc^{3})', $
      ytitle='Log(Number of Redshifts)', position=pos

    for ii = 0, n_elements(these)-1 do begin
       if strmatch(specz[these[ii]].survey,'*DESI*') then begin
          xyouts, alog10(specz[these[ii]].vol)+xoff[ii], alog10(specz[these[ii]].nz)+yoff[ii], $
            specz[these[ii]].survey, charsize=1.7, align=0.5
          plots, alog10(specz[these[ii]].vol), alog10(specz[these[ii]].nz), $
            psym=symcat(15), symsize=symsize[ii]*1.0, color=cgcolor('dodger blue')
          plots, alog10(specz[these[ii]].vol), alog10(specz[these[ii]].nz), $
            psym=symcat(6,thick=8.0), symsize=symsize[ii]*1.4, color=cgcolor('navy')
       endif else begin
          xyouts, alog10(specz[these[ii]].vol)+xoff[ii], alog10(specz[these[ii]].nz)+yoff[ii], $
            specz[these[ii]].survey, charsize=1.4, align=0.5
          plots, alog10(specz[these[ii]].vol), alog10(specz[these[ii]].nz), $
            psym=symcat(16), symsize=symsize[ii], color=cgcolor('grey')
          plots, alog10(specz[these[ii]].vol), alog10(specz[these[ii]].nz), $
            psym=symcat(9,thick=3.0), symsize=symsize[ii], color=cgcolor('black')
       endelse
    endfor

    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

    
stop    
    




    info = rsex(outpath+'survey_parameters.txt')
    nsurvey = n_elements(info)
    
    info = struct_addtags(info,replicate({xs_nz: 1.0, ys_nz: 1.0, $
      xs_depth: 1.0, ys_depth: 1.0, vol: 0.0},nsurvey))
    for ii = 0, nsurvey-1 do begin
       if info[ii].zmin gt 0 then begin
          info[ii].vol = jhnvol(info[ii].zmin,info[ii].zmax)*$
            info[ii].specz_area*3600.0/1D9 ; [Gpc^3]
       endif
       case strtrim(info[ii].survey,2) of
          'PRIMUS': begin
             info[ii].xs_nz = 0.4
             info[ii].ys_nz = 1.5
          end
          'AGES': begin
             info[ii].xs_nz = 2.9
             info[ii].ys_nz = 0.85
          end             
          'DEEP2': begin
             info[ii].xs_nz = 0.4
             info[ii].ys_nz = 1.5
          end             
          'SDSS-I/II': begin
             info[ii].xs_nz = 1.0
             info[ii].ys_nz = 1.5
          end             
          'zCOSMOS': begin
             info[ii].xs_nz = 0.25
             info[ii].ys_nz = 1.15
          end             
          'GAMA': begin
             info[ii].xs_nz = 3.0
             info[ii].ys_nz = 0.8
          end             
          '2dFGRS': begin
             info[ii].xs_nz = 4.0
             info[ii].ys_nz = 1.0
          end             
          'BOSS': begin
             info[ii].xs_nz = 1.0
             info[ii].ys_nz = 0.5
          end             
          'WiggleZ': begin
             info[ii].xs_nz = 4.0
             info[ii].ys_nz = 0.8
          end             
          'DES': begin
             info[ii].xs_depth = 1.0
             info[ii].ys_depth = 2.0
          end             
          'DESI-ELGs': begin
             info[ii].xs_nz = 0.15
             info[ii].ys_nz = 1.2
          end             
          'DESI-LRGs': begin
             info[ii].xs_nz = 0.15
             info[ii].ys_nz = 1.2
          end             
          'DESI-QSOs': begin
             info[ii].xs_nz = 1.7
             info[ii].ys_nz = 0.4
          end             
          else: 
       endcase
    endfor
;   struct_print, info

return
end
    


;; ---------------------------------------------------------------------------
;; collect info on various surveys
;    nsurvey = 17
;    info = replicate({survey: '', nz: 0.0, specz_area: 0.0, phot_area: 0.0, $
;      zminmax: [0.0,0.0], vol: 0.0, xs: 1.0, ys: 1.0, zmed: 0.0, $
;      specz_rlim: 0.0, phot_rlim: 0.0, phot_ilim: 0.0},nsurvey)
;
;    info[0].survey = 'PRIMUS'
;    info[0].xs = 0.4
;    info[0].ys = 1.5
;    info[0].nz = 1.4E5
;    info[0].specz_area = 10.0
;    info[0].zmed = 0.6
;    info[0].zminmax = [0.2,1.0]
;
;    info[1].survey = 'AGES'
;    info[1].xs = 2.9
;    info[1].ys = 0.85
;    info[1].nz = 1.4E4
;    info[1].specz_area = 7.8
;    info[1].zmed = 0.3
;    info[1].zminmax = [0.01,0.8]
;
;    info[2].survey = 'DEEP2'
;    info[2].xs = 0.4
;    info[2].ys = 1.5
;    info[2].nz = 4E4
;    info[2].specz_area = 3.0
;    info[2].zmed = 1.0
;    info[2].zminmax = [0.7,1.5]
;
;; last release of SDSS I/II: DR7 http://www.sdss.org/dr7/    
;    info[3].survey = 'SDSS I/II' ; Main
;    info[3].xs = 1.0
;    info[3].ys = 1.5
;    info[3].nz = 929555.0 ; just galaxies
;    info[3].zmed = 0.1
;    info[3].zminmax = [0.01,0.25]
;    info[3].specz_area = 9380 ; 7966.0
;    info[3].specz_rlim = 17.7
;    info[3].phot_area = 11663.0
;    info[3].phot_rlim = 22.2
;    info[3].phot_ilim = 21.3
;
;    info[4].survey = 'TKRS'
;    info[4].xs = 1.0
;    info[4].ys = 2.0
;    info[4].nz = 1440.0
;    info[4].zmed = 0.7
;    info[4].specz_area = (10.0*16.0)/3600.0
;    info[4].zminmax = [0.0,1.5]
;
;    info[5].survey = 'CNOC2'
;    info[5].xs = 1.0
;    info[5].ys = 0.35
;    info[5].nz = 2000.0
;    info[5].specz_area = 1400.0/3600.0
;    info[5].zmed = 0.3
;    info[5].zminmax = [0.12,0.55]
;
;;   info[6].survey = 'zCOSMOS'
;    info[6].survey = 'zCOSMOS-!cBright'
;    info[6].xs = 0.25
;    info[6].ys = 1.15
;    info[6].nz = 2E4
;    info[6].specz_area = 1.7
;    info[6].zmed = 0.6
;    info[6].zminmax = [0.1,1.2]
;
;    info[7].survey = 'zCOSMOS-!cDeep'
;    info[7].xs = 1.0
;    info[7].ys = 0.5
;    info[7].nz = 1E4
;    info[7].specz_area = 1.0
;    info[7].zmed = 2.1
;    info[7].zminmax = [1.4,3.0]
;
;    info[8].survey = 'COMBO17'
;    info[8].xs = 0.7
;    info[8].ys = 1.4
;    info[8].nz = 2.5E4
;    info[8].specz_area = 0.78
;    info[8].zmed = 0.6
;    info[8].zminmax = [0.2,1.1]
;
;    info[9].survey = '2dFGRS'
;    info[9].xs = 4.0
;    info[9].ys = 1.0
;    info[9].nz = 2.5E5
;    info[9].specz_area = 2000.0
;    info[9].zmed = 0.1
;    info[9].zminmax = [0.0,0.22]
;
;    info[10].survey = 'BOSS'
;    info[10].xs = 1.0
;    info[10].ys = 0.5
;    info[10].nz = 1.5E6
;    info[10].specz_area = 1E4
;    info[10].zmed = 0.3
;    info[10].zminmax = [0.1,0.7]
;
;    info[11].survey = 'DESI-ELGs'
;    info[11].xs = 0.15
;    info[11].ys = 1.2
;    info[11].nz = 23.8E6
;    info[11].specz_area = 14000.0
;    info[11].zmed = 1.1
;    info[11].zminmax = [0.6,1.6]
;
;    info[12].survey = 'DESI-LRGs'
;    info[12].xs = 0.15
;    info[12].ys = 1.2
;    info[12].nz = 4.2E6
;    info[12].specz_area = 14000.0
;    info[12].zmed = 0.7
;    info[12].zminmax = [0.4,1.0]
;
;    info[13].survey = 'DESI-QSOs'
;    info[13].xs = 1.7
;    info[13].ys = 0.4
;    info[13].nz = 1.4E6
;    info[13].specz_area = 14000.0
;    info[13].zmed = 1.5
;    info[13].zminmax = [0.9,2.2]
;
;    info[14].survey = 'GAMA'
;    info[14].xs = 3.0
;    info[14].ys = 0.8
;    info[14].nz = 79599/2.0*3 ; scaling to full 3-year survey (Baldry+10)
;    info[14].specz_area = 143 ; deg^2
;    info[14].zmed = 0.2
;    info[14].zminmax = [0.01,0.5]
;
;; http://wigglez.swin.edu.au/site/data.html; "There are 7 fields with
;; average area of 140 square degrees totaling ~1000 square
;; degrees. The 7 fields overlap with GALEX, SDSS, and RCS2 in 00, 01,
;; 03, 09, 11, 15, 22 hr fields, shown below."
;    info[15].survey = 'WiggleZ'
;    info[15].xs = 4.0
;    info[15].ys = 0.8
;    info[15].nz = 238770.0
;    info[15].zmed = 0.6
;    info[15].zminmax = [0.2,1.0]
;    info[15].specz_area = 1000.0
;    info[15].phot_rlim = 22.5
;
;; http://www.darkenergysurvey.org/survey/    
;    info[16].survey = 'DES' ; photometric
;    info[16].xs = 4.0
;    info[16].ys = 0.8
;    info[16].phot_area = 5000.0 
;    info[16].phot_rlim = 24.8
;    info[16].phot_ilim = 24.0

