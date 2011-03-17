pro ages_plot_bad_redshifts
; jm06jan20uofa - given the potentially wrongly classified galaxies
;                 from AGES_FLAG_BAD_REDSHIFTS, use this routine
;                 to look at all the relevant spectra and reclassify
;                 them

    analysis_path = ages_path(/analysis)
    info = mrdfits(analysis_path+'ages_flag_bad_redshifts.fits.gz',1,/silent)
    
; now check each of the identified spectra and flag misclassified ones

    nflag = n_elements(info)
    flag = lindgen(nflag)
;   flag = where(info.z lt 1.0,nflag)

;   for iflag = 73L, nflag-1L do begin
    for iflag = 0L, nflag-1L do begin

       title = string(info[flag[iflag]].pass,format='(I3.3)')+'/'+$
         string(info[flag[iflag]].aper,format='(I3.3)')+' z = '+$
         strtrim(string(info[flag[iflag]].z,format='(F12.4)'),2)+' '+$
         strtrim(info[flag[iflag]].class,2)

       struct_print, struct_trimtags(info[flag[iflag]],except=['WAVE','FLUX']), /no_head
       plot, info[flag[iflag]].wave, smooth(info[flag[iflag]].flux,3), ps=10, xsty=3, ysty=3, $
         xthick=2.0, ythick=2.0, charsize=1.8, charthick=2.0, title=title, $
         xrange=[min(info[flag[iflag]].wave),8500]
       djs_oplot, 1215*[1,1]*(1.0+info[flag[iflag]].z), !y.crange, line=1, color='red', thick=2
       djs_oplot, 1550*[1,1]*(1.0+info[flag[iflag]].z), !y.crange, line=1, color='red', thick=2
       djs_oplot, 2800*[1,1]*(1.0+info[flag[iflag]].z), !y.crange, line=1, color='red', thick=2
       djs_oplot, 3727*[1,1]*(1.0+info[flag[iflag]].z), !y.crange, line=1, color='blue', thick=2
       djs_oplot, 5007*[1,1]*(1.0+info[flag[iflag]].z), !y.crange, line=1, color='green', thick=2
       cc = get_kbrd(1)

    endfor
    
stop
    
return
end
    
