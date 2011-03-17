pro alpha_arctests, testarclist=testarclist, check_arclines=check_arclines
; jm09jan06nyu - excised from original code
    
;   datapath = getenv('DEEP2_ALPHA_DIR')+'/ut080414/'
    linelistpath = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/'
    thislinlist = linelistpath+'mike_thar_alpha_custom.lst'
    mikelinlist = linelistpath+'mike_thar.lst'

    setup = 1L
    if keyword_set(blue) then side = 1L else side = 2L

    ipsig = [[0,50,10],[50,200,15]]
    ifsig = [[0,50,5],[50,200,10]]

    rr = mike_allarc_sngl('Raw/r0001.fits',setup,side,/clobber,ifsig=ifsig,ipsig=ipsig)
    x_tweakarc, 'Arcs/Fits/mr0001_fit.idl', 42, linlist=mikelinlist, qafil='test.ps', $
      '/Users/ioannis/local/idl/MIKE/pro/Arcs/templ_arc_2x2R.idl', $
      ostr_fil=mike_getfil('ordr_str',setup,SIDE=side,/name)
   
; ---------------------------------------------------------------------------
; some tests to decide which line-list to use; note that I had to make
; some hard-copies of arc r0032 so that the script could be run with
; the various line-lists, thusly:
;    r0032 --> mike_thar.lst
;    r9032 --> mike_thar_alpha_custom.lst
;    r9132 --> mike_thar_alpha_murphy.lst
;    r9232 --> mike_thar_alpha_custom_shift.lst

    if keyword_set(testarclist) then begin
       doallarcs = 1L
       if keyword_set(doallarcs) then begin
          arcs = where(mike.type EQ 'ARC' AND mike.flg_anly NE 0 $
            AND mike.setup EQ setup AND mike.side EQ side,narc)
          arc_fil = 'Arcs/Arc_m'+strtrim(mike[arcs].img_root,2)
          for iarc = 0L, narc-1L do begin
             res = mike_fitarc_work(arc_fil[iarc],setup,side,$
               linlist=thislinlist,clobber=clobber,/noextrap)
          endfor
       endif else begin
          res1 = mike_fitarc_work('Arcs/Arc_mr0032.fits',setup,side,$
            linlist=linelistpath+'mike_thar_alpha_custom.lst',$
            clobber=clobber,/noextrap) ;,iordr=55,fordr=52);,inter=0
          res1 = mike_allarc_sngl('Raw/r0032.fits',setup,side,$
            linlist=linelistpath+'mike_thar.lst',$
            clobber=clobber)
          res2 = mike_allarc_sngl('Raw/r9032.fits',setup,side,$
            linlist=linelistpath+'mike_thar_alpha_custom.lst',$
            clobber=clobber)
          res3 = mike_allarc_sngl('Raw/r9132.fits',setup,side,$
            linlist=linelistpath+'mike_thar_alpha_murphy.lst',$
            clobber=clobber)
          res4 = mike_allarc_sngl('Raw/r9232.fits',setup,side,$
            linlist=linelistpath+'mike_thar_alpha_custom_shift.lst',$
            clobber=clobber)
       endelse
    endif

; ---------------------------------------------------------------------------
; write out QA files for the arcs    
    
    if keyword_set(check_arclines) then begin
       if (file_test('check_arclines',/dir) eq 0L) then begin
          splog, 'Please make check_arclines directory'
          return
       endif
       ordr_str = mike_getfil('ordr_str',setup,side=side)
       doallarcs = 1L
       if keyword_set(doallarcs) then begin
          arcs = where(mike.type EQ 'ARC' AND mike.flg_anly NE 0 $
            AND mike.setup EQ setup AND mike.side EQ side,narc)
          thesearcs = strtrim(mike[arcs].img_root,2)
          linlist = replicate(thislinlist,narc)
;         linlist = linelistpath+replicate('mike_thar.lst',narc)
       endif else begin
          thesearcs = 'r'+['0032','9032','9132','9232']+'.fits'
          narc = n_elements(thesearcs)
          match, mike.img_root, thesearcs, arcs, m2
          linlist = linelistpath+'mike_thar'+$
            ['','_alpha_custom','_alpha_murphy',$
            '_alpha_custom_shift']+'.lst'
       endelse
       for iarc = 0L, narc-1L do begin
; generate the arc name and the IDL save set name from the raw arcfile
          arc_root = 'Raw/'+thesearcs[iarc]
          arc_fil = 'Arcs/Arc_m'+thesearcs[iarc]
          arc_saveset = 'Arcs/Fits/m'+repstr(thesearcs[iarc],'.fits','_fit.idl')
;         arc_root = strcompress(mike[arcs[iarc]].rootpth+mike[arcs[iarc]].img_root,/remove)
;         arc_fil = mike_getfil('arc_fil',subfil=arc_root,/name)
;         arc_saveset = mike_getfil('arc_fit',subfil=arc_fil,/name)
; restore the IDL save set and read the line-list used in X_ARCFIT
          restore, arc_saveset
          x_arclist, linlist[iarc], lines
          lines = lines[sort(lines.wave)]
; open the output file and write out (see X_ARCFIT)
          outfil = 'check_arclines/'+'check_'+repstr($
            file_basename(arc_root),'.fits','.dat')
          splog, 'Writing '+outfil
          openw, lun, outfil, /get_lun
          printf, lun, '# '+file_basename(arc_fil)
          printf, lun, '# order, pixel, true lambda, fitted lambda, '+$
            'flag [1=good; 0=rejected]'
          for ii = 0L, n_elements(ordr_str)-1L do begin
             if rejstr[ii].ngdf EQ 0L then continue
             these = lindgen(rejstr[ii].ngdf) ; lines used
             gdfit = rejstr[ii].gdfpt[these]
             if (rejstr[ii].nrej NE 0L) then $
               rejpt = rejstr[ii].rejpt[0:rejstr[ii].nrej-1] else $
               rejpt = -1L
             fit = 10^x_calcfit(double(rejstr[ii].gdfpx[these]),FITSTR=all_arcfit[ii])
             goodbad = these*0B+1B ; 1=good, 0=bad
             if (rejpt[0] ne -1L) then goodbad[rejpt] = 0B
             for jj = 0L, n_elements(gdfit)-1L do $
               printf, lun, ordr_str[ii].order, rejstr[ii].gdfpx[these[jj]], $
               lines[gdfit[jj]].wave, fit[jj], goodbad[jj]
          endfor
          free_lun, lun
       endfor
       spawn, 'mkdir -p /tmp/junk ; cp -f QA/Arcs01/qa_arcfit_mr????.ps.gz '+$
         'check_arclines/check_r????.dat /tmp/junk ; pushd /tmp/junk ; '+$
         'tar czvf '+datapath+'alpha_linetests_08jun11.tar.gz *.ps.gz *.dat ; popd', /sh
;      spawn, 'mkdir -p /tmp/junk ; cp -f QA/Arcs01/qa_arcfit_mr[0,9]?32.ps.gz '+$
;        'check_arclines/check_r[0,9]?32.dat /tmp/junk ; pushd /tmp/junk ; '+$
;        'tar czvf '+datapath+'alpha_linetests_08jun10.tar.gz *.ps.gz *.dat ; popd', /sh
    endif
    
return
end
    
