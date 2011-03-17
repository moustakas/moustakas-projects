pro oplot_gcrgb, imags, color, infobase, datapath, trgbmag, halo=halo, core=core
; jm00aug6ucb
; overplot the globular cluster rgb's on the color-magnitude diagram

        plot, color, imags, ps=3, color=7, xsty=1, ysty=1, syms=1.2, $
          xtit = 'V-I', ytit='I', xr=[-1,4], yr=[max(imags),min(imags)], $ ; yr = [27.2,20.6], $
          title = infobase.truename+ ' CMD'

        okay = 'N'
        read, okay, prompt='Write the CMD to postscript (Y/[N])? ' & print
        if strupcase(okay) eq 'Y' then begin
            path = trgb_datapath()
            path = path[3]+'CMDs/' ; plot subdirectory
            if keyword_set(halo) then ps_open, path+infobase.object+'_halo_cmd', /ps_fonts else $
              if keyword_set(core) then ps_open, path+infobase.object+'_core_cmd', /ps_fonts else $
              ps_open, path+infobase.object+'_cmd', /ps_fonts
            device, /times
            plot, color, imags, ps=3, color=7, xsty=1, ysty=1, syms=1.2, $
              xtit = 'V-I', ytit='I', xr=[-1,4], yr=[max(imags),min(imags)], $ ; yr = [27.2,20.6], $
              title = infobase.truename+ ' CMD'
            ps_close
        endif

        print, 'As a diagnostic for finding the TRGB, click on the apparent TRGB '
        print, 'magnitude to overplot two globular cluster isochrones, one of M15'
        print, '(yellow, metal poor) and one of 47 Tuc (red, metal rich): '
        repeat begin
            cursor, x, trgbmag, /down, /data 
            mbutton = !mouse.button
        endrep until (mbutton eq 1L or mbutton eq 2L)
        oplot, [!x.crange[0],!x.crange[1]], [trgbmag,trgbmag], line=2, color=3, thick=2.5
        xyouts, [-0.3,-0.3], [trgbmag-0.2,trgbmag-0.2], strn(trgbmag,format='(F5.2)')+' mag', $
          charthick=2, charsize=1.2, color=1
        
	trgb_dacosta, gcrgb     ; globular cluster rgb
        trgbcolor = -4.06	; absolute TRGB color

        nozero = where(gcrgb[0].imag ne float(0))
        
        oplot, gcrgb[0].color, gcrgb[0].imag[nozero]+(trgbmag-trgbcolor), color=5, line=0, thick=3
        oplot, gcrgb[1].color, gcrgb[1].imag[nozero]+(trgbmag-trgbcolor), color=2, line=0, thick=3

return
end


pro trgb_colorcut, objname, hst=hst, halo=halo, core=core
;+
; NAME:
;	TRGB_COLORCUT
;
; PURPOSE:
;	Display an object's color-magnitude diagram and isolate
;	particular stellar populations interactively.
;
; CALLING SEQUENCE:
;	trgb_colorcut, objname
;
; INPUTS:
;	objname	: galaxy name (string)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;	The user can choose to write a file called
;	objname_starlist_ccut.dat of the photometry to be kept.
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;	COLORTABLE1, TRGB_READATA, TRGB_DATAPATH()
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 27, UCB
;
;-

        !x.ticklen=0.03
        p_thick = !p.thick & chars = !p.charsize
        charth = !p.charthick & xth = !x.thick & yth = !y.thick
        !p.thick = 2 & !p.charsize=1.8 & !p.charthick=2 & !x.thick=2 & !y.thick=2

        colortable1
        window, 0, xs=500, ys=500

        trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core
        imags = data[3,*] - infobase.a_i   ; apply extinction correction
        color = data[7,*] - infobase.e_v_i ; color excess correction
        nstars = n_elements(data[0,*])

        allgood = 0L
        repeat begin

            oplot_gcrgb, imags, color, infobase, datapath, trgbmag, halo=halo, core=core
            xyouts, [0.25,0.25], [0.83,0.83], strn(nstars)+' Stars', /normal, $
              charthick=2, charsize=1.5, color=1
        
            print & print, 'Click on the blue side of the color cut: '
            repeat begin
                cursor, colorstart, y, /down, /data 
                mbutton = !mouse.button
            endrep until (mbutton eq 1L or mbutton eq 2L)
            oplot, [colorstart,colorstart], [!y.crange[0],!y.crange[1]], line=2, color=2, thick=2.5
            xyouts, [colorstart-0.2,colorstart-0.2], [!y.crange[0]-0.3,!y.crange[0]-0.3], $
              strn(colorstart,format='(F5.2)')+' mag', charthick=2, charsize=1.2, color=1, align=1.0
            
            print & print, 'Now click on the red side of the color cut: '
            repeat begin
                cursor, colorend, y, /down, /data 
                mbutton = !mouse.button
            endrep until (mbutton eq 1L or mbutton eq 2L)
            oplot, [colorend,colorend], [!y.crange[0],!y.crange[1]], line=2, color=2, thick=2.5
            xyouts, [colorend+0.2,colorend+0.2], [!y.crange[0]-0.3,!y.crange[0]-0.3], $
              strn(colorend,format='(F5.2)')+' mag', charthick=2, charsize=1.2, color=1, align=0.0
            
            rgbstars = where((color gt colorstart) and (color lt colorend),rgbcount)
            xyouts, [0.25,0.25], [0.78,0.78], strn(rgbcount)+' Stars selected', /normal, $
              charthick=2, charsize=1.5, color=1

            okay = '' & print
            read, okay, prompt = 'Press ENTER to continue or type any key to redraw the color lines: '
            if okay eq '' then allgood = 1L

        endrep until allgood

; write the photometry to a file

        ccutdata = data[*,rgbstars]

        okayphot = 'N' & print
        read, okayphot, prompt='Write the photometry to a file? (Y/[N])? '
        if strupcase(okayphot) eq 'Y' then begin

            if keyword_set(hst) then begin
                if keyword_set(halo) then ccutfile = datapath+'/'+objname+'_starlist_halo_ccut.dat' else $
                  if keyword_set(core) then ccutfile = datapath+'/'+objname+'_starlist_core_ccut.dat' else $
                  ccutfile = datapath+'/'+objname+'_starlist_ccut.dat'
            endif else begin
                if keyword_set(halo) then ccutfile = datapath+'/'+objname+'_IVmags_halo_ccut.dat' else $
                  if keyword_set(core) then ccutfile = datapath+'/'+objname+'_IVmags_core_ccut.dat' else $
                  ccutfile = datapath+'/'+objname+'_IVmags_ccut.dat'
            endelse
                
            print
            print, 'Writing '+ccutfile+' . . . '

            openw, lun1, ccutfile, /get_lun
            printf, lun1, '#  ID   Xcenter   Ycenter    I       Ierr      V       Verr     V-I'
            printf, lun1, ' '
            for k = 0L, rgbcount-1L do printf, lun1, ccutdata[*,k], format = '(1x,I6,7F9.3)'
            free_lun, lun1        

        endif

; generate a text file of results

        path = trgb_datapath()
        respath = path[4]          ; results subdirectory
        if keyword_set(halo) then resfile = respath+objname+'_halo_ccut.txt' else $
          if keyword_set(core) then resfile = respath+objname+'_core_ccut.txt' else $
          resfile = respath+objname+'_ccut.txt'

        okay = 'N' & print
        read, okay, prompt='Write the color-cut values to a data file (Y/[N])? '
        if strupcase(okay) eq 'Y' then begin
            
            print & print, 'Writing '+resfile+' . . .'

            openw, lun1, resfile, /get_lun
            printf, lun1, '# Object         : ', infobase.truename
            printf, lun1, '# Date           : ', systime()
            printf, lun1, '# -------------------------------------------'
            if keyword_set(halo) then printf, lun1, '# Keywords Used  : Halo'
            if keyword_set(core) then printf, lun1, '# Keywords Used  : Core'
            if ((not keyword_set(halo)) and (not keyword_set(core))) then $
              printf, lun1, '# Keywords Used  : None'
            if strupcase(okayphot) eq 'Y' then printf, lun1, '# Photometry File: ', strn(ccutfile)
            printf, lun1, '# Blue Cut       : ', strn(colorstart,format='(F5.2)')+' mag'
            printf, lun1, '# Red Cut        : ', strn(colorend,format='(F5.2)')+' mag'
            printf, lun1, '# Apparent TRGB  : ', strn(trgbmag,format='(F5.2)')+' mag'
            printf, lun1, '# Stars          : ', strn(nstars)
            printf, lun1, '# Good Stars     : ', strn(rgbcount)
            printf, lun1, '# Rejected Stars : ', strn(nstars-rgbcount)

            free_lun, lun1

        endif
            
        print
            
        !x.ticklen=0
        !p.thick = p_thick & !p.charsize=chars & !p.charthick=charth
        !x.thick = xth & !y.thick=yth

return
end

