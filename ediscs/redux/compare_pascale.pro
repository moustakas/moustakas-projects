pro compare_pascale, postscript=postscript
; jm04oct28uofa
; compare my flux-calibrated spectra with Pascale's

    rootpath = ediscs_path(/d2002)

    bigcluster = ['cl1040-1155','cl1054-1146','cl1054-1245','cl1216-1201','cl1232-1250']
    ncluster = n_elements(bigcluster)

    psname = 'compare_pascale.ps'
    if keyword_set(postscript) then begin
       dfpsplot, rootpath+psname, /color, /square
       postthick = 5.0
    endif else postthick = 2.0

    for icluster = 0L, ncluster-1L do begin

       cluster = bigcluster[icluster]
       
       datapath = rootpath+cluster+'/'
       fluxedpath = datapath+'fluxed/'

       pushd, fluxedpath
       melist = file_search('*.fits',count=mecount)
       pascalemelist = repstr(repstr(melist,'t_b_',''),'.ms','')
       popd

       pascalepath = ediscs_path()+'MPE_2002/MARS2002/FluxCalibrated/'

       pushd, pascalepath
       shelist2 = file_search(cluster+'_mask2/*.fits',count=shecount)
       shelist3 = file_search(cluster+'_mask3/*.fits',count=shecount)
       shelist4 = file_search(cluster+'_mask4/*.fits',count=shecount)
       bigshelist = [shelist2,shelist3,shelist4]
       nshelist = n_elements(bigshelist)
       shelist = strarr(nshelist)
       for j = 0L, nshelist-1L do shelist[j] = strmid(bigshelist[j],$
         strpos(bigshelist[j],'/')+1,strlen(bigshelist[j]))
       popd

       pagemaker, nx=1, ny=2, position=pos, /normal, yspace=0.0, $
         xmargin=1.2
       
       for i = 0L, mecount-1L do begin

          match = where(strtrim(shelist,2) eq strtrim(pascalemelist[i],2),nmatch)
          if (nmatch eq 1L) then begin

             me = rd1dspec(melist[i],datapath=fluxedpath)
             she = rd1dspec(bigshelist[match],datapath=pascalepath)

             title = repstr(repstr(shelist[match],'_',' '),'.fits','')

             djs_iterstat, me.spec, mask=memask, sigrej=2.0
             djs_iterstat, she.spec, mask=shemask, sigrej=2.0

             mew = where(memask)
             shew = where(shemask)
             
             yrange = 1D17*[0.0,max(me[mew].spec)>max(she[shew].spec)]*[1.0,1.2]
;            yrange = 1D17*[min(me[mew].spec<she[shew].spec),max(me[mew].spec)>max(she[shew].spec)]*[1.0,1.2]
             
             djs_plot, me.wave, 1D17*me.spec, ps=10, xsty=3, ysty=3, color=djs_icolor('grey'), $
               title=title, charsize=1.5, charthick=postthick, yrange=yrange, $
               xthick=postthick, ythick=postthick, thick=2, position=pos[*,0], $
               xtickname=replicate(' ',10), ytitle=textoidl('f_{\lambda} (10^{-17} '+flam_units()+')')
             djs_oplot, she.wave, 1D17*she.spec, ps=10, color='green'

             legend, ['Moustakas','Jablonka'], /left, /top, box=0, $
               charsize=2.0, charthick=postthick, line=[0,0], thick=postthick, $
               color=djs_icolor(['grey','green'])

             djs_plot, me.wave, me.spec/she.spec, ps=10, xsty=3, ysty=3, $
               position=pos[*,1], /noerase, color='red', charsize=1.5, $
               charthick=postthick, xthick=postthick, ythick=postthick, thick=2, $
               xtitle='Wavelength [\AA]', ytitle='Ratio'
             if not keyword_set(postscript) then cc = get_kbrd(1)

             icleanup, me
             icleanup, she
             
          endif
          
       endfor

    endfor
       
    if keyword_set(postscript) then begin
       dfpsclose
       spawn, ['gzip -f '+rootpath+psname], /sh
    endif

stop    
    
return
end
    
