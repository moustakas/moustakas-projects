pro primus_2dqaplots, date=date1, test=test, html=html
; jm08dec01nyu - all the QA plots from a given mask into a single file

    html = 1
    
; wget --user=primus --password=unicron -r -nd -nH -np --accept ".pdf" http://sdss.physics.nyu.edu/primus/redux/0017/uberqaplots/    
    
;   gs -dNOPAUSE -sDEVICE=pdfwrite -sOutputFile=new.pdf -dBATCH XXXX.pdf XXXX.pdf <XXXX.ps XXXX.ps etc.>
    
    rerun = '0017'
    basedir = getenv('PRIMUS_DATA')+'/redux/'+rerun+'/'
    outdir = basedir+'uberqaplots/'
    html_path = basedir+'www_qaplots/'
    if (file_test(outdir,/dir) eq 0) then spawn, 'mkdir '+outdir, /sh

;   date1 = 'ut071206'

    if (n_elements(date1) eq 0L) then begin
       date = file_basename(file_search(basedir+'ut??????'))
    endif else date = date1
    ndate = n_elements(date)
    if (ndate gt 1L) then begin
       for idate = 0L, ndate-1L do primus_2dqaplots, date=date[idate], test=test
       return
    endif

    for ii = 0L, ndate-1L do begin
       thisdir = basedir+date[ii]+'/'
       qaplotsdir = thisdir+'qaplots/'
       allmask = file_basename(file_search(thisdir+'*_summary.txt',count=nmask))
;      for jj = 2L, 2L do begin
       for jj = 0L, nmask-1L do begin
          thismask = strmid(allmask[jj],0,strpos(allmask[jj],'_summary.txt'))
; read the summary file to figure out the exposures
          summary = djs_readlines(thisdir+allmask[jj])
          summary = summary[where(strmatch(summary,'*Extraction Statistics*'))+1:$
            where(strmatch(summary,'*Signal/Noise Statistics*'))-1]
          summary = summary[where(strcompress(summary,/remove) ne '')]
;         niceprint, summary
          openw, lun, '/tmp/junk', /get_lun
          for mm = 0L, n_elements(summary)-1L do $
            printf, lun, summary[mm]
          free_lun, lun
          readcol, '/tmp/junk', ccd, skipline=1, format='A', /silent
          nccd = n_elements(ccd)
;         niceprint, ccd
          ccd_start = repstr(ccd[0],'ccd','')
          ccd_end = repstr(ccd[nccd-1],'ccd','')
; build the QA plots list
;         qaplots = file_search(qaplotsdir+repstr(maskplots,'thismask',thismask),count=nqa)
;         qaplots = qaplotsdir+repstr(maskplots,'thismask',thismask)
          if ((date[ii] eq 'ut061223') and (thismask eq 'vvds0138')) or $
            ((date[ii] eq 'ut061224') and (thismask eq 'vvds0139')) then begin
             snr_qaplots = ['snr/'+thismask+'_qa_snr_Combo-17.ps.gz','snr/'+thismask+'_qa_snr_VVDS.ps.gz']
             snr_ccdplots = ['snr/'+ccd+'_qa_snr_Combo-17.ps.gz',$
               'snr/'+ccd+'_qa_snr_VVDS.ps.gz']
          endif else begin
             snr_qaplots = 'snr/'+thismask+'-qaplot.ps.gz'
             snr_ccdplots = 'snr/'+ccd+'-qaplot.ps.gz'
          endelse
          wavecheck_ccdplots = 'wavecal/'+ccd+'_sky_wavecheck.ps.gz'
          qaplots = qaplotsdir+[$
            'findslits/'+thismask+'_findslit_qa.ps.gz',$
            snr_qaplots,$
            snr_ccdplots,$
            ''+thismask+'_losses.ps',$
            'extract/'+thismask+'_halo.ps',$
            'extract/'+thismask+'_skypsf.ps',$
            'wavecal/'+thismask+'_coadd_sky_wavecheck.ps.gz',$
            wavecheck_ccdplots,$
            thismask+'_qa_ccd'+ccd_start+'-'+ccd_end+'_extract.ps.gz',$
            thismask+'_qa_ccd'+ccd_start+'-'+ccd_end+'_skysubtract.ps.gz',$
            thismask+'_emptyslits.ps',$
            thismask+'_fstar_qa.ps'$
            ]
          nqa = n_elements(qaplots)
;         niceprint, qaplots
; uncompress all .gz files to /tmp
          gz = where(strmatch(qaplots,'*.ps.gz*'),ngz)
          if (ngz ne 0L) then begin
;            for kk = 0L, 0L do begin
             for kk = 0L, ngz-1L do begin
                spawn, '/bin/cp -f '+qaplots[gz[kk]]+' /tmp/'+file_basename(qaplots[gz[kk]]), /sh
                qaplots[gz[kk]] = '/tmp/'+file_basename(qaplots[gz[kk]])
                spawn, 'gunzip -f '+qaplots[gz[kk]], /sh
             endfor
          endif
          qaplots = repstr(qaplots,'.ps.gz','.ps')
; check whether all the relevant plots exist and build a summary file
          for qq = 0L, nqa-1L do begin
             info1 = struct_trimtags(file_info(qaplots[qq]),select=['name','size'])
             info1 = create_struct(info1,{npage: 0})
             if (info1.size gt 0) then begin
                spawn, "grep '^%%Pages' "+info1.name+" | awk '{ print $2 }' | grep -v '(atend)'", npage, /sh
;               print, qaplots[qq], npage
                info1.npage = npage
             endif
             info1.name = file_basename(info1.name)
             if (qq eq 0L) then info = info1 else info = [info,info1]
          endfor
;         struct_print, info
;         niceprint, qaplots
          outfile = outdir+'uberqaplot_'+thismask+'_'+date[ii]+'.pdf'
          txtfile = repstr(outfile,'.pdf','.txt')
;         outfile = qaplotsdir+'uberqaplot_'+thismask+'.pdf'

          htmlbase = thismask+'_'+date[ii]
          if (file_test(html_path+htmlbase,/dir) eq 0L) then spawn, $
            'mkdir -p '+html_path+htmlbase, /sh
          if keyword_set(html) then spawn, 'cp -f '+strjoin(qaplots,' ')+' '+$
            html_path+htmlbase+'/', /sh

; build the relevant PNG files          
          im_ps2html, htmlbase, html_path=html_path, /psfind, /only_png
          im_ps2html, htmlbase, html_path=html_path, /pngfind
          
stop
          splog, 'Writing '+outfile
          if (not keyword_set(test)) then begin
             spawn, 'gs -q -dNOPAUSE -sDEVICE=pdfwrite '+$
               '-sOutputFile='+outfile+' -dBATCH '+strjoin(qaplots,' '), /sh
             struct_print, info, file=txtfile
          endif
       endfor
    endfor
    
return
end
