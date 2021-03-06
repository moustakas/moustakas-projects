function get_title, qq
    words = strsplit(qq,' ',/extract)
    nw = n_elements(words)
    if nw lt 5 then tit = qq else $
      tit = strjoin([words[0:nw/2-1],'!c',words[nw/2:nw-1]],' ')
return, tit    
end

function remap_answers, str, possible=possible, npos=npos, $
  hard=hard, agree=agree, never=never
    if keyword_set(agree) then possible = ['Strongly agree',$
      'Agree','Neutral','Disagree','Strongly disagree']
    if keyword_set(hard) then possible = ['Too hard','Somewhat hard',$
      'About right','Somewhat easy','Too easy']
    if keyword_set(never) then possible = ['Never','Less than half the time',$
      'About half the time','More than half the time','Always']
    npos = n_elements(possible)
    out = fltarr(n_elements(str))
    for pp = 0, npos-1 do begin
       this = where(strlowcase(strcompress(str,/remove)) eq $
         strlowcase(strcompress(possible[pp],/remove)))
       if this[0] ne -1 then out[this] = pp
    endfor 
return, out
end

pro phys390_f12_evaluations
; jm12oct20siena

    path = '~/Dropbox/teaching/12fall/phys390/evaluations/'
    
    semester = 'Fall 2012'
    class = 'Physics 390 - Introductory Astrophysics I'
    
; read the evaluation form downloaded from GoogleDocs
    data = read_evaluation(path+'phys390_evaluation1.csv',questions=qq)
    nstudent = n_elements(data[0,*])

; get rid of the timestamp
    keep = where(strmatch(qq,'*timestamp*',/fold) eq 0)
    qq = qq[keep]
    data = data[keep,*]

    bin = 0.5
    psfile = path+'phys390_evaluation1.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, $
      ymargin=[1.0,1.1], charsize=2.0, xmargin=[1.1,0.2], $
      width=7.2

; never...always questions
    indx = 5
    for ii = 0, n_elements(indx)-1 do begin
       num = remap_answers(data[indx[ii],*],possible=possible,npos=npos,/never)
       im_plothist, num, bin=bin, xx, yy, /noplot
       yrange = [0,max(yy)*1.2]
       im_plothist, num-bin/2.0, bin=bin, xtickinterval=1.0, $
         xtickname=strtrim(repstr(possible,' ','!c'),2), $
         xrange=[-0.5,npos-0.5], yrange=yrange, xsty=1, ysty=1, $
         /fill, position=pos, ytitle='Number'
       xyouts, (pos[2]-pos[0])/2.0+pos[0], pos[3]+(pos[3]-pos[1])*0.1, $
         get_title(qq[indx[ii]]), align=0.5, /normal, charsize=1.8
    endfor

; hard..easy questions
    indx = 4
    for ii = 0, n_elements(indx)-1 do begin
       num = remap_answers(data[indx[ii],*],possible=possible,npos=npos,/hard)
       im_plothist, num, bin=bin, xx, yy, /noplot
       yrange = [0,max(yy)*1.2]
       im_plothist, num-bin/2.0, bin=bin, xtickinterval=1.0, $
         xtickname=strtrim(repstr(possible,' ','!c'),2), $
         xrange=[-0.5,npos-0.5], yrange=yrange, xsty=1, ysty=1, $
         /fill, position=pos, ytitle='Number'
       xyouts, (pos[2]-pos[0])/2.0+pos[0], pos[3]+(pos[3]-pos[1])*0.1, $
         get_title(qq[indx[ii]]), align=0.5, /normal, charsize=1.8
    endfor
    
; agree...disagree questions
    indx = [0,1,2,6]
    for ii = 0, n_elements(indx)-1 do begin
       num = remap_answers(data[indx[ii],*],possible=possible,npos=npos,/agree)
       im_plothist, num, bin=bin, xx, yy, /noplot
       yrange = [0,max(yy)*1.2]
       im_plothist, num-bin/2.0, bin=bin, xtickinterval=1.0, $
         xtickname=strtrim(repstr(possible,' ','!c'),2), $
         xrange=[-0.5,npos-0.5], yrange=yrange, xsty=1, ysty=1, $
         /fill, position=pos, ytitle='Number'
       xyouts, (pos[2]-pos[0])/2.0+pos[0], pos[3]+(pos[3]-pos[1])*0.1, $
         get_title(qq[indx[ii]]), align=0.5, /normal, charsize=1.8
    endfor
    
    im_plotconfig, psfile=psfile, /pdf, /psclose, /pskeep

; push all the paragraph responses into a latex document
    indx = [3,7,8]
    texfile = path+'phys390_evaluation1_long.tex'
    openw, lun, texfile, /get_lun
    printf, lun, '\documentclass[12pt,preprint]{aastex}'
    printf, lun, '\usepackage{mdwlist}'
    printf, lun, '\pagestyle{plain}'
    printf, lun, '\begin{document}'
    for ii = 0, n_elements(indx)-1 do begin
;      printf, lun, '\vspace*{-2mm}'
       printf, lun, '\noindent {\large {\bf '+strtrim(qq[indx[ii]],2)+'}}'
;      printf, lun, '\vspace*{1mm}'
       printf, lun, '\begin{itemize*}'
       for jj = 0, nstudent-1 do begin
          txt = strcompress(data[indx[ii],jj],/remove)
          if (txt ne '') and (txt ne '...') then printf, lun, $
            '\item{'+strtrim(repstr(data[indx[ii],jj],'"',''),2)+'}'
          printf, lun, '\vspace*{2mm}'
       endfor
       printf, lun, '\end{itemize*}'
       printf, lun, '\vspace*{4mm}'
;      printf, lun, '\newpage'
    endfor
    printf, lun, '\end{document}'
    free_lun, lun
    spawn, 'latex '+texfile, /sh
    spawn, 'dvips -Ppdf -N0 -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex',''), /sh
    spawn, 'ps2pdf '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','.pdf'), /sh
    rmfile, repstr(texfile,'.tex','.aux')
    rmfile, repstr(texfile,'.tex','.log')
    rmfile, repstr(texfile,'.tex','.dvi')
    rmfile, repstr(texfile,'.tex','.ps')
    
;; build a montage of the results
;    spawn, 'convert '+psfile+' '+repstr(psfile,'.ps','.png'), /sh
;    pngfile = file_search(path+'phys130_'+suffix+'_evaluation1-*.png',count=npage)
;
;    nrow = 2
;    ncol = 2
;    size = 500
;    cmd = 'montage -bordercolor black -borderwidth 1 '+ $
;      '-tile '+strtrim(nrow,2)+'x'+strtrim(nrow,2)+' -geometry +0+0 '+$
;      '-quality 100 -resize '+size+'x'+size+' '+strjoin(rc3[these],' ')+$
;      ' '+outfile

; clean up    
    rmfile, psfile
;   spawn, '/bin/rm '+path+'phys130_'+suffix+'_evaluation1-*.png'
    
return
end
