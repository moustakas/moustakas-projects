function get_title, qq
    words = strsplit(qq,' ',/extract)
    nw = n_elements(words)
    if nw lt 5 then tit = qq else $
      tit = strjoin([words[0:nw/2-1],'!c',words[nw/2:nw-1]],' ')
return, tit    
end

function remap_answers, str, possible=possible, npos=npos, $
  hard=hard, agree=agree, never=never, hours=hours
    if keyword_set(agree) then possible = ['Strongly agree',$
      'Agree','Neutral','Disagree','Strongly disagree']
    if keyword_set(hard) then possible = ['Too hard','Somewhat hard',$
      'About right','Somewhat easy','Too easy']
    if keyword_set(never) then possible = ['Never','Once or twice',$
      'Several times','Every week']
    if keyword_set(hours) then possible = ['0-2','2-5','5-10',$
      'More than 10']
    npos = n_elements(possible)
    out = fltarr(n_elements(str))
    for pp = 0, npos-1 do begin
       this = where(strlowcase(strcompress(str,/remove)) eq $
         strlowcase(strcompress(possible[pp],/remove)))
       if this[0] ne -1 then out[this] = pp
    endfor 
return, out
end

pro phys130_f13_evaluations
; jm13oct04siena

    path = getenv('TEACHING_DIR')+'/Phys130/130-F13/evaluations/'
    
    semester = 'Fall 2013'
    class = 'Physics 130 - General Physics I'
    
; read the evaluation form downloaded from GoogleDocs
    suffix = '06'
    data = read_evaluation(path+'phys130-'+suffix+'-evaluation1.csv',questions=qq)
    nstudent = n_elements(data[0,*])

; get rid of the timestamp and "Sample Question 2"
    keep = where(strmatch(qq,'*timestamp*',/fold) eq 0)
    qq = qq[keep]
    data = data[keep,*]

    bin = 0.5
    psfile = path+'phys130-'+suffix+'-evaluation1.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, $
      ymargin=[1.0,1.1], charsize=2.0, xmargin=[1.1,0.2], $
      width=7.2

; agree...disagree questions - questions 0-7
    indx = lindgen(8)
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
    
; hard..easy questions - questions 10-12
    indx = [10,11,12]
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
    
; number of hours questions - question 8,9
    indx = [8,9]
    for ii = 0, n_elements(indx)-1 do begin
       num = remap_answers(data[indx[ii],*],possible=possible,npos=npos,/hours)
       im_plothist, num, bin=bin, xx, yy, /noplot
       yrange = [0,max(yy)*1.2]
       im_plothist, num-bin/2.0, bin=bin, xtickinterval=1.0, $
         xtickname=strtrim(repstr(possible,' ','!c'),2), $
         xrange=[-0.5,npos-0.5], yrange=yrange, xsty=1, ysty=1, $
         /fill, position=pos, ytitle='Number'
       xyouts, (pos[2]-pos[0])/2.0+pos[0], pos[3]+(pos[3]-pos[1])*0.1, $
         get_title(qq[indx[ii]]), align=0.5, /normal, charsize=1.8
    endfor
    
; never...every week questions - questions 13,14
    indx = [13,14]
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

    im_plotconfig, psfile=psfile, /pdf, /psclose, /pskeep

; push all the paragraph responses into a latex document; questions
; 15-17
    indx = [15,16,17]
    texfile = path+'phys130-'+suffix+'-evaluation1-long.tex'
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
          if (txt ne '') and (txt ne '...') then begin
             printf, lun, '\item{'+strtrim(repstr(repstr(data[indx[ii],jj],'"',''),'&','\&'),2)+'}'
             printf, lun, '\vspace*{2mm}'
          endif
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
