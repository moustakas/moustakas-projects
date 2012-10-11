pro process_grades, data, assign=assign, allassign=allassign, $
  weight=weight, class=class, semester=semester, lab=lab
; jm12oct05siena - process the final grades

    nstudent = n_elements(data)
    
    cumupoints = weight*0       ; cumulative number of points
    cumupossible = weight*0     ; cumulative number of points possible
    isdata = fix(weight*0)

; loop on each student and then on all the possible assignments
    for ss = 0, nstudent-1 do begin
       tmpfile = 'tmp_student.txt'
       openw, lun, tmpfile, /get_lun
       
       printf, lun, 'Grade Summary for '+strupcase(strtrim(data[ss].first_name,2))+' '+$
         strupcase(strtrim(data[ss].last_name,2))
       printf, lun, class
       printf, lun, semester
       printf, lun, 'Current Date: '+im_today()
       printf, lun, ' '
       printf, lun, '                       ******'
       printf, lun, '   The following material is CONFIDENTIAL and intended'
       printf, lun, '       exclusively for the student named above.'
       printf, lun, '                       ******'
       printf, lun, ' '
       if keyword_set(lab) then begin
          printf, lun, 'Listed below are your individual lab scores.  The three'
          printf, lun, 'numbers listed are the number of points earned, the number'
          printf, lun, 'of possible points, and the percent score.  If any of your'
          printf, lun, 'scores are zero then you must make up that lab in order to'
          printf, lun, 'pass the class.  If you identify any errors please contact'
          printf, lun, 'me immediately.'
          printf, lun, '   -Prof. Moustakas (jmoustakas@siena.edu)'
       endif else begin
          printf, lun, 'Listed below is each component of your grade (e.g., Homework),'
          printf, lun, 'and the score you received on each assignment or examination.'
          printf, lun, 'The three numbers listed are the number of points earned, the'
          printf, lun, 'number of possible points, and the percent score.  If you '
          printf, lun, 'identify any errors please contact me immediately.'
          printf, lun, '   -Prof. Moustakas (jmoustakas@siena.edu)'
       endelse
       printf, lun, ' '

       for ii = 0, n_elements(allassign)-1 do begin
          printf, lun, strupcase(allassign[ii])+' (Weight='+string(100*weight[ii],format='(I0)')+'%)'

; has this been assigned yet (or has it happened yet, e.g., final?)
          yes = where(strlowcase(allassign[ii]) eq strlowcase(strtrim(assign,2)))
          if yes[0] ne -1 then begin
             isdata[ii] = 1
             
; get the relevant structure indices
             points = data[ss].(tag_indx(data,strupcase(repstr(allassign[ii],' ','_'))))
             possible = data[ss].(tag_indx(data,strupcase(repstr(allassign[ii]+'_points',' ','_'))))
             details = data[ss].(tag_indx(data,strupcase(repstr(allassign[ii]+'_details',' ','_'))))
             dates = data[ss].(tag_indx(data,strupcase(repstr(allassign[ii]+'_date',' ','_'))))

; detail each assignment...
             nassign = n_elements(points)
             for nn = 0, nassign-1 do begin
                if possible[nn] gt 0 then begin ; normal assignments
                   perc = 100.0*points[nn]/possible[nn]
                   printf, lun, '  '+details[nn]+' due '+dates[nn]
                   printf, lun, '     '+strtrim(string(points[nn],format='(F12.1)'),2)+$
                     ', '+strtrim(string(possible[nn],format='(F12.1)'),2)+$
                     ', '+strtrim(string(perc,format='(F12.2)'),2)+'%'
                endif
                if possible[nn] eq 0 and points[nn] gt 0.0 then begin ; extra credit
                   printf, lun, '  Extra Credit: '+details[nn]+' on '+dates[nn]
                   printf, lun, '     Homework points earned: '+$
                     strtrim(string(points[nn],format='(F12.1)'),2)
                endif
             endfor

; ...now add it up             
             totpoints = total(points)
             totpossible = total(possible)
             totperc = 100.0*totpoints/totpossible
             
             printf, lun, '  TOTAL:'
             printf, lun, '     '+strtrim(string(totpoints,format='(F12.1)'),2)+$
               ', '+strtrim(string(totpossible,format='(F12.1)'),2)+$
               ', '+strtrim(string(totperc,format='(F12.2)'),2)+'%'

             cumupoints[ii] = totpoints
             cumupossible[ii] = totpossible
          endif else begin
             printf, lun, '   No scores yet.'

; assume the student gets 90%
             isdata[ii] = -1
             cumupoints[ii] = 90.0
             cumupossible[ii] = 100.0
          endelse
       printf, lun, ' '
    endfor 

; print or predict the final grade
       final = 100*im_weighted_mean(cumupoints/cumupossible,weights=weight)
       if total(isdata eq -1) gt 0 then begin
          printf, lun, 'Some of your scores in the class have not been tabulated yet.'
          printf, lun, 'However, if we ASSUME that you will earn 90% of the possible'
          printf, lun, 'points on ALL the remaining assignments or examinations, then'
          printf, lun, 'your final percent grade in the class will be '+$
            strtrim(string(final,format='(F12.2)'),2)+'%.'
          printf, lun, ' '
          printf, lun, 'Keep studying!!'
       endif else begin
          if keyword_set(lab) then begin
             printf, lun, 'Your current percent grade in lab is '+$
               strtrim(string(final,format='(F12.2)'),2)+'%.'
          endif else begin
             printf, lun, 'Your final percent grade in the class is '+$
               strtrim(string(final,format='(F12.2)'),2)+'%.'
          endelse
       endelse
; close and email the file
       free_lun, lun
       spawn, 'cat '+tmpfile
       rmfile, tmpfile
    endfor 

return
end
    
