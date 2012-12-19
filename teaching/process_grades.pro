pro process_grades, data, assign=assign, allassign=allassign, $
  weight=weight, class=class, alldata=alldata, semester=semester, $
  lab=lab, test=test, sendit=sendit, final=final
; jm12oct05siena - process the final grades

    nstudent = n_elements(data)
    
    cumupoints = weight*0       ; cumulative number of points
    cumupossible = weight*0     ; cumulative number of points possible
    isdata = fix(weight*0)

; build the output data structure
    alldata = struct_trimtags(data,select=['last_name','first_name',$
      'email','major','year'])
    create_struct, junk, '', [allassign,'current_grade','final_grade','letter_grade'], $
      [replicate('F',n_elements(allassign)),'F','F','A']
    alldata = struct_addtags(alldata,replicate(junk,nstudent))
    finaltags = tag_names(alldata)
    
; loop on each student and then on all the possible assignments
    for ss = 0, nstudent-1 do begin
       tmpfile = 'tmp_student.txt'
       openw, lun, tmpfile, /get_lun

       name = strupcase(strtrim(data[ss].first_name,2))+' '+$
         strupcase(strtrim(data[ss].last_name,2))
       printf, lun, 'Grade Summary for '+name
       printf, lun, class
       printf, lun, semester
       printf, lun, 'Current Date: '+im_today()
       printf, lun, ' '
       printf, lun, '                           ******'
       printf, lun, '   The following material is CONFIDENTIAL and intended'
       printf, lun, '       exclusively for the student named above.'
       printf, lun, '                           ******'
       printf, lun, ' '
       if keyword_set(lab) then begin
          printf, lun, 'Listed below are your individual lab scores.  The three'
          printf, lun, 'numbers listed are the number of points earned, the number'
          printf, lun, 'of possible points, and the percent score.  If any of your'
          printf, lun, 'scores are zero then you must make up that lab in order to'
          printf, lun, 'pass the class.  If you identify any errors please contact'
          printf, lun, 'me immediately.'
;         printf, lun, '   -Prof. Moustakas (john.moustakas@gmail.com)'
       endif else begin
          printf, lun, 'Listed below is each component of your grade (e.g., Homework),'
          printf, lun, 'and the score you have received on each assignment or examination.'
          printf, lun, 'The three numbers listed are the number of points earned, the'
          printf, lun, 'number of possible points, and the percent score.  If you '
          printf, lun, 'identify any errors please contact me immediately.'
;         printf, lun, '   -Prof. Moustakas (john.moustakas@gmail.com)'
       endelse
       printf, lun, ' '

       for ii = 0, n_elements(allassign)-1 do begin
          printf, lun, strupcase(allassign[ii])+' (Weight='+$
            string(100*weight[ii],format='(I0)')+'%)'

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
                   printf, lun, '  '+details[nn]+': '+dates[nn]
                   printf, lun, '     '+strtrim(string(points[nn],format='(F12.1)'),2)+$
                     ', '+strtrim(string(possible[nn],format='(F12.1)'),2)+$
                     ', '+strtrim(string(perc,format='(F12.2)'),2)+'%'
                endif
                if possible[nn] eq 0 and points[nn] gt 0.0 then begin ; extra credit
                   printf, lun, '  Extra Credit: '+details[nn]+': '+dates[nn]
                   if keyword_set(lab) then printf, lun, '     Points earned: '+$
                     strtrim(string(points[nn],format='(F12.1)'),2) else $
                       printf, lun, '     Homework points earned: '+$
                     strtrim(string(points[nn],format='(F12.1)'),2)
                endif
             endfor

; ...now add it up             
             totpoints = total(points)
             totpossible = total(possible)
             totperc = 100.0*totpoints/totpossible

             if nassign gt 1 then begin
                printf, lun, '  TOTAL:'
                printf, lun, '     '+strtrim(string(totpoints,format='(F12.1)'),2)+$
                  ', '+strtrim(string(totpossible,format='(F12.1)'),2)+$
                  ', '+strtrim(string(totperc,format='(F12.2)'),2)+'%'
             endif

             cumupoints[ii] = totpoints
             cumupossible[ii] = totpossible

; pack it into the output structure
             indx = tag_indx(alldata,idl_validname(allassign[ii],/convert_all))
             alldata[ss].(indx) = totpoints/totpossible
          endif else begin
             printf, lun, '   No scores yet.'

; assume the student gets 90%
             isdata[ii] = -1
             cumupoints[ii] = 90.0
             cumupossible[ii] = 100.0
          endelse
       printf, lun, ' '
    endfor 

; print the current grade and predict the final grade, unless
; we're at the end of the semester
       these = where(isdata eq 1)
       current = 100*im_weighted_mean(cumupoints[these]/$
         cumupossible[these],weights=weight[these])
       indx = tag_indx(alldata,'current_grade')
       alldata[ss].(indx) = current
       grade = 100*im_weighted_mean(cumupoints/cumupossible,weights=weight)
       if total(isdata eq -1) gt 0 then begin
          printf, lun, 'Some of your scores in the class have not been tabulated yet.'
          printf, lun, 'However, based on the fraction of your grade that you have'
          printf, lun, 'earned so far, your current percent grade in the class is '+$
            strtrim(string(current,format='(F12.2)'),2)+'%.'
          printf, lun, ' '
          printf, lun, 'If we ASSUME that you will earn 90% of the possible remaining'
          printf, lun, 'points on ALL the remaining assignments or examinations, then'
          printf, lun, 'your final percent grade in the class will be '+$
            strtrim(string(grade,format='(F12.2)'),2)+'%.'
          printf, lun, ' '
          printf, lun, 'Keep studying!!'
       endif else begin
          if keyword_set(lab) then begin
             if keyword_set(final) then begin
                indx = tag_indx(alldata,'final_grade')
                alldata[ss].(indx) = grade
                printf, lun, 'All labs have been completed!  Your final'
                printf, lun, 'percent grade in the course is: '+$
                  strtrim(string(grade,format='(F12.2)'),2)+'%'
             endif else begin
                printf, lun, 'Your current percent grade in lab is '+$
                  strtrim(string(grade,format='(F12.2)'),2)+'%.'
             endelse
          endif else begin
             indx = tag_indx(alldata,'final_grade')
             alldata[ss].(indx) = grade

             indx = tag_indx(alldata,'letter_grade')
             iindx = tag_indx(data,'letter_grade')
             alldata[ss].(indx) = data[ss].(iindx)
             
             printf, lun, 'Your final percent grade in the class is: '+$
               strtrim(string(grade,format='(F12.2)'),2)+'%'
             printf, lun, ' '
             printf, lun, 'Your final letter grade in the class is: '+$
               strtrim(alldata[ss].(indx),2)
          endelse
       endelse
       printf, lun, ' '
       printf, lun, '   -Prof. Moustakas (john.moustakas@gmail.com)'
; close and email the file
       free_lun, lun
       if keyword_set(test) then begin
          junkfile = 'junk.txt'
          openw, lun, junkfile, /get_lun
          printf, lun, 'Grade Summary for '+name
          printf, lun, class
          printf, lun, semester
          printf, lun, 'Current Date: '+im_today()
          printf, lun, ' '
          printf, lun, 'This email is intended for '+name+'.  If this '
          printf, lun, 'email has been sent to the wrong individual then'
          printf, lun, 'please let me know as soon as possible.'
          printf, lun, '   -Prof. Moustakas (john.moustakas@gmail.com)'
          free_lun, lun
;         data[ss].email = 'kineta.beach@gmail.com'
          if keyword_set(sendit) then begin
             if ss eq 0 then begin
                splog, 'About to send test grade summary emails - are you sure?!?'
                cc = get_kbrd(1)
                if strupcase(cc) ne 'Y' then stop
             endif
             spawn, 'send_an_email_to_many_people_using_gmail.py '+$
               '--email '+strtrim(data[ss].email,2)+' --subject '+$
               '"Test of automated grade summary script" --text-body '+$
               junkfile, /sh
          endif else spawn, 'cat '+junkfile, /sh
          rmfile, junkfile
;         wait, 5
       endif else begin          
          if keyword_set(sendit) then begin
             if ss eq 0 then begin
                splog, 'About to send actual grade summary emails - are you sure?!?'
                cc = get_kbrd(1)
                if strupcase(cc) ne 'Y' then stop
             endif
             spawn, 'send_an_email_to_many_people_using_gmail.py '+$
               '--email '+strtrim(data[ss].email,2)+' --subject '+$
               '"'+im_today()+' Grade Summary" --text-body '+$
               tmpfile, /sh
          endif else spawn, 'cat '+tmpfile, /sh
       endelse
       rmfile, tmpfile
;      wait, 5
    endfor ; close student loop

return
end
    
