pro phys130_f13_grades, alldata, test=test, sendit=sendit, final=final
; jm13aug25siena - parse the grades for my two sections of this class 

    date = '13dec15' ; update this
    semester = 'Fall 2013'
    class = 'Physics 130 - General Physics I'

    section = ['06','24']
    for ii = 0, 1 do begin
       path = getenv('TEACHING_DIR')+'/Phys130/130-F13/grades/'+section[ii]+'/'

; read the grade spreadsheet downloaded from GoogleDocs
       gradefile = path+'phys130-'+section[ii]+'-'+date+'.csv'
       data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
       allassign = ['Homework','Lab','Quizzes','Final Exam']
       weight = [0.25,0.15,0.40,0.20]
       droplowest = [0,0,1,0]

;      keep = where(strmatch(data.last_name,'*Lavera*') eq 0)
       keep = where(data.final_exam gt 0.0)
       data = data[keep]

       process_grades, data, assign=assign, allassign=allassign, $
         weight=weight, droplowest=droplowest, class=class, $
         semester=semester, test=test, sendit=sendit, alldata=alldata, $
         final=final
       srt = sort(alldata.current_grade)
       struct_print, alldata[srt]
;      struct_print, alldata

       print
       print
       if keyword_set(sendit) eq 0 then cc = get_kbrd(1)
    endfor

return
end
