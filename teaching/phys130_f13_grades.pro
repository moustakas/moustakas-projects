pro phys130_f13_grades, alldata, test=test, sendit=sendit, final=final
; jm13aug25siena - parse the grades for my two sections of this class 

    date = '' ; update this
    semester = 'Fall 2013'
    class = 'Physics 130 - General Physics I'

    section = ['06','24']
    for ii = 0, 1 do begin
       path = getenv('TEACHING_DIR')+'/130-F13/grades/'+section[ii]+'/'

; read the grade spreadsheet downloaded from GoogleDocs
       gradefile = path+'phys130-'+section[ii]+'-grades-'+date+'.csv'
       data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
       allassign = ['Homework','Lab','Midterm 1','Midterm 2','Midterm 3','Final Exam']
       weight = [0.2,0.15,0.15,0.15,0.15,0.2]
;      keep = where(strmatch(data.last_name,'*Ippo*'))
;      keep = where(data.final_exam gt 0.0)
;      data = data[keep]
    
       process_grades, data, assign=assign, allassign=allassign, $
         weight=weight, class=class, semester=semester, test=test, $
         sendit=sendit, alldata=alldata, final=final
       srt = sort(alldata.current_grade)
       struct_print, alldata[srt]
    endfor

return
end
