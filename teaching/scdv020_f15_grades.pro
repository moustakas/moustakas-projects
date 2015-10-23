pro scdv020_f15_grades, alldata, test=test, sendit=sendit, final=final
; jm15sep29 - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/020-F15/grades/'
    
    date = '15oct12' ; update this
    semester = 'Fall 2015'
    class = 'SCDV 020 - Intro to Engineering'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'scdv020_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework']
    weight = [1.0]
    droplowest = [0]

;   keep = where(strmatch(data.last_name,'*Young*'))
;   keep = where(data.final_exam gt 0.0)
;   data = data[keep]

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
