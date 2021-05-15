#!/usr/local/bin/python3


import mysql.connector
import jinja2
import cgi

# This line tells the template loader where to search for template files
templateLoader = jinja2.FileSystemLoader( searchpath="./templates" )

# This creates your environment and loads a specific template
env = jinja2.Environment(loader=templateLoader)
template = env.get_template('history_template.html')

form = cgi.FieldStorage()
Sequence = str(form.getfirst('Sequence'))
Reverse = str(form.getfirst('Reverse'))
Length = str(form.getfirst('Length'))
GC = str(form.getfirst('GC'))
tm = str(form.getfirst('tm'))





outputs=list()

def main():
    conn = mysql.connector.connect(user='hzhan141', password='huiminzhang',
                                   host='localhost', database='hzhan141_chado')
    

    curs = conn.cursor()
    qry=  """
        SELECT * from history
        WHERE Sequence like %s
        OR Reverse like %s
       OR Length like %s
        OR GC like %s
        OR Melting_temperature like %s;
        """  
        
    curs.execute(qry, (Sequence , Reverse , Length,  GC ,tm ))
 
  
    
    for row in curs:
        outputs.append({'seq':row[0], 'rev':row[1], 'len':row[2], 'gc':row[3], 'TM':row[4]})
     
    curs.close()
    conn.close()

if __name__ == '__main__':
    main()

print("Content-type: text/html\n\n")
print(template.render(outputs=outputs))
#print ("%s" % output)

