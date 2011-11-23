import os
import sys

def parse_line(line):
    line = line.split(',')
    info = line[1:] # take out the time sig at the beginning of each line
    return info
 
def create_line(info):
    new_line = ""
    for point in info:
        num_spaces = 11 - len(point)
        new_line += " "*num_spaces + str(point)+"\t"
    return new_line

def run():
    file_to_open = str(sys.argv[1])
    file_to_save = str(sys.argv[2])

    file(file_to_save,'wb').write('')
    f_new = open(file_to_save,'ab')
    counter = 0
    for line in open(file_to_open):
        if counter == 0:
            counter += 1
            continue
        info = parse_line(line)
        new_line = create_line(info)
        print new_line
        f_new.write(new_line)
        counter += 1

    f_new.close()
run()
