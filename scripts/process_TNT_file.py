import os
import re
import csv


# Open the TNT file 
#After download, there's an error with the "Martinez" name. I first fixed it in notepad++ and then ran the script
#with open('mbank_X25860_7-9-2022_542.TNT') as file:
with open('./../../mbank_X25860_date.tnt') as file:
    text = file.read()


#Keep only the relevant part of the TNT file
stripped_begin = re.search(r"\n&\[cont\]\n",text)
stripped_end = re.search(r"\n\n\n;\n\n\ncnames",text)
stripped = text[stripped_begin.end():stripped_end.start()]



#Prepare the text for the CSV file
while '  ' in stripped:
    stripped = re.sub(' {2,}',' ',stripped)


stripped = stripped.replace(',','')

stripped = stripped.replace(' ',',')

stripped = stripped.replace('?','')

print(stripped)
#Replace ranges with the average value
while True:
    range_obj = re.search(r'[^,]+-[^,]+', stripped)
    if range_obj == None:
        break
    range = range_obj.group()
    range = range.split('-')
    range[0] = float(range[0])
    range[1] = float(range[1])
    new_value = str(round(sum(range)/2,1))
    stripped = stripped[:range_obj.span()[0]] + new_value + stripped[range_obj.span()[1]:]
    #print('replaced', range_obj.group(), 'with', new_value)

sep_strip = stripped.splitlines()

#hominins = re.findall(r'\n([^,]+)',stripped)
#stripped = re.sub(r'\n[^,]+,',r'\n',stripped)

#print(stripped)

#Write the CSV file


#with open('mbank_X25860_7-9-2022_541_character_list.txt','r',encoding= 'utf-8') as file:
with open('./../../mbank_X25860_date_character_list.txt','r',encoding= 'utf-8') as file:
    text = file.read()

#Get the character list
con_chr = re.search(r'235.*',text).group()
short_con_chr = re.findall(r'\d\d\d.\s([^@]*?).?\s\s\s',con_chr)
short_con_chr.insert(0,'specimen')
short_con_chr.append(re.search(r'555.\s([^@]*?)$',text).group(1))


con_chr = con_chr.split('   ')
con_chr.insert(0,'specimen')


with open('./../../raw_craniometric_data_straight_from_tnt.csv', 'w', encoding="utf-8",newline= '') as f:
    writer = csv.writer(f)
    writer.writerow(short_con_chr)
    for row in sep_strip:
        print(row.split(','))
        writer.writerow(row.split(','),)


print('done')
