import pandas as pd
import numpy as np
import re

# THIS SCRIPT ASSUMES ONLY LETTERS USED ARE y AND m 
# CONVERT TO THIS BEFORE USING SCRIPT!

# BE CAREFUL ALL - SYMBOLS ARE ACTUAL - (not –)

# This script first converts manually-curated age/date/sex to consistant age formatting
# and adds to the rest of the metadata.
# Then it goes through ALL metadata to put into ranges.

# See below for list of formats allowed
# Note that range 0-1 is put as <1 because I think this is what is meant.
# Note that range 0-18 is put as <18 because I think this is what is meant
#       (x-18 is always recorded as <18 in age_range 2 and 3)


#Converts months to years. Pass ints
def month_to_year(mon):
    return round(float(mon)/12, 2)

#should only be passed single ages! not ranges!
#(pass range ages individually)
def convert_to_year(dat):
    #0.5 or 15 - basically, already a year!
    if len(re.findall('[a-z]', dat)) == 0:
        new_yr = float(dat)
    #1y5m
    elif len(re.findall('[a-z]', dat)) == 2:
        parts = re.split('[a-z]', dat)
        new_yr = float(parts[0]) + month_to_year(int(parts[1]))
    #4y
    elif re.search('y', dat):
        new_yr = float(dat.replace("y",""))
    #14m
    elif re.search('m', dat):
        new_yr = round(float(dat.replace("m", ""))/12, 2)
    #15d
    elif re.search('d', dat):
        new_yr = round(float(dat.replace("d", ""))/365, 2)
    else:
        print("ERROR! Problem with age ", dat)
    return new_yr


if __name__ == '__main__':
    import argparse

    parser = parser = argparse.ArgumentParser(description='parse metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ages-in', help="input ages/dates/sex file")           
    parser.add_argument('--meta-in', help="input meta file")
    parser.add_argument('--meta-out-ages', help="output metadata just with ages added")
    parser.add_argument('--meta-out', help="final output metadata with age categories too")
    args = parser.parse_args()

    #################################################################
    ##      Part 1
    #################################################################

    ages = pd.read_csv(args.ages_in, sep='\t', index_col=False)
    meta = pd.read_csv(args.meta_in, sep='\t', index_col=False)

    #create new column for AFM (& reorder)
    meta = meta.reindex(columns = ['strain', 'accession', 'date', 'sex', 'age', 
        'symptom', 'country', 'collab_country', 'region', 'host',
        'subgenogroup', 'Lab-ID', 'orig_strain', 'seq-len', 'date_added', 'genbank-host', 'moltype',
        'virus', 'authors', 'title', 'url', 'paper_url'])


    for i, row in ages.iterrows():
        #only do this if the accession is in the metadata
        if row.accession in meta.accession.values:
            
            if not pd.isnull(row.date):
                meta.loc[meta.accession == row.accession, 'date'] = row.date.strip()
            if not pd.isnull(row.sex):
                meta.loc[meta.accession == row.accession, 'sex'] = row.sex.strip()
            if not pd.isnull(row.symptom):
                meta.loc[meta.accession == row.accession, 'symptom'] = row.symptom.strip()

            #0.8    years   leave
            #77     years   leave
            #4y     years   remove y
            #14m    months  convert to years
            #15d    days    convert to years
            #1y2m   both    convert to years
            #<11y   range   convert to range (0-11)
            #>79y   range   convert to range (79-200)
            #29-78y range   convert to range (29-78)
            #1-18   range   convert to range (1-18)
            #.5-18y range   convert to range (.5-18)
            #2m-10y range   convert to years and range

            #need to process ages!
            if not pd.isnull(row.age):
                #is a range with -
                if re.search('-', row.age):
                    parts = re.split('-', row.age)
                    new_age = str(convert_to_year(parts[0])) + '-' + str(convert_to_year(parts[1]))
                #is range with <
                elif re.search('<', row.age):
                    new_age = '0-' + str(convert_to_year(row.age.replace("<", "")))
                #is range with >
                elif re.search('>', row.age):
                    new_age = str(convert_to_year(row.age.replace(">", ""))) + '-200'
                #is single number
                else:
                    new_age = str(convert_to_year(row.age))

                meta.loc[meta.accession == row.accession, 'age'] = new_age

    meta.to_csv(args.meta_out_ages, sep='\t', index=False)


    #################################################################
    ##      Part 2
    #################################################################

    #now put them into categories as well.

    #raw_age    all ages in whatever format after parsing.
    #age        exact ages (no ranges)
    #age_range1 <1yr 1-5yr  6-17yr  18-64yr >=65yr
    #age_range2 <18yr  18-65yr >=65yr
    #age_range3 <18yr   >=18yr

    #rename age column raw_age
    meta.rename(columns={'age':'raw_age'}, inplace=True)

    #create new columns for new ages (& reorder)
    meta = meta.reindex(columns = ['strain', 'accession', 'date', 'sex', 'raw_age',
        'symptom', 'age', 'age_range1', 'age_range2', 'age_range3', 'country', 'region', 'collab_country',
        'host',
        'subgenogroup', 'Lab-ID', 'orig_strain', 'seq-len', 'date_added', 'genbank-host', 'moltype',
        'virus', 'authors', 'title', 'url', 'paper_url'])

    #create temp holders
    age = []
    age_range1 = []
    age_range2 = []
    age_range3 = []

    for raw_age in meta.raw_age:
        if pd.isnull(raw_age):
            age.append('')
            age_range1.append('')
            age_range2.append('')
            age_range3.append('')
        else:
            #is a range with -  (only try search if its a string)
            if type(raw_age)==str and re.search('-', raw_age):
                age.append('')
                parts = np.array([float(x) for x in re.split('-', raw_age)])
                #special treatment 0-1
                if parts[0] == 0 and parts[1] == 1:
                    age_range1.append('<1y')
                    age_range2.append('<18y')
                    age_range3.append('<18y')
                #special treatment 0-18
                elif parts[0] == 0 and parts[1] == 18:
                    age_range1.append('')
                    age_range2.append('<18y')
                    age_range3.append('<18y')
                elif all(parts < 1):
                    age_range1.append('<1y')
                    age_range2.append('<18y')
                    age_range3.append('<18y')
                elif all(parts < 6):
                    if all(parts >=1) and all(parts < 6):
                        age_range1.append('1-5y')
                    else:
                        age_range1.append('')
                    age_range2.append('<18y')
                    age_range3.append('<18y')
                elif all(parts < 18):
                    if all(parts >=6) and all(parts < 18):
                        age_range1.append('6-17y')
                    else:
                        age_range1.append('')
                    age_range2.append('<18y')
                    age_range3.append('<18y')
                #special treatment xx-18
                elif parts[0] < 18 and parts[1] == 18:
                    age_range1.append('')
                    age_range2.append('<18y')
                    age_range3.append('<18y')
                elif all(parts < 65):
                    if all(parts >=18) and all(parts < 65):
                        age_range1.append('18-64y')
                        age_range2.append('18-64y')
                    else:
                        age_range1.append('')
                        age_range2.append('')
                    if all(parts >= 18):
                        age_range3.append('>=18y')
                    else:
                        age_range3.append('')
                else:
                    if all(parts >= 65):
                        age_range1.append('>=65y')
                        age_range2.append('>=65y')
                    else:
                        age_range1.append('')
                        age_range2.append('')
                    if all(parts >= 18):
                        age_range3.append('>=18y')
                    else:
                        age_range3.append('')
            #is not a range
            else:
                rage = float(raw_age)
                age.append(rage)
                if rage < 1:
                    age_range1.append('<1y')
                    age_range2.append('<18y')
                    age_range3.append('<18y')
                elif rage < 6:
                    age_range1.append('1-5y')
                    age_range2.append('<18y')
                    age_range3.append('<18y')
                elif rage < 18:
                    age_range1.append('6-17y')
                    age_range2.append('<18y')
                    age_range3.append('<18y')
                elif rage < 65:
                    age_range1.append('18-64y')
                    age_range2.append('18-64y')
                    age_range3.append('>=18y')
                else:
                    age_range1.append('>=65y')
                    age_range2.append('>=65y')
                    age_range3.append('>=18y')

    meta.loc[:, 'age'] = age
    meta.loc[:, 'age_range1'] = age_range1
    meta.loc[:, 'age_range2'] = age_range2
    meta.loc[:, 'age_range3'] = age_range3

    meta.to_csv(args.meta_out, sep='\t', index=False)