# -*- coding: utf-8 -*-
import dateutil
'''
A parser object to work with our custom field configs and uncurated fields

'''



def parse_sdf_record(headers, obj, destination_field, mol,fielderrors):
    custom_fields = {}

    for hdr in headers:
        try:
            value = unicode(mol.GetProp(hdr)).strip()
            custom_fields[hdr] = value
            test_specific_parse_errors(hdr, value, fielderrors)

        except KeyError:
            pass
            #custom_fields[hdr]   = u""
                            
    setattr(obj,destination_field, custom_fields)





def parse_pandas_record(headers, obj, destination_field, row,fielderrors, headerswithdata):
    custom_fields = {}
    
    for hdr in headers:
        if unicode(row[hdr]) == u"nan":
            pass
            # custom_fields[hdr] = ""
        else:
            value = unicode(row[hdr]).strip()
            if value:
                headerswithdata.add(hdr)
                custom_fields[hdr] = unicode(value)
                test_specific_parse_errors(hdr, value, fielderrors)
    #Set excel fields as uncurated
    setattr(obj,destination_field, custom_fields)





def test_specific_parse_errors(hdr, value, fielderrors):
    if hdr not in fielderrors["stringdate"]:
        try:
            curated_value = dateutil.parser.parse(value).strftime("%Y-%m-%d")
        except:
            fielderrors["stringdate"].add(hdr)
    if hdr not in fielderrors["number"]:      
        try:
            curated_value = float(value)
            if hdr not in fielderrors["integer"]: 
                if "." in value:
                    fielderrors["integer"].add(hdr)
        except:
            fielderrors["integer"].add(hdr)
            fielderrors["number"].add(hdr)
