# -*- coding: utf-8 -*-

'''
A parser object to work with our custom field configs and uncurated fields

'''



def parse_sdf_record(headers, obj, destination_field, mol):
    custom_fields = {}
    for hdr in headers:
        try:
            custom_fields[hdr] = unicode(mol.GetProp(hdr))
        except KeyError:
            custom_fields[hdr]   = u""
                            
    setattr(obj,destination_field, custom_fields)



def parse_pandas_record(headers, obj, destination_field):
    custom_fields = {}
    for hdr in headers:
        if unicode(row[hdr]) == u"nan":
            custom_fields[hdr] = ""
        else:
            custom_fields[hdr] = unicode(row[hdr]) 
    #Set excel fields as uncurated
    setattr(obj,destination_field, custom_fields)