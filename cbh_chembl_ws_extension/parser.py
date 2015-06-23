# -*- coding: utf-8 -*-
import dateutil
import jsonpatch
'''
A parser object to work with our custom field configs and uncurated fields

'''

class CovertDateOperation(jsonpatch.PatchOperation):
    """Ensures that a data point is formatted correctly for a date field"""

    def apply(self, obj):
        try:
            from_ptr = jsonpatch.JsonPointer(self.location)
        except KeyError as ex:
            raise jsonpatch.InvalidJsonPatch(
                "The operation does not contain a 'from' member")

        subobj, part = from_ptr.to_last(obj)
        try:
            value = subobj[part]
        except (KeyError, IndexError) as ex:
            raise jsonpatch.JsonPatchConflict(str(ex))

        if isinstance(subobj, jsonpatch.MutableMapping) and \
                self.pointer.contains(from_ptr):
            raise jsonpatch.JsonPatchConflict('Cannot move values into its own children')

        obj = jsonpatch.RemoveOperation({
            'op': 'remove',
            'path': self.location
        }).apply(obj)

        obj = jsonpatch.AddOperation({
            'op': 'add',
            'path': self.location,
            'value': dateutil.parser.parse(value).strftime("%Y-%m-%d")
        }).apply(obj)

        return obj






class MyJsonPatch(jsonpatch.JsonPatch):
    def __init__(self, patch):
        instance = super(MyJsonPatch, self).__init__(patch)
        self.operations["convertdate"] = CovertDateOperation




def apply_json_patch(dictdata, patch):
    mjp = MyJsonPatch(patch)
    try:
        mjp.apply(dictdata, in_place=True)
    except jsonpatch.JsonPatchConflict:
        pass






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
