# -*- coding: utf-8 -*-

import codecs
import csv
import cStringIO

from django.conf import settings
from django.http import Http404, HttpResponse
from django.template import Context
from django.template.loader import get_template
from django.utils import timezone
from tastypie.serializers import Serializer

import re
import json
import xlsxwriter
import os
import xlrd
import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools

import pybel

import copy

def get_field_name_from_key(key):
    return key.replace(u"__space__", u" ")

def get_key_from_field_name(name):
    return name.replace(u" ", u"__space__")






def flatten_dict(d, base=None):
    """Converts a dictionary of dictionaries or lists into a simple
    dictionary.

    For example the following dictionary

    foobar = {'key1': 'value1',
              'key2': {'skey1': 'svalue1'},
              'key3': ['svalue2', 'svalue3']}

    gets converted to

    foobar = {'key1': 'value1',
              'key2.skey1': 'svalue1',
              'key3.0': 'svalue2',
              'key3.1': 'svalue3'}

    """
    new_dict = {}
    for key, value in d.iteritems():
        if isinstance(value, dict):
            new_base = ''
            if base:
                new_base = '%s.' % base
            new_base += key
            new_dict.update(flatten_dict(value, base=new_base))
        elif isinstance(value, list):
            new_base = ''
            if base:
                new_base += '%s.' % base
            new_base += '%s' % key
            i = 0
            for item in value:
                new_base_index = new_base + '.%d' % i
                if isinstance(item, dict):
                    new_dict.update(flatten_dict(item, base=new_base_index))
                else:
                    new_dict.update({new_base_index: item})
                i += 1
        elif base:
            new_dict.update({'%s.%s' % (base, key): value})
        else:
            new_dict.update({key: value})

    return new_dict


class CSVUnicodeWriter(object):
    """
    A CSV writer which will write rows to CSV file "f",
    which is encoded in the given encoding.

    Original code from http://docs.python.org/library/csv.html#csv-examples
    Altered by Giorgos Logiotatidis <giorgos@mozilla.com>

    """

    def __init__(self, f, dialect=csv.excel, encoding='utf-8', **kwds):
        # Redirect output to a queue
        self.queue = cStringIO.StringIO()
        self.writer = csv.writer(self.queue, dialect=dialect, **kwds)
        self.stream = f
        self.encoder = codecs.getincrementalencoder(encoding)()

    def writerow(self, row):
        newrow = []
        for s in row:
            newrow.append(unicode(s))

        self.writer.writerow([s.encode('utf-8') for s in newrow])
        # Fetch UTF-8 output from the queue ...
        data = self.queue.getvalue()
        data = data.decode('utf-8')
        # ... and reencode it into the target encoding
        data = self.encoder.encode(data)
        # write to the target stream
        self.stream.write(data)
        # empty queue
        self.queue.truncate(0)

    def writerows(self, rows):
        for row in rows:
            self.writerow(row)


class CSVSerializer(Serializer):
    """Extend tastypie's serializer to export to CSV format."""
    # formats = ['json', 'jsonp', 'xml', 'yaml', 'html', 'csv']
    # content_types = {'json': 'application/json',
    #                  'jsonp': 'text/javascript',
    #                  'xml': 'application/xml',
    #                  'yaml': 'text/yaml',
    #                  'html': 'text/html',
    #                  'csv': 'text/csv'}

    def to_csv(self, data, options=None):
        """Convert data to CSV."""
        options = options or {}
        data = self.to_simple(data, options)
        raw_data = cStringIO.StringIO()

        writer = CSVUnicodeWriter(raw_data, delimiter=';', quotechar='"',
                                  quoting=csv.QUOTE_MINIMAL)

        for category in data:
            if category == 'objects' and len(data[category]) > 0:
                items = []
                available_keys = []
                for item in data[category]:
                    flatitem = flatten_dict(item)
                    items.append(flatitem)
                    available_keys += [key for key in flatitem.keys()
                                       if key not in available_keys]

                available_keys = sorted(available_keys)
                writer.writerow(available_keys)

                for item in data[category]:
                    flatitem = flatten_dict(item)
                    writer.writerow(map(lambda x: flatitem.get(x),
                                        available_keys))

        raw_data.seek(0)
        return raw_data

class XLSSerializer(Serializer):
    # formats = ['json', 'jsonp', 'xml', 'yaml', 'html', 'csv', 'xls']
    # content_types = {'json': 'application/json',
    #                  'jsonp': 'text/javascript',
    #                  'xml': 'application/xml',
    #                  'yaml': 'text/yaml',
    #                  'html': 'text/html',
    #                  'csv': 'text/csv',
    #                  'xls': 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'}

    def to_xls(self, data, options=None):
        '''write excel file here'''
        output = cStringIO.StringIO()
        
        #make a pandas dataframe from the data here
        #then export as xls or to xlsxwriter
        data = self.to_simple(data, {})
        exp_json = json.loads(data.get('export',[]))
        ordered_fields = [ 'UOx ID', 'SMILES', 'Added By',  'Std InChi', 'Mol Weight', 'alogp'  ]
        headers = data.get('headers', {})
        ordered_fields += headers["custom_fields"]
        ordered_fields += headers["uncurated_fields"]
        df = pd.DataFrame(exp_json)


        df.fillna('', inplace=True)

        cols = df.columns.tolist()
        #now for the list we have in the order we have it, move the columns by name
        #this way you end up with your core fields at the start and custom fields at the end.
        
        for idx, item in enumerate(ordered_fields):
            cols.insert(idx, cols.pop(cols.index(item)))
            


        #reindex the dataframe
        df = df.ix[:, cols]
        widths = []
        for col in df.columns.tolist():
            col = str(col)
            titlewidth = len(col)
            try:
                w = df[col].astype(unicode).str.len().max()
                if w > titlewidth:
                    widths.append(int(w*1.2))
                else:
                    widths.append(int(titlewidth* 1.2))
            except:
                widths.append(int(titlewidth* 1.2))


        writer = pd.ExcelWriter('temp.xlsx', engine='xlsxwriter')
        writer.book.filename = output
        df.to_excel(writer, sheet_name='Sheet1', index=False)
        workbook = writer.book
        format = workbook.add_format()
        worksheet = writer.sheets['Sheet1']
        format.set_text_wrap()
        #make the UOx ID and SMILES columns bigger
        #BUG - can't set column format until pandas 0.16
        #https://github.com/pydata/pandas/issues/9167
        for index, width in enumerate(widths):
            if width > 150:
                width = 150
            elif width < 15:
                width = 15
            worksheet.set_column(index ,index , width)
        writer.save()
        
        return output.getvalue()



SDF_TEMPLATE = ">  <{name}>\n{value}\n\n"



class SDFSerializer(Serializer):
    '''For exporting query sets as SD/Mol files'''
    formats = ['json', 'jsonp', 'xml', 'yaml', 'html', 'csv', 'xls', 'sdf']
    content_types = {'json': 'application/json',
                     'jsonp': 'text/javascript',
                     'xml': 'application/xml',
                     'yaml': 'text/yaml',
                     'html': 'text/html',
                     'csv': 'text/csv',
                     'xls': 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                     'sdf': 'chemical/x-mdl-sdfile'}

    def to_sdf(self, data, options=None):
        '''Convert to SDF'''
        mols = []
        index = 0
        options = options or {}
        exp_json = json.loads(data.get('export',[]))
        df = pd.DataFrame(exp_json)
        df.fillna('', inplace=True)
        cols = df.columns.tolist()
        #now for the list we have in the order we have it, move the columns by name
        #this way you end up with your core fields at the start and custom fields at the end.
        ordered_fields = [ 'UOx ID', 'SMILES', 'Known Drug', 'Added By', 'MedChem Friendly', 'Std InChi', 'Mol Weight', 'alogp'  ]
        for idx, item in enumerate(ordered_fields):
            cols.insert(idx, cols.pop(cols.index(item)))
        #reindex the dataframe
        df = df.ix[:, cols]
        #pull data back out of dataframe to put into rdkit tools

        row_iterator = df.iterrows()
        headers = list(df)

        mol_strings = []
        for index, row in row_iterator:
            #Simple string based SDF formatter
            mol_str = row['ctab'].replace("RDKit          2D\n", "Generated by ChemBio Hub ChemReg http://chembiohub.ox.ac.uk/chemreg\n").split("END")[0]
            properties = []
            for field in headers:
                if field !="ctab":
                    properties.append(
                        SDF_TEMPLATE.format(
                            **{"name" :str(field), "value" : str(row[field])}
                            )
                        )
            mol_strings.append("".join([mol_str,"END\n"] + properties))
        return "$$$$\n".join(mol_strings) + "$$$$"
        




class iCalSerializer(Serializer):
    """Extend tastypie's serializer to export to iCal format."""
    formats = ['json', 'jsonp', 'xml', 'yaml', 'html', 'ical']
    content_types = {'json': 'application/json',
                     'jsonp': 'text/javascript',
                     'xml': 'application/xml',
                     'yaml': 'text/yaml',
                     'html': 'text/html',
                     'ical': 'text/calendar'}

    def to_ical(self, data, options=None):
        """Convert data to iCal."""
        options = options or {}

        if 'error' in data:
            raise Http404

        if isinstance(data, dict) and 'objects' in data:
            events = [event.obj for event in data['objects']]
        else:
            events = [data.obj]

        date_now = timezone.now()
        ical = get_template('multi_event_ical_template.ics')

        return ical.render(Context({'events': events, 'date_now': date_now,
                                    'host': settings.SITE_URL}))

class CamelCaseJSONSerializer(Serializer):
    # formats = ['json']
    # content_types = {
    #     'json': 'application/json',
    # }

    def to_json(self, data, options=None):
        # Changes underscore_separated names to camelCase names to go from python convention to javacsript convention
        data = self.to_simple(data, options)

        def underscoreToCamel(match):
            return match.group()[0] + match.group()[2].upper()

        def camelize(data):
            if isinstance(data, dict):
                new_dict = {}
                for key, value in data.items():
                    new_key = re.sub(r"[a-z]_[a-z]", underscoreToCamel, key)
                    if new_key in ["customFields", "uncuratedFields"]:
                        
                        for k, v in value.iteritems():
                            if isinstance(v, basestring):
                                if  v.startswith("[") and v.endswith("]"):
                                    try:
                                        value[k] = json.loads(v)
                                        continue
                                    except:
                                        pass
                                elif "." in v:
                                    try:
                                        value[k] = float(v)
                                        continue
                                    except:
                                        pass
                                else:
                                    try:
                                        value[k] = int(v)
                                        continue
                                    except:
                                        pass
                                value[k] = v
                        new_dict[new_key] = value
                    else:
                        new_dict[new_key] = camelize(value)
                return new_dict
            if isinstance(data, (list, tuple)):
                for i in range(len(data)):
                    data[i] = camelize(data[i])
                return data
            return data

        camelized_data = camelize(data)
        for key, value in camelized_data.iteritems():
            try:
                dictd = json.loads(value)
                if isinstance(dictd, dict):
                    camelized_data[key] = dictd
            except:
                pass
        return json.dumps(camelized_data, sort_keys=True)

    def from_json(self, content):
        # Changes camelCase names to underscore_separated names to go from javascript convention to python convention
        data = json.loads(content)

        def camelToUnderscore(match):
            return match.group()[0] + "_" + match.group()[1].lower()

        def underscorize(data):
            if isinstance(data, dict):
                new_dict = {}
                for key, value in data.items():
                    new_key = re.sub(r"[a-z][A-Z]", camelToUnderscore, key)
                    if new_key in ["custom_fields", "uncurated_fields", "compoundstats", "batchstats", "warnings", "properties"]:
                        new_dict[new_key] = value

                    else:
                        new_dict[new_key] = underscorize(value)
                return new_dict
            if isinstance(data, (list, tuple)):
                for i in range(len(data)):
                    data[i] = underscorize(data[i])
                return data
            return data

        underscored_data = underscorize(data)

        return underscored_data


class CBHCompoundBatchSerializer(CamelCaseJSONSerializer,CSVSerializer,XLSSerializer,SDFSerializer):
    pass

def convert_query(data):
    if isinstance(data, dict):
        new_dict = {}
        for key, value in data.items():
            new_key = get_key_from_field_name(key)
            new_key = new_key.replace("uncuratedFields","uncurated_fields")
            new_key = new_key.replace("customFields","custom_fields")

            new_dict[new_key] = convert_query(value)
        return new_dict
    if isinstance(data, (list, tuple)):
        for i in range(len(data)):
            data[i] = convert_query(data[i])
        return data
    return data

class CBHCompoundBatchElasticSearchSerializer(Serializer):
    formats = ['json']
    content_types = {
        'json': 'application/json',
    }

    def convert_query(self, es_request):

        def camelToUnderscore(match):
            return match.group()[0] + "_" + match.group()[1].lower()

        es_request["filter"] = convert_query(es_request["filter"])
        es_request["sort"] = convert_query(es_request["sort"])
        newsort = []
        for item in es_request["sort"]:
            newItem = {}
            for sort, direction in item.items():
                if ("." in sort):
                    newItem[sort + ".raw"] = direction
                else:
                    new_key = re.sub(r"[a-z][A-Z]", camelToUnderscore, sort)
                    if new_key != "id":
                        new_key += ".raw"
                    newItem[new_key] = direction
            newsort.append(newItem)
        es_request["sort"] = newsort



    def handle_data_from_django_hstore(self, value):
        '''Hstore passes data in the wrong format'''
        for k, v in value.iteritems():
            if isinstance(v, basestring):
                if  v.startswith("[") and v.endswith("]"):
                    try:
                        value[k] = json.loads(v)
                        continue
                    except:
                        pass
                # elif "." in v:
                #     try:
                #         value[k] = float(v)
                #         continue
                #     except:
                #         pass
                # else:
                #     try:
                #         value[k] = int(v)
                #         continuecompound_stats
                #     except:
                #         pass
                value[k] = unicode(v)

    def to_es_ready_data(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        data['custom_field_list'] = []
        self.handle_data_from_django_hstore(data["custom_fields"])

        for key, value in data["custom_fields"].items():
            if type(value) == list:
                for val in value:
                    if val:
                        val = val.replace(u"\n|\r", " ")
                        data['custom_field_list'].append({'name':key, 'value':val, 'searchable_name': key.split(" ")[0].lower(), 'aggregation': '%s|%s' % (key, val) })
            else:
                if value:
                    value = value.replace(u"\n|\r", " ")
                    data['custom_field_list'].append({'name':key, 'value':value, 'searchable_name': key.split(" ")[0].lower(), 'aggregation': '%s|%s' % (key, value) })
                
        data['custom_field_list'].append({'name': "Project", 'value':data['project'], 'searchable_name': 'project', 'aggregation': '%s|%s' % ('Project', data['project']) })
        data['custom_field_list'].append({'name': "Upload Id", 'value':data['multiple_batch_id'], 'searchable_name': 'upload', 'aggregation': '%s|%d' % ('Upload', data['multiple_batch_id']) })
        

        for key, value in data.items():
            if key in ["custom_fields", "uncurated_fields"]:
                if options and options.get("underscorize", False):
                    data[key] = self.underscorize_fields(value)           
                self.handle_data_from_django_hstore( value)
        return data



    def to_es_ready_non_chemical_data(self, data, options=None):
        options = options or {}
        newdata = {}
        non_chem_fields = ['batch_number', 
                            'blinded_batch_id', 
                            'chemblId', 
                            'created', 
                            'created_by', 
                            'editable_by', 
                            'id', 
                            'modified',
                            'multiple_batch_id',
                            'project',
                            'timestamp',
                            'uncurated_fields',
                            'custom_fields']
        data = self.to_simple(data, options)
        #pull out the non-chemical fields for the main search - 
        #we will send the hits to the backend for any subsequent structure searching
        # for item in non_chem_fields:
        #     if data[item]:
        #         newdata[item] = data[item]
        newdata = copy.deepcopy(data)
        newdata['custom_field_list'] = []
        self.handle_data_from_django_hstore(data["custom_fields"])

        for key, value in data["custom_fields"].items():
            if type(value) == list:
                for val in value:
                    if val:
                        v = val.replace('\n', ' ').replace('\r', '')
                        agg = '%s|%s' % (key, v)
                        newdata['custom_field_list'].append({'name':key, 'value':v, 'searchable_name': key.split(" ")[0].lower(), 'aggregation': agg})
            else:
                if value:
                    v = value.replace('\n', ' ').replace('\r', '')
                    agg = '%s|%s' % (key, v)

                    newdata['custom_field_list'].append({'name':key, 'value':v, 'searchable_name': key.split(" ")[0].lower(), 'aggregation': agg })
                
        newdata['custom_field_list'].append({'name': "Project", 'value':newdata['project'], 'searchable_name': 'project', 'aggregation': '%s|%s' % ('Project', newdata['project']) })
        newdata['custom_field_list'].append({'name': "Upload Id", 'value':newdata['multiple_batch_id'], 'searchable_name': 'upload', 'aggregation': '%s|%d' % ('Upload', newdata['multiple_batch_id']) })
        
        for key, value in data.items():
            if key in ["custom_fields", "uncurated_fields"]:
                if options and options.get("underscorize", False):
                    newdata[key] = self.underscorize_fields(value)           
                self.handle_data_from_django_hstore( value)
        return newdata
    


    def to_python_ready_data(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        data.pop("custom_field_list", False)
        for key, value in data.items():
            if key in ["custom_fields", "uncurated_fields"]:
                data[key] = self.deunderscorize_fields(value)
        return data


    def to_json(self, data, options=None):
        self.to_es_ready_data( data, options=options) 
        return json.dumps(data, sort_keys=True)

    def underscorize_fields(self,dictionary):
        return {
                    get_key_from_field_name(key):value 
                    for key, value in dictionary.items()
                }



    def deunderscorize_fields(self,dictionary):
        return {
                    get_field_name_from_key(key):value 
                    for key, value in dictionary.items()
                }









