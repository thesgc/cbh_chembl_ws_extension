from tastypie.resources import ModelResource, ALL, ALL_WITH_RELATIONS
from django.conf import settings
from chembl_webservices.base import ChEMBLApiBase
from tastypie.utils import trailing_slash
from django.conf.urls import *
from django.core.exceptions import ObjectDoesNotExist
from tastypie.authorization import Authorization
from tastypie import http
from django.http import HttpResponse
import base64
import time
from collections import OrderedDict
from tastypie.resources import ModelResource, Resource


try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import Draw
except ImportError:
    Chem = None
    Draw = None
    AllChem = None

try:
    from rdkit.Chem.Draw import DrawingOptions
except ImportError:
    DrawingOptions = None

try:
    import indigo
    import indigo_renderer
except ImportError:
    indigo = None
    indigo_renderer = None

from chembl_webservices.cache import ChemblCache

from tastypie.exceptions import BadRequest
from chembl_core_db.chemicalValidators import validateSmiles, validateChemblId, validateStandardInchiKey
from tastypie.utils.mime import build_content_type
from tastypie.exceptions import ImmediateHttpResponse
from django.db.utils import DatabaseError
from django.db import transaction
from django.db import connection
from chembl_beaker.beaker.core_apps.rasterImages.impl import _ctab2image
try:
    from chembl_compatibility.models import MoleculeDictionary
    from chembl_compatibility.models import CompoundMols
    from chembl_compatibility.models import MoleculeHierarchy
except ImportError:
    from chembl_core_model.models import MoleculeDictionary
    from chembl_core_model.models import CompoundMols
    from chembl_core_model.models import MoleculeHierarchy

try:
    DEFAULT_SYNONYM_SEPARATOR = settings.DEFAULT_COMPOUND_SEPARATOR
except AttributeError:
    DEFAULT_SYNONYM_SEPARATOR = ','

try:
    WS_DEBUG = settings.WS_DEBUG
except AttributeError:
    WS_DEBUG = False

from chembl_webservices.compounds import CompoundsResource
from chembl_webservices.base import ChEMBLApiSerializer
from cbh_chembl_ws_extension.base import CBHApiBase
from cbh_chembl_ws_extension.serializers import CBHCompoundBatchSerializer

from tastypie.utils import dict_strip_unicode_keys
from tastypie.serializers import Serializer
from django.core.serializers.json import DjangoJSONEncoder
from tastypie import fields, utils
from cbh_chembl_model_extension.models import CBHCompoundBatch, CBHCompoundMultipleBatch
from tastypie.authentication import SessionAuthentication
import json
from tastypie.paginator import Paginator
from chembl_beaker.beaker.core_apps.conversions.impl import _smiles2ctab, _apply

from flowjs.models import FlowFile
import xlrd
import pandas as pd
import numpy as np
import urllib




class CBHCompoundsReadResource(CBHApiBase, CompoundsResource):

#-----------------------------------------------------------------------------------------------------------------------

    class Meta:
        resource_name = 'compounds'
        authorization = Authorization()
        include_resource_uri = False
        paginator_class = None
        serializer = ChEMBLApiSerializer('compound')
        allowed_methods = ['get']
        default_format = 'application/xml'
        cache = ChemblCache()

#-----------------------------------------------------------------------------------------------------------------------




    def base_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/stdinchikey/(?P<stdinchikey>\w[\w-]*)%s$" % (
                self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'),
                name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/stdinchikey/(?P<stdinchikey>\w[\w-]*)\.(?P<format>json|xml)%s$" % (
                self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'),
                name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/smiles/?(?P<smiles>[\S]*)\.(?P<format>json|xml)%s$" % (
                self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_list'),
                name="api_dispatch_list"),
            url(r"^(?P<resource_name>%s)/smiles/?(?P<smiles>[\S]*)%s$" % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'),
                name="api_dispatch_list"),
            url(r"^(?P<resource_name>%s)/substructure/?(?P<smiles>[\S]*)\.(?P<format>json|xml)%s$" % (
                self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_sim'),
                name="api_get_substructure"),
            url(r"^(?P<resource_name>%s)/substructure/?(?P<smiles>[\S]*)%s$" % (
                self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_sim'),
                name="api_get_substructure"),
            url(r"^(?P<resource_name>%s)/similarity/(?P<smiles>[\S]*)/(?P<simscore>\d+)\.(?P<format>json|xml)%s$" % (
                self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_sim'),
                name="api_get_similarity"),
            url(r"^(?P<resource_name>%s)/similarity/(?P<smiles>[\S]*)/(?P<simscore>\d+)%s$" % (
                self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_sim'),
                name="api_get_similarity"),
            url(r"^(?P<resource_name>%s)/similarity.(?P<format>json|xml)%s$" % (
                self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_sim'),
                name="api_get_similarity"),
            url(r"^(?P<resource_name>%s)/similarity%s$" % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_sim'),
                name="api_get_similarity"),
            url(r"^(?P<resource_name>%s)/(?P<chemblid>\w[\w-]*)\.(?P<format>json|xml)%s$" % (
                self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'),
                name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<chemblid>\w[\w-]*)%s$" % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'),
                name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<chemblid>\w[\w-]*)/image%s$" % (
                self._meta.resource_name, trailing_slash()), self.wrap_view('get_cached_image'),
                name="api_get_image"),
            # url(r"^(?P<resource_name>%s)$" % self._meta.resource_name, self.wrap_view('dispatch_compounds'),
            #     name="api_dispatch_compounds"),
            url(r"^(?P<resource_name>%s)$" % self._meta.resource_name, self.wrap_view('post_list'),
                name="api_post_list"),
            url(r"^(?P<resource_name>%s)\.(?P<format>json|xml|csv|xls|sdf)$" % self._meta.resource_name,
                self.wrap_view('dispatch_compounds'), name="api_dispatch_compounds"),
        ]







class CBHCompoundBatchResource(ModelResource):

    class Meta:
        filtering = {
            "std_ctab": ALL_WITH_RELATIONS,
            "ctab": ALL,
            "multiple_batch_id": ALL_WITH_RELATIONS,
        }
        always_return_data = True
        prefix = "related_molregno"
        fieldnames = [('chembl_id', 'chemblId'),
                  ('pref_name', 'preferredCompoundName'),
                  ('max_phase', 'knownDrug'),
                  ('compoundproperties.med_chem_friendly', 'medChemFriendly'),
                  ('compoundproperties.ro3_pass', 'passesRuleOfThree'),
                  ('compoundproperties.full_molformula', 'molecularFormula'),
                  ('compoundstructures.canonical_smiles', 'smiles'),
                  ('compoundstructures.standard_inchi_key', 'stdInChiKey'),
                  ('compoundproperties.molecular_species', 'species'),
                  ('compoundproperties.num_ro5_violations', 'numRo5Violations'),
                  ('compoundproperties.rtb', 'rotatableBonds'),
                  ('compoundproperties.mw_freebase', 'molecularWeight'),
                  ('compoundproperties.alogp', 'alogp'),
                  ('compoundproperties.acd_logp', 'acdLogp'),
                  ('compoundproperties.acd_logd', 'acdLogd'),
                  ('compoundproperties.acd_most_apka', 'acdAcidicPka'),
                  ('compoundproperties.acd_most_bpka', 'acdBasicPka')]
        csv_fieldnames = [('chembl_id', 'UOX ID'),
                  ('pref_name', 'Preferred Name'),
                  ('max_phase', 'Known Drug'),
                  ('compoundproperties.med_chem_friendly', 'MedChem Friendly'),
                  ('compoundproperties.ro3_pass', 'passesRuleOfThree'),
                  ('compoundproperties.full_molformula', 'Mol Formula'),
                  ('compoundstructures.canonical_smiles', 'SMILES'),
                  ('compoundstructures.standard_inchi_key', 'Std InChiKey'),
                  ('compoundproperties.num_ro5_violations', 'Rule of 5 violations'),
                  ('compoundproperties.rtb', 'Rotatable Bonds'),
                  ('compoundproperties.mw_freebase', 'Mol Weight')]
        
        queryset = CBHCompoundBatch.objects.all()
        resource_name = 'cbh_compound_batches'
        authorization = Authorization()
        include_resource_uri = False
        serializer = CBHCompoundBatchSerializer()
        allowed_methods = ['get', 'post', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        paginator_class = Paginator


    def post_validate(self, request, **kwargs):
        """Runs the validation for a single or small set of molecules"""
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)

        updated_bundle = self.obj_build(bundle, dict_strip_unicode_keys(deserialized))
        bundle.obj.validate()
        dictdata = bundle.obj.__dict__
        dictdata.pop("_state")
        updated_bundle = self.build_bundle(obj=bundle.obj, data=dictdata)
        return self.create_response(request, updated_bundle, response_class=http.HttpAccepted)

    def get_project_custom_field_names(self, request, **kwargs):
        # deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'text/plain'))
        
        # deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(request=request)

        fields = CBHCompoundBatch.objects.get_all_keys()
        bundle.data['field_names'] =[{'name': item, 'count': 1, 'last_used': ''} for item in fields]        
        return self.create_response(request, bundle, response_class=http.HttpAccepted)


        #return HttpResponse("{ 'field_names': [ {'name': 'test1', 'count': 1, 'last_used': ''}, {'name': 'test2', 'count': 1, 'last_used': ''} ] }")


    def save_related(self, bundle):
        #bundle.obj.created_by = request.user.
        bundle.obj.generate_structure_and_dictionary()
        

    def full_hydrate(self, bundle):
        '''As the object is created we run the validate code on it'''
        bundle = super(CBHCompoundBatchResource, self).full_hydrate(bundle)
        bundle.obj.validate()
        return bundle


    def obj_build(self, bundle, kwargs):
        """
        A ORM-specific implementation of ``obj_create``.
        """
        bundle.obj = self._meta.object_class()
        for key, value in kwargs.items():
            setattr(bundle.obj, key, value)
        setattr(bundle.obj, "id", -1)
        
        return bundle

    def prepend_urls(self):
        return [
        url(r"^(?P<resource_name>%s)/validate/$" % self._meta.resource_name,
                self.wrap_view('post_validate'), name="api_validate_compound_batch"),
        url(r"^(?P<resource_name>%s)/validate_list/$" % self._meta.resource_name,
                self.wrap_view('post_validate_list'), name="api_validate_compound_list"),
        url(r"^(?P<resource_name>%s)/existing/$" % self._meta.resource_name,
                self.wrap_view('get_project_custom_field_names'), name="api_batch_existing_fields"),
                

        url(r"^(?P<resource_name>%s)/multi_batch_save/$" % self._meta.resource_name,
                self.wrap_view('multi_batch_save'), name="multi_batch_save"),
        url(r"^(?P<resource_name>%s)/multi_batch_custom_fields/$" % self._meta.resource_name,
                self.wrap_view('multi_batch_custom_fields'), name="multi_batch_custom_fields"),
        url(r"^(?P<resource_name>%s)/validate_files/$" % self._meta.resource_name,
                self.wrap_view('post_validate_files'), name="api_compound_validate_files"),
        url(r"^(?P<resource_name>%s)/export_file/$" % self._meta.resource_name,
                self.wrap_view('export_file'), name="api_compound_export_file"),
        url(r"^(?P<resource_name>%s)/smiles2svg/(?P<structure>\w[\w-]*)/$" % 
                self._meta.resource_name, self.wrap_view('get_image_from_pipe'),name="smiles2svg"),
        ]

    def multi_batch_save(self, request, **kwargs):
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        id = bundle.data["current_batch"]
        batches = CBHCompoundMultipleBatch.objects.get(pk=id).uploaded_data
        bundle.data["saved"] = 0
        bundle.data["errors"] = []
        print "this many batches"
        print len(batches)
        for batch in batches:


                batch.save(validate=False)
                batch.generate_structure_and_dictionary()
                bundle.data["saved"] += 1
  #          except Exception , e:
   #             bundle.data["errors"] += e

        return self.create_response(request, bundle, response_class=http.HttpCreated)


    def multi_batch_custom_fields(self, request, **kwargs):
        '''Save custom fields from the mapping section when adding ID/SMILES list'''
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)

        id = bundle.data["current_batch"]

        mb = CBHCompoundMultipleBatch.objects.get(pk=id)
        for b in mb.uploaded_data:
            b.custom_fields = bundle.data["custom_fields"]
        mb.save()

        return self.create_response(request, bundle, response_class=http.HttpAccepted)



    def validate_multi_batch(self,multi_batch, bundle, request):
        total = len(multi_batch.uploaded_data)
        bundle.data["objects"] = {"pains" :[], "changed" : [], "errors" :[]}
        for batch in multi_batch.uploaded_data:

            batch = batch.__dict__
            batch.pop("_state")
            if batch["warnings"]["pains_count"] != "0":
                bundle.data["objects"]["pains"].append(batch)
            if batch["errors"] != {}:
                bundle.data["objects"]["errors"].append(batch)
                total = total - 1

            # if batch["warnings"]["hasChanged"].lower() == "true":
            #     bundle.data["objects"]["changed"].append(batch)  

        bundle.data["objects"]["total"] = total

        bundle.data["current_batch"] = multi_batch.pk
        return self.create_response(request, bundle, response_class=http.HttpAccepted)



    def post_validate_list(self, request, **kwargs):
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        type = bundle.data.get("type",None).lower()
        objects =  bundle.data.get("objects", [])
       
        batches = []
        if type == "smiles":
            batches = [CBHCompoundBatch.objects.from_rd_mol(Chem.MolFromSmiles(obj), smiles=obj) for obj in objects ]

        elif type == "inchi":
            mols = _apply(objects,Chem.MolFromInchi)
            batches = [CBHCompoundBatch.objects.from_rd_mol(mol, smiles=Chem.MolToSmiles(mol)) for mol in mols]
        # for b in batches:
        #     b["created_by"] = request.user.username

        multiple_batch = CBHCompoundMultipleBatch.objects.create()
        for b in batches:
            b.multiple_batch_id = multiple_batch.pk

        multiple_batch.uploaded_data=batches
        multiple_batch.save()
        return self.validate_multi_batch(multiple_batch, bundle, request)


    def post_validate_files(self, request, **kwargs):
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)

        file_name = bundle.data['file_name']
        mappings = bundle.data['mapping']
        correct_file = FlowFile.objects.filter(original_filename=file_name)[0]

        batches = []
        headers = []

        if (correct_file.extension == ".sdf"):
            #read in the file
            suppl = Chem.ForwardSDMolSupplier(correct_file.file)

            #read the headers from the first molecule
            for mol in suppl:
                if mol is None: continue
                if not headers: 
                    headers = list(mol.GetPropNames())
                b = CBHCompoundBatch.objects.from_rd_mol(mol)
                custom_fields = {}
                for hdr in headers:
                    if hdr in mappings["ignored_fields"]:
                        continue
                    elif hdr in mappings["new_fields"]:
                        custom_fields[hdr] = mol.GetProp(hdr) 
                    else:
                        for key, mapping in mappings["remapped_fields"].iteritems():
                            if hdr in mapping:
                                custom_fields[key] = mol.GetProp(hdr)

                b.custom_fields = custom_fields
                batches.append(b)
                

        elif(correct_file.extension in (".xls", ".xlsx")):
            #we need to know which column contains structural info - this needs to be defined on the mapping page and passed here
            #read in the specified structural column
            structure_col = bundle.data["struc_col"]
            df = pd.read_excel(correct_file.file)
            #read the smiles string value out of here, when we know which column it is.
            row_iterator = df.iterrows()
            headers = list(df)
            for index, row in row_iterator:
                smiles_str = row[structure_col]
                b = CBHCompoundBatch.objects.from_rd_mol(Chem.MolFromSmiles(smiles_str), smiles=smiles_str)
                #work out custom fields from mapping object
                #new_fields, remapped_fields, ignored_fields
                
                custom_fields = {}
                for hdr in headers:
                    if hdr in mappings["ignored_fields"]:
                        continue
                    elif hdr in mappings["new_fields"]:
                        custom_fields[ hdr] = row[hdr] 
                    else:
                        for key, mapping in mappings["remapped_fields"].iteritems():
                            if hdr in mapping:
                                custom_fields[key] = row[hdr]

                b.custom_fields = custom_fields
                batches.append(b)


        multiple_batch = CBHCompoundMultipleBatch.objects.create()
        for b in batches:
            b.multiple_batch_id = multiple_batch.pk

        multiple_batch.uploaded_data=batches
        multiple_batch.save()
        return self.validate_multi_batch(multiple_batch, bundle, request)

    def alter_list_data_to_serialize(self, request, data):
        '''use the request type to determine which fields should be limited for file download,
           add extra fields if needed (eg images) and enumerate the custom fields into the 
           rest of the calculated fields'''
        # for names in self.Meta.fieldnames:
        #     data[names[1]] = deepgetattr(data, names[0], None)
        #objs = json.loads(data["objects"])
        # for key, val in json.loads(data["objects"]).iteritems():
        #     print(val)
        print(self.determine_format(request))

        if(self.determine_format(request) == ('application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' or 'chemical/x-mdl-sdfile') ):
            fields_to_keep = {'chemblId':'UOx ID',
                              'canonical_smiles':'SMILES',
                              'knownDrug':'Known Drug',
                              'medChemFriendly':'MedChem Friendly',
                              'standard_inchi':'Std InChi',
                              'rtb':'Rotatable Bonds',
                              'molecularWeight':'Mol Weight',
                              'molecularFormula':'Mol Formula',
                              'acdLogp': 'alogp',
                              'custom_fields':'custom_fields',}

            for index, b in enumerate(data["objects"]):
                #print(b.data['standard_inchi'])
                #remove items which are not listed as being kept
                new_data = {}
                for k, v in b.data.iteritems():
                    for name, display_name in fields_to_keep.iteritems():
                        if k == name:
                            #b.data[display_name] = v
                            #del(b.data[k])
                            new_data[display_name] = v



                for field, value in b.data['custom_fields'].iteritems():
                    new_data[field] = value
                #now remove custom_fields
                del(new_data['custom_fields'])
                b.data = new_data


        return data

    def create_response(self, request, data, response_class=HttpResponse, **response_kwargs):
        """
        Extracts the common "which-format/serialize/return-response" cycle.
        Mostly a useful shortcut/hook.
        """

        desired_format = self.determine_format(request)
        serialized = self.serialize(request, data, desired_format)
        rc = response_class(content=serialized, content_type=build_content_type(desired_format), **response_kwargs)

        if(desired_format == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'):
            rc['Content-Disposition'] = 'attachment; filename=export.xlsx'
        elif(desired_format == 'chemical/x-mdl-sdfile'):
            rc['Content-Disposition'] = 'attachment: filename=export.sdf'
        return rc


    def dehydrate(self, bundle):
        
        try:
            data = bundle.obj.related_molregno
            for names in self.Meta.fieldnames:
                bundle.data[names[1]] = deepgetattr(data, names[0], None)

            mynames = ["editable_by","viewable_by", "warnings", "properties", "custom_fields", "errors"]
            for name in mynames:
                bundle.data[name] = json.loads(bundle.data[name]) 
        except:
            pass
    
        return bundle

    def get_list(self, request, **kwargs):
        """
        Returns a serialized list of resources.
        Calls ``obj_get_list`` to provide the data, then handles that result
        set and serializes it.
        Should return a HttpResponse (200 OK).
        """
        # TODO: Uncached for now. Invalidation that works for everyone may be
        #       impossible.
        base_bundle = self.build_bundle(request=request)
        objects = self.obj_get_list(bundle=base_bundle, **self.remove_api_resource_names(kwargs))
        sorted_objects = self.apply_sorting(objects, options=request.GET)

        paginator = self._meta.paginator_class(request.GET, sorted_objects, resource_uri=self.get_resource_uri(), limit=self._meta.limit, max_limit=self._meta.max_limit, collection_name=self._meta.collection_name)
        to_be_serialized = paginator.page()

        # Dehydrate the bundles in preparation for serialization.
        bundles = []

        for obj in to_be_serialized[self._meta.collection_name]:
            bundle = self.build_bundle(obj=obj, request=request)
            bundles.append(self.full_dehydrate(bundle, for_list=True))

        to_be_serialized[self._meta.collection_name] = bundles
        to_be_serialized = self.alter_list_data_to_serialize(request, to_be_serialized)
        return self.create_response(request, to_be_serialized)




    def get_object_list(self, request):
        return super(CBHCompoundBatchResource, self).get_object_list(request).select_related("related_molregno", "related_molregno__compound_properties")





def deepgetattr(obj, attr, ex):
    """Recurses through an attribute chain to get the ultimate value."""
    try:
        return reduce(getattr, attr.split('.'), obj)

    except:
        return ex










class CBHCompoundBatchUpload(ModelResource):

    class Meta:
        always_return_data = True
        fieldnames = [('chembl_id', 'chemblId')],
        queryset = FlowFile.objects.all()
        resource_name = 'cbh_batch_upload'
        authorization = Authorization()
        include_resource_uri = False
        allowed_methods = ['get', 'post', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()

    def prepend_urls(self):
        return [
        url(r"^(?P<resource_name>%s)/headers/$" % self._meta.resource_name,
                self.wrap_view('return_headers'), name="api_compound_batch_headers"),
        ]

    def return_headers(self, request, **kwargs):
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)

        request_json = bundle.data

        file_name = request_json['file_name']
        correct_file = self.get_object_list(request).filter(original_filename=file_name)[0]
        headers = []
        header_json = { }
        #get this into a datastructure if excel

        #or just use rdkit if SD file
        if (correct_file.extension == ".sdf"):
            #read in the file
            suppl = Chem.ForwardSDMolSupplier(correct_file.file)

            #read the headers from the first molecule

            for mol in suppl:
                if mol is None: continue
                if not headers: 
                    headers = list(mol.GetPropNames())
                    break

        elif(correct_file.extension in (".xls", ".xlsx")):
            #read in excel file, use pandas to read the headers
            df = pd.read_excel(correct_file.file)
            headers = list(df)

        #this converts to json in preparation to be added to the response
        bundle.data["headers"] = headers

            #send back

        return self.create_response(request, bundle, response_class=http.HttpAccepted)

    

