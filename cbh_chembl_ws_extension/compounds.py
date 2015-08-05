from tastypie.resources import ModelResource, ALL, ALL_WITH_RELATIONS
from django.conf import settings
from django.conf.urls import *
from django.core.exceptions import ObjectDoesNotExist
from tastypie.authorization import Authorization
from tastypie import http
from django.http import HttpResponse
import base64
import time
from collections import OrderedDict
from tastypie.resources import ModelResource, Resource
from itertools import chain
from pybel import readfile , readstring
import re
import shortuuid
from dateutil.parser import parse
import copy
import elasticsearch_client
from difflib import SequenceMatcher as fuzzymatch
import re
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


from tastypie.exceptions import BadRequest
from chembl_core_db.chemicalValidators import validateSmiles, validateChemblId, validateStandardInchiKey



from tastypie.exceptions import ImmediateHttpResponse
from django.db.utils import DatabaseError
from django.db import transaction
from django.db import connection

try:
    from chembl_compatibility.models import MoleculeDictionary
    from chembl_compatibility.models import CompoundMols
    from chembl_compatibility.models import MoleculeHierarchy
except ImportError:
    from chembl_core_model.models import MoleculeDictionary
    from chembl_core_model.models import CompoundMols
    from chembl_core_model.models import MoleculeHierarchy

try:
    WS_DEBUG = settings.WS_DEBUG
except AttributeError:
    WS_DEBUG = False

from cbh_chembl_ws_extension.authorization import ProjectAuthorization
from cbh_chembl_ws_extension.projects import ProjectResource
from cbh_chembl_ws_extension.serializers import CBHCompoundBatchSerializer, CBHCompoundBatchElasticSearchSerializer, get_key_from_field_name
from chembl_business_model.models import CompoundStructures
#from cbh_chembl_ws_extension.base import NBResource
from tastypie.utils import dict_strip_unicode_keys
from tastypie.serializers import Serializer
from tastypie import fields, utils
from cbh_chembl_model_extension.models import CBHCompoundBatch, CBHCompoundMultipleBatch, Project, PinnedCustomField
from tastypie.authentication import SessionAuthentication
import json
from tastypie.paginator import Paginator

from flowjs.models import FlowFile
import xlrd
import pandas as pd
import numpy as np
import urllib
from tastypie.validation import Validation

from django.db.models import Max, Q

from tastypie.serializers import Serializer

from tastypie.authentication import SessionAuthentication

import  chemdraw_reaction

from django.contrib.auth import get_user_model

from rdkit.Chem.AllChem import Compute2DCoords

from django.db.models import Prefetch

import dateutil.parser

            
from cbh_chembl_ws_extension.parser import parse_pandas_record, parse_sdf_record, apply_json_patch



# from tastypie.utils.mime import build_content_type
def build_content_type(format, encoding='utf-8'):
    """
    Appends character encoding to the provided format if not already present.
    """
    if 'charset' in format:
        return format

    return "%s; charset=%s" % (format, encoding)




class CBHCompoundBatchResource(ModelResource):
    project = fields.ForeignKey(ProjectResource, 'project', blank=False, null=False)
    substructure_smarts = ""
    class Meta:
        filtering = {
            "std_ctab": ALL_WITH_RELATIONS,
            "ctab": ALL,
            "multiple_batch_id": ALL_WITH_RELATIONS,
            "project": ALL_WITH_RELATIONS,
            "with_substructure": ALL_WITH_RELATIONS,
            "similar_to": ALL_WITH_RELATIONS,
            "flexmatch": ALL_WITH_RELATIONS,
            "created": ['gte','lte'],
            "created_by": ALL_WITH_RELATIONS,
            "excludes": ALL_WITH_RELATIONS,
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
                  ('compoundproperties.acd_most_bpka', 'acdBasicPka'),
                  ('compoundproperties.full_molformula', 'fullMolformula')]
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
        fields_to_keep = {'chemblId':'UOx ID',
                              'id':'Batch ID',
                              'canonical_smiles':'SMILES',
                              'created_by': 'Added By',
                              'knownDrug':'Known Drug',
                              'medChemFriendly':'MedChem Friendly',
                              'standard_inchi':'Std InChi',
                              'molecularWeight':'Mol Weight',
                              'molecularFormula':'Mol Formula',
                              'acdLogp': 'alogp',
                              'custom_fields':'custom_fields',}
        ordrered_ftk = OrderedDict([('chemblId','UOx ID'),
                                     ('canonical_smiles','SMILES'),
                                     ('knownDrug','Known Drug'),
                                     ('medChemFriendly','MedChem Friendly'),
                                     ('standard_inchi','Std InChi'),
                                     ('molecularWeight','Mol Weight'),
                                     ('acdLogp', 'alogp')])
        
        
        queryset = CBHCompoundBatch.objects.all()
        resource_name = 'cbh_compound_batches'
        authorization = ProjectAuthorization()
        include_resource_uri = False
        serializer = CBHCompoundBatchSerializer()
        allowed_methods = ['get', 'post', 'put', 'patch']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        paginator_class = Paginator
        #elasticsearch config items
        es_index_name = "chemreg_chemical_index"

    

    def apply_filters(self, request, applicable_filters):
        """
        An ORM-specific implementation of ``apply_filters``.
        The default simply applies the ``applicable_filters`` as ``**kwargs``,
        but should make it possible to do more advanced things.
        """
        self.substructure_smarts = ""
        ws = request.GET.get("with_substructure", None)
        st = request.GET.get("similar_to", None)
        fm = request.GET.get("flexmatch", None)
        #this is the similarity index for fingerprint-like searching
        fp = request.GET.get("fpValue", None)
        cms = None
        pids = self._meta.authorization.project_ids(request)
        #get filters which indicate blanks should be ignored
        if ws:
            smiles = self.convert_mol_string(ws)
            cms = CompoundMols.objects.with_substructure(smiles)
        elif st:
            smiles = self.convert_mol_string(st)
            if fp == None:
                cms = CompoundMols.objects.similar_to(smiles,90)
            else:
                cms = CompoundMols.objects.similar_to(smiles,fp)
        elif fm:
            smiles = self.convert_mol_string(fm)
            cms = CompoundMols.objects.flexmatch(smiles)

        eb = request.GET.get("excludeBlanks")
        if eb:
            print("exclude blanks") 
        #else:
        #    cms = CompoundMols.objects.all()
        #To be generalised
        cust = request.GET.get("search_custom_fields__kv_any", None)
        #initialise this with project ids, sets up the correct AND initialisation with custom fields
        cust_queries = Q(project_id__in=set(pids))


        if cust:
            #loop through the custom fields
            #split on pipe (|)
            #put items which are from the same custom field into an OR Q query
            #https://docs.djangoproject.com/en/1.7/topics/db/queries/#complex-lookups-with-q-objects

            cfields = json.loads(cust)
            grouped_fields = {}
            for cfield in cfields:
                cfield_parts = cfield.split("|")
                if grouped_fields.has_key(cfield_parts[0]):
                    grouped_fields[cfield_parts[0]].append(cfield)
                else:
                    grouped_fields[cfield_parts[0]] = [cfield]
            
            #grouped_fields = json.dumps(grouped_fields)
            for key,val in grouped_fields.iteritems():
                print val
                field_specific_queries = [Q(custom_fields__kv_single=value) for value in val]
                #initialise with the first object
                inner_queries = field_specific_queries.pop()
                for item in field_specific_queries:
                    #OR the subqueries from the same custom field column
                    inner_queries |= item
                #AND the sets of custom field queries
                cust_queries &= inner_queries

            #applicable_filters["custom_fields__kv_any"] = cust
        if cms != None:
            #run the sql for pulling in new compounds into compound_mols
            indexed = CBHCompoundBatch.objects.index_new_compounds()
            applicable_filters["related_molregno_id__in"] = cms.values_list("molecule_id", flat=True)

        if request.GET.get("related_molregno__chembl__chembl_id__in", None):
            applicable_filters["related_molregno__chembl__chembl_id__in"] = request.GET.get("related_molregno__chembl__chembl_id__in").split(",")
        
        dateend = applicable_filters.get("created__lte", None)
        if dateend:
            applicable_filters["created__lte"] += " 23:59:59"
        dataset = self.get_object_list(request).filter(**applicable_filters).filter(cust_queries)
        func_group = request.GET.get("functional_group", None)
        
        if func_group:
            funccms = CompoundMols.objects.with_substructure(func_group)
            dataset = dataset.filter(related_molregno_id__in=funccms.values_list("molecule_id", flat=True))

        return dataset.order_by("-id")
    
    

    def get_chembl_ids(self, request, **kwargs):
        '''Get a single list of pinned fields for the project previously listed all custom fields in the DB but this was unwieldy'''
        bundle = self.build_bundle(request=request)
        pids = self._meta.authorization.project_ids(request)
        filters = {"project__id__in" : pids}
        prefix = request.GET.get("chembl_id__chembl_id__startswith", None).upper()
        desired_format = self.determine_format(request)

        if(prefix):
            filters["chembl_id__chembl_id__startswith"] = prefix
            uox_ids = list(MoleculeDictionary.objects.filter(**filters).values_list("chembl_id", flat=True)[0:20])
            bundle.data = [{"value" :uox, "label" : uox} for uox in uox_ids]
            serialized = json.dumps(bundle.data)
        else:
            serialized = "[]"

       
        rc = HttpResponse(content=serialized, content_type=build_content_type(desired_format), )

        return rc

    def get_elasticsearch_ids(self, request, **kwargs):
        bundle = self.build_bundle(request=request)
        pids = self._meta.authorization.project_ids(request)
        filters = {"project__id__in" : pids}
        prefix = request.GET.get("chembl_id__chembl_id__startswith", None).upper()
        desired_format = self.determine_format(request)

        if(prefix):
            #filters["chembl_id__chembl_id__startswith"] = prefix
            #uox_ids = list(MoleculeDictionary.objects.filter(**filters).values_list("chembl_id", flat=True)[0:20])
            uox_ids = list(elasticsearch_client.get_autocomplete(pids, prefix, 'chemblId'))
            bundle.data = [{"value" :uox, "label" : uox} for uox in uox_ids]
            serialized = json.dumps(bundle.data)
        else:
            serialized = "[]"

       
        rc = HttpResponse(content=serialized, content_type=build_content_type(desired_format), )

        return rc


    def get_elasticsearch_autocomplete(self, request, **kwargs):
        
        bundle = self.build_bundle(request=request)
        pids = self._meta.authorization.project_ids(request)
        prefix = request.GET.get("custom__field__startswith", None)
        custom_field = request.GET.get("custom_field", None)
        desired_format = self.determine_format(request)
        #send these project ids to the elasticsearch query?
            #filters["search_custom_fields__kv_any"] = prefix
        uox_ids = list(elasticsearch_client.get_autocomplete(pids, 
                                                             prefix, 
                                                             'custom_field_list.aggregation', 
                                                             custom_fields=True, 
                                                             single_field=custom_field))
        #bundle.data = ["%s|%s" % (uox, uox) for uox in uox_ids]
        bundle.data = [{"value" :uox, "label" : self.labelify_aggregate(uox, single_field=custom_field)} for uox in uox_ids]
        serialized = json.dumps(bundle.data)
       
       
        rc = HttpResponse(content=serialized, content_type=build_content_type(desired_format), )
        
        return rc

    def labelify_aggregate(self, agg, single_field=None):
        #remove the pipe and add square brackets to the custom field type
        splits = agg.split("|")
        label = agg
        #if this is for a single field or custom field, we don't need to show the facet name
        if(len(splits) > 1 and single_field==None):
            label = '%s: %s' % (splits[0], splits[1])
        elif(len(splits) > 1):
            label = splits[1]
        return label

    def reindex_elasticsearch(self, request, **kwargs):
        
        desired_format = self.determine_format(request)
        batches = self.get_object_list(request)
        #we only want to store certain fields in the search index
        batch_dicts = self.batches_to_es_ready(batches, request, non_chem_data_only=True)
        #reindex compound data
        index_name = elasticsearch_client.get_main_index_name()
        es_reindex = elasticsearch_client.create_temporary_index(batch_dicts, request, index_name)

        return HttpResponse(content=json.dumps({"data" :es_reindex}), content_type=build_content_type(desired_format) )

    def reindex_compound(self, request, **kwargs):
        #call this when we need to re-index a compound record which has had fields edited
        desired_format = self.determine_format(request)
        deserialized = self.deserialize(request, request.body, format=desired_format)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        id = bundle.data.get("id", None)
        dataset = self.get_object_list(request).filter(id=id)
        batch_dicts = self.batches_to_es_ready(dataset, request, non_chem_data_only=True)
        es_reindex = elasticsearch_client.reindex_compound(batch_dicts[0], id)

        return HttpResponse(content=json.dumps({"data" :es_reindex}), content_type=build_content_type(desired_format) )


    def convert_mol_string(self, strn):
        #commit
        try:
            mol = Chem.MolFromMolBlock(strn)
            smiles = Chem.MolToSmiles(mol)
        except:
            smiles = strn
        self.substructure_smarts = smiles 
        return smiles


    def match_list_to_moleculedictionaries(self, batch, project, structure_type="MOL"):
        if structure_type == "MOL":
            structure_key = batch.standard_inchi_key
        else:
            raise NotImplemented
        proj_data = MoleculeDictionary.objects.by_project_and_natural_key(structure_type,
                                                                    structure_key,  
                                                                    project.pk)
        same_project = proj_data.values("molregno", "chembl", "created_by")
        
        pub = MoleculeDictionary.objects.by_natural_key_public_except_project(
                                            structure_type, 
                                            structure_key, 
                                           project.pk).order_by("-insert_date")
        different_project_but_public = pub.values("molregno", "chembl", "created_by")
        all_items = list(same_project) + list(different_project_but_public)
        linkedproject = 0
        linkedpublic = 0
        new = 0

        for item in all_items:
            item["tobelinked"] = False

        if same_project.count() > 0:
            #Increment the forced registration number compared to what is already in the database as this can then be used to force the registration of the molecule
            forced_reg_no = proj_data.aggregate(Max('forced_reg_index'))["forced_reg_index__max"] + 1
            linkedproject += 1
            batch.warnings["forced_reg_no"] = forced_reg_no
        elif different_project_but_public.count() > 0:
            linkedpublic += 1        
        else:
            new += 1
        batch.warnings["linkedpublic"] = linkedpublic
        batch.warnings["linkedproject"] = linkedproject
        if len(all_items) > 0:
            all_items[0]["tobelinked"] = True
        batch.warnings["linkable_molecules"] = all_items


    def post_validate(self, request, **kwargs):
        """Runs the validation for a single or small set of molecules"""
        #self.authorized_update_detail(self.get_object_list(bundle.request), bundle)



        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(self.get_object_list(bundle.request), bundle)

        #updated_bundle = self.obj_build(bundle, dict_strip_unicode_keys(deserialized))
        try:
            m = Chem.MolFromMolBlock(bundle.data["ctab"])                
            if not m:
                raise Exception("valancy_or_other_error")
        except:
            raise Exception("valancy_or_other_error")

        obj = CBHCompoundBatch.objects.from_rd_mol(m, orig_ctab=bundle.data["ctab"], )
        bundle.obj = obj
        bundle.obj.validate()
        self.match_list_to_moleculedictionaries(bundle.obj,bundle.data["project"] )
        dictdata = bundle.obj.__dict__
        dictdata.pop("_state")
        
        updated_bundle = self.build_bundle(obj=bundle.obj, data=dictdata)
        return self.create_response(request, updated_bundle, response_class=http.HttpAccepted)



    def save_related(self, bundle):

        bundle.obj.generate_structure_and_dictionary()

 #   def save_related(self, bundle):
 #       #For a single molecule, set the chirality
 #       if bundle.obj.multiple_batch_id:
 #           bundle.obj.generate_structure_and_dictionary(chirality=bundle.data["stereo_selected"]["name"])
 #       else:
 #           bundle.obj.generate_structure_and_dictionary()
        
    def alter_deserialized_list_data(self, request, deserialized):
        proj = Project.objects.get(project_key=deserialized["project_key"])
        deserialized["project"] = proj
        return deserialized

    def alter_deserialized_detail_data(self, request, deserialized):
        '''A project may be necessary for create statements'''
        if deserialized["project_key"]:
            proj = Project.objects.get(project_key=deserialized["project_key"])
            deserialized["project"] = proj
        return deserialized

    def full_hydrate(self, bundle):
        '''As the object is created we run the validate code on it'''
        bundle = super(CBHCompoundBatchResource, self).full_hydrate(bundle)
        if not bundle.obj.id:
            bundle.obj.created_by=bundle.request.user.username
            # bundle.obj.validate()
            # self.match_list_to_moleculedictionaries(bundle.obj,bundle.data["project"] )
        return bundle



    def prepend_urls(self):
        return [
        url(r"^(?P<resource_name>%s)/delete_index/$" % self._meta.resource_name,
            self.wrap_view('delete_index'), name="delete_index"),
        url(r"^(?P<resource_name>%s)/update_temp_batches/$" % self._meta.resource_name,
            self.wrap_view('update_temp_batches'), name="update_temp_batches"),
        url(r"^(?P<resource_name>%s)/get_part_processed_multiple_batch/$" % self._meta.resource_name,
            self.wrap_view('get_part_processed_multiple_batch'), name="api_get_part_processed_multiple_batch"),
        url(r"^(?P<resource_name>%s)/get_list_elasticsearch/$" % self._meta.resource_name,
            self.wrap_view('get_list_elasticsearch'), name="api_get_list_elasticsearch"),
        url(r"^(?P<resource_name>%s)/get_chembl_ids/$" % self._meta.resource_name,
            self.wrap_view('get_chembl_ids'), name="api_get_chembl_ids"),
        url(r"^(?P<resource_name>%s)/get_elasticsearch_ids/$" % self._meta.resource_name,
            self.wrap_view('get_elasticsearch_ids'), name="api_get_elasticsearch_ids"),
        url(r"^(?P<resource_name>%s)/reindex_elasticsearch/$" % self._meta.resource_name,
            self.wrap_view('reindex_elasticsearch'), name="api_compounds_reindex_elasticsearch"),
        url(r"^(?P<resource_name>%s)/reindex_compound/$" % self._meta.resource_name,
            self.wrap_view('reindex_compound'), name="api_reindex_compound"),
        url(r"^(?P<resource_name>%s)/get_elasticsearch_autocomplete/$" % self._meta.resource_name,
            self.wrap_view('get_elasticsearch_autocomplete'), name="api_get_elasticsearch_autocomplete"),
        url(r"^(?P<resource_name>%s)/validate/$" % self._meta.resource_name,
                self.wrap_view('post_validate'), name="api_validate_compound_batch"),
        url(r"^(?P<resource_name>%s)/validate_list/$" % self._meta.resource_name,
                self.wrap_view('post_validate_list'), name="api_validate_compound_list"),
        url(r"^(?P<resource_name>%s)/multi_batch_save/$" % self._meta.resource_name,
                self.wrap_view('multi_batch_save'), name="multi_batch_save"),
        url(r"^(?P<resource_name>%s)/multi_batch_custom_fields/$" % self._meta.resource_name,
                self.wrap_view('multi_batch_custom_fields'), name="multi_batch_custom_fields"),
        url(r"^(?P<resource_name>%s)/validate_files/$" % self._meta.resource_name,
                self.wrap_view('post_validate_files'), name="api_compound_validate_files"),
        url(r"^(?P<resource_name>%s)/export_file/$" % self._meta.resource_name,
                self.wrap_view('export_file'), name="api_compound_export_file"),]
       



    def patch_dict(self, dictdata, headers):
        for header in headers:
            json_patches = copy.copy(header.get("operations", False))

            if json_patches:
                apply_json_patch(dictdata, json_patches)



    def multi_batch_save(self, request, **kwargs):
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        id = bundle.data["multiplebatch"]
        mb = CBHCompoundMultipleBatch.objects.get(pk=id)
        project = Serializer().serialize(mb)
        limit = 500
        offset = 0
        batches = []
        hasMoreData = True
        while hasMoreData :
            
            bundles = self.get_cached_temporary_batch_data(mb.id,  {"limit":limit, "offset": offset}, request)
            for d in bundles["objects"]:
                self.patch_dict(d, copy.deepcopy(bundle.data["headers"]))
            set_of_batches = self.get_cached_temporary_batches( bundles ,request,  bundledata=bundle.data)
            batches.extend(set_of_batches["objects"])
            offset += limit
            if len(set_of_batches["objects"]) == 0:
                hasMoreData = None

        mb.saved=True
        mb.created_by = request.user.username
        mb.save()

        bundle.data["saved"] = 0
        bundle.data["ignored"] = 0
        to_be_saved = []
        for batch in batches: 
            if(batch.obj.properties.get("action", "") == "New Batch"):   
                batch.obj.id = None    
                batch.obj.generate_structure_and_dictionary()
                batch.multi_batch_id = id
                bundle.data["saved"] += 1
                to_be_saved.append(batch.obj)

        elasticsearch_client.delete_index(elasticsearch_client.get_temp_index_name(request, mb.id))
        batch_dicts = self.batches_to_es_ready(to_be_saved, request)
        index_name=elasticsearch_client.get_main_index_name()
        elasticsearch_client.create_temporary_index(batch_dicts, request, index_name)
        #this needs to be the main elasticsearch compound index
        #and should update any existing records in there? That might be in another method

        return self.create_response(request, bundle, response_class=http.HttpCreated)


    def delete_index(self, request, **kwargs):
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        id = bundle.data["multiplebatch"]
        mb = CBHCompoundMultipleBatch.objects.get(pk=id)
        elasticsearch_client.delete_index(elasticsearch_client.get_temp_index_name(request, mb.id))
        return self.create_response(request, bundle, response_class=http.HttpAccepted)



    def convert_custom_field(self, uncurated_value, field_schema):
        curated_value = uncurated_value
        if field_schema.get("format", "") == "yyyy-mm-dd":
            if uncurated_value:
                
                curated_value = dateutil.parser.parse(uncurated_value).strftime("%Y-%m-%d")

        elif (field_schema.get("field_type", "") == PinnedCustomField.UISELECTTAGS):
            
            curated_value = json.dumps(uncurated_value.split(","))
        elif (field_schema.get("field_type", "") == PinnedCustomField.INTEGER):
            curated_value = int(uncurated_value)     
        elif (field_schema.get("field_type", "") in [ PinnedCustomField.NUMBER, PinnedCustomField.PERCENTAGE]):
            curated_value = float(uncurated_value)    
        return curated_value


    def update_temp_batches(self, request, **kwargs):
        '''change the structure column for an excel file'''
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        multi_batch_id = bundle.data["multiplebatch"]    
        es_serializer = CBHCompoundBatchElasticSearchSerializer()

        es_ready_updates = [es_serializer.to_es_ready_data(dictdata, 
            options={"underscorize": True}) for dictdata in bundle.data["objects"]]
        index_name=elasticsearch_client.get_temp_index_name(request, multi_batch_id)
        elasticsearch_client.create_temporary_index(es_ready_updates, request, index_name)
        elasticsearch_client.get_action_totals(index_name, bundle.data)
        return self.create_response(request, bundle, response_class=http.HttpAccepted)






    def multi_batch_custom_fields(self, request, **kwargs):
        '''change the structure column for an excel file'''
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        id = bundle.data["multiplebatch"]
        headers = bundle.data["headers"]
        # structure_col = bundle.data.get("structure_col", None)

        mb = CBHCompoundMultipleBatch.objects.get(pk=id)
        processSmiles = False
        # if structure_col and structure_col != mb.uploaded_data.get("structure_col", ""):
        #     processSmiles =  True
        index_name=elasticsearch_client.get_temp_index_name(request, mb.id)
        elasticsearch_client.get_action_totals(index_name, bundle.data)
        mb.uploaded_data = bundle.data
        mb.save()

        return self.create_response(request, bundle, response_class=http.HttpAccepted)



    def validate_multi_batch(self,multi_batch, bundle, request, batches):
        batches_not_errors = [batch for batch in batches if batch and not batch.warnings.get("parseError", None) and not batch.warnings.get("smilesParseError", None)]
        
        if (len(batches) ==0):
            raise BadRequest("no_data")
        for b in batches_not_errors:
            b.properties["action"] = "New Batch"

        batches_with_structures = [batch for batch in batches_not_errors if not batch.blinded_batch_id]
        blinded_data =  [batch for batch in batches_not_errors if batch.blinded_batch_id]
        sdfstrings = [batch.ctab for batch in batches_with_structures]
        sdf = "\n".join(sdfstrings)
        
        filename = "/tmp/" + shortuuid.ShortUUID().random()
        text_file = open(filename, "w")
        text_file.write(sdf)
        text_file.close()

        import subprocess
        from subprocess import PIPE, Popen
        p = Popen([settings.INCHI_BINARIES_LOCATION['1.02'], "-STDIO",  filename], stdout=PIPE, stderr=PIPE)
        a = p.communicate()
        inchis = {}
        
        inchiparts = a[0].split("\nStructure:")
        
        for i, inch in enumerate(inchiparts):
            parts = inch.split("\n")
            if len(parts) == 1:
                continue
            ints =[s for s in parts[0].split() if s.isdigit()]
            part = "".join(ints)
            inchis[part] = parts[1]



        if not bundle.data.get("fileerrors"):
            bundle.data["fileerrors"] = []

        new_uploaded_data = []
        already_found = set([])
        duplicates = set([])


        for i, batch in enumerate(batches_with_structures):
            batch.standard_inchi = inchis[str(i+1)]
            batch.validate(temp_props=False)  
            if batch.standard_inchi_key in already_found:
                #setting this in case we change it later
                duplicates.add(batch.standard_inchi_key)
            else:
                already_found.add(batch.standard_inchi_key)
                
            new_uploaded_data.append(batch)

        already_in_db = MoleculeDictionary.objects.filter(project=bundle.data["project"] , structure_type="MOL", structure_key__in=already_found).values_list("structure_key", flat=True)
        already_in_db = set(already_in_db)
        

        bundle.data["new"] = 0
        new_data = set([])
        duplicate_overlaps = set([])
        duplicate_new = set([])

        for batch in batches_with_structures:
            if batch.standard_inchi_key in duplicates:
                batch.warnings["duplicate"] = True
            if batch.standard_inchi_key in already_in_db:
                batch.warnings["overlap"] = True
                if batch.standard_inchi_key in duplicates:
                    batch.warnings["duplicate"] = True
                    duplicate_overlaps.add(batch.standard_inchi_key)
            else:
                batch.warnings["new"] = True
                
                new_data.add(batch.standard_inchi_key)
                if batch.standard_inchi_key in duplicates:
                    batch.warnings["duplicate"] = True
                    duplicate_new.add(batch.standard_inchi_key)
            
        for batch in batches_with_structures:
            if batch.warnings.get("noStructure") == True:
                del batch.warnings["noStructure"]
        for batch in blinded_data:
            batch.warnings["noStructure"] = True


        bundle.data["batchstats"] = {}
        bundle.data["batchstats"]["withstructure"] = len(batches_with_structures)
        bundle.data["batchstats"]["parseErrors"] = len(batches) - len(batches_not_errors) + len([b for b in batches_not_errors if b.warnings.get("parseError", False) == "true"])
        bundle.data["batchstats"]["withoutstructure"] = len(blinded_data)
        bundle.data["batchstats"]["total"] = len(batches)

        bundle.data["compoundstats"] = {}
        bundle.data["compoundstats"]["total"] = len(already_in_db) + len(new_data)
        bundle.data["compoundstats"]["overlaps"] = len(already_in_db)
        bundle.data["compoundstats"]["new"] = len(new_data)
        bundle.data["compoundstats"]["duplicateoverlaps"] = len(duplicate_overlaps)
        bundle.data["compoundstats"]["duplicatenew"] = len(duplicate_new)
        bundle.data["multiplebatch"] = multi_batch.pk


        fifty_batches_for_first_page = self.set_cached_temporary_batches(batches, multi_batch.id, request)

        multi_batch.uploaded_data = bundle.data
        multi_batch.save()
        bundle.data["objects"] = fifty_batches_for_first_page

        index_name=elasticsearch_client.get_temp_index_name(request, multi_batch.id)
        elasticsearch_client.get_action_totals(index_name, bundle.data)

        return self.create_response(request, bundle, response_class=http.HttpAccepted)


    def get_part_processed_multiple_batch(self, request, **kwargs):
        """
        Get the part processed data from elasticsearch and the stats about the
        multiple batch
        """
        # TODO: Uncached for now. Invalidation that works for everyone may be
        #       impossible.
        bundle = self.build_bundle(request=request)
        
        # self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        id = request.GET.get("current_batch")
        mb = CBHCompoundMultipleBatch.objects.get(pk=id)
        
        to_be_serialized = mb.uploaded_data

        to_be_serialized = self.get_cached_temporary_batch_data( id, request.GET, request, bundledata = to_be_serialized)   
        index_name=elasticsearch_client.get_temp_index_name(request, id)
        elasticsearch_client.get_action_totals(index_name, to_be_serialized)     
        # to_be_serialized = self.alter_list_data_to_serialize(request, to_be_serialized)
        return self.create_response(request, to_be_serialized)
 




    def post_validate_list(self, request, **kwargs):
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(self.get_object_list(bundle.request), bundle)

        smilesdata =  bundle.data.get("smilesdata", "")
        objects = smilesdata.splitlines(True)
        batches = []
        #first assume smiles

        allmols = [(obj, Chem.MolFromSmiles(str(obj))) for obj in objects]
        #Next test for inchi
        for m in allmols:
            if m[1] is None:
                inchimol = Chem.MolFromInchi(str(m[0].encode('ascii','ignore')))
                if inchimol is not None:
                    m = (Chem.MolToSmiles( inchimol), inchimol) 


        for m in allmols:
            if(m[1]):
                Compute2DCoords(m[1])
        batches = []
        multiple_batch = CBHCompoundMultipleBatch.objects.create(project=bundle.data["project"])

        for mol2 in allmols:
            if mol2[1]:
                b = CBHCompoundBatch.objects.from_rd_mol(mol2[1], smiles=mol2[0], project=bundle.data["project"], reDraw=True)
            else:
                b = CBHCompoundBatch.objects.blinded(project=bundle.data["project"])
                b.warnings["smilesParseError"] = "true"
                b.properties["action"] = "Ignore"
                b.original_smiles = mol2[0]
            b.multiple_batch_id = multiple_batch.pk
            b.created_by = bundle.request.user.username
            batches.append(b)
  

                
        bundle.data["current_batch"] = multiple_batch.pk
        bundle.data["headers"] = []
        return self.validate_multi_batch(multiple_batch, bundle, request, batches)



    def post_validate_files(self, request, **kwargs):
        
        automapped_structure = False
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        file_name = bundle.data['file_name']
        session_key = request.COOKIES[settings.SESSION_COOKIE_NAME]
        correct_file = FlowFile.objects.get(identifier="%s-%s" % (session_key, file_name))
        batches = []
        headers = []
        errors = []
        fielderrors = {}
        fielderrors["stringdate"] = set([])
        fielderrors["number"] = set([])
        fielderrors["integer"] = set([])
        structure_col = bundle.data.get("struccol","")
        
        if (".cdx" in correct_file.extension ):
            mols = [mol for mol in readfile( str(correct_file.extension[1:]), str(correct_file.file.name), )]
            rxn = None
            index = 0
            if correct_file.extension == '.cdxml':
                #Look for a stoichiometry table in the reaction file
                rxn = chemdraw_reaction.parse( str(correct_file.file.name))
                headers = ["%Completion", 
                            "%Yield", 
                            "Expected Moles", 
                            "Product Moles", 
                            "Expected Mass", 
                            "Product Mass", 
                            "MW", 
                            "role",
                            "Purity", 
                            "Limit Moles", 
                            # "Formula", 
                            "Equivalents", 
                            "Measured Mass"]
            
            for pybelmol in mols:
                molfile = pybelmol.write("mdl")

                if molfile.strip() and molfile != "*":
                    rd_mol = Chem.MolFromMolBlock(molfile, sanitize=False)
                    '''
                    Pybel can read some hypervalent molecules that RDKit cannot read
                    Therefore currently these molecules are outputted as images and sent back to the front end
                    https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg04466.html
                    '''
                    if rd_mol:
                        smiles = Chem.MolToSmiles(rd_mol)

                        if smiles.strip():
                            try:
                                b = CBHCompoundBatch.objects.from_rd_mol(rd_mol, orig_ctab=molfile, smiles=smiles, project=bundle.data["project"])
                            except Exception, e:
                                b =None
                                index = index -1
                            if b:
                                if rxn:
                                    #Here we set the uncurated fields equal to the reaction data extracted from Chemdraw
                                    b.uncurated_fields = rxn.get(pybelmol.title, {})
                                batches.append(b)
                            else:
                                errors.append({"index" : index+1, "image" : pybelmol.write("svg"), "message" : "Unable to produce inchi from this molecule"})

                    else:
                        errors.append({"index" : index+1, "image" : pybelmol.write("svg"), "message" : "Invalid valency or other error parsing this molecule"})
                    index += 1
                    

        else: 
            if (correct_file.extension == ".sdf"):
                #read in the file
                suppl = Chem.ForwardSDMolSupplier(correct_file.file)
                mols = [mo for mo in suppl]
                if(len(mols) > 1000):
                    raise BadRequest("file_too_large")
                #read the headers from the first molecule
                
                headers = get_all_sdf_headers(correct_file.file.name)
                data = correct_file.file.read()
                data = data.replace("\r\n","\n").replace("\r","\n")
                ctabs = data.split("$$$$")
                for index, mol in enumerate(mols):
                    if mol is None: 
                        b=None
                        errors.append({"index" : index+1, "message" : "Invalid valency or other error parsing this molecule"})
                    else:
                        orig_data = ctabs[index]
                        # try:
                        try:
                            b = CBHCompoundBatch.objects.from_rd_mol(mol, orig_ctab=orig_data, project=bundle.data["project"])
                        except Exception, e:
                            errors.append({"index" : index+1,  "message" : str(e)})
                            b =None
                        if b and dict(b.uncurated_fields) == {}:
                            #Only rebuild the uncurated fields if this has not been done before
                            parse_sdf_record(headers, b, "uncurated_fields", mol, fielderrors)
                            
                    batches.append(b)
               

            elif(correct_file.extension in (".xls", ".xlsx")):
                #we need to know which column contains structural info - this needs to be defined on the mapping page and passed here
                #read in the specified structural column
                
                headerswithdata = set([])
                
                df = None
                try:
                    df = pd.read_excel(correct_file.file)
                except IndexError:
                    raise BadRequest("no_headers")

                if len(df.index) > 1000:
                    raise BadRequest("file_too_large")
                #read the smiles string value out of here, when we know which column it is.
                row_iterator = df.iterrows()
                headers = list(df)
                headers = [h.replace(".","__") for h in headers]
                df.columns = headers

                for index, row in row_iterator:
                    #Only automap on the first attempt at mapping the smiles column
                    if not structure_col and not bundle.data.get("headers", None):
                        max_score = 0 

                        for header in headers:
                            #fuzzy matching for smiles - this should also match things like "canonical_smiles"
                            hdr = re.sub('[^0-9a-zA-Z]+', ' ', header)
                            for h in hdr.split(" "):
                                h = h.strip()
                                if h:
                                    score = fuzzymatch(a="smiles", b=h.lower()).ratio()
                                    if score > max_score and score > 0.9:
                                        structure_col = header
                                        max_score = score
                                        automapped_structure = True

                    if structure_col:
                        smiles_str = row[structure_col]
                        try:
                            struc = Chem.MolFromSmiles(smiles_str)
                            if struc:
                                Compute2DCoords(struc)
                                b = CBHCompoundBatch.objects.from_rd_mol(struc, smiles=smiles_str, project=bundle.data["project"], reDraw=True)
                                b.blinded_batch_id = None
                            else:
                                raise Exception("Smiles not processed")
                                    
                        except Exception, e:
                            b = CBHCompoundBatch.objects.blinded(project=bundle.data["project"])
                            b.original_smiles = smiles_str
                            b.warnings["smilesParseError"] = "true"
                            b.properties["action"] = "Ignore"
                    else:
                        b = CBHCompoundBatch.objects.blinded(project=bundle.data["project"])                    

                    if b: 
                        if dict(b.uncurated_fields) == {}:
                        #Only rebuild the uncurated fields if this has not been done before
                            parse_pandas_record(headers, b, "uncurated_fields", row, fielderrors, headerswithdata)
                    else:
                        errors.append({"index" : index+1, "message" : "Invalid valency or other error parsing this identifier",  "SMILES": smiles_str})
                    batches.append(b)
                headers = [hdr for hdr in headers if hdr in headerswithdata]
            else:
                raise BadRequest("file_format_error")

        multiple_batch = CBHCompoundMultipleBatch.objects.create(project = bundle.data["project"])
        for b in batches:
            if b:
                b.multiple_batch_id = multiple_batch.pk
                b.created_by = bundle.request.user.username
           
        bundle.data["fileerrors"] = errors
        bundle.data["automapped"] = 0
        cfr = ProjectResource()
        schemaform = cfr.get_schema_form(bundle.data["project"].custom_field_config,"" )
        if not bundle.data.get("headers", None):
            bundle.data["headers"] = []
            for header in headers:
                copyto = ""
                automapped = False
                operations = []
                if header == structure_col:  
                    copyto = "SMILES for chemical structures" 
                    if automapped_structure:
                        automapped = True
                else:
                    form = copy.deepcopy(schemaform["form"])
                    copyto = ""
                    max_score = 0

                    for form_item in form:
                        score = fuzzymatch(a=form_item["key"].lower(), b=header.lower()).ratio()
                        if score > max_score and score > 0.9:
                            matched_item = form_item
                            copyto = matched_item["key"]
                            automapped = True 
                            if(matched_item["field_type"]=="uiselecttags"):
                                operations.append({"op": "split", "path": "/uncurated_fields/" + header})
                                operations.append({"op": "move", "path": "/custom_fields/" + matched_item["key"] , "from" : "/uncurated_fields/" + header })
                            else:
                                operation = {"op": "move", "path": "/custom_fields/" + matched_item["key"]  , "from" : "/uncurated_fields/" + header }
                                operations.append(operation);
                                if(matched_item.get("format", "")=="date"):
                                    operations.append({"op": "convertdate", "path": "/custom_fields/" + matched_item["key"]})

                bundle.data["headers"].append({
                                                    "name": header,
                                                    "automapped": automapped, 
                                                    "copyto": copyto,
                                                    "operations" : operations,
                                                    "fieldErrors" : { 
                                                        "stringdate": header in fielderrors["stringdate"],
                                                        "integer": header in fielderrors["integer"],
                                                        "number": header in fielderrors["number"]
                                                    }
                                            })

        return self.validate_multi_batch(multiple_batch, bundle, request, batches)


    def alter_list_data_to_serialize(self, request, data):
        '''use the request type to determine which fields should be limited for file download,
           add extra fields if needed (eg images) and enumerate the custom fields into the 
           rest of the calculated fields'''
        if self.substructure_smarts:
            for index, b in enumerate(data["objects"]):
                ctab = b.data["properties"]["substructureMatch"] = self.substructure_smarts

        if(self.determine_format(request) == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' or request.GET.get("format") == "sdf" or self.determine_format(request) == 'chemical/x-mdl-sdfile' ):
            
            ordered_cust_fields = PinnedCustomField.objects.filter(custom_field_config__project__project_key__in = request.GET.get("project__project_key__in","").split(",")).order_by("custom_field_config__project__id","position").values("name", "field_type")
            seen = set()
            deduplicated_cfs = [x for x in ordered_cust_fields if x["name"] not in seen and not seen.add(x["name"])]
            df_data = []
            ordered_cust_fields = []
            projects = set([]) 
            uncurated_field_names = set()
            for index, b in enumerate(data["objects"]):
                #remove items which are not listed as being kept
                new_data = {}
                projects.add(b.obj.project_id)
                for k, v in b.data.iteritems():
                    for name, display_name in self.Meta.fields_to_keep.iteritems():
                        if k == name:
                            new_data[display_name] = v
                #we need sd format exported results to retain stereochemistry - use mol instaed of smiles
                if(self.determine_format(request) == 'chemical/x-mdl-sdfile' or request.GET.get("format") == "sdf"):
                    new_data['ctab'] = b.data['ctab']
                #dummy
                #not every row has a value for every custom field
                for item in deduplicated_cfs:
                    cf_value = b.data["custom_fields"].get(item["name"], "")
                    if item["field_type"] == PinnedCustomField.UISELECTTAGS:
                        if isinstance(cf_value, basestring):
                            try:
                                cf_value = json.loads(cf_value).join(",")
                            except:
                                pass
                        elif isinstance(cf_value, list):
                            cf_value = cf_value.join(",")
                    new_data[item["name"]] = cf_value 
                    
                #now remove custom_fields
                del(new_data['custom_fields'])
                for field, value in b.data['uncurated_fields'].iteritems():
                    new_data[field] = value
                    uncurated_field_names.add(field)
                    
                # #now remove custom_fields
                # del(new_data['uncurated_fields'])
                b.data = new_data
                df_data.append(new_data)               
            df = pd.DataFrame(df_data)
            data['export'] = df.to_json() 
            data["headers"] = {
                "cbh" : [display_name for name, display_name in self.Meta.fields_to_keep.iteritems()],
                "custom_fields" : [cf["name"] for cf in deduplicated_cfs] ,
                "uncurated_fields" : sorted(list(uncurated_field_names))
            }
            

        return data

    def create_response(self, request, data, response_class=HttpResponse, **response_kwargs):
        """
        Extracts the common "which-format/serialize/return-response" cycle.
        Mostly a useful shortcut/hook.
        """

        desired_format = self.determine_format(request)
        serialized = self.serialize(request, data, desired_format, options=self.Meta.ordrered_ftk)
        rc = response_class(content=serialized, content_type=build_content_type(desired_format), **response_kwargs)

        if(desired_format == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'):
            rc['Content-Disposition'] = 'attachment; filename=export.xlsx'
        elif(desired_format == 'chemical/x-mdl-sdfile'):
            rc['Content-Disposition'] = 'attachment; filename=export.sdf'
        return rc


    def dehydrate(self, bundle):
        
        #try:
        data = bundle.obj.related_molregno
        user = None
        
        for names in self.Meta.fieldnames:
            if not bundle.obj.blinded_batch_id:
                try:
                    bundle.data[names[1]] = deepgetattr(data, names[0], None)
                except(AttributeError):   
                    bundle.data[names[1]] = ""
            else:
                if names[1] == "chemblId":
                    bundle.data[names[1]] = bundle.obj.blinded_batch_id
                else:
                    bundle.data[names[1]] = ""
        if bundle.obj.created_by:
          #user = User.objects.get(username=bundle.obj.created_by)
            User = get_user_model()
            try:
                user = User.objects.get(username=bundle.obj.created_by)
            except ObjectDoesNotExist:
                pass


        mynames = ["editable_by","uncurated_fields", "warnings", "properties", "custom_fields", "errors"]
        for name in mynames:
            bundle.data[name] = json.loads(bundle.data[name])
        #bundle.data["created_by"] = user.__dict__ 
        if user != None:
            if user.first_name:
                bundle.data["created_by"] = "%s %s" % (user.first_name, user.last_name)
            else:
                bundle.data["created_by"] = user.username
        else:
            bundle.data["created_by"] = ""
        bundle.data["timestamp"] = str(bundle.data["created"])[0:10]
        #except:
        #    pass
    
        return bundle

    def batches_to_es_ready(self, batches, request, non_chem_data_only=None):
        batch_dicts = []
        index = 1
        es_serializer = CBHCompoundBatchElasticSearchSerializer()
        for batch in batches:
            if batch:
                if(not batch.id):
                    batch.id = index
                bun = self.build_bundle(obj=batch, request=request)

                bun = self.full_dehydrate(bun, for_list=True)
                if bun.data["project"] == '':
                    bun.data["project"] = '/%s/cbh_projects/%d' % (settings.WEBSERVICES_NAME, bun.obj.project_id)
                if non_chem_data_only:
                    ready = es_serializer.to_es_ready_non_chemical_data(bun.data, options={"underscorize": True})
                else:
                    ready = es_serializer.to_es_ready_data(bun.data , options={"underscorize": True})
                batch_dicts.append(ready)
            else:
                #preserve the line number of the batch that could not be processed
                batch_dicts.append(
                    {"id": index, 
                    "warnings" : {"parseError": "true"
                                    },
                    "properties" : {"action": "Ignore"},
                    "project": str(self.project)
                    }
                    )
            index += 1
        return batch_dicts

    def set_cached_temporary_batches(self, batches, multi_batch_id, request):
        es_serializer = CBHCompoundBatchElasticSearchSerializer()
        batch_dicts = self.batches_to_es_ready(batches, request)
        index_name=elasticsearch_client.get_temp_index_name(request, multi_batch_id)
        elasticsearch_client.create_temporary_index(batch_dicts, request, index_name)
        #Now get rid of my ES preparation again
        return [es_serializer.to_python_ready_data(d) for d in batch_dicts[0:10]]


    def get_cached_temporary_batches(self, bundles, request, bundledata={}):
        fakeobj = CBHCompoundBatch(project = bundledata["project"], id=10)
        bun = self.build_bundle(obj=fakeobj)
        bun = self.full_dehydrate(bun)
        proj = bun.data["project"]
        data = []
        for datum in bundles["objects"]:
            datum = self.build_bundle(data=datum, request=request)
            datum.data["project"] = proj
            datum = self.full_hydrate(datum)
            data.append(datum)
        bundles["objects"] = data
        return bundles




    def get_list_elasticsearch(self, request, **kwargs):
        """
        Returns a serialized list of resources.
        Calls ``obj_get_list`` to provide the data, then handles that result
        set and serializes it.
        Should return a HttpResponse (200 OK).
        """
        # TODO: Uncached for now. Invalidation that works for everyone may be
        #       impossible.
        base_bundle = self.build_bundle(request=request)
        must_list = []
            #If there is a substructure query then we use the non-elasticsearch implementation to pull back the required fields
        objects = self.obj_get_list(bundle=base_bundle, **self.remove_api_resource_names(kwargs))
        ids = objects.values_list("id", flat=True)[0:100000]
        
        get_data = request.GET
        blanks_filter = json.loads(get_data.get("showBlanks", "[]"))
        blanks_filter = [get_key_from_field_name(field) for field in blanks_filter]
        nonblanks_filter = json.loads(get_data.get("showNonBlanks", "[]"))
        nonblanks_filter = [get_key_from_field_name(field) for field in nonblanks_filter]
        archived = request.GET.get("archived", None)
        #we need to alter the query to look for blanks or limit to non-blanks if specified
        query = {"ids" : {
                                "type" : "batches",
                                "values" : [str(i) for i in ids]
                            }

                }

        #modify the query to include
        blanks_queries = []
        if blanks_filter:
            #blanks_queries = []
            for blank in blanks_filter:
                
                blanks_queries.append({"bool" :
                                            {"should" :[
                                              {"term": {blank + ".raw": ""} },
                                               {"missing": {"field": blank}}
                                               ]
                                             }
                                            })

        if archived == "true" :
            blanks_queries.append({"term": {"properties.archived" : "true"}})
        else:
            blanks_queries.append({"bool" : 
                                            {"should" :[{"term": {"properties.archived" : "false"}},
                                            {"missing": {"field": "properties.archived"}}]}
                                })
        nonblanks_queries = []
        if nonblanks_filter:
            
            for nonblank in nonblanks_filter:

                nonblanks_queries.append({"bool" :
                                            {"should" :[
                                              {"term": {nonblank + ".raw": ""} },
                                               {"missing": {"field": nonblank}}
                                               ]
                                             }
                                            })



        modified_query = {
                        "bool": {
                            "must": [query] + blanks_queries,
                            "must_not": nonblanks_queries,
                        },
            }

        
        es_request = {
            "version" : True,
            "from" : get_data.get("offset", 0),
            "size" : get_data.get("limit", 50),
            "filter" : modified_query,
            "sort" : json.loads(get_data.get("sorts",'[{"id": {"order": "desc"}}]'))
        }
        index = elasticsearch_client.get_main_index_name()
        es_serializer = CBHCompoundBatchElasticSearchSerializer()
        es_serializer.convert_query(es_request)
        bundledata = elasticsearch_client.get(index, es_request, {})
        bundledata["objects"] = [es_serializer.to_python_ready_data(d) for d in bundledata["objects"]]

        # paginator = self._meta.paginator_class(request.GET, sorted_objects, resource_uri=self.get_resource_uri(), limit=self._meta.limit, max_limit=self._meta.max_limit, collection_name=self._meta.collection_name)
        # to_be_serialized = paginator.page()

        # # Dehydrate the bundles in preparation for serialization.
        # bundles = []

        # for obj in to_be_serialized[self._meta.collection_name]:
        #     bundle = self.build_bundle(obj=obj, request=request)
        #     bundles.append(self.full_dehydrate(bundle, for_list=True))

        return self.create_response(request, bundledata)






    def get_cached_temporary_batch_data(self, multi_batch_id, get_data,request, bundledata={}):
        es_request = {
            "from" : get_data.get("offset", 0),
            "size" : get_data.get("limit", 50),
            "filter" : json.loads(get_data.get("query", '{ "match_all" : {}}')),
            "sort" : json.loads(get_data.get("sorts",'[{"id": {"order": "asc"}}]'))
        }
        index = elasticsearch_client.get_temp_index_name(request, multi_batch_id)
        es_serializer = CBHCompoundBatchElasticSearchSerializer()
        es_serializer.convert_query(es_request)
        bundledata = elasticsearch_client.get(index, es_request, bundledata)
        bundledata["objects"] = [es_serializer.to_python_ready_data(d) for d in bundledata["objects"]]
        return bundledata



    




    def get_object_list(self, request):
        return super(CBHCompoundBatchResource, self).get_object_list(request).prefetch_related(Prefetch( "related_molregno__compoundproperties")).prefetch_related(Prefetch( "project"))





def deepgetattr(obj, attr, ex):
    """Recurses through an attribute chain to get the ultimate value."""
    try:
        return reduce(getattr, attr.split('.'), obj)

    except:
        return ex






class CBHCompoundMultipleBatchResource(ModelResource):
    #comp_batch = fields.ForeignKey(CBHCompoundBatchResource, 'cbh_compound_batches', blank=False, null=False)
    #batches = fields.ToManyField(CBHCompoundBatchResource, 'batches', full=True)
    project = fields.ForeignKey(ProjectResource, 'project', blank=False, null=False, full=True)
    class Meta:
        filtering = {
            "created_by": ALL_WITH_RELATIONS,
            "project": ALL_WITH_RELATIONS,
        }
        always_return_data = True
        queryset = CBHCompoundMultipleBatch.objects.all()
        resource_name = 'cbh_multiple_batches'
        authorization = ProjectAuthorization()
        include_resource_uri = False
        allowed_methods = ['get']
        default_format = 'application/json'
        authentication = SessionAuthentication()  

    def apply_filters(self, request, applicable_filters):
        pids = self._meta.authorization.project_ids(request)
        dataset = self.get_object_list(request).filter(**applicable_filters).filter(project_id__in=set(pids))
        return dataset.order_by("-created")

    def get_object_list(self, request):
        return super(CBHCompoundMultipleBatchResource, self).get_object_list(request).defer('uploaded_data').prefetch_related(Prefetch( "project"))



class CBHCompoundBatchUpload(ModelResource):

    class Meta:
        excludes = ['uploaded_data']

        always_return_data = True
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
        session_key = request.COOKIES[settings.SESSION_COOKIE_NAME]
        correct_file = self.get_object_list(request).get(identifier="%s-%s" % (session_key, file_name))
        
        header_json = { }

        if (correct_file.extension == ".sdf"):
            #read in the file
            headers = get_all_sdf_headers(correct_file.file.name)


        elif(correct_file.extension in (".xls", ".xlsx")):
            #read in excel file, use pandas to read the headers
            df = pd.read_excel(correct_file.file)
            headers = list(df)

        #this converts to json in preparation to be added to the response
        bundle.data["headers"] = list(set(headers))
            #send back
            #we should allow SD file uploads with no meta data
        if (len(headers) == 0 and correct_file.extension in (".xls", ".xlsx") ):
            raise BadRequest("no_headers")
        return self.create_response(request, bundle, response_class=http.HttpAccepted)


def get_all_sdf_headers(filename):

    from subprocess import Popen, PIPE
    from shlex import split

    p1 = Popen(split('/bin/grep "^>" %s' % filename), stdout=PIPE)
    p2 = Popen(split('/usr/bin/cut -d "<" -f2'), stdin=p1.stdout, stdout=PIPE)
    p3 = Popen(split('/usr/bin/cut -d ">" -f1'), stdin=p2.stdout, stdout=PIPE)
    p4 = Popen(split('/usr/bin/sort'), stdin=p3.stdout, stdout=PIPE)
    p5 = Popen(split('/usr/bin/uniq'), stdin=p4.stdout, stdout=PIPE)
    out = p5.communicate()
    return [i for i in out[0].split("\n") if i]
