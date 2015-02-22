from tastypie.resources import ModelResource, ALL, ALL_WITH_RELATIONS
from django.conf import settings
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
from itertools import chain
from pybel import readfile , readstring

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

from cbh_chembl_ws_extension.authorization import ProjectAuthorization
from cbh_chembl_ws_extension.projects import ProjectResource
from chembl_webservices.compounds import CompoundsResource
from chembl_webservices.base import ChEMBLApiSerializer
from cbh_chembl_ws_extension.serializers import CBHCompoundBatchSerializer

from tastypie.utils import dict_strip_unicode_keys
from tastypie.serializers import Serializer
from django.core.serializers.json import DjangoJSONEncoder
from tastypie import fields, utils
from cbh_chembl_model_extension.models import CBHCompoundBatch, CBHCompoundMultipleBatch, Project, PinnedCustomField
from tastypie.authentication import SessionAuthentication
import json
from tastypie.paginator import Paginator
from chembl_beaker.beaker.core_apps.conversions.impl import _smiles2ctab, _apply

from flowjs.models import FlowFile
import xlrd
import pandas as pd
import numpy as np
import urllib
from tastypie.validation import Validation

from django.db.models import Max

from tastypie.serializers import Serializer

from tastypie.authentication import SessionAuthentication



# class MoleculeValidation(Validation):
#     def is_valid(self, bundle, request=None):
#         if not bundle.data:
#             return {'No Data': 'invalid_call_to_api, no data'}

#         errors = {}

#         structure_type = bundle.data["structure_type"]
#         structure_key = bundle.data["structure_key"]
#         project = bundle.data["project"]
#         chirality = bundle.data["chirality"]
#         forced_reg_index = bundle.data.get("forced_reg_index",0)
#         forced_reg_reason = bundle.data.get("forced_reg_reason","")

#         public = bundle.data.get("public", False)
#         molregno = bundle.data.get("molregno", None)
#         #Test uniqueness

#         if not molregno:
#             public_or_project_moldicts = MoleculeDictionary.objects.filter(
#                                             (Q(structure_type=structure_type)
#                                            Q(structure_key=structure_key),
#                                            Q(chirality=chirality), 
#                                            Q(public=True)) |
#                                             (Q(structure_type=structure_type)
#                                            Q(structure_key=structure_key),
#                                            Q(chirality=chirality), 
#                                            Q(project_id=project.id)
#                                                 ))

            
#             if public_moldicts.count() > 0:
#                 max_reg = moldicts.aggregate(Max('forced_reg_index'))["forced_reg_index__max"]
#                 if forced_reg_index > 0:
#                     if forced_reg_index <= max_reg:
#                         #Somohow the forced registration index is too low, the client must re submit
#                         errors["incorrect_force_reg_index"] ="Forced Registration unsuccesful, try with higher index"
#                         errors["next_forced_reg_index"] = max_reg + 

#                     if len(forced_reg_reason) < 3:
#                         errors["invalid_force_reg_reason"] ="Forced registration reason should be longer than 3 characters"
#                 else:
#                     errors["already_registered_as_public"] ="Already registered for this project or publically"
#                     errors["next_forced_reg_index"] = max_reg + 1
#         return errors

            
        

            



            


#         return errors


# class MoleculeDictionaryResource(ModelResource):
#     project = fields.ForeignKey(ProjectResource, 'project', blank=True, null=True)
#     class Meta:    
#         queryset = MoleculeDictionary.objects.all()
#         resource_name = 'molecule_dictionaries'
#         authorization = ProjectAuthorization()
#         include_resource_uri = False
#         allowed_methods = ['get', 'post', 'put']
#         default_format = 'application/json'
#         authentication = SessionAuthentication()
#         paginator_class = Paginator






class CBHCompoundBatchResource(ModelResource):
    project = fields.ForeignKey(ProjectResource, 'project', blank=False, null=False)
    class Meta:
        filtering = {
            "std_ctab": ALL_WITH_RELATIONS,
            "ctab": ALL,
            "multiple_batch_id": ALL_WITH_RELATIONS,
            "project": ALL_WITH_RELATIONS,
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
        authorization = ProjectAuthorization()
        include_resource_uri = False
        serializer = CBHCompoundBatchSerializer()
        allowed_methods = ['get', 'post', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        paginator_class = Paginator


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
        print batch.warnings

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
        batch.warnings["new"] = new
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

        updated_bundle = self.obj_build(bundle, dict_strip_unicode_keys(deserialized))
        bundle.obj.validate()
        self.match_list_to_moleculedictionaries(bundle.obj,bundle.data["project"] )
        dictdata = bundle.obj.__dict__
        dictdata.pop("_state")
        
        updated_bundle = self.build_bundle(obj=bundle.obj, data=dictdata)
        return self.create_response(request, updated_bundle, response_class=http.HttpAccepted)

    def get_project_custom_field_names(self, request, **kwargs):
        '''Combine the pinned fields for the project with the most frequently used fields on the project into a single list'''
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'text/plain'))  
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)

        pinned_fields = list(PinnedCustomField.objects.filter(custom_field_config__project=bundle.data["project"]).values())
        print pinned_fields
        secondwhere = True
        project_id = bundle.data["project"].id
        if len(pinned_fields) > 0:
            secondwhere = "key not in %s" % json.dumps([pf["name"] for pf in pinned_fields]).replace("\"","'").replace("[","(").replace("]",")")
        fields = CBHCompoundBatch.objects.get_all_keys(where="project_id = %d" % project_id, 
                                                        secondwhere=secondwhere)
        project_fields = [{'name': item[0], 'count': item[1]} for item in fields]
        pinned_fields.extend(project_fields)
        bundle.data['field_names'] = pinned_fields
        return self.create_response(request, bundle, response_class=http.HttpAccepted)


    def save_related(self, bundle):
        #bundle.obj.created_by = request.user.

        bundle.obj.generate_structure_and_dictionary()

 #   def save_related(self, bundle):
 #       #For a single molecule, set the chirality
 #       if bundle.obj.multiple_batch_id:
 #           bundle.obj.generate_structure_and_dictionary(chirality=bundle.data["stereo_selected"]["name"])
 #       else:
 #           bundle.obj.generate_structure_and_dictionary()
        

    def alter_deserialized_detail_data(self, request, deserialized):
        proj = Project.objects.get(project_key=deserialized["project_key"])
        deserialized["project"] = proj
        return deserialized

    def full_hydrate(self, bundle):
        '''As the object is created we run the validate code on it'''
        bundle = super(CBHCompoundBatchResource, self).full_hydrate(bundle)
        bundle.obj.validate()
        self.match_list_to_moleculedictionaries(bundle.obj,bundle.data["project"] )
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
        if bundle.obj.pk:
            self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        id = bundle.data["current_batch"]
        batches = CBHCompoundMultipleBatch.objects.get(pk=id).uploaded_data
        bundle.data["saved"] = 0
        bundle.data["errors"] = []

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
        if bundle.obj.pk:
            self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        id = bundle.data["current_batch"]
        mappings = bundle.data.get('mapping', [])

        mb = CBHCompoundMultipleBatch.objects.get(pk=id)
        headers = None
        for b in mb.uploaded_data:
            if mappings:
                if not headers:
                    headers = b.custom_fields.keys()
                custom_fields = {}
                for hdr in headers:
                    if hdr in mappings["ignored_fields"]:
                        continue
                    elif hdr in mappings["new_fields"]:
                        custom_fields[hdr] = b.custom_fields[hdr]
                    else:
                        for key, mapping in mappings["remapped_fields"].iteritems():
                            if hdr in mapping:
                                custom_fields[key] = b.custom_fields[hdr]
                b.custom_fields = custom_fields
            else:
                b.custom_fields = bundle.data["custom_fields"]
        mb.save()
        #Might not be needed
        return self.validate_multi_batch(mb, bundle, request)



    def validate_multi_batch(self,multi_batch, bundle, request):
        total = len(multi_batch.uploaded_data)
        bundle.data["objects"] = []
        bundle.data["new"] = 0
        bundle.data["linkedpublic"] = 0
        bundle.data["linkedproject"] = 0
        bundle.data["errors"] = 0
        
        for batch in multi_batch.uploaded_data:
            self.match_list_to_moleculedictionaries(batch,bundle.data["project"] )
            for key in ["new", "linkedproject", "linkedpublic"]:
                bundle.data[key] += int(batch.warnings[key])
            b = batch.__dict__
            bundle.data["objects"].append(b)
            #catch all molecules that should have some choices associated with them
            if b["errors"] != {}:
                bundle.data["errors"] += 1
                total = total - 1


                
            #elif b["warnings"]["pains_count"] != "0":
                
        multi_batch.save()

        bundle.data["total"] = total

        bundle.data["current_batch"] = multi_batch.pk

        return self.create_response(request, bundle, response_class=http.HttpAccepted)



    def post_validate_list(self, request, **kwargs):
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        type = bundle.data.get("type",None).lower()
        objects =  bundle.data.get("objects", [])
       
        batches = []
        if type == "smiles":
            batches = [CBHCompoundBatch.objects.from_rd_mol(Chem.MolFromSmiles(obj), smiles=obj, project=bundle.data["project"]) for obj in objects ]

        elif type == "inchi":
            mols = _apply(objects,Chem.MolFromInchi)
            batches = [CBHCompoundBatch.objects.from_rd_mol(mol, smiles=Chem.MolToSmiles(mol), project=bundle.data["project"]) for mol in mols]
        # for b in batches:
        #     b["created_by"] = request.user.username

        multiple_batch = CBHCompoundMultipleBatch.objects.create()
        for b in batches:
            b.multiple_batch_id = multiple_batch.pk

        multiple_batch.uploaded_data=batches
        multiple_batch.save()
        bundle.data["current_batch"] = multiple_batch.pk
        return self.validate_multi_batch(multiple_batch, bundle, request)



    def post_validate_files(self, request, **kwargs):
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        file_name = bundle.data['file_name']
        session_key = request.COOKIES[settings.SESSION_COOKIE_NAME]
        correct_file = FlowFile.objects.get(identifier="%s-%s" % (session_key, file_name))

        batches = []
        headers = []
        if (".cdx" in correct_file.extension ):
            mols = [mol.write("smi").split("\t")[0] for mol in readfile( str(correct_file.extension[1:]), str(correct_file.file.name), )]
            for smiles in mols:

                if smiles.strip() and smiles != "*":
                        b = CBHCompoundBatch.objects.from_rd_mol(Chem.MolFromSmiles(smiles), smiles=smiles, project=bundle.data["project"])
                        batches.append(b)
        else: 
            if (correct_file.extension == ".sdf"):
                #read in the file
                suppl = Chem.ForwardSDMolSupplier(correct_file.file)

                #read the headers from the first molecule
                for mol in suppl:
                    if mol is None: continue
                    if not headers: 
                        headers = list(mol.GetPropNames())


                    b = CBHCompoundBatch.objects.from_rd_mol(mol, smiles=Chem.MolToSmiles(mol), project=bundle.data["project"])
                    custom_fields = {}
                    for hdr in headers:

                        custom_fields[hdr] = mol.GetProp(hdr) 
                        
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
                    b = CBHCompoundBatch.objects.from_rd_mol(Chem.MolFromSmiles(smiles_str), smiles=smiles_str, project=bundle.data["project"])
                    #work out custom fields from mapping object
                    #new_fields, remapped_fields, ignored_fields
                    
                    custom_fields = {}
                    for hdr in headers:
                        custom_fields[ hdr] = row[hdr] 
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
        if len(headers) == 0:
            return BadRequest("no_headers")
        return self.create_response(request, bundle, response_class=http.HttpAccepted)

    

