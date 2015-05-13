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
from tastypie.utils.mime import build_content_type
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
from cbh_chembl_ws_extension.serializers import CBHCompoundBatchSerializer
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

from django.db.models import Max

from tastypie.serializers import Serializer

from tastypie.authentication import SessionAuthentication

import  chemdraw_reaction

from django.contrib.auth import get_user_model

from rdkit.Chem.AllChem import Compute2DCoords

from django.db.models import Prefetch


            



class CBHCompoundBatchResource(ModelResource):
    project = fields.ForeignKey(ProjectResource, 'project', blank=False, null=False)

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

    def apply_filters(self, request, applicable_filters):
        """
        An ORM-specific implementation of ``apply_filters``.
        The default simply applies the ``applicable_filters`` as ``**kwargs``,
        but should make it possible to do more advanced things.
        """
        ws = request.GET.get("with_substructure", None)
        st = request.GET.get("similar_to", None)
        fm = request.GET.get("flexmatch", None)
        #this is the similarity index for fingerprint-like searching
        fp = request.GET.get("fpValue", None)
        cms = None

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
        #else:
        #    cms = CompoundMols.objects.all()
        #To be generalised
        cust = request.GET.get("search_custom_fields__kv_any", None)
        if cust:
            applicable_filters["custom_fields__kv_any"] = cust
        if cms != None:
            #run the sql for pulling in new compounds into compound_mols
            indexed = CBHCompoundBatch.objects.index_new_compounds()
            applicable_filters["related_molregno_id__in"] = cms.values_list("molecule_id", flat=True)

        pids = self._meta.authorization.project_ids(request)
        return self.get_object_list(request).filter(**applicable_filters).filter(project_id__in=set(pids)).order_by("-created")
    
    

   


    def convert_mol_string(self, strn):
        #commit
        try:
            mol = Chem.MolFromMolBlock(strn)
            smiles = Chem.MolToSmiles(mol)
        except:
            smiles = strn
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
        bundle.obj.std_ctab = bundle.obj.ctab
        bundle.obj.validate()
        self.match_list_to_moleculedictionaries(bundle.obj,bundle.data["project"] )
        dictdata = bundle.obj.__dict__
        dictdata.pop("_state")
        
        updated_bundle = self.build_bundle(obj=bundle.obj, data=dictdata)
        return self.create_response(request, updated_bundle, response_class=http.HttpAccepted)


    def get_project_custom_field_names(self, request, **kwargs):
        '''Get a single list of pinned fields for the project previously listed all custom fields in the DB but this was unwieldy'''
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'text/plain'))  
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)

        pinned_fields = list(PinnedCustomField.objects.filter(custom_field_config__project=bundle.data["project"]).values())
        bundle.data['field_names'] = pinned_fields
        return self.create_response(request, bundle, response_class=http.HttpAccepted)


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
                self.wrap_view('export_file'), name="api_compound_export_file"),]
       

    def multi_batch_save(self, request, **kwargs):
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        id = bundle.data["current_batch"]
        mb = CBHCompoundMultipleBatch.objects.get(pk=id)
        batches = mb.uploaded_data
        mb.saved=True
        mb.created_by = request.user.username
        mb.save()

        bundle.data["saved"] = 0
        bundle.data["errors"] = []


        for batch in batches:

                
                batch.generate_structure_and_dictionary()
                #batch.save(validate=False)
                bundle.data["saved"] += 1
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
                #If there is a mappings object this data has come from the SD upload
                if not headers:
                    headers = b.uncurated_fields.keys()
                custom_fields = {}
                uncurated_fields = {}
                for hdr in headers:
                    if hdr in mappings["ignored_fields"]:
                        #ignored fields are removed
                        b.uncurated_fields.pop(hdr, False)
                    elif hdr in mappings["new_fields"]:
                        #New fields are left as uncurated for now
                        continue
                    else:
                        for key, mapping in mappings["remapped_fields"].iteritems():
                            #Remapped fields are added to the curated fields as only curated data will be present
                            if hdr in mapping:
                                custom_fields[key] = b.uncurated_fields.pop(hdr, "")
                b.custom_fields = custom_fields
            else:
                #Only curated fields can be added via the standard UI
                b.custom_fields = bundle.data["custom_fields"]

        mb.save()
        return self.create_response(request, bundle, response_class=http.HttpAccepted)



    def validate_multi_batch(self,multi_batch, bundle, request):
        t = time.time()
        sdfstrings = [batch.ctab for batch in multi_batch.uploaded_data]
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
            

        

        total = len(multi_batch.uploaded_data)
        bundle.data["objects"] = []
        bundle.data["new"] = 0
        bundle.data["linkedpublic"] = 0
        bundle.data["linkedproject"] = 0
        if not bundle.data.get("fileerrors"):
            bundle.data["fileerrors"] = []
        bundle.data["errors"] = len(bundle.data.get("fileerrors", []))
        bundle.data["dupes"] = 0
        new_uploaded_data = []
        already_found = set([])
        for i, batch in enumerate(multi_batch.uploaded_data):
            batch.standard_inchi = inchis[str(i+1)]
            batch.validate(temp_props=False)
            if batch.__dict__["errors"] != {}:
                bundle.data["errors"] += 1
                bundle.data["fileerrors"].append({"index" : i+1, "error": "Inchi Parse Error"})
            else:
                batch_key = batch.get_uk()
                if batch_key in already_found:
                    #setting this in case we change it later
                    batch.properties["dupe"] = True
                    bundle.data["dupes"] += 1
                else:
                    already_found.add(batch_key)
                    # self.match_list_to_moleculedictionaries(batch,bundle.data["project"] )
                    # for key in ["new", "linkedproject", "linkedpublic"]:
                    #     bundle.data[key] += int(batch.warnings[key])
                new_uploaded_data.append(batch)

        total_already_in_db = MoleculeDictionary.objects.filter(structure_type="MOL", structure_key__in=already_found).count()
        bundle.data["linkedproject"] = total_already_in_db
        bundle.data["new"] = len(multi_batch.uploaded_data) - total_already_in_db - bundle.data["dupes"]
        multi_batch.uploaded_data = new_uploaded_data
        multi_batch.save()
        print t -time.time()
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
            allmols = [(obj, Chem.MolFromSmiles(str(obj))) for obj in objects]

        elif type == "inchi":
            allmols = [(Chem.MolToSmiles( Chem.MolFromInchi(str(obj.encode('ascii','ignore')))), Chem.MolFromInchi(obj.encode('ascii','ignore'))) for obj in objects]

        errors = [{"index" : i+1, "error": "Parse Error"} for i, mol in enumerate(allmols) if not mol[1]]
        bundle.data["fileerrors"] = errors

        mols = [ mol1 for mol1 in allmols if mol1[1]]
        for m in mols:
            Compute2DCoords(m[1])
        batches = [CBHCompoundBatch.objects.from_rd_mol(mol2[1], smiles=mol2[0], project=bundle.data["project"], reDraw=True) for mol2 in mols]
    

        multiple_batch = CBHCompoundMultipleBatch.objects.create(project=bundle.data["project"])
        for b in batches:
            b.multiple_batch_id = multiple_batch.pk
            b.created_by = bundle.request.user.username

        multiple_batch.uploaded_data=batches
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
        errors = []
        
        if (".cdx" in correct_file.extension ):
            mols = [mol for mol in readfile( str(correct_file.extension[1:]), str(correct_file.file.name), )]
            rxn = None
            index = 0
            if correct_file.extension == '.cdxml':
                #Look for a stoichiometry table in the reaction file
                rxn = chemdraw_reaction.parse( str(correct_file.file.name))
            
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
                                    b.uncurated_fields = rxn[index]
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
                #read the headers from the first molecule
                
                headers = get_all_sdf_headers(correct_file.file.name)
                data = correct_file.file.read()
                data = data.replace("\r\n","\n").replace("\r","\n")
                ctabs = data.split("$$$$")
                t = time.time()
                for index, mol in enumerate(mols):
                    if mol is None: 
                        errors.append({"index" : index+1, "message" : "Invalid valency or other error parsing this molecule"})
                    else:
                        orig_data = ctabs[index]
                        # try:
                        try:
                            b = CBHCompoundBatch.objects.from_rd_mol(mol, orig_ctab=orig_data, project=bundle.data["project"])
                        except Exception, e:
                            errors.append({"index" : index+1,  "message" : str(e)})
                            b =None
                        if b:
                            custom_fields = {}
                            for hdr in headers:
                                try:
                                    custom_fields[hdr] = mol.GetProp(hdr) 
                                except KeyError:
                                    custom_fields[hdr]   = ""
                            
                            b.uncurated_fields = custom_fields
                            batches.append(b)
               

                print t -time.time()
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
                    struc = Chem.MolFromSmiles(smiles_str)
                    if struc:
                        Compute2DCoords(struc)
                        try:
                            b = CBHCompoundBatch.objects.from_rd_mol(struc, smiles=smiles_str, project=bundle.data["project"], reDraw=True)
                        except Exception, e:
                            errors.append({"index" : index+1,  "message" : str(e)})
                            b =None
                        #work out custom fields from mapping object
                        #new_fields, remapped_fields, ignored_fields
                        if b:
                            custom_fields = {}
                            for hdr in headers:
                                custom_fields[ hdr] = row[hdr] 
                            #Set excel fields as uncurated
                            b.uncurated_fields = custom_fields
                            batches.append(b)
                    else:
                        errors.append({"index" : index+1, "message" : "Invalid valency or other error parsing this identifier",  "SMILES": smiles_str})
            else:
                raise BadRequest("Invalid File Format")

        multiple_batch = CBHCompoundMultipleBatch.objects.create(project = bundle.data["project"])
        for b in batches:
            b.multiple_batch_id = multiple_batch.pk
            b.created_by = bundle.request.user.username
           
        bundle.data["fileerrors"] = errors

        multiple_batch.uploaded_data=batches
        # multiple_batch.save()


        return self.validate_multi_batch(multiple_batch, bundle, request)







    def alter_list_data_to_serialize(self, request, data):
        '''use the request type to determine which fields should be limited for file download,
           add extra fields if needed (eg images) and enumerate the custom fields into the 
           rest of the calculated fields'''
        if(self.determine_format(request) == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' or request.GET.get("format") == "sdf" or self.determine_format(request) == 'chemical/x-mdl-sdfile' ):
            
            df_data = []
            ordered_cust_fields = []
            keys_list = []
            for index, b in enumerate(data["objects"]):
                #remove items which are not listed as being kept
                new_data = {}
                for k, v in b.data.iteritems():
                    for name, display_name in self.Meta.fields_to_keep.iteritems():
                        if k == name:
                            new_data[display_name] = v
                #we need sd format exported results to retain stereochemistry - use mol instaed of smiles
                if(self.determine_format(request) == 'chemical/x-mdl-sdfile' or request.GET.get("format") == "sdf"):
                    new_data['ctab'] = b.data['ctab']
                #dummy
                #not every row has a value for every custom field
                for field, value in b.data['custom_fields'].iteritems():
                    new_data[field] = value
                    keys_list.append(field)
                    
                #now remove custom_fields
                del(new_data['custom_fields'])
                # for field, value in b.data['uncurated_fields'].iteritems():
                #     new_data[field] = value
                #     keys_list.append(field)
                    
                # #now remove custom_fields
                # del(new_data['uncurated_fields'])
                b.data = new_data
                df_data.append(new_data)               
            df = pd.DataFrame(df_data)
            data['export'] = df.to_json()            
            
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
        
        if bundle.obj.created_by:
          #user = User.objects.get(username=bundle.obj.created_by)
            User = get_user_model()
            user = User.objects.get(username=bundle.obj.created_by)
        for names in self.Meta.fieldnames:
            try:
                bundle.data[names[1]] = deepgetattr(data, names[0], None)
            except(AttributeError):
                bundle.data[names[1]] = ""

        mynames = ["editable_by","uncurated_fields", "warnings", "properties", "custom_fields", "errors"]
        for name in mynames:
            bundle.data[name] = json.loads(bundle.data[name])
        #bundle.data["created_by"] = user.__dict__ 
        if user != None:
            bundle.data["created_by"] = user.username
        else:
            bundle.data["created_by"] = ""
        #except:
        #    pass
    
        return bundle

    # def get_list(self, request, **kwargs):
    #     """
    #     Returns a serialized list of resources.
    #     Calls ``obj_get_list`` to provide the data, then handles that result
    #     set and serializes it.
    #     Should return a HttpResponse (200 OK).
    #     """
    #     # TODO: Uncached for now. Invalidation that works for everyone may be
    #     #       impossible.
    #     base_bundle = self.build_bundle(request=request)

    #     ########
    #     #TODO - convert your molfile into smiles or smarts for substructure
    #     #you should interject here, pull the molfile out of the kwargs and do your rdkit conversion etc
    #     #then put back into the request
    #     t1 = time.time()
    #     objects = self.obj_get_list(bundle=base_bundle, **self.remove_api_resource_names(kwargs))
    #     sorted_objects = self.apply_sorting(objects, options=request.GET)
    #     t2 = time.time()
    #     print t2-t1
    #     paginator = self._meta.paginator_class(request.GET, sorted_objects, resource_uri=self.get_resource_uri(), limit=self._meta.limit, max_limit=self._meta.max_limit, collection_name=self._meta.collection_name)
    #     to_be_serialized = paginator.page()
    #     t3 = time.time()
    #     # Dehydrate the bundles in preparation for serialization.
    #     bundles = []

    #     print t3-t2
    #     for obj in to_be_serialized[self._meta.collection_name]:
    #         bundle = self.build_bundle(obj=obj, request=request)
    #         bundles.append(self.full_dehydrate(bundle, for_list=True))
    #     t4 = time.time()

    #     print t4-t3
    #     to_be_serialized[self._meta.collection_name] = bundles
    #     to_be_serialized = self.alter_list_data_to_serialize(request, to_be_serialized)
    #     t5 = time.time()
    #     print t5-t4
    #     return self.create_response(request, to_be_serialized)




    def get_object_list(self, request):
        return super(CBHCompoundBatchResource, self).get_object_list(request).prefetch_related(Prefetch( "related_molregno__compoundproperties"))





def deepgetattr(obj, attr, ex):
    """Recurses through an attribute chain to get the ultimate value."""
    try:
        return reduce(getattr, attr.split('.'), obj)

    except:
        return ex






class CBHCompoundMultipleBatchResource(ModelResource):
    class Meta:
        always_return_data = True
        queryset = FlowFile.objects.all()
        resource_name = 'cbh_batch_upload'
        authorization = Authorization()
        include_resource_uri = False
        allowed_methods = ['get']
        default_format = 'application/json'
        authentication = SessionAuthentication()    



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
            raise BadRequest("No Headers")
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
