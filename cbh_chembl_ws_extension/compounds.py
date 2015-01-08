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
from tastypie.resources import Resource, ModelResource

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
from cbh_chembl_ws_extension.base import CBHApiBase, CamelCaseJSONSerializer
from tastypie.utils import dict_strip_unicode_keys
from tastypie.serializers import Serializer
from django.core.serializers.json import DjangoJSONEncoder
from tastypie import fields, utils
from cbh_chembl_model_extension.models import CBHCompoundBatch
from tastypie.authentication import SessionAuthentication
import json
from tastypie.paginator import Paginator

from flowjs.models import FlowFile

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
            url(r"^(?P<resource_name>%s)\.(?P<format>json|xml)$" % self._meta.resource_name,
                self.wrap_view('dispatch_compounds'), name="api_dispatch_compounds"),
        ]






class CBHCompoundBatchResource(ModelResource):

    class Meta:
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
        queryset = CBHCompoundBatch.objects.all()
        resource_name = 'cbh_compound_batches'
        authorization = Authorization()
        include_resource_uri = False
        serializer = CamelCaseJSONSerializer()
        allowed_methods = ['get', 'post', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        paginator_class = Paginator


    def post_list_validate(self, request, **kwargs):
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
        return HttpResponse("{ field_names: [ {'name': 'test1', 'count': 1, 'last_used': ''}, {'name': 'test2', 'count': 1, 'last_used': ''} ] }")


    def save_related(self, bundle):
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
                self.wrap_view('post_list_validate'), name="api_validate_compound_batch"),
        ]


    def dehydrate(self, bundle):
        data = bundle.obj.related_molregno
        for names in self.Meta.fieldnames:
            
            bundle.data[names[1]] = deepgetattr(data, names[0], None)
        return bundle


    def get_object_list(self, request):
        return super(CBHCompoundBatchResource, self).get_object_list(request).select_related("related_molregno", "related_molregno__compound_properties")



def deepgetattr(obj, attr, ex):
    """Recurses through an attribute chain to get the ultimate value."""
    try:
        return reduce(getattr, attr.split('.'), obj)

    except:
        print attr
        return ex










class CBHCompoundBatchUpload(ModelResource):

    class Meta:
        always_return_data = True
        #fieldnames = [('files', 'files')]
        queryset = FlowFile.objects.all()
        resource_name = 'cbh_batch_upload'
        authorization = Authorization()
        include_resource_uri = False
        # serializer = CamelCaseJSONSerializer()
        allowed_methods = ['get', 'post', 'put']
        default_format = 'application/json'
        #default_format = 'text/html'
        authentication = SessionAuthentication()
        # paginator_class = Paginator

    def prepend_urls(self):
        return [
        url(r"^(?P<resource_name>%s)/headers/$" % self._meta.resource_name,
                self.wrap_view('dummy_response'), name="api_compound_batch_headers"),
        ]

    def dummy_reponse(self, request, **kwargs):
        #obj = self.queryset[0].name
        
        return self.create_response(request, "obj", response_class=http.HttpAccepted)


    # def get_object_list(self, request):
    #     return super(CBHCompoundBatchUpload, self).get_object_list(request)