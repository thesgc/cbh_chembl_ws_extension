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
from cbh_chembl_ws_extension.base import CBHApiBase
from tastypie.utils import dict_strip_unicode_keys



#-----------------------------------------------------------------------------------------------------------------------






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

    def post_list(self, request, **kwargs):
        """
        Creates a new resource/object with the provided data.
        Calls ``obj_create`` with the provided data and returns a response
        with the new resource's location.
        If a new resource is created, return ``HttpCreated`` (201 Created).
        If ``Meta.always_return_data = True``, there will be a populated body
        of serialized data.
        """
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        updated_bundle = self.obj_create(bundle, **self.remove_api_resource_names(kwargs))
        location = self.get_resource_uri(updated_bundle)

        if not self._meta.always_return_data:
            return http.HttpCreated(location=location)
        else:
            updated_bundle = self.full_dehydrate(updated_bundle)
            updated_bundle = self.alter_detail_data_to_serialize(request, updated_bundle)
            return self.create_response(request, updated_bundle, response_class=http.HttpCreated, location=location)



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






