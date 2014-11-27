from tastypie.resources import Resource
from tastypie.serializers import Serializer
from tastypie.serializers import XML_ENCODING
from tastypie.api import Api
from tastypie.exceptions import ImmediateHttpResponse
from django.core.exceptions import ObjectDoesNotExist, MultipleObjectsReturned
from tastypie import http
from chembl_core_db.utils import plural
from tastypie.exceptions import UnsupportedFormat
from tastypie.exceptions import BadRequest
from StringIO import StringIO
from django.core.exceptions import ImproperlyConfigured
import mimeparse
from tastypie.utils.mime import build_content_type
import simplejson
from django.utils import six
from django.http import HttpResponse
from django.http import HttpResponseRedirect

import time
import logging
import django
import urlparse
from django.conf import settings
from django.http import HttpResponseNotFound
from tastypie.exceptions import NotFound

# If ``csrf_exempt`` isn't present, stub it.
try:
    from django.views.decorators.csrf import csrf_exempt
except ImportError:
    def csrf_exempt(func):
        return func

try:
    import defusedxml.lxml as lxml
    from defusedxml.common import DefusedXmlException
    from defusedxml.lxml import parse as parse_xml
    from lxml.etree import Element, tostring, LxmlError, XMLParser
except ImportError:
    lxml = None

try:
    TOP_LEVEL_PAGE = settings.TASTYPIE_TOP_LEVEL_PAGE
except AttributeError:
    TOP_LEVEL_PAGE = 'https://www.ebi.ac.uk/chembl/ws'

try:
    WS_DEBUG = settings.WS_DEBUG
except AttributeError:
    WS_DEBUG = False

from chembl_webservices.base import ChEMBLApiBase


class CBHApiBase(ChEMBLApiBase):


    def __init__(self):
        self.log = logging.getLogger(__name__)
        super(Resource, self).__init__()

    def cached_obj_get_list(self, request=None, **kwargs):
        """
        A version of ``obj_get_list`` that uses the cache as a means to get
        commonly-accessed data faster.
        """
        cache_key = self.generate_cache_key('list', **kwargs)
        get_failed = False
        in_cache = True

        try:
            obj_list = self._meta.cache.get(cache_key)
        except Exception:
            obj_list = None
            get_failed = True
            self.log.error('Caching get exception', exc_info=True, extra={'request': request,})

        if obj_list is None:
            in_cache = False
            obj_list = self.obj_get_list(request=request, **kwargs)
            if not get_failed:
                try:
                    self._meta.cache.set(cache_key, obj_list)
                except Exception:
                    self.log.error('Caching set exception', exc_info=True, extra={'request': request,})

        return obj_list, in_cache

#-----------------------------------------------------------------------------------------------------------------------

    def cached_obj_get(self, request=None, **kwargs):
        """
        A version of ``obj_get`` that uses the cache as a means to get
        commonly-accessed data faster.
        """
        cache_key = self.generate_cache_key('detail', **kwargs)
        get_failed = False
        in_cache = True

        try:
            bundle = self._meta.cache.get(cache_key)
        except Exception:
            bundle = None
            get_failed = True
            self.log.error('Caching get exception', exc_info=True, extra={'request': request,})

        if bundle is None:
            in_cache = False
            bundle = self.obj_get(request=request, **kwargs)
            if not get_failed:
                try:
                    self._meta.cache.set(cache_key, bundle)
                except Exception:
                    self.log.error('Caching set exception', exc_info=True, extra={'request': request,})

        return bundle, in_cache

#-----------------------------------------------------------------------------------------------------------------------
