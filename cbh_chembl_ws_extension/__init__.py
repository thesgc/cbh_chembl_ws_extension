
try:
    __version__ = __import__('pkg_resources').get_distribution('chembl_webservices').version
except Exception as e:
    __version__ = 'development'

from chembl_webservices.base import ChEMBLApi
# from chembl_webservices.status import *
# from chembl_webservices.compounds import *
# from chembl_webservices.assays import *
# from chembl_webservices.targets import *
# from chembl_webservices.bioactivities import *
# from chembl_webservices.drugs import *
from cbh_chembl_ws_extension.compounds import *
from cbh_chembl_ws_extension.projects import *

from cbh_chembl_ws_extension.base import *
from django.conf import settings

DEFAULT_API_NAME='chemblws'

try:
    api_name = settings.WEBSERVICES_NAME
except AttributeError:
    api_name = DEFAULT_API_NAME

api = ChEMBLApi(api_name=api_name)

api.register(CBHCompoundBatchResource())
api.register(CBHCompoundBatchUpload())
api.register(UserResource())
api.register(ProjectResource())
