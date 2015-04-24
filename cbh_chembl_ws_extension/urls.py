
from django.conf.urls import patterns, url, include
from cbh_chembl_ws_extension.base import Login, Logout
from chembl_core_db.utils import DirectTemplateView
from django.conf import settings
from flowjs import urls as flow
from django.contrib import admin



from tastypie.api import Api

from cbh_chembl_ws_extension.compounds import *
from cbh_chembl_ws_extension.projects import *

from cbh_chembl_ws_extension.base import *
from django.conf import settings
from django_webauth.views import LoginView
DEFAULT_API_NAME='chemblws'

try:
    api_name = settings.WEBSERVICES_NAME
except AttributeError:
    api_name = DEFAULT_API_NAME



api = Api(api_name=api_name)

api.register(CBHCompoundBatchResource())
api.register(CBHCompoundBatchUpload())
api.register(UserResource())
api.register(ProjectResource())
api.register(CustomFieldConfigResource())
api.register(PinnedCustomFieldResource())

admin.autodiscover()


from django.contrib.auth.decorators import login_required

urlpatterns = patterns('',
    url(r'^%s/webauth' % api_name.split("/")[0] ,LoginView.as_view(), name="login1"),
    url(r'^%s/login' % api_name.split("/")[0] ,Login.as_view(), name="login"),
    url(r'^%s/logout' % api_name ,Logout.as_view(), name="logout"),
    url(r'^%s/flow/' % api_name, include(flow)), #adding this to allow configured upload URL within django-flowjs
    url(r'^%s/admin/' % api_name, include(admin.site.urls)),
	url(r'^grappelli/', include('grappelli.urls')),
    url(r'^%s/$' % api_name.split("/")[0], login_required(Index.as_view()) )

)

urlpatterns += api.urls

if "django_webauth" in settings.INSTALLED_APPS:
    urlpatterns += patterns('',
                            url(r'^%s/webauth' % api_name.split("/")[0], include('django_webauth.urls', 'webauth')),
                            )