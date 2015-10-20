
from django.conf.urls import patterns, url, include
from cbh_core_ws.resources import Login, Logout
from django.conf import settings
from flowjs import urls as flow
from django.contrib import admin


from tastypie.api import Api

from cbh_chembl_ws_extension.compounds import *
from cbh_chembl_ws_extension.projects import *

from cbh_core_ws.resources import *
from django.conf import settings
DEFAULT_API_NAME = 'chemblws'

try:
    api_name = settings.WEBSERVICES_NAME
except AttributeError:
    api_name = DEFAULT_API_NAME


from cbh_chembl_id_generator.resources import CBHPluginResource
api = Api(api_name=api_name)

api.register(CBHCompoundBatchResource())
api.register(CBHCompoundMultipleBatchResource())
api.register(CBHCompoundBatchUpload())
api.register(UserResource())
api.register(ChemregProjectResource())
api.register(SkinningResource())
api.register(CBHPluginResource())

admin.autodiscover()


from django.contrib.auth.decorators import login_required

urlpatterns = patterns('',
                       url(r'^%s/login' % api_name.split("/")
                           [0], Login.as_view(), name="login"),
                       url(r'^%s/logout' %
                           api_name, Logout.as_view(), name="logout"),
                       # adding this to allow configured upload URL within
                       # django-flowjs
                       url(r'^%s/flow/' % api_name, include(flow)),
                       url(r'^%s/admin/' % api_name, include(admin.site.urls)),
                       url(r'^grappelli/', include('grappelli.urls')),
                       url(r'^%s/$' % api_name.split("/")
                           [0], login_required(Index.as_view()))

                       )

urlpatterns += api.urls

if "django_webauth" in settings.INSTALLED_APPS:
    from django_webauth.views import LoginView, LogoutView

    urlpatterns += patterns('',
                            url(r'^%s/webauth' % api_name.split("/")
                                [0], LoginView.as_view(), name="webauthlogin"),
                            url(r'^%s/webauthlogout' % api_name.split("/")
                                [0], LogoutView.as_view(), name="webauthlogout"),
                            )
