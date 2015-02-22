from django.conf.urls import patterns, url, include
from cbh_chembl_ws_extension import api as webservices
from cbh_chembl_ws_extension import Login
from chembl_webservices import __version__ as ws_version
from chembl_webservices import api_name
from chembl_core_db.utils import DirectTemplateView
from django.conf import settings
from flowjs import urls as flow
from django.contrib import admin

admin.autodiscover()


urlpatterns = patterns('',
    url(r'^%s/login' % api_name ,Login.as_view(), name="login"),
    url(r'^%s/flow/' % api_name, include(flow)), #adding this to allow configured upload URL within django-flowjs
    url(r'^%s/admin/' % api_name, include(admin.site.urls)),
	url(r'^grappelli/', include('grappelli.urls')),
)
urlpatterns += webservices.urls